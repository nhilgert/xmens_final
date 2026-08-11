#!/usr/bin/env python3
"""Create staging-only fib/dec and Dataverse artifacts from a staged final atlas."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil

import scanpy as sc
import yaml


sm = globals()["snakemake"]
species = sm.params.species
final_path = Path(sm.input.final)
vocab_path = Path(sm.input.vocab)
fibdec_path = Path(sm.output.fibdec)
deposit_all = Path(sm.output.deposit_all)
deposit_fibdec = Path(sm.output.deposit_fibdec)
manifest_path = Path(sm.output.manifest)
vocab = yaml.safe_load(vocab_path.read_text())
expected = {"human": 47533, "mouse": 43197}[species]
FIBDEC = {"decidual_cell", "endometrial_fibroblast"}


def refuse_existing(*paths):
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite staged artifacts: {existing}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_copy(source: Path, destination: Path):
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    shutil.copyfile(source, tmp)
    os.replace(tmp, destination)


refuse_existing(fibdec_path, deposit_all, deposit_fibdec, manifest_path)
adata = sc.read_h5ad(final_path)
for column in ("cell_type", "cell_subtype"):
    if column not in adata.obs:
        raise ValueError(f"{final_path}: missing {column}")
if "counts" not in adata.layers:
    raise ValueError(f"{final_path}: missing counts layer")
adata.uns["annotation_schema_version"] = vocab["schema_version"]
adata.uns["display_names"] = vocab["display_names"]

mask = adata.obs.cell_type.astype(str).isin(FIBDEC).values
fibdec = adata[mask].copy()
actual_types = set(fibdec.obs.cell_type.astype(str))
status = "PASS"
suspects = []
if actual_types != FIBDEC:
    status = "SUSPECT"
    suspects.append(f"fib/dec types {sorted(actual_types)} != {sorted(FIBDEC)}")
if fibdec.n_obs != expected:
    status = "SUSPECT"
    suspects.append(f"fib/dec count {fibdec.n_obs} != {expected}")

fibdec_path.parent.mkdir(parents=True, exist_ok=True)
tmp_fibdec = fibdec_path.with_name(f".{fibdec_path.name}.{os.getpid()}.tmp.h5ad")
fibdec.write_h5ad(tmp_fibdec)
os.replace(tmp_fibdec, fibdec_path)

# Deposits are independent byte copies so their provenance is unambiguous and
# exact equality can be checked by SHA-256.
atomic_copy(final_path, deposit_all)
atomic_copy(fibdec_path, deposit_fibdec)

manifest = {
    "species": species,
    "status": status,
    "suspects": suspects,
    "schema_version": vocab["schema_version"],
    "artifacts": {
        "final": {"path": str(final_path), "sha256": sha256(final_path)},
        "fibdec": {"path": str(fibdec_path), "shape": list(fibdec.shape),
                   "sha256": sha256(fibdec_path)},
        "deposit_all": {"path": str(deposit_all), "sha256": sha256(deposit_all)},
        "deposit_fibdec": {"path": str(deposit_fibdec),
                           "sha256": sha256(deposit_fibdec)},
    },
}
manifest["checks"] = {
    "deposit_all_matches_final": (
        manifest["artifacts"]["deposit_all"]["sha256"] ==
        manifest["artifacts"]["final"]["sha256"]),
    "deposit_fibdec_matches_fibdec": (
        manifest["artifacts"]["deposit_fibdec"]["sha256"] ==
        manifest["artifacts"]["fibdec"]["sha256"]),
    "fibdec_expected": expected,
    "fibdec_actual": fibdec.n_obs,
}
manifest_path.parent.mkdir(parents=True, exist_ok=True)
tmp_manifest = manifest_path.with_name(f".{manifest_path.name}.{os.getpid()}.tmp")
tmp_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
os.replace(tmp_manifest, manifest_path)
print(json.dumps(manifest, indent=2, sort_keys=True))
