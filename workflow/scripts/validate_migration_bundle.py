#!/usr/bin/env python3
"""Validate one species' complete staged annotation-migration bundle."""
import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml


sm = globals()["snakemake"]
species = sm.params.species
paths = {name: Path(getattr(sm.input, name)) for name in
         ("refined", "final", "fibdec", "deposit_all", "deposit_fibdec")}
vocab = yaml.safe_load(Path(sm.input.vocab).read_text())
report_path = Path(sm.output.report)
counts_path = Path(sm.output.counts)
FIBDEC = {"decidual_cell", "endometrial_fibroblast"}
FIBDEC_N = {"human": 47533, "mouse": 43197}
LEGACY_OBS = {"coarse_refined", "fine_refined", "subtype", "annotation_level",
              "cell_type_annotations", "cell_subtype_annotations",
              "cell_type_subannotations"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def strings_digest(values) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode())
        digest.update(b"\0")
    return digest.hexdigest()


def h5_numeric_tree_digest(path: Path, key: str) -> str:
    digest = hashlib.sha256()
    with h5py.File(path, "r") as handle:
        root = handle[key]
        datasets = []
        root.visititems(lambda name, obj: datasets.append((name, obj))
                        if isinstance(obj, h5py.Dataset) else None)
        for name, dataset in sorted(datasets):
            digest.update(name.encode())
            digest.update(str(dataset.shape).encode())
            digest.update(str(dataset.dtype).encode())
            if dataset.ndim == 0:
                digest.update(np.asarray(dataset[()]).tobytes())
            else:
                step = max(1, min(dataset.shape[0], 1_000_000))
                for start in range(0, dataset.shape[0], step):
                    digest.update(np.asarray(dataset[start:start + step]).tobytes())
    return digest.hexdigest()


def gate(rows, check, expected, actual):
    rows.append({"species": species, "check": check, "expected": expected,
                 "actual": actual, "delta": actual - expected,
                 "status": "PASS" if actual == expected else "SUSPECT"})


failures = []
suspects = []
gates = []
objects = {name: ad.read_h5ad(path, backed="r") for name, path in paths.items()}
refined, final, fibdec = objects["refined"], objects["final"], objects["fibdec"]

for name, obj in objects.items():
    if obj.obs_names.has_duplicates:
        failures.append(f"{name}: duplicate obs_names")
    if obj.var_names.has_duplicates:
        failures.append(f"{name}: duplicate var_names")
    for column in ("cell_type", "cell_subtype"):
        if column not in obj.obs:
            failures.append(f"{name}: missing {column}")
    legacy = sorted(LEGACY_OBS & set(obj.obs.columns))
    if legacy:
        failures.append(f"{name}: legacy annotation columns remain {legacy}")
    if "counts" not in obj.layers:
        failures.append(f"{name}: missing counts layer")
    if {"cell_type", "cell_subtype"} <= set(obj.obs.columns):
        cell_types = obj.obs.cell_type.astype(str)
        cell_subtypes = obj.obs.cell_subtype.astype(str)
        allowed_types = set(vocab["expected_coarse_counts"][species])
        invalid_types = sorted(set(cell_types) - allowed_types)
        if invalid_types:
            failures.append(f"{name}: invalid {species} cell_type values {invalid_types}")
        subtype_parent = {
            subtype: parent
            for parent, subtypes in vocab["allowed"].items()
            for subtype in (subtypes or [])
        }
        invalid_subtypes = sorted(set(cell_subtypes) - {""} - set(subtype_parent))
        if invalid_subtypes:
            failures.append(f"{name}: invalid cell_subtype values {invalid_subtypes}")
        invalid_pairs = sorted({
            (parent, subtype)
            for parent, subtype in zip(cell_types, cell_subtypes)
            if subtype and subtype_parent.get(subtype) != parent
        })
        if invalid_pairs:
            failures.append(f"{name}: invalid subtype-parent pairs {invalid_pairs}")

# A malformed object must still yield a machine-readable FAIL result rather
# than raising later while accessing a missing annotation/counts field.
critical = any(
    column not in obj.obs or "counts" not in obj.layers
    for obj in objects.values()
    for column in ("cell_type", "cell_subtype")
)
if critical:
    payload = {
        "species": species, "schema_version": vocab["schema_version"],
        "status": "FAIL", "failures": failures, "suspects": [], "gates": [],
        "paths": {name: str(path) for name, path in paths.items()},
        "shapes": {name: [obj.n_obs, obj.n_vars] for name, obj in objects.items()},
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    counts_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=["species", "check", "expected", "actual", "delta", "status"]).to_csv(
        counts_path, sep="\t", index=False)
    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(0)

if not refined.obs_names.equals(final.obs_names):
    failures.append("refined/final obs_names or order differ")
if not refined.var_names.equals(final.var_names):
    failures.append("refined/final var_names or order differ")
for column in ("cell_type", "cell_subtype", "origin_lineage"):
    if column in refined.obs and column in final.obs:
        if not refined.obs[column].astype(str).equals(final.obs[column].astype(str)):
            failures.append(f"refined/final {column} differs")

expected_types = vocab["expected_coarse_counts"][species]
type_counts = final.obs.cell_type.astype(str).value_counts() if "cell_type" in final.obs else pd.Series()
for label, expected in expected_types.items():
    gate(gates, f"cell_type:{label}", int(expected), int(type_counts.get(label, 0)))
gate(gates, "total_cells", int(vocab["expected_coarse_counts"]["totals"][species]),
     final.n_obs)
gate(gates, "fibdec_total", FIBDEC_N[species], fibdec.n_obs)
if "cell_type" in fibdec.obs and set(fibdec.obs.cell_type.astype(str)) != FIBDEC:
    failures.append(f"fibdec types are {sorted(set(fibdec.obs.cell_type.astype(str)))}")
if "cell_type" in final.obs:
    final_fibdec_mask = final.obs.cell_type.astype(str).isin(FIBDEC).values
    expected_fibdec_ids = final.obs_names[final_fibdec_mask]
    if not fibdec.obs_names.equals(expected_fibdec_ids):
        failures.append("fibdec obs_names/order differs from exact final fib/dec mask")
    if not fibdec.var_names.equals(final.var_names):
        failures.append("fibdec var_names/order differs from final")

if species == "human" and "cell_subtype" in final.obs:
    sub = final.obs.cell_subtype.astype(str)
    gate(gates, "endothelial_blood", 828, int((sub == "endothelial_blood").sum()))
    gate(gates, "endothelial_lymph", 506, int((sub == "endothelial_lymph").sum()))
    gate(gates, "perivascular_blank", 471,
         int(((final.obs.cell_type.astype(str) == "perivascular_cell") & (sub == "")).sum()))
    gate(gates, "nonciliated_epithelial_blank", 10608,
         int(((final.obs.cell_type.astype(str) == "epithelial_cell") & (sub == "")).sum()))
    gate(gates, "smooth_muscle_blank", 2084,
         int(((final.obs.cell_type.astype(str) == "smooth_muscle_cell") & (sub == "")).sum()))
else:
    sub = final.obs.cell_subtype.astype(str) if "cell_subtype" in final.obs else pd.Series()
    gate(gates, "endothelial_blood", 3487, int((sub == "endothelial_blood").sum()))
    gate(gates, "endothelial_lymph", 0, int((sub == "endothelial_lymph").sum()))
    if "cell_type" in final.obs:
        for label, expected in (("epithelial_cell", 6265),
                                ("smooth_muscle_cell", 1443),
                                ("mesothelial_cell", 550)):
            gate(gates, f"{label}_blank", expected,
                 int(((final.obs.cell_type.astype(str) == label) & (sub == "")).sum()))

for name in ("refined", "final", "fibdec", "deposit_all", "deposit_fibdec"):
    obj = objects[name]
    if obj.uns.get("annotation_schema_version") != vocab["schema_version"]:
        failures.append(f"{name}: wrong annotation_schema_version")
    if obj.uns.get("display_names") != vocab["display_names"]:
        failures.append(f"{name}: display_names mismatch")
    if "cell_subtype" in obj.obs:
        subtype = obj.obs.cell_subtype.astype(str)
        if obj.obs.cell_subtype.isna().any() or set(subtype.str.lower()) & {"nan", "none", "na", "<na>"}:
            failures.append(f"{name}: invalid cell_subtype missing values")

if "X_scVI" not in final.obsm or final.obsm["X_scVI"].shape != (final.n_obs, 30):
    failures.append("final: missing/wrong X_scVI")
elif not np.isfinite(np.asarray(final.obsm["X_scVI"])).all():
    failures.append("final: nonfinite X_scVI")
if "X_umap_scVI" not in final.obsm or final.obsm["X_umap_scVI"].shape != (final.n_obs, 2):
    failures.append("final: missing/wrong X_umap_scVI")
elif not np.isfinite(np.asarray(final.obsm["X_umap_scVI"])).all():
    failures.append("final: nonfinite X_umap_scVI")

digests = {
    "refined_obs": strings_digest(refined.obs_names),
    "final_obs": strings_digest(final.obs_names),
    "refined_var": strings_digest(refined.var_names),
    "final_var": strings_digest(final.var_names),
    "refined_counts": h5_numeric_tree_digest(paths["refined"], "layers/counts"),
    "final_counts": h5_numeric_tree_digest(paths["final"], "layers/counts"),
}
if digests["refined_counts"] != digests["final_counts"]:
    failures.append("refined/final counts-layer digest differs")

# Exact subset provenance: compare every fib/dec counts row to the same final row.
if "cell_type" in final.obs and "counts" in final.layers and "counts" in fibdec.layers:
    selected = np.flatnonzero(final.obs.cell_type.astype(str).isin(FIBDEC).values)
    for start in range(0, len(selected), 1000):
        rows = selected[start:start + 1000]
        left = final.layers["counts"][rows, :]
        right = fibdec.layers["counts"][start:start + len(rows), :]
        if hasattr(left, "toarray"):
            differs = (left != right).nnz
        else:
            differs = int(np.count_nonzero(np.asarray(left) != np.asarray(right)))
        if differs:
            failures.append(f"fibdec counts differ from final mask in rows {start}:{start+len(rows)}")
            break

file_hashes = {name: sha256(path) for name, path in paths.items()}
if file_hashes["final"] != file_hashes["deposit_all"]:
    failures.append("all-cell deposit is not byte-identical to staged final")
if file_hashes["fibdec"] != file_hashes["deposit_fibdec"]:
    failures.append("fib/dec deposit is not byte-identical to staged fib/dec")

suspects.extend([row for row in gates if row["status"] == "SUSPECT"])
status = "FAIL" if failures else ("SUSPECT" if suspects else "PASS")
payload = {
    "species": species, "schema_version": vocab["schema_version"], "status": status,
    "failures": failures, "suspects": suspects, "gates": gates,
    "paths": {name: str(path) for name, path in paths.items()},
    "shapes": {name: [obj.n_obs, obj.n_vars] for name, obj in objects.items()},
    "digests": digests, "sha256": file_hashes,
}

report_path.parent.mkdir(parents=True, exist_ok=True)
tmp_report = report_path.with_name(f".{report_path.name}.{os.getpid()}.tmp")
tmp_report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
os.replace(tmp_report, report_path)
counts_path.parent.mkdir(parents=True, exist_ok=True)
tmp_counts = counts_path.with_name(f".{counts_path.name}.{os.getpid()}.tmp")
pd.DataFrame(gates).to_csv(tmp_counts, sep="\t", index=False)
os.replace(tmp_counts, counts_path)
print(json.dumps(payload, indent=2, sort_keys=True))
