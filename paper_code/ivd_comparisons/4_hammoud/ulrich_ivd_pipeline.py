#!/usr/bin/env python3
"""Objects -> pseudobulk DE -> support tables, with Ulrich et al. 2024 as the human atlas.

Successor to the HECA-F release of cagri_ivd_end_to_end_20260815. What changes:
  * the human arm is Ulrich 2024 (20 donors, whole-cell) instead of HECA F;
  * both species are fitted with one engine, PyDESeq2, with a fixed block term
    (human: ~ donor + state; mouse: ~ timepoint + group). No limma, no R.
What does not change: the mouse pseudobulk preparation and every frozen check on it,
the arm-aware eligibility rule, the support rule, the consensus-123 -> Ensembl 116
mapping, the 1:1 orthology resource, the Wang lists and their matching, and the
decidual-zone outputs. Those come from the upstream run_pipeline.py, imported
read-only and hash-verified; the upstream repository is never written to.

Run modes
  --record   fit everything, write the release, and write acceptance_observed.yaml with
             every count that config_ulrich.yaml pins. No acceptance checks.
  (default)  same, but every pinned count must reproduce exactly or the run fails.

Usage
  python ulrich_ivd_pipeline.py --config config_ulrich.yaml [--record] [--n-cpus 8]
"""
from __future__ import annotations

import argparse
import atexit
import hashlib
import importlib.util
import json
import platform
import shutil
import sys
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import yaml
from anndata.io import read_elem

HERE = Path(__file__).resolve().parent


def fail(message: str) -> None:
    raise SystemExit(f"FAILED: {message}")


def sha256_file(path: Path, chunk: int = 16 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def load_upstream(config: dict[str, Any]):
    """Import the upstream run_pipeline.py after verifying its hash."""
    spec = config["upstream"]
    path = Path(spec["repo"]) / "run_pipeline.py"
    if sha256_file(path) != spec["run_pipeline_sha256"]:
        fail(f"Upstream run_pipeline.py changed: {path}")
    module_spec = importlib.util.spec_from_file_location("upstream_run_pipeline", path)
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module


# ----------------------------------------------------------------------------------------
# counts

def stream_pseudobulk(h5ad: Path, matrix_path: str, assignments: np.ndarray, n_groups: int,
                      source_to_canonical: np.ndarray, row_chunk: int = 4096) -> tuple[np.ndarray, np.ndarray]:
    """Sum CSR rows into group x canonical-gene counts (upstream stream_canonical_pseudobulk,
    generalized to any CSR group path, e.g. raw/X)."""
    out = np.zeros((n_groups, int(source_to_canonical.max()) + 1), dtype=np.int64)
    n_cells = np.zeros(n_groups, dtype=np.int64)
    kept_sum = 0
    with h5py.File(h5ad, "r") as handle:
        grp = handle[matrix_path]
        enc = grp.attrs.get("encoding-type")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        if enc != "csr_matrix":
            fail(f"/{matrix_path} must be CSR, found {enc}")
        n_rows, n_cols = (int(x) for x in grp.attrs["shape"])
        if len(assignments) != n_rows or len(source_to_canonical) != n_cols:
            fail("Assignment / feature-map dimensions differ from the matrix")
        indptr, indices_ds, data_ds = grp["indptr"], grp["indices"], grp["data"]
        for start in range(0, n_rows, row_chunk):
            stop = min(n_rows, start + row_chunk)
            local = assignments[start:stop]
            if not np.any(local >= 0):
                continue
            ptr = np.asarray(indptr[start: stop + 1], dtype=np.int64)
            off = int(ptr[0]); ptr = ptr - off
            idx = np.asarray(indices_ds[off: off + int(ptr[-1])], dtype=np.int64)
            raw = np.asarray(data_ds[off: off + int(ptr[-1])])
            if raw.dtype.kind == "f" and (not np.isfinite(raw).all() or not np.equal(raw, np.floor(raw)).all()):
                fail("Counts contain non-finite or fractional values")
            if np.any(raw < 0):
                fail("Counts contain negative values")
            data = raw.astype(np.int64, copy=False)
            for r, target in enumerate(local):
                if target < 0:
                    continue
                lo, hi = int(ptr[r]), int(ptr[r + 1])
                cidx = source_to_canonical[idx[lo:hi]]
                keep = cidx >= 0
                if np.any(keep):
                    vals = data[lo:hi][keep]
                    np.add.at(out[target], cidx[keep], vals)
                    kept_sum += int(vals.sum(dtype=np.int64))
                n_cells[target] += 1
    if int(out.sum(dtype=np.int64)) != kept_sum:
        fail("Pseudobulk count conservation failed")
    return out, n_cells


def arm_aware_eligible(counts: pd.DataFrame, is_numerator: np.ndarray, cpm_threshold: float) -> pd.Series:
    """Upstream rule: CPM >= threshold in the majority of the pseudobulks of at least one arm."""
    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    num_major = (int(is_numerator.sum()) + 1) // 2
    den_major = (int((~is_numerator).sum()) + 1) // 2
    return (cpm.loc[is_numerator].ge(cpm_threshold).sum(axis=0).ge(num_major)
            | cpm.loc[~is_numerator].ge(cpm_threshold).sum(axis=0).ge(den_major))


# ----------------------------------------------------------------------------------------
# human: Ulrich 2024

def prepare_ulrich_base(RP, obj: Path, obj_spec: dict, manifest: Path, human: dict, workdir: Path):
    """One pseudobulk per donor x author label, for every label any contrast uses."""
    fields = [human["label_field"], human["donor_field"], *human.get("donor_level_fields", []), human.get("tissue_field", "")]
    fields = [f for f in fields if f]
    with h5py.File(obj, "r") as handle:
        obs = read_elem(handle["obs"])[fields]
        var_names = pd.Index(read_elem(handle[obj_spec["var_path"]]).index.astype(str), name="source_feature_id")
        shape = tuple(int(x) for x in handle[obj_spec["matrix"]].attrs["shape"])
    if shape != tuple(int(x) for x in obj_spec["shape"]):
        fail(f"Ulrich matrix shape changed: {shape}")
    source_to_canonical, features = RP.read_feature_manifest(manifest, var_names)
    labels = sorted({lab for con in human["contrasts"] for lab in (con["labels"] if con.get("kind") == "between_donor" else con["decidual"] + con["reference"])})
    label = obs[human["label_field"]].astype(str)
    donor = obs[human["donor_field"]].astype(str)
    missing = sorted(set(labels) - set(label))
    if missing:
        fail(f"Labels absent from the object: {missing}")
    for f in human.get("donor_level_fields", []):
        if obs.groupby(donor, observed=True)[f].nunique().gt(1).any():
            fail(f"{f} is not constant within donor")
    selected = label.isin(labels).to_numpy()
    keys = sorted(set(zip(donor[selected], label[selected])))
    key_index = {k: i for i, k in enumerate(keys)}
    assignments = np.full(len(obs), -1, dtype=np.int64)
    for pos in np.flatnonzero(selected):
        assignments[pos] = key_index[(donor.iloc[pos], label.iloc[pos])]
    matrix, n_cells = stream_pseudobulk(obj, obj_spec["matrix"], assignments, len(keys), source_to_canonical)
    rows = []
    for i, (d, lab) in enumerate(keys):
        mask = (donor == d) & (label == lab)
        row = {"base_id": f"donor={d}__label={lab}", "donor": d, "label": lab, "n_cells": int(n_cells[i])}
        for f in human.get("donor_level_fields", []):
            row[f] = str(obs.loc[mask, f].iloc[0])
        if human.get("tissue_field"):
            for t, n in obs.loc[mask, human["tissue_field"]].astype(str).value_counts().items():
                row[f"n_cells_{t}"] = int(n)
        rows.append(row)
    base_meta = pd.DataFrame(rows).fillna(0)
    base_counts = pd.DataFrame(matrix, columns=features["ensembl_gene_id"].astype(str), index=base_meta["base_id"])
    if int(base_meta["n_cells"].sum()) != int(selected.sum()):
        fail("Base pseudobulks do not account for every selected cell")
    out = base_counts.copy(); out.insert(0, "base_id", out.index)
    RP.write_tsv(out.reset_index(drop=True), workdir / "ulrich_base_pseudobulk_counts.tsv")
    RP.write_tsv(base_meta, workdir / "ulrich_base_pseudobulk_metadata.tsv")
    return features, base_counts, base_meta


def contrast_pseudobulks(base_counts: pd.DataFrame, base_meta: pd.DataFrame, con: dict, floor: int):
    """Sum base pseudobulks into donor x {decidual, reference}; keep donors with both arms at the floor."""
    arm = pd.Series(pd.NA, index=base_meta.index, dtype="string")
    arm[base_meta["label"].isin(con["decidual"])] = "decidual"
    arm[base_meta["label"].isin(con["reference"])] = "reference"
    bm = base_meta.assign(arm=arm).dropna(subset=["arm"])
    flt = con.get("donor_filter")                      # optional: keep donors whose donor-level field is in `keep` (e.g. luteal-phase donors only)
    if flt:
        if flt["field"] not in bm.columns:
            fail(f"donor_filter field {flt['field']} is not a donor-level field")
        bm = bm[bm[flt["field"]].astype(str).isin([str(v) for v in flt["keep"]])]
        if bm.empty:
            fail(f"donor_filter {flt} leaves no cells")
    meta_rows, mats = [], []
    tissue_cols = [c for c in base_meta.columns if c.startswith("n_cells_")]
    for (d, a), g in bm.groupby(["donor", "arm"], sort=True):
        row = {"pseudobulk_id": f"donor={d}__state={a}", "donor": d, "state": a, "n_cells": int(g.n_cells.sum()),
               "labels": ";".join(f"{l}:{n}" for l, n in zip(g.label, g.n_cells))}
        for c in base_meta.columns:
            if c not in ("base_id", "donor", "label", "n_cells") and not c.startswith("n_cells_"):
                row[c] = g[c].iloc[0]
        for c in tissue_cols:
            row[c] = int(g[c].sum())
        meta_rows.append(row)
        mats.append(base_counts.loc[g.base_id].sum(axis=0).to_numpy())
    meta = pd.DataFrame(meta_rows)
    counts = pd.DataFrame(np.vstack(mats), columns=base_counts.columns, index=meta.index).astype(np.int64)
    keep = meta["n_cells"].ge(floor)
    meta, counts = meta[keep], counts[keep]
    both = meta.groupby("donor")["state"].nunique().eq(2)
    keep = meta["donor"].isin(both[both].index)
    meta, counts = meta[keep].reset_index(drop=True), counts[keep].reset_index(drop=True)
    meta["library_size"] = counts.sum(axis=1).to_numpy()
    return meta, counts


def donor_pseudobulks(base_counts: pd.DataFrame, base_meta: pd.DataFrame, con: dict, floor: int):
    """One pseudobulk per donor from the listed labels; factor and covariate are donor-level fields.
    Donors whose factor value is neither the numerator nor the denominator are excluded."""
    bm = base_meta[base_meta["label"].isin(con["labels"])]
    rows, mats = [], []
    for d, g in bm.groupby("donor", sort=True):
        row = {"pseudobulk_id": f"donor={d}", "donor": d, "n_cells": int(g.n_cells.sum()), "labels": ";".join(f"{l}:{n}" for l, n in zip(g.label, g.n_cells))}
        for c in base_meta.columns:
            if c not in ("base_id", "donor", "label", "n_cells") and not c.startswith("n_cells_"):
                row[c] = g[c].iloc[0]
        rows.append(row); mats.append(base_counts.loc[g.base_id].sum(axis=0).to_numpy())
    meta = pd.DataFrame(rows); counts = pd.DataFrame(np.vstack(mats), columns=base_counts.columns, index=meta.index).astype(np.int64)
    keep = meta["n_cells"].ge(floor) & meta[con["factor_field"]].astype(str).isin([con["numerator"], con["denominator"]])
    meta, counts = meta[keep].reset_index(drop=True), counts[keep].reset_index(drop=True)
    meta["state"] = meta[con["factor_field"]].astype(str)
    meta["block"] = meta[con["covariate_field"]].astype(str).str.replace(r"\W", "_", regex=True)
    meta["library_size"] = counts.sum(axis=1).to_numpy()
    return meta, counts


# ----------------------------------------------------------------------------------------
# one engine for both species

def fit_deseq2(counts: pd.DataFrame, meta: pd.DataFrame, *, block: str, factor: str, numerator: str,
               denominator: str, n_cpus: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    """PyDESeq2 Wald test of factor numerator vs denominator under ~ block + factor.
    Returns per-gene statistics (unshrunk and apeglm-shrunk LFC, numerator over denominator) and fit diagnostics."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    md = pd.DataFrame({block: meta[block].astype(str).to_numpy(), factor: meta[factor].astype(str).to_numpy()},
                      index=meta["pseudobulk_id"].astype(str).to_numpy())
    cts = counts.copy(); cts.index = md.index
    design = f"~ {block} + {factor}"
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(counts=cts, metadata=md, design=design, refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    dm = dds.obsm["design_matrix"]
    rank = int(np.linalg.matrix_rank(dm.to_numpy(dtype=float)))
    if rank != dm.shape[1]:
        fail(f"Design {design} is rank-deficient: rank {rank} < {dm.shape[1]} columns")
    # independent_filter=False: the prespecified arm-aware eligibility rule is the only gene filter, so padj is BH over
    # exactly the eligible genes (what padj meant in the limma release). Cook's outlier filtering stays at the DESeq2 default.
    stats = DeseqStats(dds, contrast=[factor, numerator, denominator], independent_filter=False, inference=inference, quiet=True)
    stats.summary()
    res = stats.results_df.copy()
    coef = [c for c in dm.columns if c.startswith(f"{factor}[T.")]
    if len(coef) != 1:
        fail(f"Expected one {factor} coefficient, found {coef}")
    sign = 1.0 if coef[0] == f"{factor}[T.{numerator}]" else -1.0
    if coef[0] not in (f"{factor}[T.{numerator}]", f"{factor}[T.{denominator}]"):
        fail(f"Unexpected coefficient {coef[0]}")
    stats.lfc_shrink(coeff=coef[0])
    shrunk = stats.results_df["log2FoldChange"] * sign
    # the contrast result is numerator-over-denominator by construction; the shrunk coefficient must agree in sign
    # wherever both are clearly nonzero
    clear = res["log2FoldChange"].abs().gt(1) & shrunk.abs().gt(0.5)
    if clear.any() and (np.sign(res.loc[clear, "log2FoldChange"]) != np.sign(shrunk[clear])).mean() > 0.01:
        fail("Shrunk and unshrunk fold changes disagree in sign; coefficient orientation is wrong")
    out = pd.DataFrame({
        "ensembl_gene_id": res.index.astype(str),
        "baseMean": res["baseMean"].to_numpy(),
        "log2FoldChange_unshrunk": res["log2FoldChange"].to_numpy(),
        "log2FoldChange_shrunk": shrunk.reindex(res.index).to_numpy(),
        "lfcSE": res["lfcSE"].to_numpy(),
        "stat": res["stat"].to_numpy(),
        "pvalue": res["pvalue"].to_numpy(),
        "padj": res["padj"].to_numpy(),
    })
    out["estimable"] = np.isfinite(out["log2FoldChange_unshrunk"]) & np.isfinite(out["pvalue"])
    sf = pd.Series(np.asarray(dds.obsm["size_factors"] if "size_factors" in dds.obsm else dds.obs["size_factors"]), index=md.index)
    diag = {"design": design, "design_columns": list(dm.columns), "design_rank": rank,
            "residual_df": int(dm.shape[0] - rank), "coefficient": coef[0], "n_pseudobulks": int(len(md)),
            "n_genes_fitted": int(len(out)), "size_factor_min": float(sf.min()), "size_factor_median": float(sf.median()),
            "size_factor_max": float(sf.max()), "n_padj_na": int(out["padj"].isna().sum()),
            "n_not_estimable": int((~out["estimable"]).sum())}
    return out, diag


def assemble_result(carrier_features: pd.DataFrame, counts: pd.DataFrame, eligible: pd.Series, fitted: pd.DataFrame,
                    meta: pd.DataFrame, factor: str, numerator: str, *, source_id: str, positive_meaning: str,
                    calls: dict) -> pd.DataFrame:
    carrier = pd.concat([carrier_features.reset_index(drop=True), pd.DataFrame({
        "carrier": True, "production_eligible": eligible.to_numpy(),
        "total_count": counts.sum(axis=0).to_numpy(dtype=np.int64),
        "nonzero_pseudobulks": counts.gt(0).sum(axis=0).to_numpy(dtype=np.int64)})], axis=1)
    if fitted["ensembl_gene_id"].tolist() != carrier.loc[eligible.to_numpy(), "ensembl_gene_id"].astype(str).tolist():
        fail(f"{source_id}: fitted gene order differs from the eligible input")
    result = carrier.merge(fitted, on="ensembl_gene_id", how="left", validate="one_to_one")
    if (~result["production_eligible"] & result["log2FoldChange_unshrunk"].notna()).any():
        fail(f"{source_id}: statistics returned for ineligible genes")
    result["estimable"] = result["production_eligible"] & result["estimable"].astype("boolean").fillna(False).astype(bool)
    stat_cols = ["log2FoldChange_unshrunk", "log2FoldChange_shrunk", "lfcSE", "stat", "pvalue", "padj"]
    result.loc[~result["estimable"], stat_cols] = np.nan
    result["padj_for_call"] = result["padj"].where(np.isfinite(result["padj"]), 1.0)
    result["significant_fdr05_lfc1"] = (result["estimable"] & result["padj_for_call"].lt(float(calls["padj_lt"]))
                                        & result["log2FoldChange_unshrunk"].abs().gt(float(calls["abs_unshrunk_lfc_gt"])))
    result["source_id"] = source_id
    result["positive_lfc_meaning"] = positive_meaning
    result["analysis_method"] = "pydeseq2_wald_fixed_block"
    result["n_pseudobulks_numerator"] = int((meta[factor].astype(str) == numerator).sum())
    result["n_pseudobulks_denominator"] = int((meta[factor].astype(str) != numerator).sum())
    return result


# ----------------------------------------------------------------------------------------
# release bookkeeping

def artifact_inventory(root: Path):
    files = sorted(p for p in root.rglob("*") if p.is_file() and p.name != "release_manifest.json")
    entries, lines = [], []
    for p in files:
        rel, dg, sz = p.relative_to(root).as_posix(), sha256_file(p), p.stat().st_size
        entries.append({"relative_path": rel, "sha256": dg, "size_bytes": sz})
        lines.append(f"{rel}\t{dg}\t{sz}\n")
    return entries, hashlib.sha256("".join(lines).encode()).hexdigest()


def write_manifest(root: Path, config: dict, objects: dict[str, Path], execution_hashes: dict[str, str]) -> None:
    entries, tree = artifact_inventory(root)
    manifest = {
        "schema_version": "1.0", "release_status": "complete_immutable", "analysis_id": config["analysis_id"],
        "created_utc": pd.Timestamp.utcnow().isoformat(),
        "artifact_tree_algorithm": "sha256(sorted(relative_path\\tsha256\\tsize_bytes\\n)); manifest excluded",
        "artifact_tree_sha256": tree, **{f"{k}_sha256": v for k, v in execution_hashes.items()},
        "object_hashes_verified": True,
        "objects": {k: {"path": str(p), "sha256": config["objects"][k]["sha256"], "hash_verified": True} for k, p in objects.items()},
        "resource_inputs": config["resources"], "files": entries}
    (root / "release_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    again, tree2 = artifact_inventory(root)
    if again != entries or tree2 != tree:
        fail("Release inventory changed while writing the manifest")


def check_equal(checks: list, label: str, observed: Any, expected: Any, record: bool) -> None:
    ok = expected is not None and observed == expected
    if not record and not ok:
        fail(f"{label} changed: expected {expected}, observed {observed}")
    checks.append({"check": label, "observed": json.dumps(observed, default=str), "expected": json.dumps(expected, default=str),
                   "pass": bool(ok) if not record else "recorded"})


# ----------------------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=HERE / "config_ulrich.yaml")
    ap.add_argument("--output-dir", type=Path, default=None)
    ap.add_argument("--n-cpus", type=int, default=8)
    ap.add_argument("--record", action="store_true", help="write observed acceptance values instead of checking them")
    args = ap.parse_args()

    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    if config.get("schema_version") != "1.0" or config.get("analysis_id") != "cagri_ivd_ulrich_deseq2_20261019":
        fail("Unexpected config schema/analysis_id")
    if platform.python_version() != str(config["runtime"]["python"]):
        fail(f"Expected Python {config['runtime']['python']}, observed {platform.python_version()}")
    import pydeseq2
    if pydeseq2.__version__ != str(config["runtime"]["pydeseq2"]):
        fail(f"Expected PyDESeq2 {config['runtime']['pydeseq2']}, observed {pydeseq2.__version__}")
    RP = load_upstream(config)
    record = bool(args.record)
    acc = config.get("acceptance") or {}

    execution_inputs = {"config": config_path, "pipeline": Path(__file__).resolve(),
                        "upstream_pipeline": Path(config["upstream"]["repo"]) / "run_pipeline.py"}
    exec_hash = {k: sha256_file(p) for k, p in execution_inputs.items()}
    res = {k: RP.resolve_resource(HERE, spec) for k, spec in config["resources"].items()}
    ulrich_obj = Path(config["objects"]["ulrich_cells"]["path"])
    xmens_obj = Path(config["objects"]["xmens_mouse"]["path"])
    for key, path in (("ulrich_cells", ulrich_obj), ("xmens_mouse", xmens_obj)):
        if sha256_file(path) != config["objects"][key]["sha256"]:
            fail(f"Object hash mismatch: {path}")
    RP.verify_object(xmens_obj, config["objects"]["xmens_mouse"], check_hash=False)

    output = (args.output_dir or Path(config["output_directory"])).resolve()
    if record:
        output = output.with_name(output.name + "_record")
    staging = output.with_name(output.name + ".inprogress")
    if output.exists() or staging.exists():
        fail(f"Refusing existing output/staging path: {output}")
    staging.mkdir(parents=True)
    atexit.register(lambda: staging.exists() and shutil.rmtree(staging, ignore_errors=True))
    inter = staging / "intermediate"; inter.mkdir()
    checks: list[dict[str, Any]] = []
    calls = config["support_calls"]
    observed_acc: dict[str, Any] = {"human": {}, "mouse": {}, "status_partitions": {}}

    # ---- human
    human = config["human"]
    features_h, base_counts, base_meta = prepare_ulrich_base(RP, ulrich_obj, config["objects"]["ulrich_cells"], res["ulrich_feature_manifest"], human, inter)
    human_results: dict[str, pd.DataFrame] = {}
    human_diag: dict[str, Any] = {}
    for con in human["contrasts"]:
        cid = con["id"]
        between = con.get("kind") == "between_donor"
        if between:
            meta, counts = donor_pseudobulks(base_counts, base_meta, con, int(human["cell_floor"]))
            num, den, block = con["numerator"], con["denominator"], "block"
            if meta.groupby("block")["state"].nunique().lt(2).all():
                fail(f"{cid}: factor and covariate are confounded")
        else:
            meta, counts = contrast_pseudobulks(base_counts, base_meta, con, int(con.get("cell_floor", human["cell_floor"])))   # optional per-contrast floor
            num, den, block = "decidual", "reference", "donor"
        is_num = meta["state"].eq(num).to_numpy()
        eligible = arm_aware_eligible(counts, is_num, float(human["eligibility_cpm_threshold"]))
        fitted, diag = fit_deseq2(counts.loc[:, eligible.to_numpy()], meta, block=block, factor="state",
                                  numerator=num, denominator=den, n_cpus=args.n_cpus)
        result = assemble_result(features_h, counts, eligible, fitted, meta, "state", num,
                                 source_id=f"ulrich_{cid}", positive_meaning=f"{num}_minus_{den}", calls=calls)
        obs_c = {"pseudobulks": int(len(meta)), "donors": int(meta.donor.nunique()), "eligible_genes": int(eligible.sum()),
                 "count_sum": int(counts.to_numpy().sum()), "strict_calls": int(result["significant_fdr05_lfc1"].sum()),
                 "residual_df": diag["residual_df"]}
        for k, v in obs_c.items():
            check_equal(checks, f"human_{cid}_{k}", v, (acc.get("human", {}).get(cid) or {}).get(k), record)
        observed_acc["human"][cid] = obs_c
        human_diag[cid] = {**diag, **obs_c, "label": con["label"], "kind": "between_donor" if between else "within_donor",
                           "numerator": con["labels"] if between else con["decidual"], "denominator": con["denominator"] if between else con["reference"],
                           "numerator_value": num, "denominator_value": den, "unit": "donors" if between else "donor pairs"}
        out_meta = meta.copy(); RP.write_tsv(out_meta, inter / f"ulrich_{cid}_metadata.tsv")
        out_c = counts.copy(); out_c.insert(0, "pseudobulk_id", meta["pseudobulk_id"]); RP.write_tsv(out_c, inter / f"ulrich_{cid}_counts.tsv")
        RP.write_csv(result, staging / f"ulrich_{cid}_all_genes.csv")
        result.to_parquet(staging / f"ulrich_{cid}_all_genes.parquet", index=False)
        human_results[cid] = result

    # ---- mouse: upstream preparation (all its frozen cell-count checks), PyDESeq2 fit
    xm = config["xmens_broad"]
    x_carrier, x_counts, x_meta, x_eligible = RP.prepare_xmens_broad(xmens_obj, config["objects"]["xmens_mouse"], res["xmens_feature_manifest"], xm, inter)
    x_meta = x_meta.assign(group=x_meta["group"].astype(str), timepoint=x_meta["timepoint"].astype(str))
    x_fit, x_diag = fit_deseq2(x_counts.loc[:, x_eligible.to_numpy()].set_axis(x_carrier.loc[x_eligible.to_numpy(), "ensembl_gene_id"].astype(str), axis=1),
                               x_meta, block="timepoint", factor="group", numerator="decidual", denominator="fibroblast", n_cpus=args.n_cpus)
    x_features = x_carrier[["ensembl_gene_id", "gene_symbol", "biotype", "source_feature_ids", "source_feature_n"]]
    xmens = assemble_result(x_features, x_counts, x_eligible, x_fit, x_meta, "group", "decidual",
                            source_id="mouse_all_decidual_vs_fibroblast_timepoint_blocked", positive_meaning="decidual_minus_fibroblast", calls=calls)
    obs_m = {"eligible_genes": int(x_eligible.sum()), "strict_calls": int(xmens["significant_fdr05_lfc1"].sum()), "residual_df": x_diag["residual_df"]}
    for k, v in obs_m.items():
        check_equal(checks, f"mouse_{k}", v, (acc.get("mouse") or {}).get(k), record)
    observed_acc["mouse"] = obs_m
    if int(x_eligible.sum()) != int(xm["expected_eligible_genes"]):
        fail("X-Mens eligible gene count differs from the upstream release")      # same rule, same counts
    RP.write_csv(xmens, staging / "xmens_broad_all_genes.csv"); xmens.to_parquet(staging / "xmens_broad_all_genes.parquet", index=False)

    # ---- support tables (upstream builders)
    primary = human["primary"]
    support_h = {}
    for cid, result in human_results.items():
        s = RP.build_heca_support(res["consensus123"], res["consensus123_human_mapping"], result, calls)
        s["source_id"] = f"ulrich_{cid}"
        support_h[cid] = s
        RP.write_csv(s, staging / ("consensus123_ulrich_support_rows.csv" if cid == primary else f"consensus123_ulrich_{cid}_support_rows.csv"))
    ortho = RP.build_consensus_ortholog_support(res["consensus123"], res["consensus123_mouse_orthology"], xmens, calls, 114)
    wang = RP.build_wang_support(res["wang_ivd"], res["wang_oil"], res["wang_pregnancy"], xmens, calls)
    RP.write_csv(ortho, staging / "consensus_ortholog_xmens_broad_support_rows.csv")
    # Wang -> X-Mens matching is fit-independent except for candidate ranking by estimability: matched/unmatched/alias
    # counts and the Nrg1 collision must equal the upstream release; estimable and strict-called counts are recorded.
    wm = {"wang_matched": int(wang["xmens_match_candidate_n"].gt(0).sum()), "wang_unmatched": int(wang["xmens_match_candidate_n"].eq(0).sum()),
          "wang_alias_route": int(wang["xmens_match_route"].eq("source_feature_id_alias").sum())}
    for k, v in wm.items():
        if v != int(xm[f"expected_{k}_n"]):
            fail(f"{k} changed from the upstream release: {v} vs {xm[f'expected_{k}_n']}")
        checks.append({"check": k, "observed": v, "expected": int(xm[f"expected_{k}_n"]), "pass": True})
    nrg1 = wang.set_index("wang_ivd_gene").loc["Nrg1"]
    if int(nrg1["xmens_match_candidate_n"]) != 2 or str(nrg1["xmens_selected_ensembl_gene_id"]) != "ENSMUSG00000062991":
        fail("Nrg1 collision resolution changed")
    for k, v in {"wang_estimable": int(wang["xmens_estimable"].map(RP.bool_value).sum()),
                 "wang_strict_called": int(wang["xmens_called_strict_lfc1"].map(RP.bool_value).sum())}.items():
        check_equal(checks, f"mouse_{k}", v, (acc.get("mouse") or {}).get(k), record)
        observed_acc["mouse"][k] = v
    RP.write_csv(wang, staging / "wang_ivd_support_rows.csv")

    parts = {f"consensus123_vs_ulrich_{cid}": (s["support_status"], ["supported", "reversed", "tested_not_called", "ineligible_or_unavailable"]) for cid, s in support_h.items()}
    parts.update({
        "consensus_ortholog_vs_xmens_broad": (ortho["support_status"], ["supported", "reversed", "tested_not_called", "ineligible_or_unavailable"]),
        "wang_ivd_vs_oil": (wang["wang_oil_status"], ["supported", "reversed", "not_in_published_deg_list", "ineligible_or_unmatched"]),
        "wang_ivd_vs_pregnancy": (wang["wang_pregnancy_status"], ["supported", "reversed", "not_in_published_deg_list", "ineligible_or_unmatched"]),
        "wang_ivd_vs_xmens_broad": (wang["xmens_status"], ["supported", "reversed", "tested_not_called", "ineligible_or_unmatched"]),
    })
    for name, (series, cats) in parts.items():
        observed = RP.normalized_partition(series, cats)
        if set(series.astype(str)) - set(cats):
            fail(f"{name}: unexpected status values")
        check_equal(checks, f"partition_{name}", observed, (acc.get("status_partitions") or {}).get(name), record)
        observed_acc["status_partitions"][name] = observed
    # direction-only Wang comparisons do not depend on any fit: they must equal the upstream release
    up_parts = config["upstream"]["unchanged_partitions"]
    for name in ("wang_ivd_vs_oil", "wang_ivd_vs_pregnancy"):
        if observed_acc["status_partitions"][name] != up_parts[name]:
            fail(f"{name} differs from the upstream release; it should not depend on the fit")

    summary_rows = [RP.comparison_summary_row(f"consensus123_vs_ulrich_{cid}", f"123 human consensus IVD genes vs Ulrich {human_diag[cid]['label']}",
                                              s["support_status"], "ineligible_or_unavailable") for cid, s in support_h.items()]
    summary_rows += [
        RP.comparison_summary_row("wang_ivd_vs_oil", "Wang IVD genes vs Wang oil/artificial in vivo", wang["wang_oil_status"], "ineligible_or_unmatched"),
        RP.comparison_summary_row("wang_ivd_vs_pregnancy", "Wang IVD genes vs Wang pregnancy in vivo", wang["wang_pregnancy_status"], "ineligible_or_unmatched"),
        RP.comparison_summary_row("wang_ivd_vs_xmens_broad", "Wang IVD genes vs broad X-Mens", wang["xmens_status"], "ineligible_or_unmatched"),
        RP.comparison_summary_row("consensus_ortholog_vs_xmens_broad", "Consensus IVD genes with a strict 1:1 mouse ortholog vs broad X-Mens", ortho["support_status"], "ineligible_or_unavailable"),
    ]
    summary = pd.DataFrame(summary_rows)
    summary.insert(1, "is_primary_human", summary["comparison_id"].eq(f"consensus123_vs_ulrich_{primary}"))
    RP.write_csv(summary, staging / "ivd_support_summary.csv")
    contrasts = pd.DataFrame([{"contrast_id": cid, "primary": cid == primary, **{k: (";".join(v) if isinstance(v, list) else v) for k, v in d.items()
                                                                                 if k not in ("design_columns",)}} for cid, d in human_diag.items()])
    RP.write_csv(contrasts, staging / "ulrich_contrasts.csv")

    # ---- decidual-zone outputs (upstream, unchanged)
    cell_table, _, _ = RP.build_zone_outputs(xmens_obj, config["objects"]["xmens_mouse"], res["zone_mapping"], res["zone_markers"], config["zones"], staging)
    RP.write_zone_derivative(xmens_obj, cell_table, res["zone_mapping"], res["zone_markers"], staging / "final_mouse_with_putative_decidual_zone.h5ad")

    RP.write_csv(pd.DataFrame(checks), staging / "validation_checks.csv")
    if {k: sha256_file(p) for k, p in execution_inputs.items()} != exec_hash:
        fail("Config or code changed during the run")
    run_meta = {"analysis_id": config["analysis_id"], "created_utc": pd.Timestamp.utcnow().isoformat(), "mode": "record" if record else "verify",
                "python_executable": sys.executable, "python_version": platform.python_version(), "packages": RP.package_versions(),
                "n_cpus": args.n_cpus, "execution_sha256": exec_hash, "primary_human_contrast": primary,
                "human_fits": human_diag, "mouse_fit": {**x_diag, **obs_m},
                "support_rule": "eligible and estimable and padj<0.05 and abs(unshrunk log2FC)>1; supported requires same direction"}
    (staging / "run_metadata.json").write_text(json.dumps(run_meta, indent=2, sort_keys=True, default=str) + "\n")
    if record:
        (staging / "acceptance_observed.yaml").write_text(yaml.safe_dump(observed_acc, sort_keys=False))
    for key, path in (("ulrich_cells", ulrich_obj), ("xmens_mouse", xmens_obj)):
        if sha256_file(path) != config["objects"][key]["sha256"]:
            fail(f"Source object changed during the run: {path}")
    write_manifest(staging, config, {"ulrich_cells": ulrich_obj, "xmens_mouse": xmens_obj}, exec_hash)
    staging.rename(output)
    print(json.dumps({"status": "complete_immutable", "mode": run_meta["mode"], "output": str(output),
                      "summary": summary[["comparison_id", "supported_n", "reversed_n", "denominator_n"]].to_dict("records")}, indent=1))


if __name__ == "__main__":
    main()
