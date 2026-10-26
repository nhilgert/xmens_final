#!/usr/bin/env python
"""Per-nucleus plot table for the decidual-zone figures (rule zone_plot_table, workflow/rules/16_zone_modules.smk).

One CSV holds everything fig_zone_modules.py needs, so every figure can be redrawn from it without the h5ad:
    barcode (index), sample, timepoint, cell_type, cell_subtype, x_um_qc, y_um_qc    obs of the release object
    embedding_outlier                                                               barcode listed in --flag-barcodes
    auc_<MODULE> (11 columns)                                                       the run's zone_aucell.csv, unchanged
Rows are the scored nuclei (the run's timepoints) in object order. x_um_qc / y_um_qc are NaN for nuclei that were never
placed or failed spatial QC; they still count in the heatmaps.

A sidecar <table>.provenance.json records the inputs (md5/sha256) and copies the scoring metadata (window, seed,
normalization, scored timepoints, gene sets) from the run's provenance.json, so the CSV and its sidecar describe
themselves. Written 2026-10-24. Existing outputs are never overwritten.
"""
import argparse, datetime, hashlib, json, os, sys
import numpy as np, pandas as pd, h5py
from anndata.io import read_elem

OBS = ["sample", "timepoint", "cell_type", "cell_subtype", "x_um_qc", "y_um_qc"]
SCORING_KEYS = ["auc_threshold", "rank_cutoff", "seed", "normalization", "timepoints", "n_scored", "n_vars", "genesets", "versions"]


def md5(path):
    h = hashlib.md5(); f = open(path, "rb")
    for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--obj", required=True)                  # release h5ad (obs is read, nothing else)
    ap.add_argument("--scores", required=True)               # the run's zone_aucell.csv
    ap.add_argument("--scores-provenance", required=True)    # the run's provenance.json
    ap.add_argument("--flag-barcodes", required=True)        # CSV with a 'barcode' column (embedding outliers)
    ap.add_argument("--out", required=True)                  # <name>.csv; <name>.provenance.json is written next to it
    a = ap.parse_args()
    assert a.out.endswith(".csv"), a.out
    side = a.out[:-len(".csv")] + ".provenance.json"
    clash = [f for f in (a.out, side) if os.path.exists(f)]
    if clash:
        sys.exit(f"refusing to overwrite {clash}")

    S = pd.read_csv(a.scores, index_col=0); S.index = S.index.astype(str)
    assert S.index.is_unique and all(c.startswith("auc_") for c in S.columns), "unexpected score table"
    with h5py.File(a.obj, "r") as h:
        idx = pd.Index(np.asarray(read_elem(h["obs/_index"])).astype(str))
        obs = pd.DataFrame({c: np.asarray(read_elem(h["obs/" + c])) for c in OBS}, index=idx)
    dtypes = {c: str(obs[c].dtype) for c in OBS}
    for c in ("sample", "timepoint", "cell_type", "cell_subtype"):
        obs[c] = obs[c].astype(str)
    missing = S.index.difference(obs.index)
    assert missing.empty, f"{len(missing)} scored barcodes are not in the object"
    order = obs.index[obs.index.isin(S.index)]             # object order
    assert len(order) == len(S)

    flags = set(pd.read_csv(a.flag_barcodes)["barcode"].astype(str))
    T = obs.loc[order].copy()
    T["embedding_outlier"] = T.index.isin(flags)
    T = T.join(S.loc[order])
    T.index.name = "barcode"
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    T.to_csv(a.out)

    # read-back check: the scores in the table are bit-identical to the score table
    R = pd.read_csv(a.out, index_col=0)
    assert np.array_equal(R[S.columns].to_numpy(), S.loc[order].to_numpy()), "scores changed on the way through the CSV"

    P = json.load(open(a.scores_provenance))
    placed = T[["x_um_qc", "y_um_qc"]].notna().all(1)
    prov = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "script": os.path.abspath(__file__), "script_md5": md5(os.path.abspath(__file__)),
        "columns": list(T.reset_index().columns), "obs_dtypes_in_object": dtypes,
        "n_rows": int(len(T)), "n_rows_per_timepoint": {k: int(v) for k, v in T["timepoint"].value_counts().sort_index().items()},
        "n_placed": int(placed.sum()), "n_embedding_outlier": int(T["embedding_outlier"].sum()),
        "n_embedding_outlier_placed": int((T["embedding_outlier"] & placed).sum()),
        "inputs": {"obj": {"path": os.path.abspath(a.obj), "bytes": os.path.getsize(a.obj), "sha256": P.get("obj", {}).get("sha256")},
                   "scores": {"path": os.path.abspath(a.scores), "md5": md5(a.scores)},
                   "scores_provenance": {"path": os.path.abspath(a.scores_provenance), "md5": md5(a.scores_provenance)},
                   "flag_barcodes": {"path": os.path.abspath(a.flag_barcodes), "md5": md5(a.flag_barcodes), "n_listed": len(flags)}},
        "scoring": {k: P.get(k) for k in SCORING_KEYS},
    }
    json.dump(prov, open(side, "w"), indent=1)
    print(json.dumps({k: prov[k] for k in ("n_rows", "n_rows_per_timepoint", "n_placed", "n_embedding_outlier", "n_embedding_outlier_placed", "obs_dtypes_in_object")}), flush=True)


if __name__ == "__main__":
    main()
