"""Plot tables of the C7-fibroblast figure: human nuclei, dot plot, mouse cells with their SAMap-transferred values,
and the other positioned mouse cells (tissue background).

Called by rule c7_tables (workflow/rules/41_fig_c7_presence.smk). A port of extract_final_objects.py (FIG_c7_presence
v8, 2 Nov 2026). No nucleus is excluded.
  human        the 47,533 fibroblast/decidual nuclei: sample, timepoint, cell type, Leiden 0.5 (rule c7_leiden_labels),
               contrastive Stro-7 cosine (rule c7_stro7_cosine), C7 (log-normalized X and counts), X_umap, x_um_qc, y_um_qc
  dot plot     per Leiden cluster and gene (--genes): mean of X and % of nuclei with X > 0 (scanpy dotplot semantics)
  mouse        the fibroblast/decidual cells: sample, timepoint, cell type and subtype, X_umap, positions, and through
               the SAMap graph (connectivity-weighted mean over human neighbours, all human nuclei as sources) the share
               of neighbour weight in each Leiden cluster and the neighbours' Stro-7 cosine
  background   the other mouse cells of the all-cells object that have a position
Checks: every nucleus has a label and a cosine; the cosine equals stored obs stro7_cosine and the transferred cosine
        stored obs stro7_cosine_transferred within --tolerance; shares sum to 1; every cell is in the SAMap graph.

usage: c7_tables.py --human-obj H5AD --human-sha256 S --mouse-obj H5AD --mouse-sha256 S --mouse-all-obj H5AD
                    --mouse-all-sha256 S --joint H5AD --joint-sha256 S --labels CSV --cosine CSV.GZ --genes G,G,...
                    [--target Stro-7 --tolerance 1e-5] --human-table F --dotplot-table F --mouse-table F
                    --background-table F --checks JSON
"""
import argparse
import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

import c7_common as common

LOG = logging.getLogger("c7_tables")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for name in ("human", "mouse", "mouse-all", "joint"):
        p.add_argument(f"--{name}-obj" if name != "joint" else "--joint", required=True, type=Path)
        p.add_argument(f"--{name}-sha256", required=True)
    p.add_argument("--labels", required=True, type=Path)
    p.add_argument("--cosine", required=True, type=Path)
    p.add_argument("--genes", required=True)
    p.add_argument("--target", default="Stro-7")
    p.add_argument("--tolerance", type=float, default=1e-5)
    for name in ("human-table", "dotplot-table", "mouse-table", "background-table", "checks"):
        p.add_argument(f"--{name}", required=True, type=Path)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    checks = common.Checks(LOG.info)
    genes = a.genes.split(",")
    sha = {"human_fibdec": common.verify_sha256(a.human_obj, a.human_sha256, "human fib/dec object"),
           "mouse_fibdec": common.verify_sha256(a.mouse_obj, a.mouse_sha256, "mouse fib/dec object"),
           "mouse_allcells": common.verify_sha256(a.mouse_all_obj, a.mouse_all_sha256, "mouse all-cells object"),
           "joint": common.verify_sha256(a.joint, a.joint_sha256, "SAMap object")}

    # human
    obs = common.read_obs(a.human_obj); n = len(obs)
    umap = common.read_obsm(a.human_obj, "X_umap")
    X = common.gene_columns(a.human_obj, genes)
    C7_counts = common.gene_columns(a.human_obj, ["C7"], element="layers/counts")["C7"]
    lab = pd.read_csv(a.labels, index_col=0)["leiden_0.5"]; lab.index = lab.index.astype(str)
    cos = pd.read_csv(a.cosine, index_col=0)[a.target]; cos.index = cos.index.astype(str)
    checks.equal("every human nucleus has a Leiden label", int(obs.index.isin(lab.index).sum()), n)
    checks.equal("every human nucleus has a Stro-7 cosine", int(obs.index.isin(cos.index).sum()), n)
    d = cos.reindex(obs.index).to_numpy(float) - obs["stro7_cosine"].to_numpy(float)
    checks.within("max |Stro-7 cosine - stored stro7_cosine|", float(np.nanmax(np.abs(d))), a.tolerance)
    H = pd.DataFrame({"sample": obs["sample"].astype(str), "timepoint": obs["timepoint"].astype(str),
                      "cell_type": obs["cell_type"].astype(str),
                      "leiden_0.5": lab.reindex(obs.index).astype(int).to_numpy(),
                      "stro7_cosine": cos.reindex(obs.index).to_numpy(float),
                      "C7_log": X["C7"], "C7_counts": C7_counts, "umap_1": umap[:, 0], "umap_2": umap[:, 1],
                      "x_um_qc": obs["x_um_qc"].to_numpy(), "y_um_qc": obs["y_um_qc"].to_numpy()},
                     index=pd.Index(obs.index, name="obs_name"))
    for path in (a.human_table, a.dotplot_table, a.mouse_table, a.background_table, a.checks):
        path.parent.mkdir(parents=True, exist_ok=True)
    H.to_csv(a.human_table)

    rows = []
    for k in sorted(H["leiden_0.5"].unique()):
        m = H["leiden_0.5"].to_numpy() == k
        for gene in genes:
            x = X[gene][m]
            rows.append({"group": int(k), "human_gene": gene, "mean_expr": float(x.mean()),
                         "pct_expressing": float(100 * (x > 0).mean()), "n_cells": int(m.sum())})
    pd.DataFrame(rows).to_csv(a.dotplot_table, index=False)

    # mouse fibroblast/decidual cells
    mobs = common.read_obs(a.mouse_obj); mumap = common.read_obsm(a.mouse_obj, "X_umap")
    M = pd.DataFrame({"sample": mobs["sample"].astype(str), "timepoint": mobs["timepoint"].astype(str),
                      "cell_type": mobs["cell_type"].astype(str), "cell_subtype": mobs["cell_subtype"].astype(str),
                      "umap_1": mumap[:, 0], "umap_2": mumap[:, 1],
                      "x_um_qc": mobs["x_um_qc"].to_numpy(), "y_um_qc": mobs["y_um_qc"].to_numpy()},
                     index=pd.Index(mobs.index, name="obs_name"))

    # SAMap transfer from every human nucleus
    with h5py.File(a.joint, "r") as handle:
        jobs = read_elem(handle["obs"]); C = read_elem(handle["obsp"]["connectivities"]).tocsr()
    jobs.index = jobs.index.astype(str); N = len(jobs)
    checks.equal("SAMap graph is square on its cells", C.shape, (N, N))
    is_h = jobs.index.isin(H.index); is_m = jobs.index.isin(M.index)
    checks.equal("human nuclei in the SAMap graph", int(is_h.sum()), len(H))
    checks.equal("mouse cells in the SAMap graph", int(is_m.sum()), len(M))
    checks.equal("species of the human nuclei in the graph", set(jobs.loc[is_h, "species"].astype(str)), {"hs"})
    checks.equal("species of the mouse cells in the graph", set(jobs.loc[is_m, "species"].astype(str)), {"mm"})
    labj = H["leiden_0.5"].reindex(jobs.index).to_numpy(dtype=float)
    clusters = sorted(H["leiden_0.5"].unique())
    for k in clusters:
        ind = np.where(np.isfinite(labj), (labj == k).astype(float), np.nan)
        M[f"share_c{k}"] = pd.Series(common.transfer(C, ind, is_h), index=jobs.index).reindex(M.index).to_numpy()
    S = M[[f"share_c{k}" for k in clusters]]
    has = S.notna().all(axis=1)
    M["argmax_cluster"] = np.where(has, "c" + S.idxmax(axis=1).str.replace("share_c", ""), "none")
    M["source_frac_c3"] = float((H["leiden_0.5"] == 3).mean())
    checks.within("max |sum of shares - 1|", float((S[has].sum(axis=1) - 1).abs().max()), 1e-9)
    cosj = H["stro7_cosine"].reindex(jobs.index).to_numpy(dtype=float)
    M["stro7_cosine_transferred"] = pd.Series(common.transfer(C, cosj, is_h), index=jobs.index).reindex(M.index).to_numpy()
    stored = mobs["stro7_cosine_transferred"].to_numpy(float)
    checks.equal("mouse cells with a value: transferred vs stored", int(M["stro7_cosine_transferred"].notna().sum()),
                 int(np.isfinite(stored).sum()))
    d = M["stro7_cosine_transferred"].to_numpy(float) - stored
    checks.within("max |transferred Stro-7 cosine - stored stro7_cosine_transferred|", float(np.nanmax(np.abs(d))),
                  a.tolerance)
    M.to_csv(a.mouse_table)

    # other mouse cells with a position (tissue background)
    aobs = common.read_obs(a.mouse_all_obj, ["sample", "timepoint", "cell_type", "x_um_qc", "y_um_qc"])
    BG = aobs.loc[~aobs.index.isin(M.index) & aobs["x_um_qc"].notna()].astype(
        {"sample": str, "timepoint": str, "cell_type": str})
    BG.index.name = "obs_name"
    BG.to_csv(a.background_table)

    inputs = {"human_fibdec": a.human_obj, "mouse_fibdec": a.mouse_obj, "mouse_allcells": a.mouse_all_obj,
              "joint": a.joint, "labels": a.labels, "cosine": a.cosine}
    outputs = {"human": a.human_table, "dotplot": a.dotplot_table, "mouse": a.mouse_table,
               "background": a.background_table, "checks": a.checks}
    record = common.base_record(Path(__file__), inputs, outputs, sha)
    record["counts"] = {"human": len(H), "human_positioned": int(H["x_um_qc"].notna().sum()), "mouse": len(M),
                        "mouse_with_human_neighbour": int(has.sum()), "mouse_positioned": int(M["x_um_qc"].notna().sum()),
                        "background": len(BG),
                        "cluster_sizes": {int(k): int(v) for k, v in H["leiden_0.5"].value_counts().sort_index().items()}}
    record["checks"] = checks.rows
    common.write_json(a.checks, record)
    checks.fail_if_any("c7_tables")


if __name__ == "__main__":
    main()
