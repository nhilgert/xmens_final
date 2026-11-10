"""Review: how the C7-fibroblast cluster (Leiden cluster 3 of record) depends on the Leiden seed and on the graph.

Called by rule review_c7_leiden_seeds (workflow/rules/41_fig_c7_presence.smk; not in the default target).
Leiden (resolution 0.5, igraph, 2 iterations) for seeds 0..n-1 on (a) the object's stored neighbour graph (obsp
connectivities, 15 neighbours) and (b) the 30-neighbour graph of the recipe of record rebuilt on obsm X_scVI. Per run:
number of clusters, ARI to the labels of record, the cluster holding most cluster-3 nuclei (share of cluster 3 in it, its
purity, own cluster = both >= 90%), its C7 detection against the highest other cluster, its median Stro-7 cosine and
whether that is the highest. Also writes every run's labels.

usage: leiden_seed_sweep.py --object H5AD --object-sha256 S --labels CSV --cosine CSV.GZ --seeds N
                            --runs CSV --labels-out CSV.GZ --summary JSON
"""
import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

import c7_common as common

LOG = logging.getLogger("review_c7_leiden_seeds")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--object", required=True, type=Path)
    p.add_argument("--object-sha256", required=True)
    p.add_argument("--labels", required=True, type=Path, help="labels of rule c7_leiden_labels")
    p.add_argument("--cosine", required=True, type=Path, help="cosines of rule c7_stro7_cosine")
    p.add_argument("--target", default="Stro-7")
    p.add_argument("--seeds", type=int, default=20)
    p.add_argument("--runs", required=True, type=Path)
    p.add_argument("--labels-out", required=True, type=Path)
    p.add_argument("--summary", required=True, type=Path)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    sha = common.verify_sha256(a.object, a.object_sha256, "human fibroblast/decidual object")

    import h5py
    from anndata.io import read_elem
    obs = common.read_obs(a.object)
    with h5py.File(a.object, "r") as handle:
        g_stored = read_elem(handle["obsp"]["connectivities"]).tocsr()
    latent = common.read_obsm(a.object, "X_scVI").astype(np.float32)
    c7 = common.gene_columns(a.object, ["C7"])["C7"] > 0
    lab = pd.read_csv(a.labels, index_col=0)["leiden_0.5"]; lab.index = lab.index.astype(str)
    y = lab.reindex(obs.index).to_numpy()
    cos = pd.read_csv(a.cosine, index_col=0)[a.target]; cos.index = cos.index.astype(str)
    cos = cos.reindex(obs.index).to_numpy()
    is3 = y == 3
    A = ad.AnnData(obs=pd.DataFrame(index=obs.index)); A.obsm["X_scVI"] = latent
    sc.pp.neighbors(A, use_rep="X_scVI", n_neighbors=30)
    graphs = {"stored_15nn": g_stored, "recipe_30nn": A.obsp["connectivities"].tocsr()}

    rows, labels_out = [], pd.DataFrame(index=pd.Index(obs.index, name="obs_name"))
    for gname, G in graphs.items():
        for seed in range(a.seeds):
            B = ad.AnnData(obs=pd.DataFrame(index=obs.index))
            sc.tl.leiden(B, resolution=0.5, flavor="igraph", n_iterations=2, random_state=seed, adjacency=G, key_added="L")
            L = B.obs["L"].astype(int).to_numpy(); labels_out[f"{gname}_seed{seed}"] = L
            main = pd.Series(L[is3]).value_counts(); k3 = int(main.index[0])
            in_main = float(main.iloc[0] / is3.sum()); purity = float((y[L == k3] == 3).mean())
            med = pd.Series(cos).groupby(L).median(); c7r = pd.Series(c7).groupby(L).mean()
            rows.append({"graph": gname, "seed": seed, "n_clusters": int(len(np.unique(L))),
                         "ARI_vs_labels_of_record": float(adjusted_rand_score(y, L)),
                         "c3_main_cluster_size": int((L == k3).sum()), "c3_in_main_cluster_pct": round(100 * in_main, 1),
                         "main_cluster_is_c3_pct": round(100 * purity, 1),
                         "c3_own_cluster": bool(in_main >= 0.9 and purity >= 0.9),
                         "main_cluster_C7_pct": round(100 * float(c7r[k3]), 2),
                         "max_C7_pct_other_clusters": round(100 * float(c7r.drop(k3).max()), 2),
                         "main_cluster_cosine_median": round(float(med[k3]), 4),
                         "main_cluster_has_top_cosine": bool(med.idxmax() == k3)})
            LOG.info("%s seed %d: %s", gname, seed, rows[-1])
    R = pd.DataFrame(rows)
    for path in (a.runs, a.labels_out, a.summary):
        path.parent.mkdir(parents=True, exist_ok=True)
    R.to_csv(a.runs, index=False)
    labels_out.to_csv(a.labels_out)
    summary = R.groupby("graph").agg(runs=("seed", "size"), c3_own_cluster=("c3_own_cluster", "sum"),
                                     n_clusters_min=("n_clusters", "min"), n_clusters_max=("n_clusters", "max"),
                                     ARI_median=("ARI_vs_labels_of_record", "median"),
                                     ARI_min=("ARI_vs_labels_of_record", "min"), ARI_max=("ARI_vs_labels_of_record", "max"),
                                     c3_in_main_min=("c3_in_main_cluster_pct", "min"),
                                     top_cosine_runs=("main_cluster_has_top_cosine", "sum")).reset_index()
    record = common.base_record(Path(__file__), {"human_fibdec": a.object, "labels": a.labels, "cosine": a.cosine},
                                {"runs": a.runs, "labels": a.labels_out, "summary": a.summary}, {"human_fibdec": sha})
    record["summary"] = summary.to_dict("records")
    record["n_clusters_counts"] = {g: {int(k): int(v) for k, v in R[R.graph == g].n_clusters.value_counts().sort_index().items()}
                                   for g in R.graph.unique()}
    common.write_json(a.summary, record)


if __name__ == "__main__":
    main()
