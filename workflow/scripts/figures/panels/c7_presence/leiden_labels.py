"""Leiden clusters of the human fibroblast/decidual nuclei for the C7-fibroblast figure, recomputed with the recipe of
the labels of record (paper_code/ulrich_mapping/run_mapping.py, 30 Aug 2026).

Called by rule c7_leiden_labels (workflow/rules/41_fig_c7_presence.smk).
Recipe: sc.pp.neighbors(use_rep=X_scVI, n_neighbors=30), scanpy defaults otherwise; sc.tl.leiden(resolution=0.5,
        flavor='igraph', n_iterations=2, random_state=0).
Check:  every nucleus gets its label of record (config c7_presence.leiden.labels_of_record); the run fails otherwise.
Writes: the labels (obs_name, leiden_0.5) and a record JSON (recipe, cluster sizes, agreement, versions).

usage: leiden_labels.py --object H5AD --object-sha256 S --labels-of-record CSV --use-rep K --n-neighbors N
                        --resolution R --flavor F --n-iterations N --seed S --labels CSV --record JSON
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

LOG = logging.getLogger("c7_leiden_labels")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--object", required=True, type=Path)
    p.add_argument("--object-sha256", required=True)
    p.add_argument("--labels-of-record", required=True, type=Path)
    p.add_argument("--use-rep", default="X_scVI")
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--resolution", type=float, default=0.5)
    p.add_argument("--flavor", default="igraph")
    p.add_argument("--n-iterations", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--labels", required=True, type=Path)
    p.add_argument("--record", required=True, type=Path)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    checks = common.Checks(LOG.info)

    sha = common.verify_sha256(a.object, a.object_sha256, "human fibroblast/decidual object")
    obs = common.read_obs(a.object)
    latent = common.read_obsm(a.object, a.use_rep).astype(np.float32)
    LOG.info("%s: %d nuclei, %s %s", a.object, len(obs), a.use_rep, latent.shape)

    adata = ad.AnnData(obs=pd.DataFrame(index=obs.index))
    adata.obsm[a.use_rep] = latent
    sc.pp.neighbors(adata, use_rep=a.use_rep, n_neighbors=a.n_neighbors)
    sc.tl.leiden(adata, resolution=a.resolution, flavor=a.flavor, n_iterations=a.n_iterations,
                 random_state=a.seed, key_added="leiden_0.5")
    labels = adata.obs["leiden_0.5"].astype(int)

    record_labels = pd.read_csv(a.labels_of_record, index_col=0)["leiden_0.5"].astype(int)
    record_labels.index = record_labels.index.astype(str)
    checks.equal("every nucleus has a label of record", int(obs.index.isin(record_labels.index).sum()), len(obs))
    checks.equal("labels of record name no other nuclei", len(record_labels), len(obs))
    same = record_labels.reindex(obs.index).to_numpy() == labels.to_numpy()
    checks.equal("nuclei whose recomputed label is not their label of record", int((~same).sum()), 0)
    ari = float(adjusted_rand_score(record_labels.reindex(obs.index).to_numpy(), labels.to_numpy()))

    a.labels.parent.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame({"leiden_0.5": labels.to_numpy()}, index=pd.Index(obs.index, name="obs_name"))
    out.to_csv(a.labels)
    record = common.base_record(Path(__file__), {"human_fibdec": a.object, "labels_of_record": a.labels_of_record},
                                {"labels": a.labels, "record": a.record}, {"human_fibdec": sha})
    record["recipe"] = {"neighbors": {"use_rep": a.use_rep, "n_neighbors": a.n_neighbors, "other": "scanpy defaults"},
                        "leiden": {"resolution": a.resolution, "flavor": a.flavor, "n_iterations": a.n_iterations,
                                   "random_state": a.seed}}
    record["n_nuclei"] = int(len(obs))
    record["cluster_sizes"] = {int(k): int(v) for k, v in labels.value_counts().sort_index().items()}
    record["ari_vs_labels_of_record"] = ari
    record["checks"] = checks.rows
    common.write_json(a.record, record)
    LOG.info("cluster sizes %s; ARI vs labels of record %.6f", record["cluster_sizes"], ari)
    checks.fail_if_any("c7_leiden_labels")


if __name__ == "__main__":
    main()
