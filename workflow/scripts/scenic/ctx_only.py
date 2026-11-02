#!/usr/bin/env python3
"""
ctx_only.py
-----------
Stage 2 of the multi-seed SCENIC sweep: cisTarget pruning for ONE seed.

Takes that seed's GRNBoost2 adjacencies + the expression loom, builds co-expression
modules (modules_from_adjacencies), prunes them by motif enrichment (prune2df ->
cisTarget), converts to regulons (df2regulons), and pickles them.

Deterministic given the adjacencies — all of SCENIC's stochasticity lived in Stage 1 —
so there is no seed argument here; the seed only determines which adjacencies file we
read and which regulons file we write.

Memory-heavy (loads the feather ranking DBs), compute-light. Idempotent + atomic:
skips if the output pickle already exists, writes to a temp path and renames on success.

Usage:
    python ctx_only.py \
        --adjacencies /.../adjacencies_sketch15k/adjacencies_seed1000.tsv \
        --loom        /.../mouse_dsc_sketch15k.loom \
        --db-dir      /.../databases \
        --species     mouse \
        --out         /.../regulons_sketch15k/regulons_seed1000.pkl \
        --n-cpu       8
"""

import argparse
import logging
import os
import pickle
import sys
import traceback

import numpy as np
import pandas as pd
import loompy as lp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# Motif annotations + ranking DBs (v10 / clustered), matching run_scenic_revisions.py
DB_CONFIG = {
    "mouse": {
        "motif_anno": "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl",
        "rankings": [
            "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
            "mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        ],
    },
    "human": {
        "motif_anno": "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        "rankings": [
            "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
            "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        ],
    },
}


def parse_args():
    p = argparse.ArgumentParser(description="Single-seed cisTarget pruning (Stage 2)")
    p.add_argument("--adjacencies", required=True, help="Seed's GRNBoost2 adjacencies TSV")
    p.add_argument("--loom", required=True, help="Expression loom (same one used for GRN)")
    p.add_argument("--db-dir", required=True, help="Directory with motif annotations + feather rankings")
    p.add_argument("--species", required=True, choices=["mouse", "human"])
    p.add_argument("--out", required=True, help="Output regulons pickle (include seed in name)")
    p.add_argument("--n-cpu", type=int, default=8)
    return p.parse_args()


def resolve_db_paths(db_dir, species):
    cfg = DB_CONFIG[species]
    motif_anno = os.path.join(db_dir, cfg["motif_anno"])
    rankings = [os.path.join(db_dir, f) for f in cfg["rankings"]]
    missing = [f for f in [motif_anno] + rankings if not os.path.exists(f)]
    if missing:
        for m in missing:
            log.error(f"DB file not found: {m}")
        sys.exit(1)
    return motif_anno, rankings


def load_expr_from_loom(loom_path):
    """Loom is genes x cells; modules_from_adjacencies wants cells x genes."""
    with lp.connect(loom_path, mode="r", validate=False) as ds:
        genes = np.array(ds.ra["Gene"])
        cells = np.array(ds.ca["CellID"])
        mat = ds[:, :]  # genes x cells
    return pd.DataFrame(mat.T, index=cells, columns=genes)  # cells x genes


def main():
    args = parse_args()

    if os.path.exists(args.out) and os.path.getsize(args.out) > 0:
        log.info(f"Output exists, skipping: {args.out}")
        return

    for f in (args.adjacencies, args.loom):
        if not os.path.exists(f):
            log.error(f"Input not found: {f}")
            sys.exit(1)

    motif_anno, rankings = resolve_db_paths(args.db_dir, args.species)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    log.info(f"=== Stage 2 cisTarget | species={args.species} ===")
    log.info(f"  adjacencies: {args.adjacencies}")
    log.info(f"  loom:        {args.loom}")
    log.info(f"  motif_anno:  {motif_anno}")
    log.info(f"  rankings:    {[os.path.basename(r) for r in rankings]}")

    try:
        from pyscenic.utils import modules_from_adjacencies
        from pyscenic.prune import prune2df, df2regulons
        from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

        adj = pd.read_csv(args.adjacencies, sep="\t")
        if "transcription_factor" in adj.columns:
            adj = adj.rename(columns={"transcription_factor": "TF"})
        log.info(f"  adjacencies: {len(adj):,} edges")

        ex_mtx = load_expr_from_loom(args.loom)
        log.info(f"  expr matrix: {ex_mtx.shape[0]:,} cells x {ex_mtx.shape[1]:,} genes")

        # default rho_mask_dropouts=False (compute correlation over the full matrix —
        # the behavior the protocol prose recommends)
        modules = list(modules_from_adjacencies(adj, ex_mtx))
        log.info(f"  co-expression modules: {len(modules):,}")

        dbs = [RankingDatabase(fname=p, name=os.path.splitext(os.path.basename(p))[0])
               for p in rankings]

        df = prune2df(
            dbs,
            modules,
            motif_anno,
            client_or_address="custom_multiprocessing",
            module_chunksize=100,
            num_workers=args.n_cpu,
        )
        log.info(f"  prune2df: {len(df):,} enriched-motif rows")

        regulons = df2regulons(df)
        log.info(f"  regulons: {len(regulons):,}")

    except Exception as e:
        log.error(f"cisTarget failed: {e}")
        log.error(traceback.format_exc())
        sys.exit(1)

    if len(regulons) == 0:
        log.error("Zero regulons produced — check motif annotation / nomenclature match.")
        sys.exit(1)

    tmp_out = f"{args.out}.partial.{os.getpid()}"
    with open(tmp_out, "wb") as f:
        pickle.dump(regulons, f)
    os.replace(tmp_out, args.out)
    log.info(f"  Wrote {len(regulons):,} regulons -> {args.out}")


if __name__ == "__main__":
    main()
