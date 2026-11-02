#!/usr/bin/env python3
"""
grn_only.py
-----------
Stage 1 of the multi-seed pooled SCENIC sweep: a single seeded GRNBoost2 run.

Reads the pre-built pooled DSC loom (from build_fibroblast_inputs.py) and writes
one seeded adjacency table. cisTarget (Stage 2) and AUCell (Stage 3) are
deterministic given their inputs, so ALL of SCENIC's stochasticity lives here --
which is exactly why this is the step we sweep over seeds.

Design notes:
  * Uses arboreto_with_multiprocessing.py (NOT dask) to avoid dask-on-HPC flakiness.
  * Invokes it with the SAME Python that runs this script (sys.executable), so
    whatever environment your sbatch activates is used end to end. Keep that env
    pinned across all seeds, or the downstream consensus measures env drift
    instead of seed variance.
  * Idempotent: skips if the output adjacencies already exist and are non-empty,
    so re-submitting the array only re-runs failed seeds.
  * Atomic: writes to a temp path and renames on success, so a job killed at
    walltime never leaves a truncated file that the skip-check treats as "done".

The loom must be the canonical pySCENIC layout (genes x cells, 'Gene' row attr,
'CellID' col attr) -- which Stage 0 produces. If your arboreto build refuses loom
input, point --loom at a cells x genes CSV instead; everything else is identical.

Usage:
    python grn_only.py \
        --loom    /.../mouse_dsc_pooled.loom \
        --tf-list /.../allTFs_mm.txt \
        --out     /.../adjacencies_seed1000.tsv \
        --seed    1000 \
        --n-cpu   16
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description="Single seeded GRNBoost2 run (Stage 1)")
    p.add_argument("--loom", required=True, help="Pooled DSC loom (genes x cells) from Stage 0")
    p.add_argument("--tf-list", required=True)
    p.add_argument("--out", required=True,
                   help="Output adjacencies TSV (include the seed in the filename)")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--n-cpu", type=int, default=16)
    p.add_argument("--arboreto-script", default=None,
                   help="Path to arboreto_with_multiprocessing.py "
                        "(default: locate on PATH / next to this env's python)")
    return p.parse_args()


def find_arboreto(explicit):
    if explicit:
        if not os.path.exists(explicit):
            log.error(f"--arboreto-script not found: {explicit}")
            sys.exit(1)
        return explicit
    hit = shutil.which("arboreto_with_multiprocessing.py")
    if hit:
        return hit
    cand = os.path.join(os.path.dirname(sys.executable), "arboreto_with_multiprocessing.py")
    if os.path.exists(cand):
        return cand
    log.error("Could not locate arboreto_with_multiprocessing.py. "
              "Activate the pyscenic env or pass --arboreto-script.")
    sys.exit(1)


def main():
    args = parse_args()

    # --- idempotent skip ---------------------------------------------------
    if os.path.exists(args.out) and os.path.getsize(args.out) > 0:
        log.info(f"Output exists and is non-empty -- skipping seed {args.seed}: {args.out}")
        return

    for f in (args.loom, args.tf_list):
        if not os.path.exists(f):
            log.error(f"Input not found: {f}")
            sys.exit(1)

    arboreto = find_arboreto(args.arboreto_script)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # write to a temp path, rename atomically on success
    tmp_out = f"{args.out}.partial.{os.getpid()}"

    cmd = [
        sys.executable, "-s", arboreto,
        args.loom,
        args.tf_list,
        "--method", "grnboost2",
        "--output", tmp_out,
        "--num_workers", str(args.n_cpu),
        "--seed", str(args.seed),
    ]
    log.info(f"=== Stage 1 GRNBoost2 | seed={args.seed} | workers={args.n_cpu} ===")
    log.info(f"  loom:     {args.loom}")
    log.info(f"  tf-list:  {args.tf_list}")
    log.info(f"  arboreto: {arboreto}")
    log.info(f"  python:   {sys.executable}")
    log.info(f"  cmd: {' '.join(cmd)}")

    try:
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.stdout:
            log.info(res.stdout.strip())
        if res.stderr:
            log.info(f"[stderr] {res.stderr.strip()}")
        if res.returncode != 0:
            log.error(f"arboreto_with_multiprocessing.py exited with code {res.returncode}")
            if os.path.exists(tmp_out):
                os.remove(tmp_out)
            sys.exit(1)
    except Exception as e:
        if os.path.exists(tmp_out):
            os.remove(tmp_out)
        log.error(f"GRNBoost2 invocation failed: {e}")
        sys.exit(1)

    # --- validate before committing ---------------------------------------
    try:
        adj = pd.read_csv(tmp_out, sep="\t")
    except Exception as e:
        log.error(f"Could not read adjacencies output: {e}")
        if os.path.exists(tmp_out):
            os.remove(tmp_out)
        sys.exit(1)

    if "transcription_factor" in adj.columns:  # normalize older arboreto header
        adj = adj.rename(columns={"transcription_factor": "TF"})
        adj.to_csv(tmp_out, sep="\t", index=False)

    needed = {"TF", "target", "importance"}
    missing = needed - set(adj.columns)
    if missing:
        log.error(f"Adjacencies missing required columns {missing}; got {list(adj.columns)}")
        os.remove(tmp_out)
        sys.exit(1)
    if len(adj) == 0:
        log.error("Adjacencies file is empty -- check TF/gene nomenclature overlap "
                  "(MGI for mouse, HGNC for human).")
        os.remove(tmp_out)
        sys.exit(1)

    os.replace(tmp_out, args.out)
    log.info(f"  Wrote {len(adj):,} edges -> {args.out}")


if __name__ == "__main__":
    main()
