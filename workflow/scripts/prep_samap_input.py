#!/usr/bin/env python
"""
prep_samap_input.py — build ONE SAMap input h5ad from a refined_* pipeline object.

This is the missing edge that connects the annotation pipeline to SAMap:
    refined_{human,mouse}.h5ad
      -> subset to fib/dec on the coarse cell-type column
      -> rename coarse/fine to the SHARED SAMap schema (cell_type_annotations /
         cell_subtype_annotations) that run_samap + the stain expect
      -> put RAW COUNTS in .X (SAMap normalizes internally; feeding the pipeline's already-
         log-normalized .X double-normalizes and its min_expression=1 filter guts the matrix).
         The log matrix is preserved as layers['lognorm'].
      -> subset to one timepoint
      -> write inputs/h5ads_samap/{sp}_fib_dec_{tp}.h5ad

Both migrated species use the canonical ``cell_type`` / ``cell_subtype`` fields.

Usage:
    python prep_samap_input.py --refined data/interim/refined_human.h5ad --species hs \
        --coarse-col cell_type --fine-col cell_subtype --timepoint secretory \
        --expect-fibdec 47533 --out .../inputs/h5ads_samap/hs_fib_dec_secretory.h5ad
"""
import argparse
import anndata as ad

FIBDEC = {"decidual_cell", "endometrial_fibroblast"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--refined", required=True)
    p.add_argument("--species", choices=["hs", "mm"], required=True)
    p.add_argument("--coarse-col", required=True, help="obs col carrying the fib/dec cell_type")
    p.add_argument("--fine-col", required=True, help="obs col carrying the fine subtypes")
    p.add_argument("--timepoint", required=True)
    p.add_argument("--tp-col", default="timepoint")
    p.add_argument("--expect-fibdec", type=int, default=None,
                   help="expected whole fib/dec count (warns, does not fail, if different)")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    print(f"Loading {args.refined} ...")
    a = ad.read_h5ad(args.refined)
    print(f"  {a.shape}; using coarse={args.coarse_col!r} fine={args.fine_col!r} tp={args.tp_col!r}")
    for col in (args.coarse_col, args.fine_col, args.tp_col):
        if col not in a.obs.columns:
            raise SystemExit(f"obs column {col!r} not present. Have: {list(a.obs.columns)}")

    # --- subset to fib/dec on the coarse column -------------------------------------------
    keep = a.obs[args.coarse_col].astype(str).isin(FIBDEC).values
    a = a[keep].copy()
    got = set(a.obs[args.coarse_col].astype(str).unique())
    assert got == FIBDEC, f"after fib/dec subset the coarse labels are {got}, expected exactly {FIBDEC}"
    print(f"  fib/dec subset: {a.n_obs} cells")
    if args.expect_fibdec is not None and a.n_obs != args.expect_fibdec:
        print(f"  WARNING: fib/dec count {a.n_obs} != expected {args.expect_fibdec} "
              f"(annotation may have changed since that number was recorded)")

    # --- shared SAMap schema: the names run_samap carries into the joint obj + the stain reads
    obs = a.obs.copy()
    obs["cell_type_annotations"]    = obs[args.coarse_col].astype(str)
    obs["cell_subtype_annotations"] = obs[args.fine_col].astype(str)
    a.obs = obs

    # --- FEED RAW COUNTS to SAMap --------------------------------------------------------
    # SAMap library-normalizes + log2s internally (sam/core.preprocess_data, min_expression=1),
    # calibrated for COUNT input. refined_*.X is ALREADY log1p(normalize(counts)); passing that
    # made SAMap normalize-then-log a second time and its min_expression=1 threshold (meant for
    # counts) dropped ~half the log-scale genes. So put raw counts in .X and let SAMap do its one
    # correct normalization. Keep the pipeline's log matrix as a layer for anything that wants it.
    if "counts" not in a.layers:
        raise SystemExit("refined object has no 'counts' layer — cannot feed raw counts to SAMap. "
                         "Source counts from an upstream object before this step.")
    a.layers["lognorm"] = a.X.copy()          # preserve the pipeline's log-normalized matrix
    a.X = a.layers["counts"].copy()           # raw counts -> .X (SAMap normalizes once, correctly)
    del a.layers["counts"]
    if a.raw is not None:
        del a.raw

    # scale guard: .X must now be raw counts, not log — fail loudly if this ever regresses
    import numpy as _np
    _xs = a.X[:500]
    _xs = _xs.toarray() if hasattr(_xs, "toarray") else _np.asarray(_xs)
    assert _np.allclose(_xs, _np.round(_xs)) and _xs.max() > 50, (
        f"prep_samap_input expects RAW COUNTS in .X for SAMap, but .X looks log/normalized "
        f"(max={_xs.max():.2f}, integer={_np.allclose(_xs, _np.round(_xs))}). Refusing to write "
        f"mis-scaled SAMap input.")

    # --- subset to the one timepoint -------------------------------------------------------
    tpmask = (a.obs[args.tp_col].astype(str) == args.timepoint).values
    if int(tpmask.sum()) == 0:
        avail = sorted(a.obs[args.tp_col].astype(str).unique())
        raise SystemExit(f"timepoint {args.timepoint!r} has 0 fib/dec cells. Available: {avail}")
    a = a[tpmask].copy()
    print(f"  timepoint={args.timepoint}: {a.n_obs} cells")

    # sanity: the fine subtypes (what the mouse->human transfer votes on) survived
    n_fine = int((a.obs["cell_subtype_annotations"].astype(str).str.len() > 0).sum())
    print(f"  cells with a fine subtype: {n_fine}/{a.n_obs}")
    print(f"  carried through: obsm={list(a.obsm.keys())}  "
          f"cell_barcode in obs? {'cell_barcode' in a.obs.columns}")

    print(f"Writing {args.out} ...")
    from pathlib import Path
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(args.out)
    print("Done.")


if __name__ == "__main__":
    main()
