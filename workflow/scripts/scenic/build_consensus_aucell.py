#!/usr/bin/env python
"""
Stage 3 -- consensus regulons + AUCell, per species, for one SCENIC lane.

A CONSENSUS REGULON is one TF recovered in >= reg_thresh of the seeds (regulon-level
recurrence gate), together with its target genes scored by cross-seed occurrence.
Two ways to spend that occurrence score -- pick with --mode:

  weighted_soft  (DEFAULT, the method)
      targets = every edge recurring in >= min_occurrence seeds (drops single-run
      flukes), each WEIGHTED by occurrence/K. AUCell runs in WEIGHTED mode, so a
      50/50 target dominates and a 2/50 target is a sliver. No degree bias; no hard
      target threshold. A completeness guard drops any regulon left < min_genes.

  hard_filter  (textbook comparator, for the robustness check)
      targets = edges recurring in >= tgt_thresh of the seeds (a binary set).
      AUCell runs UNWEIGHTED. This is the FAQ recipe; tgt_thresh defaults to 0.5
      (below 0.8 on purpose, so it stops deleting diffuse lineage TFs).

Run the default to get the result; run --mode hard_filter to confirm the biology
is stable to the choice. Outputs land in consensus_{tag}_{mode}/ so they never clobber.

  PYBIN=<config scenic.python>  (the pyscenic env)
  $PYBIN build_consensus_aucell.py --species mouse --sketch 15k --n-cpu 8
  $PYBIN build_consensus_aucell.py --species mouse --sketch 15k --n-cpu 8 --mode hard_filter

The first run prints obs columns; re-run with --timepoint-col <name> for the
activity-by-timepoint table.
"""
import argparse, glob, os, pickle, sys
from collections import defaultdict, Counter

import numpy as np
import pandas as pd


def log(*a):
    print(*a, flush=True)


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True, choices=["mouse", "human"])
    ap.add_argument("--sketch", default="15k",
                    help="lane tag: 15k | 5k | full -> picks regulons_sketch{TAG} or regulons_full")
    ap.add_argument("--rev", required=True, help="SCENIC lane root (config scenic.rev): scenic_seeds/<species>/ holds the per-seed outputs")
    ap.add_argument("--h5ad", default=None,
                    help="default {rev}/scenic_seeds/{species}/{species}_dsc_pooled.h5ad")
    ap.add_argument("--layer", default=None,
                    help="adata layer for AUCell; default X. Must be counts or log-norm, NOT scaled/z-scored.")
    # --- consensus mode ---
    ap.add_argument("--mode", choices=["weighted_soft", "hard_filter"], default="weighted_soft",
                    help="weighted_soft (default, the method) | hard_filter (textbook comparator)")
    ap.add_argument("--reg-thresh", type=float, default=0.80,
                    help="regulon gate: keep the TF if recovered in >= this fraction of seeds")
    # weighted_soft knobs
    ap.add_argument("--min-occurrence", type=int, default=2,
                    help="[weighted_soft] a target must recur in >= this many seeds (drops singletons)")
    ap.add_argument("--weight-mode", choices=["occurrence", "importance", "product"],
                    default="occurrence",
                    help="[weighted_soft] AUCell weight: occurrence fraction (default), mean importance, or product")
    # hard_filter knob
    ap.add_argument("--tgt-thresh", type=float, default=0.50,
                    help="[hard_filter] keep a target if it recurs in >= this fraction of the seeds the TF appears in")
    # shared
    ap.add_argument("--min-genes", type=int, default=10,
                    help="completeness guard: drop a regulon left with fewer than this many targets")
    # AUCell (pinned for reproducibility)
    ap.add_argument("--auc-threshold", type=float, default=0.05,
                    help="top fraction of each cell's ranking the AUC integrates over")
    ap.add_argument("--aucell-seed", type=int, default=0,
                    help="seed for AUCell ranking tie-breaks (reproducibility)")
    ap.add_argument("--timepoint-col", default=None,
                    help="obs column to slice activity by; if unset/invalid the script lists obs cols and skips")
    ap.add_argument("--n-cpu", type=int, default=8)
    ap.add_argument("--out-dir", default=None,
                    help="default {rev}/scenic_seeds/{species}/consensus_{TAG}_{mode}")
    return ap.parse_args()


def lane_tag(sketch):
    return "full" if sketch == "full" else f"sketch{sketch}"


def main():
    args = parse_args()
    tag = lane_tag(args.sketch)
    root = os.path.join(args.rev, "scenic_seeds", args.species)
    reg_dir = os.path.join(root, "regulons_full" if args.sketch == "full"
                           else f"regulons_sketch{args.sketch}")

    pkls = sorted(glob.glob(os.path.join(reg_dir, "regulons_seed*.pkl")))
    if not pkls:
        sys.exit(f"No regulon pickles found in {reg_dir}")
    n_seeds = len(pkls)
    log(f"[{args.species}/{tag}] mode={args.mode} | {n_seeds} seed files in {reg_dir}")
    if args.mode == "weighted_soft":
        log(f"  reg_thresh={args.reg_thresh} min_occurrence={args.min_occurrence} "
            f"weight_mode={args.weight_mode} min_genes={args.min_genes}")
    else:
        log(f"  reg_thresh={args.reg_thresh} tgt_thresh={args.tgt_thresh} "
            f"min_genes={args.min_genes} (unweighted AUCell)")

    # --- tally across seeds -------------------------------------------------
    present = defaultdict(int)                          # name -> n seeds present
    tgt_occ = defaultdict(Counter)                      # name -> Counter(target -> n seeds)
    wsum = defaultdict(lambda: defaultdict(float))      # name -> target -> summed importance
    template = {}
    for p in pkls:
        with open(p, "rb") as fh:
            regs = pickle.load(fh)
        seen = set()
        for r in regs:
            name = r.name
            if name in seen:
                continue
            seen.add(name)
            present[name] += 1
            template.setdefault(name, r)
            for g, w in r.gene2weight.items():
                tgt_occ[name][g] += 1
                wsum[name][g] += float(w)

    # --- build consensus regulons ------------------------------------------
    import attr
    fields_cache = {}

    def evolve(t, g2w, g2o):
        cls = type(t)
        if cls not in fields_cache:
            fields_cache[cls] = {f.name for f in attr.fields(cls)}
        fn = fields_cache[cls]
        kw = {"gene2weight": g2w}
        if "gene2occurrence" in fn:
            kw["gene2occurrence"] = g2o
        if "context" in fn:
            kw["context"] = frozenset()
        if "score" in fn:
            kw["score"] = 0.0
        return attr.evolve(t, **kw)

    consensus, rows = [], []
    for name, npres in present.items():
        frac = npres / n_seeds
        keep_reg = frac >= args.reg_thresh
        occ = tgt_occ[name]
        kept = {}
        if args.mode == "weighted_soft":
            for g, c in occ.items():
                if c < args.min_occurrence:            # drop pure flukes
                    continue
                occ_frac = c / npres
                mean_imp = wsum[name][g] / c
                if args.weight_mode == "occurrence":
                    kept[g] = occ_frac
                elif args.weight_mode == "importance":
                    kept[g] = mean_imp
                else:
                    kept[g] = occ_frac * mean_imp
        else:  # hard_filter -> binary membership, AUCell unweighted
            cut = args.tgt_thresh * npres
            for g, c in occ.items():
                if c >= cut:
                    kept[g] = 1.0
        n_kept = len(kept)
        keep = keep_reg and n_kept >= args.min_genes   # completeness guard
        rows.append(dict(name=name, n_present=npres, frac_present=round(frac, 3),
                         n_targets_union=len(occ), n_targets_kept=n_kept,
                         median_occ=int(np.median([occ[g] for g in kept])) if kept else 0,
                         kept=keep))
        if keep:
            consensus.append(evolve(template[name], dict(kept), {g: occ[g] for g in kept}))

    report = pd.DataFrame(rows).sort_values(
        ["kept", "n_present", "n_targets_kept"], ascending=False)

    out_dir = args.out_dir or os.path.join(root, f"consensus_{tag}_{args.mode}")
    os.makedirs(out_dir, exist_ok=True)
    report.to_csv(os.path.join(out_dir, "recurrence_report.tsv"), sep="\t", index=False)
    with open(os.path.join(out_dir, "consensus_regulons.pkl"), "wb") as fh:
        pickle.dump(consensus, fh)

    # --- summary + histograms ----------------------------------------------
    qs = (.1, .25, .5, .75, .9)
    summ = [f"lane: {args.species}/{tag}  mode: {args.mode}",
            f"seeds: {n_seeds}",
            (f"reg_thresh={args.reg_thresh} min_occurrence={args.min_occurrence} "
             f"weight_mode={args.weight_mode} min_genes={args.min_genes}"
             if args.mode == "weighted_soft" else
             f"reg_thresh={args.reg_thresh} tgt_thresh={args.tgt_thresh} min_genes={args.min_genes}"),
            f"regulon names seen: {len(report)}",
            f"kept regulons: {int(report.kept.sum())}",
            "regulon recurrence frac quantiles: " +
            ", ".join(f"{q}:{report.frac_present.quantile(q):.2f}" for q in qs)]
    if report.kept.any():
        summ.append("targets/regulon (kept) quantiles: " +
                    ", ".join(f"{q}:{report.loc[report.kept, 'n_targets_kept'].quantile(q):.0f}" for q in qs))
    with open(os.path.join(out_dir, "recurrence_summary.txt"), "w") as fh:
        fh.write("\n".join(summ) + "\n")
    log("\n".join(summ))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].hist(report.frac_present, bins=20)
        ax[0].axvline(args.reg_thresh, c="r", ls="--")
        ax[0].set_title("regulon recurrence"); ax[0].set_xlabel("frac of seeds")
        allocc = np.concatenate([np.fromiter(tgt_occ[n].values(), dtype=int) for n in tgt_occ])
        ax[1].hist(allocc, bins=range(0, n_seeds + 2))
        cut_line = args.min_occurrence if args.mode == "weighted_soft" else args.tgt_thresh * n_seeds
        ax[1].axvline(cut_line, c="r", ls="--")
        ax[1].set_title("per-target occurrence (all edges)"); ax[1].set_xlabel("n seeds")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "recurrence_hist.png"), dpi=120)
        log("wrote recurrence_hist.png")
    except Exception as e:
        log(f"(hist skipped: {e})")

    if not consensus:
        sys.exit("No consensus regulons survived -- loosen --reg-thresh/--min-genes.")

    # --- AUCell ------------------------------------------------------------
    import anndata as ad
    import scipy.sparse as sp
    from pyscenic.aucell import aucell

    h5ad = args.h5ad or os.path.join(root, f"{args.species}_dsc_pooled.h5ad")
    log(f"loading {h5ad}")
    adata = ad.read_h5ad(h5ad)
    X = adata.layers[args.layer] if args.layer else adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    exp = pd.DataFrame(X.astype(np.float32), index=adata.obs_names, columns=adata.var_names)
    del X

    cons_genes = set().union(*[set(r.gene2weight.keys()) for r in consensus])
    miss = cons_genes - set(exp.columns)
    noweights = (args.mode == "hard_filter")
    log(f"consensus regulons: {len(consensus)} | consensus genes: {len(cons_genes)} | missing: {len(miss)}")
    log(f"AUCell ({'unweighted' if noweights else 'weighted'}, auc_threshold={args.auc_threshold}) on "
        f"{exp.shape[0]} cells x {exp.shape[1]} genes, {args.n_cpu} workers ...")
    auc_kwargs = dict(num_workers=args.n_cpu, auc_threshold=args.auc_threshold, noweights=noweights)
    try:
        auc = aucell(exp, consensus, seed=args.aucell_seed, **auc_kwargs)
    except TypeError:
        auc = aucell(exp, consensus, **auc_kwargs)
    auc.to_csv(os.path.join(out_dir, "auc_mtx.csv"))
    log(f"wrote auc_mtx.csv  ({auc.shape[0]} cells x {auc.shape[1]} regulons)")

    # --- slice by timepoint -------------------------------------------------
    tcol = args.timepoint_col
    if tcol is None or tcol not in adata.obs.columns:
        if tcol is not None:
            log(f"timepoint col '{tcol}' not in obs.")
        log("obs columns available: " + ", ".join(map(str, adata.obs.columns)))
        log("re-run with --timepoint-col <name> for activity-by-timepoint; skipping slice.")
    else:
        key = adata.obs.loc[auc.index, tcol]
        grp = auc.groupby(key.values).mean()
        grp.T.to_csv(os.path.join(out_dir, "auc_by_timepoint.tsv"), sep="\t")
        log(f"wrote auc_by_timepoint.tsv  (regulons x {grp.shape[0]} '{tcol}' levels)")

    log(f"DONE -> {out_dir}")


if __name__ == "__main__":
    main()
