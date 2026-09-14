#!/usr/bin/env python
"""Build one MI input per run for the dpt_pseudotime specificity comparison, then hand it to the
existing mi_pseudotime.py prep (which derives log1p(CP10K) from counts itself and records that in
its manifest).

RUNS
  A  dpt_unind_mouse   mouse fib/dec, timepoint == 'uninduced'          rank/N WITHIN each sample
  B  dpt_full_mouse    mouse fib/dec, all cells                         rank/N GLOBAL
  C  dpt_full_human    human fib/dec, all cells                         rank/N GLOBAL
  D  dpt_indmatch_mouse mouse fib/dec, induced, size-matched to A       rank/N WITHIN each sample

Quantile normalisation is rank/N as specified. Note for the record: MI is invariant under a monotone
transform of either variable, so ranking cannot create or destroy information -- it changes the kNN
estimator's local distances only. The reason it matters here is comparability: run A's raw dpt spans
only 0-0.425 of the global axis, so without per-arm ranking the arms would be estimated on coordinates
of different dynamic range.
"""
import sys, json, pathlib, argparse
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, warnings
warnings.filterwarnings("ignore")

X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
N = pathlib.Path("/n/netscratch/eddy_lab/Everyone/nhilgert")
MOUSE = N / "2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"
HUMAN = N / "2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"
PTK, SEED = "dpt_pseudotime", 0

SPEC = {
    "dpt_unind_mouse":    dict(src=MOUSE, arm="uninduced",  norm="within_sample"),
    "dpt_full_mouse":     dict(src=MOUSE, arm="all",        norm="global"),
    "dpt_full_human":     dict(src=HUMAN, arm="all",        norm="global"),
    "dpt_indmatch_mouse": dict(src=MOUSE, arm="induced_matched", norm="within_sample"),
}

def rank_norm(v):
    """rank / N in (0, 1]; ties averaged."""
    return pd.Series(v).rank(method="average").to_numpy() / len(v)

def build(run):
    s = SPEC[run]
    print(f"[build {run}] reading {s['src'].name}", flush=True)
    a = sc.read_h5ad(s["src"])
    assert "counts" in a.layers, f"no counts layer; layers = {list(a.layers)}"
    pt = pd.to_numeric(a.obs[PTK], errors="coerce").to_numpy(float)
    assert np.isfinite(pt).all(), f"{int((~np.isfinite(pt)).sum())} non-finite {PTK}"
    tp = a.obs["timepoint"].astype(str).to_numpy()
    sm = a.obs["sample"].astype(str).to_numpy()

    if s["arm"] == "uninduced":
        keep = tp == "uninduced"
    elif s["arm"] == "all":
        keep = np.ones(a.n_obs, bool)
    elif s["arm"] == "induced_matched":
        n_target = int((tp == "uninduced").sum())
        pool = np.flatnonzero(tp != "uninduced")
        rng = np.random.default_rng(SEED)
        keep = np.zeros(a.n_obs, bool)
        keep[np.sort(rng.choice(pool, size=n_target, replace=False))] = True
        print(f"[build {run}] size-matched draw: {n_target:,} of {len(pool):,} induced cells", flush=True)
    b = a[keep].copy()
    del a
    print(f"[build {run}] {b.n_obs:,} cells | timepoints {dict(pd.Series(tp[keep]).value_counts())}", flush=True)

    ptk = pt[keep]; smk = sm[keep]
    if s["norm"] == "global":
        b.obs["pt_norm"] = rank_norm(ptk)
        print(f"[build {run}] GLOBAL rank/N over {b.n_obs:,} cells", flush=True)
    else:
        out = np.full(b.n_obs, np.nan)
        for grp in pd.unique(smk):
            m = smk == grp
            out[m] = rank_norm(ptk[m])
            print(f"[build {run}]   within {grp}: n={int(m.sum()):,}  "
                  f"raw dpt {ptk[m].min():.4f}-{ptk[m].max():.4f} -> 0-1", flush=True)
        assert np.isfinite(out).all()
        b.obs["pt_norm"] = out
    b.obs["pt_raw"] = ptk
    from scipy.stats import spearmanr
    rho = float(spearmanr(b.obs.pt_norm, b.obs.pt_raw).statistic)
    print(f"[build {run}] spearman(pt_norm, pt_raw) = {rho:+.4f}"
          f"  ({'monotone -- global' if s['norm']=='global' else 'per-sample monotone, pooled rho < 1 by construction'})",
          flush=True)

    # keep only what prep needs: counts layer + obs + var
    out = ad.AnnData(X=b.layers["counts"].copy(), obs=b.obs[["pt_norm","pt_raw","timepoint","sample","cell_subtype"]].copy(),
                     var=b.var.copy(), layers={"counts": b.layers["counts"].copy()})
    dest = X / "results/mi_dpt" / run
    dest.mkdir(parents=True, exist_ok=True)
    h5 = dest / "input.h5ad"
    out.write_h5ad(h5)
    meta = dict(run=run, source=str(s["src"]), arm=s["arm"], normalization=s["norm"],
                pt_key_source=PTK, pt_key_used="pt_norm", n_cells=int(out.n_obs),
                n_genes_before_filter=int(out.n_vars), seed=SEED,
                spearman_ptnorm_vs_ptraw=round(rho, 4),
                note=("MI is invariant under monotone transforms of either variable; ranking is applied for "
                      "cross-arm comparability of the coordinate's dynamic range, not to change information content."))
    (dest / "build_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"[build {run}] wrote {h5} ({h5.stat().st_size/1e9:.2f} GB)", flush=True)
    return h5

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run", choices=sorted(SPEC))
    build(ap.parse_args().run)
