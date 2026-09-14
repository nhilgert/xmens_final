#!/usr/bin/env python
"""Headless real-time-informed fibroblast pseudotime — mouse or human.

moscot TemporalProblem -> CellRank RealTimeKernel -> symmetrized transitions -> DPT.

scVI latent: REUSES the per-lineage stromal fit from 04b_reintegrate
(latent_{species}_stromal.h5ad) rather than retraining. Verified 2026-08-28:
  mouse  43,197 cells, same order as the fib/dec deliverable  -> direct adopt
  human  48,004 cells, a superset (471 perivascular throwbacks) -> subset by obs_names
That latent differs substantially from the whole-tissue one on the deliverable
(max |diff| ~4.9 mouse / ~5.8 human), which is the point: 04b reselects HVGs per
lineage because the global fit blurs within-lineage structure.

Emits the same statistics the project's dac_fig notebook reports:
  Spearman(pt, day) over all cells / fibroblasts only / decidual only,
  Spearman(pt, subtype-rank), and median/mean/std/count tables per group.

Usage:
    SPECIES=mouse python run_pseudotime.py
"""
from __future__ import annotations

import json
import os
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse
from scipy.stats import spearmanr

import anndata as ad
import scanpy as sc

SPECIES = os.environ.get("SPECIES", "mouse")
VARIANT = os.environ.get("VARIANT", "v1_all_cells")
REPO = Path(os.environ.get("XMENS_REPO", "/n/lab_storage/eddy_lab/users/nhilgert/xmens_final"))
DATE = os.environ.get("DATE_TAG", "2026_08_28")
SEED = int(os.environ.get("SEED", "0"))
RETRAIN_SCVI = os.environ.get("RETRAIN_SCVI", "0") == "1"

# ---- v2 (induced-only, apoptotic-dropped) knobs --------------------------
# Defaults preserve v1 behaviour exactly: no drops, reuse 04b latent, root "far".
# DROP_SUBTYPES  : comma-list of cell_subtype values to remove (e.g. "decidual_apoptotic")
# DROP_TIMEPOINTS: comma-list of timepoint values to remove (e.g. "uninduced")
# ROOT_DIRECTION : "far" (argmax, v1 default) | "near" (argmin). Mouse needs "near"
#                  for correct orientation (rho +0.665); "far" landed the root at the
#                  late end and inverted the trajectory.
DROP_SUBTYPES = [s for s in os.environ.get("DROP_SUBTYPES", "").split(",") if s.strip()]
DROP_TIMEPOINTS = [s for s in os.environ.get("DROP_TIMEPOINTS", "").split(",") if s.strip()]
ROOT_DIRECTION = os.environ.get("ROOT_DIRECTION", "far").strip().lower()
assert ROOT_DIRECTION in ("far", "near"), f"ROOT_DIRECTION must be far|near, got {ROOT_DIRECTION!r}"

TAG = {"mouse": "XMens", "human": "Human"}[SPECIES]

# ---- species config -------------------------------------------------------
CFG = {
    "mouse": dict(
        moscot_gene_set="mouse",
        root_candidate_label="efib_proliferative",
        timepoint_to_day={
            "uninduced": 0.0, "1_dpi_early": 1.0, "1_dpi_late": 1.5,
            "2_dpi": 2.0, "3_dpi": 3.0, "4_dpi_bleeding": 4.0,
        },
        expected_n=43197,
    ),
    "human": dict(
        moscot_gene_set="human",
        root_candidate_label="efib_proliferative",   # absent -> earliest-day fallback
        timepoint_to_day={"secretory": 0.0, "menstrual": 1.0},
        expected_n=47533,
    ),
}[SPECIES]

CLUSTER_KEY = "cell_subtype"
COUNTS_LAYER = "counts"
BATCH_KEY = "sample"
N_HVG = 3000
N_DIFFMAP = 15
MOSCOT = dict(epsilon=1e-3, tau_a=0.95, tau_b=0.95)

# v2 vocabulary, ordered fibroblast -> decidual (the dac_fig SUBTYPE_ORDER, translated)
FIB_SUBTYPES = ["efib_outer", "efib_middle", "efib_middle_induced",
                "efib_subluminal", "efib_proliferative"]
DEC_SUBTYPES = ["decidual_boundary", "decidual_proliferative",
                "decidual_mature", "decidual_apoptotic"]
SUBTYPE_ORDER = FIB_SUBTYPES + DEC_SUBTYPES

OUT = Path(f"./pseudotime_{SPECIES}_{VARIANT}")
OUT.mkdir(parents=True, exist_ok=True)
report: dict = {"species": SPECIES, "variant": VARIANT, "seed": SEED}


def log(m):
    print(m, flush=True)


def sym_transition(T):
    """Symmetrize a transition matrix and re-normalize rows."""
    S = (T + T.T) / 2
    rs = S.sum(axis=1).A1
    S.data = S.data / rs[S.nonzero()[0]]
    return S


def stats_block(pt, day, subtype, label):
    """The dac_fig statistics, for one cell subset."""
    ok = np.isfinite(pt) & np.isfinite(day)
    if ok.sum() < 10:
        return {"n": int(ok.sum()), "note": "too few cells"}
    rho, p = spearmanr(pt[ok], day[ok])
    per_day = (pd.DataFrame({"day": day[ok], "pt": pt[ok]})
               .groupby("day")["pt"].agg(["median", "mean", "std", "count"]))
    log(f"\n  [{label}]  n={int(ok.sum()):,}   Spearman(pt, day) rho={rho:.4f}  p={p:.3e}")
    log("    median pt by day:")
    for d, row in per_day.iterrows():
        log(f"      day {d:>5}: median={row['median']:.4f} mean={row['mean']:.4f} "
            f"std={row['std']:.4f} n={int(row['count']):,}")
    return {
        "n": int(ok.sum()),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "per_day": {str(k): {kk: (float(vv) if pd.notna(vv) else None)
                             for kk, vv in v.items()}
                    for k, v in per_day.to_dict("index").items()},
    }


def main():
    import cellrank as cr
    import jax
    from moscot.problems.time import TemporalProblem
    from scanpy.tools._dpt import DPT

    log(f"jax backend: {jax.default_backend()}  devices: {jax.devices()}")
    report["jax_backend"] = jax.default_backend()
    if jax.default_backend() != "gpu":
        log("WARNING: JAX not on GPU — moscot will be slow")

    np.random.seed(SEED)
    sc.settings.verbosity = 2

    # ---- load ------------------------------------------------------------
    p = REPO / f"data/deliverables/{DATE}_Cevrim_{TAG}_Only_Fib_Dec_Cells.h5ad"
    adata = ad.read_h5ad(p)
    log(f"\nloaded {p.name}: {adata.n_obs:,} x {adata.n_vars:,}")
    assert COUNTS_LAYER in adata.layers, f"no {COUNTS_LAYER!r} layer: {list(adata.layers)}"

    # the full deliverable must match the expected count BEFORE any filtering,
    # so a silently-truncated input is still caught
    assert adata.n_obs == CFG["expected_n"], f"expected {CFG['expected_n']:,}, got {adata.n_obs:,}"
    report["n_cells_full"] = int(adata.n_obs)

    # ---- v2 filtering: drop apoptotic subtypes and/or uninduced timepoints ----
    if DROP_SUBTYPES or DROP_TIMEPOINTS:
        n0 = adata.n_obs
        if DROP_SUBTYPES:
            have = set(adata.obs[CLUSTER_KEY].astype(str).unique())
            missing = [s for s in DROP_SUBTYPES if s not in have]
            assert not missing, f"DROP_SUBTYPES not in {CLUSTER_KEY}: {missing} (have {sorted(have)})"
            keep = ~adata.obs[CLUSTER_KEY].astype(str).isin(DROP_SUBTYPES).values
            log(f"dropping subtypes {DROP_SUBTYPES}: -{int((~keep).sum()):,} cells")
            adata = adata[keep].copy()
        if DROP_TIMEPOINTS:
            have = set(adata.obs["timepoint"].astype(str).unique())
            missing = [s for s in DROP_TIMEPOINTS if s not in have]
            assert not missing, f"DROP_TIMEPOINTS not in timepoint: {missing} (have {sorted(have)})"
            keep = ~adata.obs["timepoint"].astype(str).isin(DROP_TIMEPOINTS).values
            log(f"dropping timepoints {DROP_TIMEPOINTS}: -{int((~keep).sum()):,} cells")
            adata = adata[keep].copy()
        log(f"after filtering: {adata.n_obs:,} cells (was {n0:,})")
        report["dropped_subtypes"] = DROP_SUBTYPES
        report["dropped_timepoints"] = DROP_TIMEPOINTS
        report["n_cells_dropped"] = int(n0 - adata.n_obs)

    report["n_cells"] = int(adata.n_obs)
    report["cell_type_counts"] = {str(k): int(v) for k, v in
                                  adata.obs["cell_type"].value_counts().items()}
    report["cell_subtype_counts"] = {str(k): int(v) for k, v in
                                     adata.obs[CLUSTER_KEY].value_counts().items()}
    log(f"cell_subtype: {report['cell_subtype_counts']}")

    # ---- day mapping -----------------------------------------------------
    tp_cat = adata.obs["timepoint"].astype("category")
    present = set(tp_cat.cat.categories)
    unmapped = sorted(present - set(CFG["timepoint_to_day"]))
    assert not unmapped, f"unmapped timepoints: {unmapped}"
    mapping = {k: v for k, v in CFG["timepoint_to_day"].items() if k in present}
    adata.obs["day"] = tp_cat.map(mapping).astype(float)
    assert adata.obs["day"].notna().all()
    adata.obs["day"] = adata.obs["day"].astype(float).astype("category")
    n_tp = len(adata.obs["day"].cat.categories)
    log(f"\nday categories: {list(adata.obs['day'].cat.categories)}  "
        f"({n_tp - 1} transition(s))")
    report["n_timepoints"] = n_tp
    report["n_transitions"] = n_tp - 1

    # ---- scVI latent ------------------------------------------------------
    # After a cell drop (apoptotic) or timepoint subset (induced-only), 04b's
    # latent is no longer valid: it was fit on the FULL stromal population, so its
    # coordinates encode variation from cells no longer present. RETRAIN_SCVI=1
    # refits scVI on exactly the filtered cells. Without a drop, reuse 04b (v1).
    filtered = bool(DROP_SUBTYPES or DROP_TIMEPOINTS)
    if filtered and not RETRAIN_SCVI:
        raise SystemExit(
            "Cells were dropped but RETRAIN_SCVI!=1. Reusing 04b's latent on a subset "
            "silently keeps embedding structure from removed cells — set RETRAIN_SCVI=1.")

    if RETRAIN_SCVI:
        import scvi
        scvi.settings.seed = SEED
        # HVGs on the FILTERED cells so the fit reselects for the induced continuum
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3",
                                    layer=COUNTS_LAYER, batch_key=BATCH_KEY,
                                    subset=False)
        log(f"HVGs (on filtered cells): {int(adata.var['highly_variable'].sum())}")
        train = adata[:, adata.var["highly_variable"]].copy()
        train.X = train.layers[COUNTS_LAYER].copy()   # scVI wants raw counts in .X
        scvi.model.SCVI.setup_anndata(train, layer=None, batch_key=BATCH_KEY)
        model = scvi.model.SCVI(train, n_latent=30)
        log(f"\ntraining scVI on {train.n_obs:,} cells x {train.n_vars} HVGs "
            f"(n_latent=30, batch_key={BATCH_KEY!r}) ...")
        model.train(max_epochs=200, early_stopping=True,
                    early_stopping_patience=20, enable_progress_bar=False)
        adata.obsm["X_scVI"] = model.get_latent_representation()
        log(f"scVI latent {adata.obsm['X_scVI'].shape} (retrained on filtered cells)")
        report["scvi_source"] = (f"retrained: {train.n_obs} cells, {train.n_vars} HVGs, "
                                 f"n_latent=30, max_epochs=200/es20")
        report["scvi_epochs_run"] = int(len(model.history.get("elbo_train", [1] * 0)) or 0)
    else:
        lat_p = REPO / f"data/interim/latent_{SPECIES}_stromal.h5ad"
        lat = ad.read_h5ad(lat_p, backed="r")
        lat_names = pd.Index(lat.obs_names.astype(str))
        dlv_names = pd.Index(adata.obs_names.astype(str))
        assert dlv_names.isin(lat_names).all(), "deliverable cells missing from stromal latent"
        idx = lat_names.get_indexer(dlv_names)
        adata.obsm["X_scVI"] = np.asarray(lat.obsm["X_scVI"])[idx]
        lat.file.close()
        log(f"\nadopted 04b stromal latent {adata.obsm['X_scVI'].shape} from {lat_p.name} "
            f"(subset+reordered by obs_names)")
        report["scvi_source"] = f"04b:{lat_p.name}"
        # HVGs still needed downstream by some report paths; harmless in v1
        sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3",
                                    layer=COUNTS_LAYER, batch_key=BATCH_KEY)
        log(f"HVGs: {int(adata.var['highly_variable'].sum())}")

    # ---- neighbours + UMAP on the latent ---------------------------------
    sc.pp.neighbors(adata, use_rep="X_scVI", random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    log(f"UMAP computed on {'retrained' if RETRAIN_SCVI else '04b'} latent -> obsm['X_umap']")

    # ---- moscot ----------------------------------------------------------
    log("\nmoscot TemporalProblem ...")
    tprob = TemporalProblem(adata)
    tprob = tprob.score_genes_for_marginals(
        gene_set_proliferation=CFG["moscot_gene_set"],
        gene_set_apoptosis=CFG["moscot_gene_set"],
    )
    tprob = tprob.prepare(time_key="day", joint_attr="X_scVI")
    tprob = tprob.solve(scale_cost="mean", **MOSCOT)
    log("  solved")

    # ---- RealTimeKernel --------------------------------------------------
    # cellrank's transition_matrix setter rejects a matrix whose rows don't sum
    # to 1 within rtol=1e-12, and the stitched moscot couplings routinely leave a
    # handful of all-zero rows (cells no coupling reached). Porting the fix from
    # realtime_informed_pseudotime_human.ipynb cells 28-30: catch the matrix
    # BEFORE the setter, repair the dead rows with a within-timepoint uniform
    # distribution, renormalize every row exactly, then assign.
    rtk = cr.kernels.RealTimeKernel.from_moscot(tprob)
    try:
        rtk = rtk.compute_transition_matrix(self_transitions="all",
                                            conn_weight=0.2, threshold="auto")
        T = rtk.transition_matrix
        log("RealTimeKernel accepted the stitched matrix directly")
        report["rtk_repair"] = {"needed": False}
    except ValueError as e:
        if "row stochastic" not in str(e):
            raise
        log(f"RealTimeKernel rejected the matrix ({e}) — applying the documented repair")
        if getattr(rtk, "_transition_matrix", None) is not None:
            T = rtk._transition_matrix
        else:
            tmap = rtk._restich_couplings(rtk.couplings, self_transitions="all",
                                          conn_weight=0.2)
            T = tmap.X

        T = csr_matrix(T).tocsr()
        rs = np.asarray(T.sum(axis=1)).ravel()
        dead = np.where(~np.isfinite(rs) | (rs == 0))[0]
        log(f"  dead rows (zero/non-finite): {len(dead)}")

        day_vals = adata.obs["day"].astype(float).to_numpy()
        if len(dead):
            T_lil = T.tolil()
            for i in dead:
                same_tp = np.where(day_vals == day_vals[i])[0]
                T_lil.rows[i] = list(same_tp)
                T_lil.data[i] = [1.0 / len(same_tp)] * len(same_tp)
            T = T_lil.tocsr()

        rs = np.asarray(T.sum(axis=1)).ravel()
        T = csr_matrix(T.multiply(1.0 / rs[:, None]))
        rs2 = np.asarray(T.sum(axis=1)).ravel()
        n_bad = int((~np.isclose(rs2, 1.0, rtol=1e-12)).sum())
        log(f"  after repair: row sums [{rs2.min():.15f}, {rs2.max():.15f}]  "
            f"non-stochastic(rtol=1e-12)={n_bad}")
        rtk.transition_matrix = T
        T = rtk.transition_matrix
        report["rtk_repair"] = {"needed": True, "n_dead_rows": int(len(dead)),
                               "n_nonstochastic_after": n_bad}

    rs = np.asarray(T.sum(axis=1)).ravel()
    log(f"transition matrix {T.shape}  row sums [{rs.min():.9f}, {rs.max():.9f}]")
    assert np.allclose(rs, 1.0, atol=1e-6)

    # ---- DPT -------------------------------------------------------------
    dpt = DPT(adata=adata, neighbors_key="neighbors")
    dpt._transitions_sym = sym_transition(T)
    dpt.compute_eigen(n_comps=N_DIFFMAP, random_state=SEED)
    adata.obsm["X_diffmap"] = dpt.eigen_basis
    adata.uns["diffmap_evals"] = dpt.eigen_values
    log(f"top eigenvalues: {np.round(adata.uns['diffmap_evals'][:8], 5)}")

    # ---- root ------------------------------------------------------------
    dc23 = adata.obsm["X_diffmap"][:, [2, 3]]
    labels = adata.obs[CLUSTER_KEY].astype(str)
    anchor = CFG["root_candidate_label"]
    if anchor in set(labels.unique()):
        root_mask = (labels == anchor).values
        root_source = f"{CLUSTER_KEY}={anchor}"
    else:
        earliest = adata.obs["day"].astype(float).min()
        root_mask = (adata.obs["day"].astype(float) == earliest).values
        root_source = f"earliest day ({earliest}) — {anchor!r} absent"
    centroid = dc23[~root_mask].mean(axis=0)
    cand = np.where(root_mask)[0]
    dists = np.linalg.norm(dc23[root_mask] - centroid, axis=1)
    # "far" (argmax) = root at the anchor cell furthest from everything else;
    # "near" (argmin) = closest to the global centroid. Mouse needs "near".
    pick = dists.argmin() if ROOT_DIRECTION == "near" else dists.argmax()
    root_idx = int(cand[pick])
    log(f"\nroot: {root_source}  pool={int(root_mask.sum()):,}  direction={ROOT_DIRECTION}  idx={root_idx}")
    report["root_direction"] = ROOT_DIRECTION
    log(f"  obs_name={adata.obs_names[root_idx]}  "
        f"timepoint={adata.obs['timepoint'].iloc[root_idx]}  "
        f"subtype={labels.iloc[root_idx]}")
    report.update(root_source=root_source, root_pool=int(root_mask.sum()),
                  root_obs_name=str(adata.obs_names[root_idx]),
                  root_timepoint=str(adata.obs["timepoint"].iloc[root_idx]),
                  root_subtype=str(labels.iloc[root_idx]),
                  anchor_available=anchor in set(labels.unique()))

    adata.uns["iroot"] = root_idx
    sc.tl.dpt(adata)
    pt = adata.obs["dpt_pseudotime"].to_numpy()
    log(f"\ndpt_pseudotime range [{np.nanmin(pt):.4f}, {np.nanmax(pt):.4f}]  "
        f"non-finite={int((~np.isfinite(pt)).sum())}")
    report["pt_min"], report["pt_max"] = float(np.nanmin(pt)), float(np.nanmax(pt))
    report["pt_non_finite"] = int((~np.isfinite(pt)).sum())

    # ---- STATISTICS (mirroring dac_fig) ---------------------------------
    day = adata.obs["day"].astype(float).to_numpy()
    sub = labels.to_numpy()
    log("\n" + "=" * 66 + "\nSTATISTICS\n" + "=" * 66)

    stats = {}
    stats["all_cells"] = stats_block(pt, day, sub, "all cells")

    fib_m = np.isin(sub, FIB_SUBTYPES)
    dec_m = np.isin(sub, DEC_SUBTYPES)
    log(f"\n  fibroblast cells: {int(fib_m.sum()):,} | decidual cells: {int(dec_m.sum()):,}")
    if fib_m.sum() >= 10:
        stats["fibroblasts_only"] = stats_block(pt[fib_m], day[fib_m], sub[fib_m], "fibroblasts only")
    if dec_m.sum() >= 10:
        stats["decidual_only"] = stats_block(pt[dec_m], day[dec_m], sub[dec_m], "decidual only")

    # Spearman(pt, subtype rank)
    rank_map = {s: i for i, s in enumerate(SUBTYPE_ORDER)}
    ranks = np.array([rank_map.get(s, np.nan) for s in sub], dtype=float)
    ok = np.isfinite(ranks) & np.isfinite(pt)
    if ok.sum() >= 10:
        rho_ct, p_ct = spearmanr(pt[ok], ranks[ok])
        log(f"\n  [subtype rank]  n={int(ok.sum()):,}  "
            f"Spearman(pt, subtype_rank) rho={rho_ct:.4f}  p={p_ct:.3e}")
        stats["subtype_rank"] = {"n": int(ok.sum()), "spearman_rho": float(rho_ct),
                                 "spearman_p": float(p_ct)}
        by_sub = (pd.DataFrame({"subtype": pd.Categorical(sub[ok], categories=SUBTYPE_ORDER,
                                                          ordered=True), "pt": pt[ok]})
                  .groupby("subtype", observed=True)["pt"]
                  .agg(["median", "mean", "std", "count"]))
        log("    median pt by subtype (ordered fib -> dec):")
        for s, row in by_sub.iterrows():
            log(f"      {str(s):26s} median={row['median']:.4f} mean={row['mean']:.4f} "
                f"std={row['std']:.4f} n={int(row['count']):,}")
        stats["by_subtype"] = {str(k): {kk: (float(vv) if pd.notna(vv) else None)
                                       for kk, vv in v.items()}
                               for k, v in by_sub.to_dict("index").items()}
    else:
        log("\n  [subtype rank] skipped — no fine subtypes present "
            "(expected for human: cell_subtype is coarse-filled)")
        stats["subtype_rank"] = {"skipped": "no fine subtypes"}

    report["stats"] = stats
    rho_all = stats["all_cells"].get("spearman_rho")
    report["orientation_ok"] = bool(rho_all is not None and rho_all > 0)
    log(f"\norientation: {'OK' if report['orientation_ok'] else 'INVERTED'} (rho={rho_all})")

    # ---- coords ----------------------------------------------------------
    has_xy = adata.obs[["x_um_dbscan", "y_um_dbscan"]].notna().all(axis=1)
    report["cells_with_coords"] = int(has_xy.sum())
    log(f"cells with spatial coords: {int(has_xy.sum()):,} / {adata.n_obs:,}")

    # ---- write -----------------------------------------------------------
    adata.obs["dpt_pseudotime_integrated"] = adata.obs["dpt_pseudotime"]
    out_h5ad = OUT / f"{DATE}_{SPECIES}_fibdec_pseudotime.h5ad"
    tmp = out_h5ad.with_suffix(".h5ad.tmp")
    adata.write_h5ad(tmp)
    os.replace(tmp, out_h5ad)
    log(f"\nwrote {out_h5ad}")

    cols = [c for c in ["sample", "timepoint", "day", "cell_type", CLUSTER_KEY,
                        "dpt_pseudotime", "x_um_dbscan", "y_um_dbscan"]
            if c in adata.obs.columns]
    adata.obs[cols].to_csv(OUT / f"{SPECIES}_pseudotime_obs.csv")
    (OUT / f"{SPECIES}_pseudotime_report.json").write_text(json.dumps(report, indent=2, default=str))
    log(f"wrote {OUT / f'{SPECIES}_pseudotime_report.json'}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        report["FAILED"] = traceback.format_exc()
        (OUT / f"{SPECIES}_pseudotime_report.json").write_text(
            json.dumps(report, indent=2, default=str))
        traceback.print_exc()
        sys.exit(1)
