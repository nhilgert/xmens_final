"""Run the Hammoud fibroblast centroid mapping and cache every result.

The notebook reads what this writes, so the expensive steps (loading the object,
clustering, the robustness battery) happen once here.
"""

import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from hammoud_fib_mapping import (
    FIB_SUBTYPES, S6C_TO_S5B, load_atlas_centroids, load_atlas_markers,
    linear_panel_matrix, group_centroids, contrastive_cosine, rank_correlation,
    correspondence, variance_panel, specificity, specific_genes,
    marker_ablation, gene_bootstrap, heldout_split, per_group_mapping,
)

HAM = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud")
SCRATCH = Path("/n/netscratch/eddy_lab/Everyone/nhilgert")
HUMAN = SCRATCH / "2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"
OUT = HAM / "fib_centroid_results"
OUT.mkdir(exist_ok=True)

RESOLUTIONS = [0.3, 0.5, 0.8, 1.2, 2.0]
N_PANEL = 3000
MIN_CLUSTER_N = 100
N_BOOT = 500
SEED = 0


def log(msg):
    print(f"[{pd.Timestamp.now():%H:%M:%S}] {msg}", flush=True)


# ------------------------------------------------------------------ 1. load
log("loading atlas centroids")
atlas = load_atlas_centroids(HAM / "pnas.2404775121.sd05.xlsx")
markers_s6c = load_atlas_markers(HAM / "pnas.2404775121.sd06.xlsx")
STRO = [c for c in atlas.columns if c.startswith("Stro-")]
log(f"  atlas: {atlas.shape[0]:,} genes x {atlas.shape[1]} subtypes")

log("loading human object")
human = sc.read_h5ad(HUMAN)
log(f"  {human.n_obs:,} nuclei x {human.n_vars:,} genes")

# X is verified log1p(CP10K); name it explicitly so intent is legible downstream.
human.layers["log1p"] = human.X.copy()

common = human.var_names.intersection(atlas.index)
log(f"  {len(common):,} genes shared with the atlas")

# Donor is not a column in this object. Sample names carry it: mck5* and mck6*
# are the two donors, and each donor contributed one cycle phase, so donor and
# phase are perfectly confounded here.
human.obs["donor"] = (
    human.obs["sample"].astype(str).str.extract(r"(mck\d+)")[0].fillna("unknown")
)
log("  donors: " + str(human.obs["donor"].value_counts().to_dict()))
log("  samples: " + str(human.obs["sample"].value_counts().to_dict()))


# ------------------------------------------------ 2. cluster the query data
log("building neighbor graph on X_scVI")
sc.pp.neighbors(human, use_rep="X_scVI", n_neighbors=30, key_added="nbr30")

X_lin_common = linear_panel_matrix(human, common, layer="log1p")
panel_fib = variance_panel(atlas, FIB_SUBTYPES, common, N_PANEL)
panel_stro = variance_panel(atlas, STRO, common, N_PANEL)

log("resolution sweep")
sweep = []
for res in RESOLUTIONS:
    key = f"leiden_{res}"
    sc.tl.leiden(human, resolution=res, neighbors_key="nbr30", key_added=key,
                 flavor="igraph", n_iterations=2, random_state=SEED)
    C, ns = group_centroids(X_lin_common, human.obs[key], common, min_n=MIN_CLUSTER_N)
    M = contrastive_cosine(C.loc[panel_fib], atlas.loc[panel_fib, FIB_SUBTYPES])
    by_q, by_a = correspondence(M, ns)
    sweep.append({
        "resolution": res,
        "n_clusters": human.obs[key].nunique(),
        "n_clusters_kept": C.shape[1],
        "n_distinct_matches": by_q["best_atlas"].nunique(),
        "n_subtypes_reciprocal": int(by_a["reciprocal"].sum()),
        "median_margin": by_a["margin"].median(),
        "smallest_cluster": int(ns.min()),
    })
    log(f"  res={res}: {sweep[-1]['n_clusters']} clusters, "
        f"{sweep[-1]['n_subtypes_reciprocal']}/6 reciprocal")

sweep = pd.DataFrame(sweep)
sweep.to_csv(OUT / "resolution_sweep.csv", index=False)

# Resolution is set deliberately rather than auto-selected. Counting reciprocal
# matches -- the obvious rule -- picks 2.0, but it does so for the wrong reason:
# the extra "match" is Stro-5 at a score of 0.08, which is noise, while the
# margins for every real match collapse (Stro-10 goes from 0.17 to 0.007 as 6
# clusters become 22). At 0.5 the three genuine correspondences are all
# reciprocal with usable margins and the figure stays readable.
#
# resolution_stability.py recomputes the whole mapping at every resolution and
# confirms the Stro-7 correspondence is reciprocal at all of them, so nothing
# here depends on this choice.
PRIMARY_RES = 0.5
RES = PRIMARY_RES
CLUSTER_KEY = f"leiden_{RES}"
log(f"primary resolution: {RES} "
    f"({int(sweep.loc[sweep.resolution == RES, 'n_clusters'].iloc[0])} clusters)")

human.obs[[CLUSTER_KEY, "sample", "donor", "timepoint"]].to_csv(
    OUT / "cluster_assignments.csv"
)


# ------------------------------------------------------ 3. primary mapping
log("primary mapping")
C_all, ns = group_centroids(X_lin_common, human.obs[CLUSTER_KEY], common,
                            min_n=MIN_CLUSTER_N)
C_all.to_parquet(OUT / "cluster_centroids.parquet")
ns.to_csv(OUT / "cluster_sizes.csv")

results = {}
for tag, cols, panel in [
    ("fib", FIB_SUBTYPES, panel_fib),
    ("stro", STRO, panel_stro),
]:
    M_cos = contrastive_cosine(C_all.loc[panel], atlas.loc[panel, cols])
    M_spr = rank_correlation(C_all.loc[panel], atlas.loc[panel, cols])
    by_q, by_a = correspondence(M_cos, ns)
    by_q_s, by_a_s = correspondence(M_spr, ns)

    M_cos.to_csv(OUT / f"similarity_{tag}_contrastive_cosine.csv")
    M_spr.to_csv(OUT / f"similarity_{tag}_spearman.csv")
    by_q.to_csv(OUT / f"correspondence_{tag}_by_cluster.csv")
    by_a.to_csv(OUT / f"correspondence_{tag}_by_subtype.csv")
    by_a_s.to_csv(OUT / f"correspondence_{tag}_by_subtype_spearman.csv")

    results[tag] = (M_cos, M_spr, by_q, by_a)
    log(f"  {tag}: {int(by_a['reciprocal'].sum())}/{len(cols)} subtypes reciprocal")

M_cos, M_spr, by_q, by_a = results["fib"]
log("\n" + by_a.round(3).to_string())


# ------------------------------------------------------- 4. replication
log("per-sample and per-phase replication")
per_sample = per_group_mapping(
    X_lin_common, human.obs, CLUSTER_KEY, "sample",
    atlas, FIB_SUBTYPES, panel_fib, common,
)
per_sample.to_csv(OUT / "replication_per_sample.csv", index=False)

per_phase = per_group_mapping(
    X_lin_common, human.obs, CLUSTER_KEY, "timepoint",
    atlas, FIB_SUBTYPES, panel_fib, common,
)
per_phase.to_csv(OUT / "replication_per_phase.csv", index=False)

consistency = (
    per_sample.groupby("subtype")["best_query"]
    .agg(modal=lambda s: s.value_counts().index[0],
         agreement=lambda s: s.value_counts(normalize=True).iloc[0],
         n_samples="size")
    .reset_index()
)
consistency.to_csv(OUT / "replication_consistency.csv", index=False)
log("\n" + consistency.round(3).to_string())


# ------------------------------------------------------- 5. robustness
log("marker ablation")
abl = marker_ablation(C_all, atlas, FIB_SUBTYPES, FIB_SUBTYPES, panel_fib)
abl.to_csv(OUT / "robustness_marker_ablation.csv", index=False)

log(f"gene bootstrap ({N_BOOT} draws)")
boot = gene_bootstrap(C_all, atlas, FIB_SUBTYPES, FIB_SUBTYPES, panel_fib,
                      n_boot=N_BOOT, seed=SEED)
boot.to_csv(OUT / "robustness_gene_bootstrap.csv", index=False)
log("\n" + boot.round(3).to_string())

log("held-out gene splits")
held = heldout_split(C_all, atlas, FIB_SUBTYPES, FIB_SUBTYPES, panel_fib,
                     n_rep=N_BOOT, seed=SEED)
held.to_csv(OUT / "robustness_heldout_split.csv", index=False)
log("\n" + held.round(3).to_string())


# ------------------------------------------------- 6. marker-level evidence
log("marker sets and per-cell scores")
marker_sets = {}
for s in FIB_SUBTYPES:
    centroid_derived = specific_genes(atlas, s, FIB_SUBTYPES, common, n=50)
    s6c_label = {v: k for k, v in S6C_TO_S5B.items()}[s]
    published = markers_s6c.query("cluster == @s6c_label").sort_values(
        "avg_log2FC", ascending=False
    )["gene"]
    published = [g for g in published if g in human.var_names][:50]
    marker_sets[s] = {
        "centroid_derived": list(centroid_derived),
        "published_s6c": published,
        "shared": sorted(set(centroid_derived) & set(published)),
    }

pd.DataFrame([
    {"subtype": s,
     "n_centroid_derived": len(v["centroid_derived"]),
     "n_published_detected": len(v["published_s6c"]),
     "n_shared": len(v["shared"]),
     "centroid_derived": ", ".join(v["centroid_derived"]),
     "published_s6c": ", ".join(v["published_s6c"]),
     "shared": ", ".join(v["shared"])}
    for s, v in marker_sets.items()
]).to_csv(OUT / "marker_sets.csv", index=False)

# Per-cell layer: a depth-corrected marker-set score, not a smoothed
# correlation. score_genes bins control genes on the log1p values.
import anndata as ad

tmp = ad.AnnData(X=human.layers["log1p"], obs=human.obs[[]].copy(),
                 var=human.var[[]].copy())
for s in FIB_SUBTYPES:
    genes = marker_sets[s]["published_s6c"] or marker_sets[s]["centroid_derived"]
    name = "score_" + s.split("_")[0]
    sc.tl.score_genes(tmp, genes, score_name=name, use_raw=False, random_state=SEED)
    human.obs[name] = tmp.obs[name].to_numpy()

score_cols = ["score_" + s.split("_")[0] for s in FIB_SUBTYPES]
human.obs[score_cols + [CLUSTER_KEY, "sample", "donor", "timepoint",
                        "x_um_dbscan", "y_um_dbscan"]].to_csv(
    OUT / "per_cell_subtype_scores.csv"
)

# Mean score per cluster: does the nominated cluster actually top the score?
score_by_cluster = (
    human.obs.groupby(CLUSTER_KEY, observed=True)[score_cols].mean()
)
score_by_cluster.to_csv(OUT / "score_by_cluster.csv")
log("\n" + score_by_cluster.round(3).to_string())

# Fraction of nuclei expressing each subtype's headline markers, per cluster.
headline = sorted({g for s in FIB_SUBTYPES for g in marker_sets[s]["published_s6c"][:8]})
headline = [g for g in headline if g in human.var_names]
import scipy.sparse as sp_

expr = sp_.csr_matrix(human[:, headline].layers["log1p"]) > 0
frac_pos = pd.DataFrame(
    np.asarray(expr.todense()), index=human.obs_names, columns=headline
).groupby(human.obs[CLUSTER_KEY].values, observed=True).mean()
frac_pos.to_csv(OUT / "headline_marker_fraction_positive.csv")

with open(OUT / "run_config.txt", "w") as fh:
    fh.write(f"resolution={RES}\ncluster_key={CLUSTER_KEY}\n"
             f"n_panel={N_PANEL}\nmin_cluster_n={MIN_CLUSTER_N}\n"
             f"n_boot={N_BOOT}\nseed={SEED}\nn_obs={human.n_obs}\n"
             f"n_shared_genes={len(common)}\n")

log(f"done -> {OUT}")
