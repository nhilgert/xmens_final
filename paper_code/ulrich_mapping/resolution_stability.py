"""Is each subtype's assignment an artefact of the clustering resolution?

The primary run picked one resolution. A reviewer will reasonably ask whether the
correspondence survives a different granularity, so this recomputes the whole
mapping at every resolution and tracks each subtype's score, margin and
reciprocity. It also checks what the nominated clusters actually express.
"""

import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

HAM = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud")
sys.path.insert(0, str(HAM / "fib_mapping"))
from hammoud_fib_mapping import (
    FIB_SUBTYPES, SHORT, load_atlas_centroids, linear_panel_matrix,
    group_centroids, contrastive_cosine, correspondence, variance_panel,
)

OUT = HAM / "fib_centroid_results"
RESOLUTIONS = [0.3, 0.5, 0.8, 1.2, 2.0]
MIN_N = 100


def log(m):
    print(f"[{pd.Timestamp.now():%H:%M:%S}] {m}", flush=True)


atlas = load_atlas_centroids(HAM / "pnas.2404775121.sd05.xlsx")
human = sc.read_h5ad("/n/netscratch/eddy_lab/Everyone/nhilgert/"
                     "2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad")
human.layers["log1p"] = human.X.copy()
common = human.var_names.intersection(atlas.index)
panel = variance_panel(atlas, FIB_SUBTYPES, common, 3000)

log("neighbors")
sc.pp.neighbors(human, use_rep="X_scVI", n_neighbors=30, key_added="nbr30")
X_lin = linear_panel_matrix(human, common, layer="log1p")

rows, labels = [], {}
for res in RESOLUTIONS:
    key = f"leiden_{res}"
    sc.tl.leiden(human, resolution=res, neighbors_key="nbr30", key_added=key,
                 flavor="igraph", n_iterations=2, random_state=0)
    labels[res] = human.obs[key].astype(str).copy()

    C, ns = group_centroids(X_lin, human.obs[key], common, min_n=MIN_N)
    M = contrastive_cosine(C.loc[panel], atlas.loc[panel, FIB_SUBTYPES])
    _, by_a = correspondence(M, ns)
    for s in FIB_SUBTYPES:
        rows.append({
            "resolution": res, "n_clusters": C.shape[1], "subtype": SHORT[s],
            "best_query": by_a.loc[s, "best_query"],
            "score": by_a.loc[s, "score"], "margin": by_a.loc[s, "margin"],
            "reciprocal": by_a.loc[s, "reciprocal"],
            "n": by_a.loc[s, "n_in_best_query"],
        })
    log(f"  res={res} done ({C.shape[1]} clusters)")

stab = pd.DataFrame(rows)
stab.to_csv(OUT / "resolution_stability.csv", index=False)

print("\n=== score of each subtype's best match, by resolution ===")
print(stab.pivot(index="subtype", columns="resolution", values="score").round(3).to_string())
print("\n=== margin ===")
print(stab.pivot(index="subtype", columns="resolution", values="margin").round(3).to_string())
print("\n=== reciprocal ===")
print(stab.pivot(index="subtype", columns="resolution", values="reciprocal").to_string())

# Does the res-2.0 Stro-7 cluster survive as a unit at coarser resolutions?
print("\n=== where the res-2.0 Stro-7 cluster (16) sits at coarser resolutions ===")
c16 = (labels[2.0] == "16").to_numpy()
for res in RESOLUTIONS[:-1]:
    ov = labels[res][c16].value_counts(normalize=True)
    print(f"  res={res}: {ov.index[0]} holds {ov.iloc[0]:.1%} of it "
          f"(that cluster is {int((labels[res]==ov.index[0]).sum()):,} cells; "
          f"cluster 16 is {int(c16.sum()):,})")

# What the nominated clusters actually express.
genes = [g for g in ["C7", "CFH", "CTSK", "CXCL12", "DCN", "LUM", "PDGFRA",
                     "IGFBP1", "PRL", "FOXO1", "WNT4", "MKI67", "TOP2A",
                     "RUNX3", "PTPRC", "ACTA2", "MYH11", "RGS5", "NOTCH3"]
         if g in human.var_names]
X = sp.csr_matrix(human[:, genes].layers["log1p"])
frac = pd.DataFrame(np.asarray((X > 0).todense()), columns=genes,
                    index=human.obs_names).groupby(labels[2.0].values).mean()
mean = pd.DataFrame(np.asarray(X.todense()), columns=genes,
                    index=human.obs_names).groupby(labels[2.0].values).mean()
frac.to_csv(OUT / "identity_marker_fraction_positive_res2.csv")
mean.to_csv(OUT / "identity_marker_mean_res2.csv")

print("\n=== fraction of nuclei positive, res-2.0 clusters (C7 cluster is 16) ===")
print(frac.round(3).to_string())

# Is C7 itself specific to cluster 16?
c7 = np.asarray((human[:, "C7"].layers["log1p"] > 0).todense()).ravel()
print(f"\nC7+ nuclei: {c7.sum():,} / {human.n_obs:,} ({c7.mean():.2%} object-wide)")
print(f"  within cluster 16: {c7[c16].mean():.2%}")
print(f"  outside cluster 16: {c7[~c16].mean():.2%}")
print(f"  enrichment: {c7[c16].mean() / max(c7[~c16].mean(), 1e-9):.1f}x")
share = c7[c16].sum() / max(c7.sum(), 1)
print(f"  cluster 16 holds {share:.1%} of all C7+ nuclei "
      f"while being {c16.mean():.1%} of nuclei")

# Spatial: does cluster 16 sit where the basalis should be?
obs = human.obs
for tp in sorted(obs["timepoint"].astype(str).unique()):
    m = (obs["timepoint"].astype(str) == tp).to_numpy()
    y = obs["y_um_dbscan"].to_numpy(float)
    ok = m & np.isfinite(y)
    yr = (y[ok] - np.nanmin(y[ok])) / (np.nanmax(y[ok]) - np.nanmin(y[ok]))
    inc = c16[ok]
    print(f"\n{tp}: cluster-16 median relative depth {np.median(yr[inc]):.2f} "
          f"vs {np.median(yr[~inc]):.2f} for the rest "
          f"(n={inc.sum():,} of {ok.sum():,})")

log("done")
