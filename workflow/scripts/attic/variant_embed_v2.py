#!/usr/bin/env python
"""
variant_embed.py — mouse re-embedding variants for the root-selection workflow.

Builds ONE variant: subset -> retrain scVI -> neighbors -> UMAP -> diffmap -> STOP.
Deliberately NO root selection and NO DPT: the object is handed to
root_selection_mouse.ipynb where Nick stars the root manually (standing rule:
roots are explicit and visually inspected, never auto-picked).

Baseline drops (all variants, = v2): decidual_apoptotic subtype, uninduced timepoint.
Variant axes:
  --drop-hiapo-mature   drop decidual_mature cells whose obs['apoptosis'] score
                        >= the 25th percentile of the ANNOTATED decidual_apoptotic
                        population's scores (threshold computed from the deliverable
                        BEFORE the apoptotic subtype is dropped — label-anchored,
                        not a magic within-mature quantile). Threshold + counts
                        recorded in the report.
  --drop-efib-outer     drop efib_outer (a priori justification: transcriptionally
                        inert, pt-proliferation rho < 0, early-day-skewed —
                        non-participating population).

Variants: A = --drop-hiapo-mature ; B = --drop-efib-outer ; C = both.

Output: {outdir}/{name}_embedded.h5ad  (+ {name}_report.json)
carries X_scVI (retrained), X_umap, X_diffmap, neighbors — root-selection-ready.
"""
import argparse
import json
import pathlib
import sys

import numpy as np
import scanpy as sc

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True)          # the mouse fib/dec deliverable
p.add_argument("--name", required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--drop-hiapo-mature", action="store_true")
p.add_argument("--drop-efib-outer", action="store_true")
p.add_argument("--keep-uninduced", action="store_true", help="varD: retain the uninduced timepoint (v2 baseline drops it)")
p.add_argument("--apo-knn", type=int, default=30, help="k for the graph-neighbourhood apoptotic fraction (computed on the deliverable X_scVI BEFORE any drop)")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--hvg", type=int, default=2000)
p.add_argument("--epochs", type=int, default=200)
args = p.parse_args()

out = pathlib.Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
rep = dict(variant=args.name, drop_hiapo_mature=args.drop_hiapo_mature,
           drop_efib_outer=args.drop_efib_outer, keep_uninduced=args.keep_uninduced, seed=args.seed)

a = sc.read_h5ad(args.h5ad)
sub = a.obs["cell_subtype"].astype(str)
tp = a.obs["timepoint"].astype(str)
print(f"loaded {a.n_obs:,}", flush=True)

# ---- apoptosis scores: NOT on the deliverable; join from the v1 all-cells ----
# pseudotime object (same 43,197 cells, scored by the notebook's AUCell pass).
V1 = ("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/notebooks/pseudotime_mouse/"
      "results/v1_all_cells/2026_08_28_mouse_fibdec_pseudotime.h5ad")
import anndata as ad_mod
v1obs = ad_mod.read_h5ad(V1, backed="r").obs
apo_s = v1obs["apoptosis"].reindex(a.obs_names)
n_missing = int(apo_s.isna().sum())
assert n_missing == 0, f"{n_missing} deliverable cells missing from v1 score join"
a.obs["apoptosis"] = apo_s.to_numpy(float)
rep["apoptosis_source"] = V1

# ---- label-anchored hi-apo threshold, computed BEFORE dropping apoptotic ----
apo = a.obs["apoptosis"].to_numpy(float)
apo_pop = apo[(sub == "decidual_apoptotic").values]
thr = float(np.quantile(apo_pop, 0.25))
rep["hiapo_threshold_q25_of_decidual_apoptotic"] = thr
rep["n_decidual_apoptotic_reference"] = int(len(apo_pop))

# ---- graph-neighbourhood apoptotic fraction (no gene set): for every cell, the fraction of its
# k nearest neighbours in the DELIVERABLE stromal latent (X_scVI, pre-drop) annotated decidual_apoptotic.
from sklearn.neighbors import NearestNeighbors
Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32)
nn = NearestNeighbors(n_neighbors=args.apo_knn + 1, n_jobs=8).fit(Z)
_, nbr = nn.kneighbors(Z)
is_apo = (sub == "decidual_apoptotic").to_numpy()
a.obs["apo_neighbor_frac"] = is_apo[nbr[:, 1:]].mean(1).astype(np.float32)
mat = (sub == "decidual_mature").to_numpy()
rep["apo_neighbor_frac"] = dict(k=args.apo_knn, latent="deliverable X_scVI (pre-drop)",
    n_mature_majority_apoptotic_neighbours=int((mat & (a.obs["apo_neighbor_frac"].to_numpy() > 0.5)).sum()),
    n_mature_any_apoptotic_neighbour=int((mat & (a.obs["apo_neighbor_frac"].to_numpy() > 0)).sum()),
    mature_frac_quantiles=[float(q) for q in np.quantile(a.obs["apo_neighbor_frac"].to_numpy()[mat], [0.5, 0.9, 0.95, 0.99])])
print("apo_neighbor_frac:", json.dumps(rep["apo_neighbor_frac"]), flush=True)

keep = np.ones(a.n_obs, bool)
keep &= ~(sub == "decidual_apoptotic").values          # baseline (v2)
if not args.keep_uninduced:
    keep &= ~(tp == "uninduced").values                # baseline (v2); varD keeps it
rep["n_uninduced_kept"] = int(((tp == "uninduced").values & keep).sum())
rep["n_after_baseline_drops"] = int(keep.sum())

if args.drop_hiapo_mature:
    hi = (sub == "decidual_mature").values & (apo >= thr)
    rep["n_hiapo_mature_dropped"] = int((hi & keep).sum())
    keep &= ~hi
if args.drop_efib_outer:
    eo = (sub == "efib_outer").values
    rep["n_efib_outer_dropped"] = int((eo & keep).sum())
    keep &= ~eo

a = a[keep].copy()
rep["n_cells_final"] = int(a.n_obs)
print(json.dumps({k: v for k, v in rep.items() if k.startswith("n_") or "thresh" in k}), flush=True)

# ---- retrain scVI on exactly these cells (counts layer) ---------------------
import scvi
assert "counts" in a.layers
sc.pp.highly_variable_genes(a, n_top_genes=args.hvg, flavor="seurat_v3", layer="counts")
av = a[:, a.var["highly_variable"]].copy()
scvi.settings.seed = args.seed
scvi.model.SCVI.setup_anndata(av, layer="counts", batch_key="sample")
m = scvi.model.SCVI(av, n_latent=30)
m.train(max_epochs=args.epochs)
a.obsm["X_scVI"] = m.get_latent_representation(av)
rep["scvi"] = dict(n_latent=30, epochs=args.epochs, hvg=args.hvg, batch_key="sample")

sc.pp.neighbors(a, use_rep="X_scVI", random_state=args.seed)
sc.tl.umap(a, random_state=args.seed)
sc.tl.diffmap(a, n_comps=15)
print("embeddings done", flush=True)

# NO root, NO DPT — root_selection notebook takes it from here.
a.uns["variant_embed"] = rep
import os
h5 = out / f"{args.name}_embedded.h5ad"
tmp = h5.with_suffix(".h5ad.tmp")
a.write_h5ad(tmp); os.replace(tmp, h5)
(out / f"{args.name}_report.json").write_text(json.dumps(rep, indent=2, default=str))
print(f"wrote {h5}", flush=True)
