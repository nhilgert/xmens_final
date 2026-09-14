#!/usr/bin/env python
"""
run_mapping_mouse.py — the mouse/XMens arm of the Hammoud centroid mapping.

Reuses hammoud_fib_mapping.py unchanged (same statistics, same robustness
checks, same do-not-undo decisions: cluster centroids not per-cell, contrastive
cosine as headline, log1p(mean(expm1)) centroids, gene-space bootstraps).

Differences from the human arm, and only these:
  1. Gene axis: atlas centroids are human symbols; the object is mouse. The
     atlas panel is translated through Ensembl 116 STRICT one2one orthologs
     (human symbol -> mouse symbol). One2many/paralog families are excluded --
     conservative; a SAMap-graph variant is the follow-up if coverage hurts.
  2. Groups: obs['fine_refined'] (eFib1/2/2i/3, eFibP, dec1/2/3, decAp) --
     no clustering step, per the handoff.
  3. Cross-species + cross-platform: the double-centering matters even more
     (species offset on top of nucleus/cell offset). Raw Spearman still
     reported for continuity.
decAp (apoptotic) is kept as a row -- informative, flagged in outputs; it is
excluded from nothing here because this is centroid correlation, not an
embedding (the apoptotic-drop rule applies to embeddings).

Outputs -> fib_centroid_results_mouse/
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud")
sys.path.insert(0, str(HERE / "fib_mapping"))
import hammoud_fib_mapping as hfm

OUT = HERE / "fib_centroid_results_mouse"
OUT.mkdir(exist_ok=True)

H5 = "/n/netscratch/eddy_lab/Everyone/nhilgert/2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"
ORTHO = "/n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv"
CLUSTER_KEY = "fine_refined"

# ---- atlas ------------------------------------------------------------------
atlas = hfm.load_atlas_centroids(HERE / "pnas.2404775121.sd05.xlsx")
cols = list(atlas.columns)
print(f"[atlas] {atlas.shape[0]} genes x {len(cols)} subtypes", flush=True)

# ---- ortholog translation of the atlas gene axis ----------------------------
comp = pd.read_csv(ORTHO, sep="\t")
one = comp[comp.homology_type == "ortholog_one2one"]
hs = one[one.species == "human"][["component_id", "gene_symbol"]].rename(columns={"gene_symbol": "hs"})
mm = one[one.species == "mouse"][["component_id", "gene_symbol"]].rename(columns={"gene_symbol": "mm"})
lut = hs.merge(mm, on="component_id")[["hs", "mm"]].drop_duplicates()
lut = lut[~lut.hs.duplicated(keep=False) & ~lut.mm.duplicated(keep=False)]  # enforce strict 1:1 both ways
atlas_mm = atlas.loc[atlas.index.intersection(lut.set_index("hs").index)]
atlas_mm.index = lut.set_index("hs").loc[atlas_mm.index, "mm"]
print(f"[ortho] atlas genes surviving strict 1:1 translation: {len(atlas_mm):,} "
      f"of {len(atlas):,}", flush=True)

# ---- mouse object ------------------------------------------------------------
a = sc.read_h5ad(H5)
assert CLUSTER_KEY in a.obs.columns, f"{CLUSTER_KEY} missing: {list(a.obs.columns)[:30]}"
rs = np.asarray(np.expm1(a.X[:100].toarray() if hasattr(a.X, "toarray") else a.X[:100]).sum(1)).ravel()
assert np.allclose(rs, 10000, rtol=0.01), f"X is not log1p(CP10K): row sums {rs[:3]}"
labels = a.obs[CLUSTER_KEY].astype(str)
print(f"[mouse] {a.n_obs:,} cells; groups: {labels.value_counts().to_dict()}", flush=True)

panel = pd.Index([g for g in atlas_mm.index if g in a.var_names])
print(f"[panel] shared genes: {len(panel):,}", flush=True)

X_lin = hfm.linear_panel_matrix(a[:, panel], panel)
C, ns = hfm.group_centroids(X_lin, labels.values, panel)

A = atlas_mm.loc[panel]
M_cos = hfm.contrastive_cosine(C, A)
M_rho = hfm.rank_correlation(C, A)
by_q, by_a = hfm.correspondence(M_cos, ns=ns)
by_q.to_csv(OUT / "TABLE_correspondence_mouse_by_query.csv")
by_a.to_csv(OUT / "TABLE_correspondence_mouse.csv")   # by-atlas view = the deliverable shape
M_cos.to_csv(OUT / "matrix_contrastive_cosine.csv")
M_rho.to_csv(OUT / "matrix_raw_spearman.csv")
print("=== by query (mouse fine_refined) ===", flush=True)
print(by_q.to_string(), flush=True)
print("=== by atlas subtype ===", flush=True)
print(by_a.to_string(), flush=True)

targets = [c for c in cols if c in hfm.FIB_SUBTYPES]
boot = hfm.gene_bootstrap(C, A, targets, cols, panel)
boot.to_csv(OUT / "gene_bootstrap_mouse.csv")
held = hfm.heldout_split(C, A, targets, cols, panel)
held.to_csv(OUT / "heldout_mouse.csv")
print(boot.to_string(), flush=True)
print(held.to_string(), flush=True)

meta = dict(object=H5, cluster_key=CLUSTER_KEY, n_cells=int(a.n_obs),
            panel_genes=len(panel), atlas_genes_translated=int(len(atlas_mm)),
            ortholog_policy="Ensembl116 strict one2one, symbol-unique both ways",
            note="cross-species arm; contrastive cosine headline, raw Spearman for continuity")
(OUT / "meta_mouse.json").write_text(json.dumps(meta, indent=2))
print("DONE", flush=True)
