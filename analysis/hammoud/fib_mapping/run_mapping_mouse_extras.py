#!/usr/bin/env python
"""
run_mapping_mouse_extras.py — the mouse-arm analyses that were never run (decision 2026-09-01,
revised same day: user asked to run them).

Writes ONLY to fib_centroid_results_mouse_extras/ (new directory). The finished directories
fib_centroid_results/ and fib_centroid_results_mouse/ are read-only inputs and are not touched.

Analyses (mirroring the human-arm skeleton where it makes sense):
  1. cluster centroids on the full ortholog panel, CACHED (never saved by the original run)
     + a provenance check: the recomputed all-39 matrix is diffed against the finished run's
       stored matrix (finished run used the netscratch copy; this uses data/deliverables)
  2. replication: mapping recomputed per section (9) and per timepoint (6)      [human 2.5]
  3. marker ablation on a discriminative variance panel                          [human 2.6]
  4. panel-matched matrices + bootstrap + held-out on the variance panel
     (tests the panel-dependence of the decidual matches; also attempts to
      reproduce the driver-JSON cosines from first principles)
  5. marker sets (published S6C + centroid-derived) translated to mouse symbols,
     detectability, dotplots                                                     [human 2.7/2.8]
  6. per-cell marker-set scores + score_by_cluster + spatial columns             [human 2.9]

Design notes, stated once:
  - groups are curated fine_refined labels; no clustering happened, so there is no
    resolution sweep here (N/A by design — the human arm swept 0.3–2.0).
  - replication uses cols = the 6 fibroblast subtypes (6-column centring, parallel to
    the human per-section replication). Mouse timepoints are an INDUCTION TIME COURSE:
    cross-timepoint disagreement can be biology; within-timepoint section pairs
    (uninduced a/b, 1_dpi_late a/b, 3_dpi a/b) are the honest replication axis.
  - every gene translation goes through Ensembl 116 strict one-to-one orthologs,
    symbol-unique both ways (same policy as run_mapping_mouse.py). Asserted panel
    floors guard the silent UPPER/Title case-collision trap.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

HERE = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud")
sys.path.insert(0, str(HERE / "fib_mapping"))
import hammoud_fib_mapping as hfm
from hammoud_fib_mapping import FIB_SUBTYPES, SHORT, S6C_TO_S5B

OUT = HERE / "fib_centroid_results_mouse_extras"
OUT.mkdir(exist_ok=True)
MOUT = HERE / "fib_centroid_results_mouse"

H5 = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/deliverables/"
          "2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad")
ORTHO = ("/n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro/orthogroups/"
         "ensembl_human_mouse_components.tsv")
CLUSTER_KEY = "fine_refined"
SEED = 0


def log(msg):
    print(f"[{pd.Timestamp.now():%H:%M:%S}] {msg}", flush=True)


# ---- 0. atlas + orthologs (identical policy to run_mapping_mouse.py) --------
log("atlas + ortholog translation")
atlas = hfm.load_atlas_centroids(HERE / "pnas.2404775121.sd05.xlsx")
markers_s6c = hfm.load_atlas_markers(HERE / "pnas.2404775121.sd06.xlsx")

comp = pd.read_csv(ORTHO, sep="\t")
one = comp[comp["homology_type"] == "ortholog_one2one"]
hs = one[one["species"] == "human"][["component_id", "gene_symbol"]].rename(columns={"gene_symbol": "hs"})
mm = one[one["species"] == "mouse"][["component_id", "gene_symbol"]].rename(columns={"gene_symbol": "mm"})
lut = hs.merge(mm, on="component_id")[["hs", "mm"]].drop_duplicates()
lut = lut[~lut["hs"].duplicated(keep=False) & ~lut["mm"].duplicated(keep=False)]
HS2MM = dict(zip(lut["hs"], lut["mm"]))
assert HS2MM.get("TSPAN6") == "Tspan6", "ortholog table sanity check failed"
log(f"  strict 1:1 pairs: {len(HS2MM):,}")

atlas_mm = atlas.loc[atlas.index.intersection(lut.set_index("hs").index)]
atlas_mm.index = lut.set_index("hs").loc[atlas_mm.index, "mm"]
log(f"  atlas genes surviving translation: {len(atlas_mm):,} of {len(atlas):,}")

# ---- 1. mouse object + cached centroids --------------------------------------
log("loading mouse object (deliverables copy)")
a = sc.read_h5ad(H5)
assert a.n_obs == 43197, f"unexpected object: {a.n_obs} cells"
rs = np.asarray(np.expm1(a.X[:100].toarray() if hasattr(a.X, "toarray") else a.X[:100]).sum(1)).ravel()
assert np.allclose(rs, 10000, rtol=0.01), f"X is not log1p(CP10K): row sums {rs[:3]}"
labels = a.obs[CLUSTER_KEY].astype(str)
log(f"  {a.n_obs:,} cells; groups: {labels.value_counts().to_dict()}")

panel = pd.Index([g for g in atlas_mm.index if g in a.var_names])
assert len(panel) > 15000, f"shared panel collapsed: {len(panel)} (case-collision trap?)"
log(f"  shared ortholog panel: {len(panel):,}")

X_lin = hfm.linear_panel_matrix(a[:, panel], panel)
C, ns = hfm.group_centroids(X_lin, labels.values, panel)
C.to_parquet(OUT / "cluster_centroids_mouse.parquet")
ns.to_csv(OUT / "cluster_sizes_mouse.csv")
log(f"  centroids cached: {C.shape[0]:,} genes x {C.shape[1]} clusters")

# provenance check: does the finished run reproduce from THIS copy of the object?
A_full = atlas_mm.loc[panel]
M_re = hfm.contrastive_cosine(C, A_full)
M_re.to_csv(OUT / "matrix_contrastive_cosine_recheck.csv")
M_orig = pd.read_csv(MOUT / "matrix_contrastive_cosine.csv", index_col=0)
maxdiff = float((M_re - M_orig.loc[M_re.index, M_re.columns]).abs().max().max())
log(f"  PROVENANCE: recomputed all-39 matrix vs finished run, max |diff| = {maxdiff:.2e}")

# ---- 2. replication per section and per timepoint ----------------------------
log("replication per sample")
per_sample = hfm.per_group_mapping(X_lin, a.obs, CLUSTER_KEY, "sample",
                                   atlas_mm, FIB_SUBTYPES, panel, panel)
per_sample.to_csv(OUT / "replication_per_sample_mouse.csv", index=False)
log("replication per timepoint")
per_tp = hfm.per_group_mapping(X_lin, a.obs, CLUSTER_KEY, "timepoint",
                               atlas_mm, FIB_SUBTYPES, panel, panel)
per_tp.to_csv(OUT / "replication_per_timepoint_mouse.csv", index=False)
consistency = (per_sample.groupby("subtype")["best_query"]
               .agg(modal=lambda s: s.value_counts().index[0],
                    agreement=lambda s: s.value_counts(normalize=True).iloc[0],
                    n_samples="size").reset_index())
consistency.to_csv(OUT / "replication_consistency_mouse.csv", index=False)
log("\n" + consistency.round(3).to_string())

# ---- 3+4. discriminative panels: ablation, bootstrap, held-out, matrices -----
# vpanelA: variance panel computed in mouse symbol space (self-consistent here)
vpanelA = hfm.variance_panel(atlas_mm, FIB_SUBTYPES, panel, 3000)
# vpanelB: variance panel computed in HUMAN space then translated — the construction
# the driver JSON most plausibly used (its n_panel=2728 < 3000 suggests translation loss)
vpanelB_h = hfm.variance_panel(atlas, FIB_SUBTYPES, atlas.index, 3000)
vpanelB = pd.Index([HS2MM[g] for g in vpanelB_h if g in HS2MM and HS2MM[g] in a.var_names])
log(f"vpanelA (mouse-space): {len(vpanelA):,} | vpanelB (human-space translated): {len(vpanelB):,}")

for tag, vp in [("vpanelA", pd.Index(vpanelA)), ("vpanelB", vpanelB)]:
    Mv = hfm.contrastive_cosine(C.loc[vp], atlas_mm.loc[vp, FIB_SUBTYPES])
    Mv.to_csv(OUT / f"matrix_contrastive_cosine_{tag}.csv")
log("panel-matched matrices written")

# compare vpanelB cosines to the driver JSON (job-log diagnostic)
DR = json.loads((MOUT / "hammoud_gene_drivers_v5.json").read_text())
MvB = pd.read_csv(OUT / "matrix_contrastive_cosine_vpanelB.csv", index_col=0)
CROSS = {"dec1": "decidual_boundary", "dec2": "decidual_proliferative", "dec3": "decidual_mature",
         "decAp": "decidual_apoptotic", "eFib1": "efib_outer", "eFib2": "efib_middle",
         "eFib2i": "efib_middle_induced", "eFib3": "efib_subluminal", "eFibP": "efib_proliferative"}
print("query          JSON pair -> cosine     vpanelB cosine (same pair)")
for k, v in DR.items():
    if not k.startswith("mouse:"):
        continue
    q, tgt = k.split(":", 1)[1].split("->", 1)
    short = {v2: k2 for k2, v2 in CROSS.items()}[q]
    print(f"{q:24s} -> {tgt[:22]:24s} {v['cosine']:+.4f}   {MvB.loc[short, tgt]:+.4f}")

abl = hfm.marker_ablation(C, atlas_mm, FIB_SUBTYPES, FIB_SUBTYPES, pd.Index(vpanelA))
abl.to_csv(OUT / "robustness_marker_ablation_mouse.csv", index=False)
boot_v = hfm.gene_bootstrap(C, atlas_mm, FIB_SUBTYPES, FIB_SUBTYPES, pd.Index(vpanelA), seed=SEED)
boot_v.to_csv(OUT / "gene_bootstrap_mouse_vpanel.csv", index=False)
held_v = hfm.heldout_split(C, atlas_mm, FIB_SUBTYPES, FIB_SUBTYPES, pd.Index(vpanelA), seed=SEED)
held_v.to_csv(OUT / "heldout_mouse_vpanel.csv", index=False)
log("\n" + boot_v.round(3).to_string())
log("\n" + held_v.round(3).to_string())

# ---- 5. marker sets, detectability, dotplots ----------------------------------
log("marker sets (translated) + detectability")
marker_sets = {}
for s in FIB_SUBTYPES:
    centroid_derived = list(hfm.specific_genes(atlas_mm, s, FIB_SUBTYPES, panel, n=50))
    s6c_label = {v: k for k, v in S6C_TO_S5B.items()}[s]
    pub_h = markers_s6c.query("cluster == @s6c_label").sort_values(
        "avg_log2FC", ascending=False)["gene"]
    published = [HS2MM[g] for g in pub_h if g in HS2MM and HS2MM[g] in a.var_names][:50]
    assert len(published) >= 15, f"{s}: published panel collapsed to {len(published)}"
    marker_sets[s] = {"centroid_derived": centroid_derived, "published": published,
                      "shared": sorted(set(centroid_derived) & set(published))}
pd.DataFrame([
    {"subtype": s, "n_centroid_derived": len(v["centroid_derived"]),
     "n_published_detected": len(v["published"]), "n_shared": len(v["shared"]),
     "centroid_derived": ", ".join(v["centroid_derived"]),
     "published": ", ".join(v["published"]), "shared": ", ".join(v["shared"])}
    for s, v in marker_sets.items()
]).to_csv(OUT / "marker_sets_mouse.csv", index=False)

a.layers["log1p"] = a.X.copy()
det_rows = []
for s, v in marker_sets.items():
    X = sp.csr_matrix(a[:, v["published"]].layers["log1p"])
    det_rows.append({"subtype": SHORT.get(s, s), "n_markers": len(v["published"]),
                     "mean_detection_rate": float(np.asarray((X > 0).mean(axis=0)).ravel().mean())})
pd.DataFrame(det_rows).to_csv(OUT / "marker_panel_detectability_mouse.csv", index=False)
log("\n" + pd.DataFrame(det_rows).round(4).to_string(index=False))

PROGRAMS = {
    "cycling (defines Stro-11)": ["MKI67", "TOP2A", "CENPF", "UBE2C", "BIRC5", "CDK1"],
    "immune (defines Stro-5)": ["PTPRC", "GNLY", "NKG7", "CCL5", "TYROBP", "TPSAB1"],
    "immediate-early / stress (defines Stro-9)": ["FOSB", "JUND", "ATF3", "NR4A1", "IER2"],
    "C7 program (defines Stro-7)": ["C7", "CFH", "CTSK", "APOD", "CCDC80", "IGFBP5"],
    "decidual (defines Stro-8/9/10)": ["IGFBP1", "PRL", "FOXO1", "WNT4"],
}
prog_rows = []
for name, gs_h in PROGRAMS.items():
    g = [HS2MM[x] for x in gs_h if x in HS2MM and HS2MM[x] in a.var_names]
    X = sp.csr_matrix(a[:, g].layers["log1p"])
    prog_rows.append({"program": name,
                      "cells_positive_for_any": float((np.asarray((X > 0).sum(axis=1)).ravel() > 0).mean()),
                      "genes_mouse": ", ".join(g), "n_translated": len(g), "n_requested": len(gs_h)})
pd.DataFrame(prog_rows).to_csv(OUT / "program_detectability_mouse.csv", index=False)
log("\n" + pd.DataFrame(prog_rows)[["program", "cells_positive_for_any", "n_translated"]].round(4).to_string(index=False))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
groups = {SHORT.get(s, s): v["published"][:6] for s, v in marker_sets.items() if v["published"]}
sc.pl.dotplot(a, groups, groupby=CLUSTER_KEY, layer="log1p", use_raw=False,
              standard_scale="var", show=False)
plt.savefig(OUT / "dotplot_mouse_published_markers.png", dpi=200, bbox_inches="tight")
plt.close("all")
log("dotplot written")

# ---- 6. per-cell scores + spatial columns -------------------------------------
log("per-cell marker-set scores")
import anndata as ad
tmp = ad.AnnData(X=a.layers["log1p"], obs=a.obs[[]].copy(), var=a.var[[]].copy())
score_cols = []
for s in FIB_SUBTYPES:
    genes = marker_sets[s]["published"] or marker_sets[s]["centroid_derived"]
    name = "score_" + s.split("_")[0]
    sc.tl.score_genes(tmp, genes, score_name=name, use_raw=False, random_state=SEED)
    a.obs[name] = tmp.obs[name].to_numpy()
    score_cols.append(name)
a.obs[score_cols + [CLUSTER_KEY, "cell_subtype", "sample", "timepoint",
                    "x_um_dbscan", "y_um_dbscan"]].to_csv(OUT / "per_cell_subtype_scores_mouse.csv")
sbc = a.obs.groupby(CLUSTER_KEY, observed=True)[score_cols].mean()
sbc.to_csv(OUT / "score_by_cluster_mouse.csv")
log("\n" + sbc.round(3).to_string())

# ---- 7. meta -------------------------------------------------------------------
meta = dict(object=str(H5), cluster_key=CLUSTER_KEY, n_cells=int(a.n_obs),
            panel_genes=int(len(panel)), atlas_genes_translated=int(len(atlas_mm)),
            vpanelA_genes=int(len(vpanelA)), vpanelB_genes=int(len(vpanelB)),
            provenance_max_abs_diff_vs_finished_run=maxdiff,
            ortholog_policy="Ensembl116 strict one2one, symbol-unique both ways",
            replication_cols="6 FIB_SUBTYPES (6-column centring)",
            seed=SEED, run_date=str(pd.Timestamp.now()))
(OUT / "meta_mouse_extras.json").write_text(json.dumps(meta, indent=2))
# read-back check (the json trap: verify by reading the file back)
assert json.loads((OUT / "meta_mouse_extras.json").read_text())["n_cells"] == 43197
log("DONE")
