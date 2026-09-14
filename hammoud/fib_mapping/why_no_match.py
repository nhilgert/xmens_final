"""Why do Stro-5, Stro-9 and Stro-11 fail to map?

Absence of a counterpart is only interesting if we can say what drives it. Each of
those three subtypes is defined by a program that our object either excludes by
construction or measures differently, and this quantifies that rather than
leaving "no counterpart" as an unexplained negative.

Also rebuilds the marker dot-plot panels with a detectability filter, since
ranking the published markers by fold-change alone surfaces generic ECM genes
that do not discriminate between our clusters.
"""

import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt

HAM = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud")
OUT = HAM / "fib_centroid_results"
sys.path.insert(0, str(HAM / "fib_mapping"))
from hammoud_fib_mapping import FIB_SUBTYPES, SHORT, load_atlas_centroids, specificity

atlas = load_atlas_centroids(HAM / "pnas.2404775121.sd05.xlsx")
human = sc.read_h5ad("/n/netscratch/eddy_lab/Everyone/nhilgert/"
                     "2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad")
human.layers["log1p"] = human.X.copy()
assign = pd.read_csv(OUT / "cluster_assignments.csv", index_col=0)
CLUSTER_KEY = [c for c in assign.columns if c.startswith("leiden")][0]
human.obs[CLUSTER_KEY] = assign[CLUSTER_KEY].reindex(human.obs_names).astype(str).values

common = human.var_names.intersection(atlas.index)
marker_sets = pd.read_csv(OUT / "marker_sets.csv")


def detect_rate(genes):
    genes = [g for g in genes if g in human.var_names]
    if not genes:
        return np.nan, []
    X = sp.csr_matrix(human[:, genes].layers["log1p"])
    return float(np.asarray((X > 0).mean(axis=0)).ravel().mean()), genes


print("=" * 78)
print("Detectability of each subtype's published marker panel in our nuclei")
print("=" * 78)
rows = []
for _, r in marker_sets.iterrows():
    pub = [g for g in str(r.published_s6c).split(", ") if g]
    rate, present = detect_rate(pub)
    rows.append({"subtype": SHORT.get(r.subtype, r.subtype),
                 "n_markers": len(pub), "n_in_object": len(present),
                 "mean_detection_rate": rate})
det = pd.DataFrame(rows)
print(det.round(4).to_string(index=False))
det.to_csv(OUT / "marker_panel_detectability.csv", index=False)

print("\n" + "=" * 78)
print("The three programs that cannot map, quantified")
print("=" * 78)

PROGRAMS = {
    "cycling (defines Stro-11)": ["MKI67", "TOP2A", "CENPF", "UBE2C", "BIRC5", "CDK1", "ASPM"],
    "immune (defines Stro-5)": ["PTPRC", "GNLY", "NKG7", "CCL5", "TYROBP", "TPSAB1", "GZMA"],
    "immediate-early / stress (defines Stro-9)": ["FOSB", "JUND", "ATF3", "NR4A1", "IER2", "SOCS3"],
    "C7 program (defines Stro-7)": ["C7", "CFH", "CTSK", "APOD", "CCDC80", "IGFBP5"],
    "decidual (defines Stro-8/9/10)": ["IGFBP1", "PRL", "FOXO1", "WNT4"],
}
for name, genes in PROGRAMS.items():
    g = [x for x in genes if x in human.var_names]
    X = sp.csr_matrix(human[:, g].layers["log1p"])
    per_gene = np.asarray((X > 0).mean(axis=0)).ravel()
    any_pos = float(np.asarray((X > 0).sum(axis=1)).ravel().mean())
    frac_any = float((np.asarray((X > 0).sum(axis=1)).ravel() > 0).mean())
    print(f"\n{name}")
    print(f"  nuclei positive for >=1 of the panel: {frac_any:.2%}")
    print("  " + ", ".join(f"{gg}={pp:.3%}" for gg, pp in zip(g, per_gene)))

print("\n" + "=" * 78)
print("Rebuilt marker panels: specific in the atlas AND detectable in our nuclei")
print("=" * 78)

X_all = sp.csr_matrix(human[:, common].layers["log1p"])
det_rate = pd.Series(np.asarray((X_all > 0).mean(axis=0)).ravel(), index=common)

panels = {}
for s in FIB_SUBTYPES:
    spec = specificity(atlas.loc[common], s, FIB_SUBTYPES)
    ok = spec[(atlas.loc[spec.index, s] > 0.5) & (det_rate.reindex(spec.index) > 0.10)]
    panels[SHORT[s]] = list(ok.head(6).index)
    print(f"{SHORT[s]:20s} {', '.join(panels[SHORT[s]]) if len(panels[SHORT[s]]) else '(no gene passes the detection filter)'}")

pd.DataFrame([{"subtype": k, "genes": ", ".join(v)} for k, v in panels.items()]).to_csv(
    OUT / "marker_panels_detectable.csv", index=False)

panels_nonempty = {k: v for k, v in panels.items() if v}
if panels_nonempty:
    sc.pl.dotplot(human, panels_nonempty, groupby=CLUSTER_KEY, layer="log1p",
                  use_raw=False, standard_scale="var", show=False)
    plt.savefig(OUT / "fig7_detectable_markers_dotplot.png", dpi=200, bbox_inches="tight")
    print(f"\nwrote {OUT / 'fig7_detectable_markers_dotplot.png'}")
