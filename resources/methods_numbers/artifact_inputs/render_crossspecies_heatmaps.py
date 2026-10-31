#!/usr/bin/env python
"""Cross-species pseudotime heatmap (human panel left, mouse panel right, same row order), Monocle smoother.

Standalone replacement for the session-kernel code that drew
conserved_profiles_heatmap_28_monocle_clip3_celltypes_notitle.svg (2026-10-20). Geometry, colours and type
sizes were read off that SVG.

Inputs
  --curves-human / --curves-mouse : raw response curves from monocle_conserved_smooth.R
                                    (TSV: gene, fit, 100 columns = evenly spaced pseudotime points)
  --genes                         : CSV with columns human, mouse; row order = heatmap order (top to bottom)
  --cells-human / --cells-mouse   : TSV with columns cell, pt_rank_norm, label (one row per cell in the run)
Transform (Monocle 2 plot_pseudotime_heatmap): log10(curve + 1), z per gene with R scale() (sd, n - 1), clip to ±zclip.
Strips: 100 equal-count bins of pt_rank_norm (stable argsort, np.array_split), stacked share of each cell type
per bin (composition) and the most common type per bin (plurality).

Outputs  <out>.svg/.pdf/.png (600 dpi), <out>_z_human.csv, <out>_z_mouse.csv, <out>_composition_{human,mouse}.csv,
         <out>_plurality.csv, <out>_provenance.json
"""
import argparse, hashlib, json, os, datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize
from matplotlib.patches import Patch
import logging
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

ap = argparse.ArgumentParser()
ap.add_argument("--curves-human", required=True); ap.add_argument("--curves-mouse", required=True)
ap.add_argument("--genes", required=True)
ap.add_argument("--cells-human", required=True); ap.add_argument("--cells-mouse", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--zclip", type=float, default=3.0)
ap.add_argument("--nbins", type=int, default=100)
ap.add_argument("--note", default="")
a = ap.parse_args()

HN = "/System/Library/Fonts/HelveticaNeue.ttc"
if os.path.exists(HN):
    fm.fontManager.addfont(HN)
plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Nimbus Sans", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "xtick.major.size": 2.0, "xtick.major.pad": 1.5})
MM = 1 / 25.4
W_MM, H_MM = 142.5, 150.4

HUMAN_TYPES = [("endometrial_fibroblast", "human eFib", "#2171b5"), ("decidual_cell", "human decidual", "#e6550d")]
MOUSE_TYPES = [("efib_outer", "mouse eFib outer", "#08306b"), ("efib_middle", "mouse eFib middle", "#2171b5"),
               ("efib_middle_induced", "mouse eFib middle (induced)", "#6baed6"),
               ("efib_subluminal", "mouse eFib subluminal", "#9ecae1"),
               ("efib_proliferative", "mouse eFib proliferative", "#41ab5d"),
               ("decidual_boundary", "mouse dec boundary", "#fdae6b"),
               ("decidual_proliferative", "mouse dec proliferative", "#e6550d"),
               ("decidual_mature", "mouse dec mature", "#a63603")]


def monocle_z(curves_tsv, genes, zclip):
    C = pd.read_csv(curves_tsv, sep="\t", index_col=0)
    C = C.drop(columns=[c for c in C.columns if c == "fit"])
    missing = [g for g in genes if g not in C.index]
    assert not missing, f"{curves_tsv} lacks {missing}"
    L = np.log10(C.loc[genes].to_numpy(float) + 1)
    Z = (L - L.mean(1, keepdims=True)) / L.std(1, ddof=1, keepdims=True)
    return pd.DataFrame(Z, index=genes), int((np.abs(Z) > zclip).sum())


def composition(cells_tsv, types, nbins):
    d = pd.read_csv(cells_tsv, sep="\t")
    labs = [t[0] for t in types]
    unknown = set(d["label"]) - set(labs)
    assert not unknown, f"{cells_tsv}: labels without a colour {unknown}"
    ix = np.argsort(d["pt_rank_norm"].to_numpy(), kind="stable")
    lab = d["label"].to_numpy()
    comp = np.vstack([pd.Series(lab[b]).value_counts(normalize=True).reindex(labs, fill_value=0).to_numpy()
                      for b in np.array_split(ix, nbins)])
    return pd.DataFrame(comp, columns=labs), len(d)


gs = pd.read_csv(a.genes)
genes_h, genes_m = gs["human"].astype(str).tolist(), gs["mouse"].astype(str).tolist()
Zh, nclip_h = monocle_z(a.curves_human, genes_h, a.zclip)
Zm, nclip_m = monocle_z(a.curves_mouse, genes_m, a.zclip)
Ch, n_h = composition(a.cells_human, HUMAN_TYPES, a.nbins)
Cm, n_m = composition(a.cells_mouse, MOUSE_TYPES, a.nbins)
nb = a.nbins; ng = len(genes_h)

fig = plt.figure(figsize=(W_MM * MM, H_MM * MM))
def axes_mm(x, y_top, w, h):              # mm from the top-left corner of the canvas
    return fig.add_axes([x / W_MM, 1 - (y_top + h) / H_MM, w / W_MM, h / H_MM])

X0 = {"human": 18.17, "mouse": 54.17}; BOX_W, BOX_H, Y_HEAT = 32.0, 120.0, 14.59
Y_PLUR, H_PLUR, Y_COMP, H_COMP = 11.89, 1.5, 6.29, 5.0
xe = np.linspace(0, 1, nb + 1); ye = np.arange(ng + 1)
norm = Normalize(-a.zclip, a.zclip); cmap = plt.get_cmap("RdBu_r")
for sp, Z, labels, comp, types in (("human", Zh, genes_h, Ch, HUMAN_TYPES), ("mouse", Zm, genes_m, Cm, MOUSE_TYPES)):
    ax = axes_mm(X0[sp], Y_HEAT, BOX_W, BOX_H)
    ax.pcolormesh(xe, ye, np.clip(Z.to_numpy(), -a.zclip, a.zclip), cmap=cmap, norm=norm, shading="flat",
                  antialiased=False, linewidth=0, edgecolors="none")
    ax.set_xlim(0, 1); ax.set_ylim(ng, 0)
    ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"], fontsize=6)
    ax.set_yticks(np.arange(ng) + 0.5); ax.set_yticklabels(labels, fontsize=6, style="italic")
    ax.tick_params(axis="y", length=0, pad=3.5)
    if sp == "mouse":
        ax.yaxis.tick_right()
    # composition strip (stacked bars, one per bin) and plurality strip
    axc = axes_mm(X0[sp], Y_COMP, BOX_W, H_COMP)
    bottom = np.zeros(nb); xc = (xe[:-1] + xe[1:]) / 2
    for lab, _, col in types:
        v = comp[lab].to_numpy()
        axc.bar(xc, v, width=1 / nb, bottom=bottom, color=col, linewidth=0, antialiased=False)
        bottom += v
    axc.set_xlim(0, 1); axc.set_ylim(0, 1); axc.axis("off")
    axc.set_title("Human" if sp == "human" else "Mouse", loc="left", fontsize=8, pad=1.5)
    axp = axes_mm(X0[sp], Y_PLUR, BOX_W, H_PLUR)
    plur = comp.to_numpy().argmax(1)
    axp.pcolormesh(xe, [0, 1], plur[None, :], cmap=ListedColormap([t[2] for t in types]),
                   norm=BoundaryNorm(np.arange(len(types) + 1) - 0.5, len(types)), shading="flat",
                   antialiased=False, linewidth=0, edgecolors="none")
    axp.set_xlim(0, 1); axp.axis("off")
    if sp == "human":
        for yy, hh, txt in ((Y_COMP, H_COMP, "composition"), (Y_PLUR, H_PLUR, "plurality")):
            fig.text((X0["human"] - 1.8) / W_MM, 1 - (yy + hh / 2) / H_MM, txt, fontsize=7, color="#888888",
                     ha="right", va="center")
fig.text((X0["human"] + X0["mouse"] + BOX_W) / 2 / W_MM, 1 - 146.0 / H_MM, "pseudotime (rank-normalized)",
         fontsize=8, ha="center", va="baseline")
# colorbar: plain rectangle, integer ticks
cax = axes_mm(102.17, 105.79, 3.0, 26.4)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, extend="neither")
cb.set_ticks(np.arange(-int(a.zclip), int(a.zclip) + 1))
cb.set_ticklabels([("−" if t < 0 else "") + str(abs(t)) for t in range(-int(a.zclip), int(a.zclip) + 1)])
cb.ax.tick_params(labelsize=6, width=0.5, length=2, pad=1.5); cb.outline.set_linewidth(0.5)
cb.set_label(f"z per gene, clipped at ±{a.zclip:g}\n(Monocle: log10, df = 3 spline)", fontsize=8)
# legend
handles = [Patch(facecolor=c, edgecolor="none", label=l) for _, l, c in HUMAN_TYPES + MOUSE_TYPES]
leg = fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(102.17 / W_MM, 1 - 20.5 / H_MM), frameon=False,
                 fontsize=7, title="cell type of each pseudotime bin:\ncomposition (stacked) and\nplurality (most common type)",
                 title_fontsize=7, alignment="left", handlelength=1.4, handleheight=1.0, borderaxespad=0, borderpad=0)

for ext in ("svg", "pdf"):
    fig.savefig(f"{a.out}.{ext}")
fig.savefig(f"{a.out}.png", dpi=600)
plt.close(fig)

Zh.set_axis(genes_h).rename_axis("gene").to_csv(f"{a.out}_z_human.csv")
Zm.set_axis(genes_m).rename_axis("gene").to_csv(f"{a.out}_z_mouse.csv")
Ch.rename_axis("bin").to_csv(f"{a.out}_composition_human.csv"); Cm.rename_axis("bin").to_csv(f"{a.out}_composition_mouse.csv")
pd.DataFrame({"col": np.arange(nb), "plurality_human": [HUMAN_TYPES[i][0] for i in Ch.to_numpy().argmax(1)],
              "plurality_mouse": [MOUSE_TYPES[i][0] for i in Cm.to_numpy().argmax(1)]}).to_csv(f"{a.out}_plurality.csv", index=False)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
prov = dict(written=datetime.datetime.now().isoformat(timespec="seconds"), note=a.note,
            inputs={k: dict(path=os.path.abspath(v), sha256=sha(v)) for k, v in
                    dict(curves_human=a.curves_human, curves_mouse=a.curves_mouse, genes=a.genes,
                         cells_human=a.cells_human, cells_mouse=a.cells_mouse).items()},
            n_cells=dict(human=n_h, mouse=n_m), n_genes=ng, zclip=a.zclip, clipped_cells=dict(human=nclip_h, mouse=nclip_m))
json.dump(prov, open(f"{a.out}_provenance.json", "w"), indent=1)
print(json.dumps(dict(out=a.out, n_cells=prov["n_cells"], clipped=prov["clipped_cells"])))
