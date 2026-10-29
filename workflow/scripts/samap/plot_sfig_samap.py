"""Supplementary figure: transfer of mouse pseudotime onto human cells through the SAMap alignment,
and permutation controls for that transfer (panels a-g, 180 mm wide).

  a-c  joint SAMap UMAP: mouse pseudotime (source), values transferred onto human cells, human pseudotime
  d    share of cross-species edge weight that links matching compartments
  e    transferred vs human pseudotime, all cells and by compartment
  f    transferred pseudotime in the human sections: observed, and one transfer after shuffling mouse
       pseudotime across all mouse cells or only within its compartment
  g    observed statistics against 1,000 shuffles of each kind (Spearman rho, Moran's I)

Display rules, as in Fig. 6 and the methods: pseudotime on one fixed 0-1 viridis scale with values as
stored (never re-ranked), and nuclei flagged embedding_outlier left out of UMAPs and spatial maps. Every
statistic (rho in e, all of g) uses all nuclei.

Called by rule samap_supp_figure (workflow/rules/10c_samap_figures.smk). Inputs:
  --joint              joint.parquet (export_samap_tables.py): species and UMAP of every joint nucleus
  --mouse-cells        mouse_cells.csv: pt_harmonized_rank (NaN for apoptotic), labels, embedding_outlier
  --human-cells        human_cells.csv: human pt_harmonized_rank, the transfer of record (pt_transferred), labels,
                       x_um_qc/y_um_qc, embedding_outlier
  --edge-share         samap_compartment_edge_share.csv (samap_edge_share.py)
  --null-replicates    null_replicates.csv and example_transfers.csv (ptnull_permutation.py)
  --example-transfers
Outputs, in --out: sfig_samap.png (300 dpi), sfig_samap_600dpi.png, sfig_samap.svg (vector text, per-panel SVGs
nested at their grid positions), panels/, sfig_samap.provenance.json. Refuses to overwrite an existing
sfig_samap.png unless --force.

Font: Helvetica Neue where installed (the version of record was drawn with it on macOS); otherwise the first of
Helvetica, Nimbus Sans, Arial, Liberation Sans and DejaVu Sans, recorded in the provenance file.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch
from PIL import Image, ImageDraw, ImageFont
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("--joint", required=True)
ap.add_argument("--mouse-cells", required=True)
ap.add_argument("--human-cells", required=True)
ap.add_argument("--edge-share", required=True)
ap.add_argument("--null-replicates", required=True)
ap.add_argument("--example-transfers", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--force", action="store_true")
args = ap.parse_args()
OUT = args.out
FIG_OUT = os.path.join(OUT, "sfig_samap.png")
if os.path.exists(FIG_OUT) and not args.force:
    sys.exit(f"REFUSE: {FIG_OUT} exists (pass --force to regenerate)")
os.makedirs(os.path.join(OUT, "panels"), exist_ok=True)


# ---- style: role-mapped size ladder (titles/labels 9 pt, legends 8 pt, ticks 7 pt) ----
PREF = ["Helvetica Neue", "Helvetica", "Nimbus Sans", "Arial", "Liberation Sans", "DejaVu Sans"]
_have = {f_.name for f_ in fm.fontManager.ttflist}
FONT = next(n for n in PREF if n in _have)
print(f"font: {FONT}" + ("" if FONT == "Helvetica Neue" else "  (Helvetica Neue not installed; fallback)"))
mpl.rcParams.update({
    "font.family": [FONT, "DejaVu Sans"],          # DejaVu supplies glyphs the main font lacks (e.g. the arrow)
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": True, "axes.spines.bottom": True,
    "axes.grid": False, "legend.frameon": False, "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.titleweight": "normal", "axes.titlelocation": "left", "axes.labelweight": "normal",
    "lines.linewidth": 1.2, "patch.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
})


def bold_face(family):
    """(path, ttc_index) of the bold face of `family`, for the PIL-drawn panel letters."""
    cands = [f_ for f_ in fm.fontManager.ttflist if f_.name == family and f_.weight in (700, "bold") and f_.style == "normal"]
    for f_ in cands:
        if f_.fname.lower().endswith(".ttc"):
            for i in range(32):
                try:
                    if ImageFont.truetype(f_.fname, 20, index=i).getname()[1] == "Bold":
                        return f_.fname, i
                except Exception:
                    break
        else:
            return f_.fname, 0
    return os.path.join(mpl.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf"), 0


LETTER_FONT, LETTER_INDEX = bold_face(FONT)


# ---- data ----
SOURCES = [args.joint, args.mouse_cells, args.human_cells, args.edge_share, args.null_replicates, args.example_transfers]

# Joint UMAP in the object's nucleus order (the seeded draw orders below depend on it).
UM = (pd.read_parquet(args.joint, columns=["species", "X_umap__0", "X_umap__1"])
      .rename(columns={"X_umap__0": "umap1", "X_umap__1": "umap2"}))
UM["species"] = UM["species"].astype(str)
assert UM.index.is_unique

MOUSE = pd.read_csv(args.mouse_cells, index_col="cell")
HUMAN = pd.read_csv(args.human_cells, index_col="cell")

# Stored global rank (average rank / N over the 42,290 non-apoptotic nuclei); not re-ranked here.
MS = MOUSE[["cell_type", "cell_subtype"]].assign(pt_rank=MOUSE["pt_harmonized_rank"])

# Human pseudotime (pt_human) and the transfer of record (pt_ranknorm), one row per human nucleus.
Jn = HUMAN[["pt_transferred", "pt_harmonized_rank", "cell_type"]].rename(
    columns={"pt_transferred": "pt_ranknorm", "pt_harmonized_rank": "pt_human"})
SPATIAL = HUMAN[["timepoint", "x_um_qc", "y_um_qc"]]
ES = pd.read_csv(args.edge_share, index_col=[0, 1])

# Permutation controls: one example transfer per design, and 1,000 replicate statistics per design.
EX = pd.read_csv(args.example_transfers).set_index("cell")
NULLS = pd.read_csv(args.null_replicates)

# The permutation run must have transferred the same pseudotime as the transfer of record shown in b and e.
_observed = EX["observed"].reindex(Jn.index)
if not np.allclose(_observed, Jn["pt_ranknorm"], atol=1e-9, equal_nan=True):
    sys.exit("REFUSE: ptnull observed transfer differs from pt_ranknorm; the inputs come from different runs")


HS_FLAGGED = set(HUMAN.index[HUMAN["embedding_outlier"].astype(bool)])
MM_FLAGGED = set(MOUSE.index[MOUSE["embedding_outlier"].astype(bool)])
outlier = ((UM.species == "hs") & UM.index.isin(HS_FLAGGED)) | ((UM.species == "mm") & UM.index.isin(MM_FLAGGED))
UMd = UM[~outlier]  # cells drawn on the UMAPs
print(f"embedding outliers left out of the UMAPs: {int((UM.species[outlier] == 'mm').sum())} mouse, "
      f"{int((UM.species[outlier] == 'hs').sum())} human")

mm = UMd[UMd.species == "mm"].join(MS[["pt_rank", "cell_type", "cell_subtype"]], how="left")
hs = UMd[UMd.species == "hs"].join(Jn[["pt_ranknorm", "pt_human", "cell_type"]], how="left")

# Colors
FIB, DEC = "#3a78b8", "#e0836b"
CTC = {"endometrial_fibroblast": FIB, "decidual_cell": DEC}
EMPTY = "#d4d4d4"
CM = mpl.colormaps["viridis"]
NORM = mpl.colors.Normalize(0, 1)  # fixed scale: values are drawn as stored
S = 0.22

# UMAP bounds, from the drawn cells
xy = UMd[["umap1", "umap2"]].values
pad = 0.02 * (xy.max(0) - xy.min(0))
XL = (xy[:, 0].min() - pad[0], xy[:, 0].max() + pad[0])
YL = (xy[:, 1].min() - pad[1], xy[:, 1].max() + pad[1])
ASP = (XL[1] - XL[0]) / (YL[1] - YL[0])

rs = lambda df, k: df.sample(frac=1, random_state=k)  # seeded draw order, so no group always sits on top
mm_s = rs(mm, 1)
hs_s = rs(hs, 3)


# ---- layout: 12-column grid, panels sized in whole columns ----
outline = {
    "claim": ("Mouse pseudotime carried across the SAMap alignment recovers human pseudotime ordering, "
              "strongly within decidual cells, and mouse pseudotime shuffled before the transfer does not."),
    "width_mm": 180, "ncol": 12, "row_heights_mm": [60, 55, 84],
    "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 4, "chart_family": "embedding", "message": "mouse source values"},
        {"letter": "b", "row": 0, "col": 4, "colspan": 4, "chart_family": "embedding", "message": "transferred values"},
        {"letter": "c", "row": 0, "col": 8, "colspan": 4, "chart_family": "embedding", "message": "human values"},
        {"letter": "d", "row": 1, "col": 0, "colspan": 4, "chart_family": "100% stacked bars", "message": "edges link matching compartments"},
        {"letter": "e", "row": 1, "col": 4, "colspan": 8, "chart_family": "scatter + binned median", "message": "order preserved by compartment"},
        {"letter": "f", "row": 2, "col": 0, "colspan": 7, "chart_family": "spatial maps", "message": "shuffles keep at most the compartment split"},
        {"letter": "g", "row": 2, "col": 7, "colspan": 5, "chart_family": "null strips", "message": "observed vs 1,000 shuffles"},
    ],
}

MMpx = 300 / 25.4
TX = 5.2  # mm from a panel's left edge to its title, leaving room for the letter


def grid_geom(outline, dpi=300, gutter_mm=4):
    mm_ = dpi / 25.4
    W = int(outline["width_mm"] * mm_)
    ncol = outline["ncol"]
    g = int(gutter_mm * mm_)
    colw = (W - g * (ncol - 1)) // ncol
    rowh = [int(h * mm_) for h in outline["row_heights_mm"]]
    row_y = [sum(rowh[:i]) + g * i for i in range(len(rowh))]
    return W, ncol, colw, rowh, row_y, g


def panel_px(outline, letter, dpi=300, gutter_mm=4):
    W, ncol, colw, rowh, row_y, g = grid_geom(outline, dpi, gutter_mm)
    p = next(q for q in outline["panels"] if q["letter"] == letter)
    cs, rs_, r = p["colspan"], p.get("rowspan", 1), p["row"]
    return colw * cs + g * (cs - 1), sum(rowh[r:r + rs_]) + g * (rs_ - 1)


def panel_xy(outline, letter, dpi=300, gutter_mm=4):
    W, ncol, colw, rowh, row_y, g = grid_geom(outline, dpi, gutter_mm)
    p = next(q for q in outline["panels"] if q["letter"] == letter)
    return p["col"] * (colw + g), row_y[p["row"]]


def pfig(L):
    """Empty figure of panel L's exact size; returns (fig, width_mm, height_mm)."""
    w, h = panel_px(outline, L)
    fig = plt.figure(figsize=(w / 300, h / 300), dpi=300)
    return fig, w / MMpx, h / MMpx


def title(fig, W, H, text):
    fig.text(TX / W, 1 - 1.0 / H, text, ha="left", va="top", fontsize=9)


# ---- panels a-c: joint UMAP ----
def umap_panel(L, title_text, draw, key=None, cbar=None):
    fig, W, H = pfig(L)
    top_mm, bot_mm, left_mm = 10.0, 4.0, 6.0
    ah = (H - top_mm - bot_mm) / H
    aw = ah * H * ASP / W
    ax = fig.add_axes([left_mm / W, bot_mm / H, aw, ah])
    ax.set_xlim(XL); ax.set_ylim(YL); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    [s_.set_visible(False) for s_ in ax.spines.values()]
    title(fig, W, H, title_text)
    draw(ax)
    dx, dy = XL[1] - XL[0], YL[1] - YL[0]
    L0 = 0.16 * dx
    kw = dict(arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#555555", mutation_scale=6,
                              shrinkA=0, shrinkB=0, capstyle="projecting", joinstyle="miter"), annotation_clip=False)
    ax.annotate("", xy=(XL[0] + L0, YL[0]), xytext=(XL[0], YL[0]), **kw)
    ax.annotate("", xy=(XL[0], YL[0] + L0), xytext=(XL[0], YL[0]), **kw)
    ax.text(XL[0] + 0.02 * dx, YL[0] - 0.012 * dy, "UMAP 1", ha="left", va="top", fontsize=7, color="#555555")
    ax.text(XL[0] - 0.03 * dx, YL[0] + 0.02 * dy, "UMAP 2", ha="right", va="bottom", rotation=90, fontsize=7, color="#555555")
    rx = (left_mm + aw * W + 3.0) / W
    if cbar:
        cax = fig.add_axes([rx, 0.42, 2.2 / W, 0.30])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CM), cax=cax, ticks=[0, 0.5, 1])
        cb.ax.set_yticklabels(["0 early", "0.5", "1 late"])
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=1.5, pad=1.5)
        cax.set_title(cbar, fontsize=8, loc="left", pad=8)
    if key:
        top = 0.33 if cbar else 1 - top_mm / H
        kax = fig.add_axes([rx - 1.2 / W, bot_mm / H, 1 - rx, top - bot_mm / H])
        kax.axis("off")
        kax.legend(
            handles=[Line2D([], [], ls="none", marker="o", ms=3.2, color=c_, label=l_) for c_, l_ in key],
            loc="upper left", frameon=False, fontsize=8, handletextpad=0.2,
            labelspacing=0.6, borderaxespad=0, handlelength=0.9
        )
    return fig


def draw_mouse_source(ax):
    ax.scatter(hs_s.umap1, hs_s.umap2, s=S, color=EMPTY, linewidths=0, rasterized=True)
    apoptotic = mm_s[mm_s.pt_rank.isna()]
    valued = mm_s[mm_s.pt_rank.notna()]
    ax.scatter(apoptotic.umap1, apoptotic.umap2, s=S, color="#b5b5b5", linewidths=0, rasterized=True)
    ax.scatter(valued.umap1, valued.umap2, s=S, c=CM(NORM(valued.pt_rank.values)), linewidths=0, rasterized=True)


def draw_human(col):
    def draw(ax):
        ax.scatter(mm_s.umap1, mm_s.umap2, s=S, color=EMPTY, linewidths=0, rasterized=True)
        valued = hs_s[hs_s[col].notna()]
        empty = hs_s[hs_s[col].isna()]
        ax.scatter(valued.umap1, valued.umap2, s=S, c=CM(NORM(valued[col].values)), linewidths=0, rasterized=True)
        if len(empty):
            ax.scatter(empty.umap1, empty.umap2, s=S * 2.5, color="#8c8c8c", linewidths=0, rasterized=True, zorder=3)
    return draw


P = {}
P["a"] = umap_panel("a", "Source: mouse pseudotime", draw_mouse_source,
    key=[(EMPTY, "human cells"), ("#b5b5b5", "mouse,\napoptotic\n(no value)")])

P["b"] = umap_panel("b", "Transferred onto human cells", draw_human("pt_ranknorm"),
    key=[(EMPTY, "mouse cells"), ("#8c8c8c", f"human, no\nmouse\nneighbor\n({int(hs.pt_ranknorm.isna().sum())} cells)")])

P["c"] = umap_panel("c", "Human pseudotime", draw_human("pt_human"),
    key=[(EMPTY, "mouse cells")], cbar="pseudotime")


# ---- panel d: share of cross-species edge weight ----
fig, W, H = pfig("d")
ax = fig.add_axes([18 / W, 11 / H, 1 - 21 / W, 1 - 28 / H])
rows = [
    ("human_to_mouse", "endometrial_fibroblast", "fibroblast"),
    ("human_to_mouse", "decidual_cell", "decidual"),
    ("mouse_to_human", "endometrial_fibroblast", "fibroblast"),
    ("mouse_to_human", "decidual_cell", "decidual"),
]
ypos = [3.6, 2.6, 1.0, 0.0]
for y, (direction, src, lab) in zip(ypos, rows):
    fshare = ES.loc[(direction, src), "endometrial_fibroblast"]
    dshare = ES.loc[(direction, src), "decidual_cell"]
    ax.barh(y, 100 * fshare, color=FIB, height=0.72, linewidth=0)
    ax.barh(y, 100 * dshare, left=100 * fshare, color=DEC, height=0.72, linewidth=0)
    match = fshare if src == "endometrial_fibroblast" else dshare
    xc = 100 * fshare / 2 if src == "endometrial_fibroblast" else 100 * (fshare + dshare / 2)
    ax.text(xc, y, f"{100 * match:.0f}%", ha="center", va="center", fontsize=8,
            color="white" if src == "endometrial_fibroblast" else "#222222")
ax.set_yticks(ypos)
ax.set_yticklabels([r[2] for r in rows])
ax.tick_params(axis="y", length=0)
ax.set_xlim(0, 100)
ax.set_xticks([0, 50, 100])
ax.set_ylim(-0.6, 4.75)
ax.spines["left"].set_visible(False)
ax.set_xlabel("share of edge weight (%)")


def direction_label(fig, ax, W, y, source, target):
    """'source -> target' with the arrow drawn as a line. Helvetica Neue has no arrow glyph, and an SVG
    that leaves the arrow to a fallback font shows a missing glyph wherever that font is not installed."""
    style = dict(xycoords=("figure fraction", "data"), ha="left", va="bottom", fontsize=8, color="#555555")
    first = ax.annotate(source, xy=(TX / W, y), **style)
    (_, y0_), (x1_, y1_) = fig.transFigure.inverted().transform(first.get_window_extent(fig.canvas.get_renderer()))
    y_mid = y0_ + 0.42 * (y1_ - y0_)  # middle of the x-height
    gap, length = 0.9 / W, 3.4 / W
    fig.add_artist(FancyArrowPatch((x1_ + gap, y_mid), (x1_ + gap + length, y_mid), transform=fig.transFigure,
                                   arrowstyle="-|>", mutation_scale=5, lw=0.7, color="#555555", shrinkA=0, shrinkB=0))
    ax.annotate(target, xy=(x1_ + 2 * gap + length, y), **style)


direction_label(fig, ax, W, 4.15, "human", "mouse")
direction_label(fig, ax, W, 1.55, "mouse", "human")
ax.legend(
    handles=[Patch(color=FIB, label="to fibroblast"), Patch(color=DEC, label="to decidual")],
    loc="upper left", bbox_to_anchor=(TX / W, 1 - 10.5 / H), bbox_transform=fig.transFigure, ncol=2, frameon=False, fontsize=8,
    handlelength=0.9, handletextpad=0.3, columnspacing=1.0, borderaxespad=0
)
title(fig, W, H, "Cross-species edges link\nmatching compartments")
P["d"] = fig


# ---- panel e: transferred vs human pseudotime, by compartment ----
fig, W, H = pfig("e")
D = Jn[Jn.pt_human.notna() & Jn.pt_ranknorm.notna()]
sets = [
    ("all human cells", D, None),
    ("decidual cells", D[D.cell_type == "decidual_cell"], DEC),
    ("fibroblasts", D[D.cell_type == "endometrial_fibroblast"], FIB),
]
x0, gap, tot = 18 / W, 3.5 / W, 1 - 22 / W
aw = (tot - 2 * gap) / 3
ay, ah = 12.0 / H, 1 - 27.0 / H
axs = [fig.add_axes([x0 + i * (aw + gap), ay, aw, ah]) for i in range(3)]
RHO = {}
for i, (ax, (ttl, g, col)) in enumerate(zip(axs, sets)):
    gg = g.sample(frac=1, random_state=4) if col is None else g
    ax.scatter(gg.pt_human, gg.pt_ranknorm, s=0.25,
               c=(gg.cell_type.map(CTC).values if col is None else col),
               alpha=0.25, linewidths=0, rasterized=True)
    if col is None:  # median-per-decile trend only on the all-cells panel; subset panels show points only
        decile = pd.qcut(g.pt_human.rank(method="first"), 10, labels=False)
        med = g.groupby(decile)[["pt_human", "pt_ranknorm"]].median()
        ax.plot(med.pt_human, med.pt_ranknorm, color="black", lw=0.9, zorder=4)
        ax.scatter(med.pt_human, med.pt_ranknorm, marker="D", s=9, color="black",
                   edgecolor="white", linewidths=0.4, zorder=5)
        ax.legend(
            handles=[Line2D([], [], color="black", lw=0.9, marker="D", ms=2.8, mec="white", mew=0.4,
                            label="median per decile")],
            loc="lower right", frameon=False, fontsize=7, handlelength=1.6, borderaxespad=0.2, handletextpad=0.4,
        )
    RHO[ttl] = spearmanr(g.pt_human, g.pt_ranknorm).statistic
    ax.text(0.04, 0.96, f"ρ = {RHO[ttl]:.2f}\nn = {len(g) / 1e3:.1f}k",
            transform=ax.transAxes, va="top", ha="left", fontsize=8)
    ax.set_title(ttl, loc="left", fontsize=8, pad=2)
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
    ax.set_xticks([0, 0.5, 1], labels=["0", "0.5", "1"]); ax.set_yticks([0, 0.5, 1], labels=["0", "0.5", "1"])
    if i:
        ax.set_yticklabels([])
axs[0].set_ylabel("transferred from mouse")
axs[1].set_xlabel("human pseudotime (rank-normalized)")
title(fig, W, H, "Transferred pseudotime preserves human ordering\nwithin decidual cells, weakly within fibroblasts")
P["e"] = fig


# ---- panel f: transferred pseudotime in the human sections, observed and shuffled ----
PHASES = {  # crop (um) and counter-clockwise rotation of each phase's section, as in Fig. 6E
    "menstrual": {"label": "menstrual", "limits": None, "rot": 90},
    "secretory": {"label": "secretory\n(decidualized)", "limits": ((400, 5500), (750, 5250)), "rot": -90},
}
MAPS = {
    "observed": "observed",
    "uniform_example": "shuffled across\nall mouse cells",
    "within_compartment_example": "shuffled within\ncompartment",
}
UM_PER_MM = 210        # one scale for every map, so the two phases are shown at the same size ratio
SECTION_PAD_UM = 60    # keeps edge dots from being clipped
MAP_DOT = 0.45         # pt^2; each section's three 10x lanes are drawn together
NO_VALUE_SPATIAL = "#b0b0b0"


def rotate(x, y, deg, cx, cy):
    t = np.deg2rad(deg)
    return (x - cx) * np.cos(t) - (y - cy) * np.sin(t) + cx, (x - cx) * np.sin(t) + (y - cy) * np.cos(t) + cy


def section_coords(phase, geom):
    """Cell ids and rotated coordinates (um) of the nuclei drawn for one phase."""
    cells = SPATIAL[(SPATIAL.timepoint == phase) & SPATIAL.x_um_qc.notna() & ~SPATIAL.index.isin(HS_FLAGGED)]
    x, y = cells.x_um_qc.to_numpy(), cells.y_um_qc.to_numpy()
    if geom["limits"] is not None:
        (x0_, x1_), (y0_, y1_) = geom["limits"]
        keep = (x >= x0_) & (x <= x1_) & (y >= y0_) & (y <= y1_)
        cells, x, y = cells[keep], x[keep], y[keep]
        cx, cy = (x0_ + x1_) / 2, (y0_ + y1_) / 2
    else:
        cx, cy = (x.min() + x.max()) / 2, (y.min() + y.max()) / 2
    xr, yr = rotate(x, y, geom["rot"], cx, cy)
    return cells.index, xr, yr


fig, W, H = pfig("f")
sections = {phase: section_coords(phase, geom) for phase, geom in PHASES.items()}
extent = {phase: ((xr.min() - SECTION_PAD_UM, xr.max() + SECTION_PAD_UM), (yr.min() - SECTION_PAD_UM, yr.max() + SECTION_PAD_UM))
          for phase, (_, xr, yr) in sections.items()}
size_mm = {phase: ((xs[1] - xs[0]) / UM_PER_MM, (ys[1] - ys[0]) / UM_PER_MM) for phase, (xs, ys) in extent.items()}
col_w = max(w for w, _ in size_mm.values())
label_w, map_gap = 9.0, 2.0
grid_w = 3 * col_w + 2 * map_gap
x_start = label_w + (W - label_w - grid_w) / 2       # centre the map grid in the space right of the labels
maps_top = H - 17.0                                   # below the two-line title and the two-line headers
DRAWN = {}
y_top = maps_top
for row, (phase, geom) in enumerate(PHASES.items()):
    ids, xr, yr = sections[phase]
    (xs, ys), (w_mm, h_mm) = extent[phase], size_mm[phase]
    y_bottom = y_top - h_mm
    DRAWN[phase] = {"nuclei": len(ids), "no_value": int(EX.loc[ids, "observed"].isna().sum())}
    for col, (column, header) in enumerate(MAPS.items()):
        cell_left = x_start + col * (col_w + map_gap)
        ax = fig.add_axes([(cell_left + (col_w - w_mm) / 2) / W, y_bottom / H, w_mm / W, h_mm / H])
        ax.set_xlim(*xs); ax.set_ylim(*ys)
        ax.axis("off")
        values = EX.loc[ids, column].to_numpy()
        missing = np.isnan(values)
        # Nuclei without a valued mouse neighbor go underneath, in grey; they are not imputed.
        ax.scatter(xr[missing], yr[missing], s=MAP_DOT, color=NO_VALUE_SPATIAL, linewidths=0, rasterized=True)
        ax.scatter(xr[~missing], yr[~missing], s=MAP_DOT, c=CM(NORM(values[~missing])), linewidths=0, rasterized=True)
        if row == 0:
            fig.text((cell_left + col_w / 2) / W, (maps_top + 1.0) / H, header, ha="center", va="bottom", fontsize=8)
    fig.text((x_start - 3.0) / W, (y_bottom + h_mm / 2) / H, geom["label"], rotation=90, ha="center", va="center", fontsize=8)
    y_top = y_bottom - map_gap

# Shared colorbar under the middle column, scale bar under the right one.
cax = fig.add_axes([(x_start + col_w + map_gap - 6.0) / W, 5.2 / H, (col_w + 12.0) / W, 1.6 / H])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CM), cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
cb.ax.set_xticklabels(["0 early", "0.5", "1 late"])
cb.outline.set_visible(False)
cb.ax.tick_params(length=1.5, pad=1.5)
cax.set_title("transferred pseudotime", fontsize=8, loc="center", pad=3)
bar_mm = 1000 / UM_PER_MM
bar_right = x_start + grid_w - 2.0
fig.add_artist(Line2D([(bar_right - bar_mm) / W, bar_right / W], [6.0 / H, 6.0 / H], color="black", lw=1.0))
fig.text((bar_right - bar_mm / 2) / W, 7.0 / H, "1 mm", ha="center", va="bottom", fontsize=7)
title(fig, W, H, "Shuffled mouse pseudotime keeps at most\nthe fibroblast/decidual split")
P["f"] = fig


# ---- panel g: observed statistics against the shuffles ----
STAT_ROWS = {
    "rho": ("Spearman ρ with\nhuman pseudotime", (-0.3, 1.0), [0, 0.5, 1.0]),
    "moran": ("Moran's I (spatial\nautocorrelation)", (-0.05, 0.8), [0, 0.4, 0.8]),
}
STAT_COLUMNS = {"all": "all cells", "decidual": "decidual", "fibroblast": "fibroblasts"}
STAT_KEYS = {
    ("rho", "all"): "rho_all", ("rho", "decidual"): "rho_decidual", ("rho", "fibroblast"): "rho_fibroblast",
    ("moran", "all"): "moran_all", ("moran", "decidual"): "moran_within_decidual",
    ("moran", "fibroblast"): "moran_within_fibroblast",
}
# The two shuffle designs differ by lightness only; they never share a mark.
NULL_COLORS = {"uniform": "#b8b8b8", "within_compartment": "#6b6b6b"}
NULL_LABELS = {
    "uniform": "1,000 shuffles across all mouse cells",
    "within_compartment": "1,000 shuffles within mouse compartment",
}
observed = NULLS[NULLS.design == "observed"].iloc[0]

fig, W, H = pfig("g")
left, right, col_gap, ah_mm, row_gap = 17.0, 1.0, 2.0, 22.0, 4.5
aw_mm = (W - left - right - 2 * col_gap) / 3
rng = np.random.default_rng(0)  # jitter only
for i, (row_key, (row_label, ylim, yticks)) in enumerate(STAT_ROWS.items()):
    top_mm = H - 13.5 - i * (ah_mm + row_gap)
    for j, (col_key, col_label) in enumerate(STAT_COLUMNS.items()):
        ax = fig.add_axes([(left + j * (aw_mm + col_gap)) / W, (top_mm - ah_mm) / H, aw_mm / W, ah_mm / H])
        stat = STAT_KEYS[(row_key, col_key)]
        for pos, design in enumerate(NULL_COLORS):
            vals = NULLS.loc[NULLS.design == design, stat].to_numpy()
            ax.scatter(pos + rng.uniform(-0.22, 0.22, vals.size), vals, s=0.8, color=NULL_COLORS[design],
                       linewidths=0, rasterized=True)
            ax.plot([pos - 0.3, pos + 0.3], [np.median(vals)] * 2, color="black", lw=0.6)
        ax.axhline(observed[stat], color="black", lw=1.1)
        ax.text(-0.45, observed[stat] + 0.015 * (ylim[1] - ylim[0]), f"{observed[stat]:.2f}",
                ha="left", va="bottom", fontsize=7)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylim(*ylim)
        ax.set_xticks([])  # design identity is carried by shade, keyed in the legend
        ax.set_yticks(yticks)
        if i == 0:
            ax.set_title(col_label, fontsize=8, pad=2)
        if j == 0:
            ax.set_ylabel(row_label, fontsize=8)
        else:
            ax.set_yticklabels([])
fig.legend(
    handles=[
        Line2D([], [], color="black", lw=1.1, label="observed transfer"),
        *[Line2D([], [], ls="none", marker="o", ms=3, color=c_, label=NULL_LABELS[d_]) for d_, c_ in NULL_COLORS.items()],
        Line2D([], [], color="black", lw=0.6, label="median of shuffles"),
    ],
    loc="lower left", bbox_to_anchor=(3.0 / W, 0.6 / H), frameon=False, fontsize=8,
    handlelength=1.4, labelspacing=0.3, borderaxespad=0,
)
title(fig, W, H, "The gradient within decidual cells\ndepends on mouse pseudotime order")
P["g"] = fig


# ---- save panels, then compose PNG and SVG ----
paths = {}
with mpl.rc_context({"savefig.bbox": "standard"}):
    for L, f in P.items():
        fn = os.path.join(OUT, "panels", f"sfig_samap_panel_{L}.png")
        f.savefig(fn, dpi=300, transparent=True)
        with mpl.rc_context({"svg.fonttype": "none"}):
            f.savefig(fn.replace(".png", ".svg"), dpi=600, transparent=True)
        assert Image.open(fn).size == panel_px(outline, L), (L, Image.open(fn).size)
        paths[L] = fn


def compose_figure(outline, panel_paths, out_path, dpi=300, gutter_mm=4, letter_font="DejaVuSans-Bold.ttf", letter_pt=11):
    W, ncol, colw, rowh, row_y, g = grid_geom(outline, dpi, gutter_mm)
    H = row_y[-1] + rowh[-1]
    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        ft = ImageFont.truetype(letter_font, int(letter_pt / 72 * dpi), index=LETTER_INDEX)
    except Exception:
        ft = ImageFont.load_default()
    for p in outline["panels"]:
        L = p["letter"]
        w, h = panel_px(outline, L, dpi, gutter_mm)
        x, y = panel_xy(outline, L, dpi, gutter_mm)
        im = Image.open(panel_paths[L]).convert("RGBA")
        if im.size != (w, h):
            im = im.resize((w, h))
        canvas.paste(im, (x, y), im)
        draw.text((x + int(1.5 / 25.4 * dpi), y + int(1 / 25.4 * dpi)), L, fill="black", font=ft)
    canvas.save(out_path)
    return out_path, (W, H)


def compose_svg(outline, panel_dir, out_svg, letter_pt=11):
    """Nest each panel SVG at its grid position (ids prefixed by the panel letter so panels cannot collide)."""
    pt = 72 / 300  # grid is computed in 300-dpi pixels
    W, ncol, colw, rowh, row_y, g = grid_geom(outline)
    page_w, page_h = W * pt, (row_y[-1] + rowh[-1]) * pt
    parts = [f'<?xml version="1.0" encoding="utf-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" '
             f'xmlns:xlink="http://www.w3.org/1999/xlink" width="{page_w:.2f}pt" height="{page_h:.2f}pt" '
             f'viewBox="0 0 {page_w:.2f} {page_h:.2f}">\n<rect width="100%" height="100%" fill="white"/>']
    for p in outline["panels"]:
        L = p["letter"]
        body = open(os.path.join(panel_dir, f"sfig_samap_panel_{L}.svg")).read()
        body = re.sub(r"<\?xml[^>]*\?>|<!DOCTYPE[^>]*>", "", body).strip()
        body = re.sub(r'id="([^"]+)"', rf'id="{L}_\1"', body)
        body = re.sub(r'(xlink:href|href)="#([^"]+)"', rf'\1="#{L}_\2"', body)
        body = re.sub(r"url\(#([^)]+)\)", rf"url(#{L}_\1)", body)
        width = float(re.search(r'<svg[^>]*\swidth="([\d.]+)pt"', body).group(1))
        height = float(re.search(r'<svg[^>]*\sheight="([\d.]+)pt"', body).group(1))
        x, y = (v * pt for v in panel_xy(outline, L))
        body = re.sub(r"<svg\b[^>]*>", (f'<svg x="{x:.3f}" y="{y:.3f}" width="{width:.3f}" height="{height:.3f}" '
                                        f'viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" '
                                        f'xmlns:xlink="http://www.w3.org/1999/xlink">'), body, count=1)
        parts.append(f'<g id="panel_{L}">{body}</g>')
        # Same place as the PNG letters: 1.5 mm in, baseline one cap height below 1 mm.
        parts.append(f'<text x="{x + 4.25:.2f}" y="{y + 11.41:.2f}" font-family="Helvetica Neue, Helvetica, Arial, sans-serif" '
                     f'font-weight="bold" font-size="{letter_pt}" fill="black">{L}</text>')
    parts.append("</svg>\n")
    open(out_svg, "w").write("\n".join(parts))


out_path, (Wf, Hf) = compose_figure(outline, paths, FIG_OUT, letter_font=LETTER_FONT)
with mpl.rc_context({"savefig.bbox": "standard"}):
    paths600 = {}
    for L, f in P.items():
        fn6 = os.path.join(OUT, "panels", f"sfig_samap_panel_{L}_600dpi.png")
        f.savefig(fn6, dpi=600, transparent=True)
        paths600[L] = fn6
compose_figure(outline, paths600, FIG_OUT.replace(".png", "_600dpi.png"), dpi=600, letter_font=LETTER_FONT)
compose_svg(outline, os.path.join(OUT, "panels"), FIG_OUT.replace(".png", ".svg"))
print(out_path, Wf, Hf, f"= {Wf / MMpx:.0f} x {Hf / MMpx:.0f} mm")


# ---- provenance sidecar (figures-suite convention) ----
def _sha(p):
    h_ = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h_.update(b)
    return h_.hexdigest()


prov = {
    "figure": os.path.abspath(FIG_OUT), "size_px": [Wf, Hf], "size_mm": [round(Wf / MMpx, 1), round(Hf / MMpx, 1)],
    "script": os.path.abspath(__file__), "script_sha256": _sha(os.path.abspath(__file__)),
    "created": datetime.datetime.now().isoformat(timespec="seconds"), "host": platform.node(),
    "font": FONT, "letter_font": [LETTER_FONT, LETTER_INDEX],
    "versions": {"python": platform.python_version(), "matplotlib": mpl.__version__, "pandas": pd.__version__, "numpy": np.__version__},
    "sources": [{"path": os.path.abspath(p), "bytes": os.path.getsize(p),
                 "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(p)).isoformat(timespec="seconds"),
                 "sha256": _sha(p)} for p in SOURCES],
    "spearman_rho": {k: round(float(v), 4) for k, v in RHO.items()},
    "embedding_outliers_not_drawn": {"mouse": int((UM.species[outlier] == "mm").sum()), "human": int((UM.species[outlier] == "hs").sum())},
    "spatial_nuclei_drawn": DRAWN,
    "permutation_observed": {k: round(float(observed[k]), 4) for k in STAT_KEYS.values()},
    "permutation_replicates": {d_: int((NULLS.design == d_).sum()) for d_ in NULL_COLORS},
    "note": "Pseudotime: pt_harmonized_rank (2026_10_22_flagged_harmonized release) and its SAMap transfer of record.",
}
json.dump(prov, open(os.path.join(OUT, "sfig_samap.provenance.json"), "w"), indent=1)
print("wrote", FIG_OUT, "+ _600dpi.png + .svg + panels/ + sfig_samap.provenance.json")
