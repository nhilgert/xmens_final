"""Decidual-zone AUCell module scores: spatial panels of one mouse section and subtype x module heatmaps.

Pipeline version (rule zone_figures, workflow/rules/16_zone_modules.smk; driver render_zone_figures.py), v16,
2026-10-24. Derived from figures/zone_modules/fig_zone_modules.py v14 (md5 e6b0ccef8e858d3ccbe6607e9289e621); the
plotting functions are unchanged.
    * v16: the ONLY data input is --plot-table, the per-nucleus CSV written by make_zone_plot_table.py (rule
      zone_plot_table): barcode, sample, timepoint, cell_type, cell_subtype, x_um_qc, y_um_qc, embedding_outlier and
      the 11 auc_<MODULE> scores. No h5ad is read, so every figure can be redrawn from that CSV and its sidecar
      <table>.provenance.json (needs numpy, pandas and matplotlib only).
    * nuclei with embedding_outlier == True are hidden from the SPATIAL panels and kept in heatmaps (rule of
      2026-10-22: embedding outliers are omitted from embedding and spatial plots and kept in every summary);
      --show-flagged draws them too.
    * the scores may cover only some timepoints (the user's 3.5 dpi decision, 2026-10-24). An all-timepoint heatmap is
      then refused (pass --heatmap-timepoint), as are spatial panels of unscored timepoints.
    * score metadata (window, seed, normalization, scored timepoints, gene sets) comes from the table's sidecar and is
      copied into every figure sidecar. The objects' stored auc_ columns (window 0.0131) are never read.

Display rotation of the spatial panels (2026-10-30): the timepoint's angle in --section-rotation
(resources/figures/mouse_section_rotation.yaml, the table every mouse spatial panel reads through
workflow/scripts/figures/section_rotation.py; rule zone_figures passes it to every job). --rotate DEG remains for
renders by hand; the two cannot be combined. The timepoint's whole placed array (all lanes, embedding outliers and
decAp included) turns as one body about the centre of its bounding box before any display filter.

Outputs per figure: <name>{,_dark}.{png,svg} + <name>{,_dark}.provenance.json, names as in v14:
    pdz_sdz_spatial_<tp>[_noapop][_rot<deg>]                 PDZ and SDZ
    zone_modules_spatial_<tp>[_noapop][_rot<deg>][_<n>zones] all 11 modules (PDZ featured), or the --modules subset
    zone_modules_heatmap_subtypes[_with_apoptotic][_<tp>][_<n>zones][_fibpooled][_p99]
    plus _<cmap> for a non-magma colour map and --suffix. Heatmap SVG is fully vector; in spatial SVG each section's
    point cloud is one embedded bitmap at --dpi.

DATA: --plot-table (see above). Its obs columns come from the mouse release object: x_um_qc, y_um_qc (microns after
    spatial QC; NaN = not placed, left out of spatial panels but counted in heatmaps), sample, timepoint (lanes of one
    timepoint are one section -> pooled), cell_type, cell_subtype. Local renders (e.g. in Helvetica Neue on a Mac) use
    the same table.

MODULES: PDZ / TDZ / SDZ = primary / transitional / secondary decidual zone (inner -> outer), VSZ (vascular sinus
    zone), MD (mesometrial decidua; top markers are uterine NK-cell genes, NK cells Li et al. describe as
    non-cytotoxic), MMy, C (conceptus), FMI (fetal-maternal interface; endothelial-rich in Li et al.'s deconvolution,
    which maps myeloid cells to C), Glan, udStr, AMMy: the 11 uterine microenvironments of Li et al. (2022) Biol
    Reprod 107:529, Table 1 (top 20 markers each; resources/signatures/aucell_modules_mouse.yaml).

CONVENTIONS (settled 2026-10-18; don't change silently):
    * AUCell window: pySCENIC default auc_threshold 0.05, seed 0 (since 2026-10-22); the sidecar names the window.
    * heatmap (--heatmap-scale p99, files *_p99*): colour = subtype MEAN / 99th pct of per-cell scores of that module
      (cells of the plotted subtypes, all timepoints), i.e. same colour = same score level as the spatial panels.
    * --pool-fib: the five eFib rows become one 'endometrial fibroblasts (pooled)' row = mean over all fibroblast cells (cell-weighted);
      the p99 ceiling is unchanged (same cells).
    * heatmap (default colmax): colour = subtype MEAN score from a true 0 to the module's largest subtype mean (per column);
      the largest mean is printed under each column so absolute levels stay readable. The exploratory
      notebook used min-max per column instead (lowest subtype always black). Subtypes are the vocabulary v4
      codes (eFib1-eFibP, dec1-dec3; a plot table from an older release is read through the vocabulary's
      legacy_aliases) and include decAp only in the *_with_apoptotic variant: that subtype (not in the exploratory
      heatmap; excluded from the pseudotime of record) has the highest mean in PDZ/TDZ/SDZ/VSZ/FMI/udStr, so
      including it rescales those columns and hides the dec3 PDZ/TDZ signal the exploratory figure showed.
    * colour = AUCell score from a TRUE 0 to the module's 99th percentile within the plotted cells, per
      module (AUCell at this threshold never approaches 1, so a fixed 0-1 axis leaves the map unused; a
      true min-max lets a few extreme cells hide the structure). Stated on the colour bar.
    * lanes of one section are pooled ("merge by timepoint"); cells without coordinates are left out
    * bright cells are drawn last (sorted by value)
    * never overwrite an existing figure without --force
    The original exploratory panel (xmens_final/analysis/figures/loose/modules_spatial_3dpi_aucell_dark.png,
    from analysis/notebooks/module_scoring_aucell.ipynb) used the 2nd-98th percentile per module, re-scored
    the 3.5 dpi cells on their own, and applied a second spatial-outlier filter; this script does none of those.
"""
import os, sys, json, argparse, datetime
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import hashlib
from pathlib import Path
import yaml

# the shared vocabulary helpers live in workflow/scripts/atlas/, the section rotation reader in workflow/scripts/figures/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "atlas"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
from apply_annotation import current_label, current_labels  # noqa: E402
import section_rotation  # noqa: E402

ORDER = ["PDZ", "TDZ", "SDZ", "VSZ", "MD", "MMy", "C", "FMI", "Glan", "udStr", "AMMy"]
TWO = [("PDZ", "primary decidual zone (PDZ)"), ("SDZ", "secondary decidual zone (SDZ)")]
SUBTYPES = ["eFib1", "eFib2", "eFib2i", "eFib3", "eFibP",
            "dec1", "dec2", "dec3"]   # rows of the exploratory heatmap (vocabulary v4 codes)
SUBTYPES_APOP = SUBTYPES + ["decAp"]   # variant; this row has the highest mean in 6 of 11 modules
N_FIB = 5   # first N_FIB rows are endometrial_fibroblast subtypes; a rule separates them from decidual
CMAP, Q = "magma", 0.99
META_GREY = "#6f6f6f"


def style():
    avail = {f.name for f in font_manager.fontManager.ttflist}
    fam = next((f for f in ("Helvetica Neue", "Helvetica", "Nimbus Sans", "Nimbus Sans L", "DejaVu Sans") if f in avail), "DejaVu Sans")
    plt.rcParams.update({"font.family": [fam, "DejaVu Sans"], "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "savefig.dpi": 300,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    return fam


def _md5(path):
    h = hashlib.md5(); f = open(path, "rb")
    for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def _src(path, **extra):
    return {"path": os.path.abspath(path), "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(path)).isoformat(timespec="seconds"),
            "bytes": os.path.getsize(path), **extra}


def load(timepoint, plot_table, vocab):
    """The per-nucleus plot table and its sidecar. Returns the spatial set g (one timepoint, placed nuclei), D (all
    rows), score metadata, sources, the scored timepoints and whether the scores cover only some timepoints."""
    auc = [f"auc_{m}" for m in ORDER]
    base = ["sample", "timepoint", "cell_type", "cell_subtype", "x_um_qc", "y_um_qc", "embedding_outlier"]
    D = pd.read_csv(plot_table, index_col=0); D.index = D.index.astype(str)
    miss = [c for c in base + auc if c not in D.columns]; assert not miss, f"missing columns in {plot_table}: {miss}"
    assert D.index.is_unique, f"duplicate barcodes in {plot_table}"
    for c in ("sample", "timepoint", "cell_type", "cell_subtype"): D[c] = D[c].astype(str)
    D["cell_subtype"] = current_labels(D["cell_subtype"], vocab)   # tables cut from older releases carry superseded names
    D["timepoint"] = current_labels(D["timepoint"], vocab, "timepoint")
    flag = D["embedding_outlier"].map({True: True, False: False, "True": True, "False": False})
    assert flag.notna().all(), "embedding_outlier must be True/False"
    D["embedding_outlier"] = flag.astype(bool)
    side = plot_table[:-len(".csv")] + ".provenance.json" if plot_table.endswith(".csv") else plot_table + ".provenance.json"
    P = json.load(open(side)).get("scoring", {}) if os.path.exists(side) else {}
    meta = {k: str(P.get(k, "?")) for k in ("auc_threshold", "rank_cutoff", "seed", "normalization", "timepoints")}
    meta["genesets"] = P.get("genesets", {}); meta["scores_file"] = os.path.abspath(plot_table)
    srcs = [_src(plot_table, md5=_md5(plot_table), role="plot table (obs + embedding_outlier + auc_<MODULE>)")]
    if os.path.exists(side): srcs.append(_src(side, md5=_md5(side), role="plot table provenance"))
    tps = sorted(D["timepoint"].unique())
    scored = P.get("timepoints")
    if isinstance(scored, list):
        scored = [current_label(t, vocab, "timepoint") for t in scored]; meta["timepoints"] = str(scored)
    if isinstance(scored, list): assert sorted(scored) == tps, f"table timepoints {tps} differ from the scored {scored}"
    partial = scored != "all"   # unknown (no sidecar) counts as partial: all-timepoint heatmaps are then refused
    assert timepoint in tps, f"timepoint {timepoint!r} is not scored; the table covers {tps}"
    g = D[(D["timepoint"] == timepoint) & D[["x_um_qc", "y_um_qc"]].notna().all(1)].copy()
    assert len(g) > 0, f"no cells with coordinates at timepoint {timepoint!r}"
    return g, D, meta, srcs, tps, partial


def rotate_section(g, timepoint, deg):
    """The placed cells of one timepoint turned by deg (counter-clockwise, y up) as one rigid body.
    The lanes of one timepoint are one section in one coordinate frame (3.5 dpi lane centres differ by < 30 um),
    so the pooled cells rotate about a single shared centre (section_rotation: the midpoint of their bounding box).
    Display only: the returned copy carries the turned coordinates, the plot table is not modified."""
    if not deg % 360: return g
    x, y = section_rotation.rotate_by_timepoint(g["x_um_qc"], g["y_um_qc"], [timepoint] * len(g), {timepoint: deg})
    return g.assign(x_um_qc=x, y_um_qc=y)


def rc(dark):
    bg, fg = ("black", "white") if dark else ("white", "black")
    return {"figure.facecolor": bg, "axes.facecolor": bg, "text.color": fg, "axes.labelcolor": fg, "xtick.color": fg,
            "ytick.color": fg, "savefig.facecolor": bg}, fg


def panel(ax, g, mod, s):
    v = g[f"auc_{mod}"].values.astype(float); hi = float(np.quantile(v, Q))
    vn = np.clip(v / hi, 0, 1) if hi > 0 else np.zeros_like(v); o = np.argsort(vn)
    ax.scatter(g.x_um_qc.values[o], g.y_um_qc.values[o], c=vn[o], cmap=CMAP, vmin=0, vmax=1, s=s, linewidths=0, rasterized=True)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.margins(0.04)
    return hi


def cbar(fig, rect, fg):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0, 1), cmap=CMAP), cax=cax, ticks=[0, 1])
    cb.ax.set_yticklabels(["0", "99th pct"]); cb.outline.set_visible(False); cb.ax.tick_params(color=fg, labelcolor=fg, length=2)
    cax.set_title("AUCell score\n(0 \u2013 99th pct,\nper module)", loc="left", fontsize=8, color=fg, pad=6)


def note(g, tp):
    lanes = sorted(g["sample"].unique()); tpl = tp.replace("_", " ")
    return f"{tpl}, n = {len(g):,} cells" + (f" ({len(lanes)} lanes of one section pooled)" if len(lanes) > 1 else "")


def fig_two(g, tp, dark):
    rcd, fg = rc(dark)
    with plt.rc_context(rcd):
        fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.6)); ceil = {}
        for ax, (m, lab) in zip(axes, TWO):
            ceil[m] = panel(ax, g, m, 1.8); ax.set_title(lab, loc="left", color=fg)
        fig.subplots_adjust(left=0.02, right=0.84, top=0.90, bottom=0.07, wspace=0.04); cbar(fig, [0.865, 0.28, 0.018, 0.40], fg)
    return fig, {"ceiling_99th_pct": {m: round(v, 4) for m, v in ceil.items()}}


def fig_six(g, tp, dark, modules):
    """Equal-size panels for a module subset (e.g. the heatmap's six), 2 rows x ceil(n/2) columns."""
    rcd, fg = rc(dark); nc = int(np.ceil(len(modules) / 2))
    with plt.rc_context(rcd):
        fig = plt.figure(figsize=(2.55 * nc + 1.2, 5.6))
        gs = matplotlib.gridspec.GridSpec(2, nc, figure=fig, left=0.01, right=1 - 1.2 / (2.55 * nc + 1.2), top=0.93, bottom=0.06, wspace=0.05, hspace=0.14)
        ceil = {}
        for i, m in enumerate(modules):
            ax = fig.add_subplot(gs[divmod(i, nc)]); ceil[m] = panel(ax, g, m, 1.4); ax.set_title(m, loc="left", color=fg)
        cbar(fig, [1 - 0.95 / (2.55 * nc + 1.2), 0.32, 0.012, 0.36], fg)
    return fig, {"modules": list(modules), "ceiling_99th_pct": {m: round(v, 4) for m, v in ceil.items()}}


def fig_all(g, tp, dark):
    rcd, fg = rc(dark)
    with plt.rc_context(rcd):
        fig = plt.figure(figsize=(11.0, 4.3))
        gs = matplotlib.gridspec.GridSpec(2, 6, figure=fig, left=0.01, right=0.90, top=0.92, bottom=0.06, wspace=0.04, hspace=0.18)
        ceil = {}
        ax = fig.add_subplot(gs[:, 0]); ceil["PDZ"] = panel(ax, g, "PDZ", 2.4); ax.set_title("PDZ", loc="left", color=fg)
        for i, m in enumerate(ORDER[1:]):
            r, c = divmod(i, 5); ax = fig.add_subplot(gs[r, c + 1]); ceil[m] = panel(ax, g, m, 1.2); ax.set_title(m, loc="left", color=fg, fontsize=8)
        cbar(fig, [0.915, 0.32, 0.012, 0.36], fg)
    return fig, {"ceiling_99th_pct": {m: round(v, 4) for m, v in ceil.items()}}


FIB_POOLED = "endometrial fibroblasts (pooled)"

def fig_heatmap(D, dark, subtypes=SUBTYPES, scale="colmax", pool_fib=False, modules=None, timepoint=None):
    if timepoint:   # one timepoint only (lanes of that section pooled); the p99 ceiling is then computed on those cells too
        D = D[D["timepoint"].astype(str) == timepoint]; assert len(D), timepoint
    ORDER_ = list(modules) if modules else ORDER
    """Mean AUCell per fib/dec subtype (rows) x module (columns); colour 0 -> column max."""
    rcd, fg = rc(dark); grey = "#bbbbbb" if dark else META_GREY
    sub = D[D["cell_subtype"].isin(subtypes)]
    missing = [s_ for s_ in subtypes if s_ not in set(sub["cell_subtype"])]
    assert not missing, f"subtypes absent from object: {missing}"
    n_fib_rows = N_FIB
    if pool_fib:   # one row = mean over ALL fibroblast cells (cell-weighted), not the mean of the five subtype means
        fibs = subtypes[:N_FIB]; assert all(f_.startswith("eFib") for f_ in fibs), fibs
        sub = sub.assign(cell_subtype=np.where(sub["cell_subtype"].isin(fibs), FIB_POOLED, sub["cell_subtype"].astype(str)))
        subtypes = [FIB_POOLED] + list(subtypes[N_FIB:]); n_fib_rows = 1
    M = sub.groupby("cell_subtype")[[f"auc_{m}" for m in ORDER_]].mean().reindex(subtypes)
    n = sub["cell_subtype"].value_counts().reindex(subtypes)
    if scale == "p99":   # 0 -> 99th pct of PER-CELL scores (cells of these subtypes, all timepoints): same colour = same score level as the spatial panels
        cmax = sub[[f"auc_{m}" for m in ORDER_]].quantile(Q)
    else:                # 0 -> largest subtype mean
        cmax = M.max(axis=0)
    S = (M / cmax.values).clip(0, 1).values
    nr, nc = S.shape; cw = 0.40                       # cell size, inches
    L, R, T, B = 1.75, 2.25, 0.45, 0.75               # margins, inches
    W, H = L + nc * cw + R, T + nr * cw + B
    with plt.rc_context(rcd):
        fig = plt.figure(figsize=(W, H))
        ax = fig.add_axes([L / W, B / H, nc * cw / W, nr * cw / H])
        # vector cells (one path per cell) so the SVG grid is editable; edgecolors="face" avoids hairline seams
        ax.pcolormesh(np.arange(nc + 1) - 0.5, np.arange(nr + 1) - 0.5, S, cmap=CMAP, vmin=0, vmax=1,
                      edgecolors="face", linewidth=0.2, rasterized=False)
        ax.set_xlim(-0.5, nc - 0.5); ax.set_ylim(nr - 0.5, -0.5)
        ax.set_xticks(range(nc)); ax.set_xticklabels(ORDER_); ax.xaxis.tick_top()
        ax.set_yticks(range(nr)); ax.set_yticklabels(subtypes)
        ax.tick_params(length=0, pad=3, colors=fg)
        for sp in ax.spines.values(): sp.set_visible(False)
        if not pool_fib:   # fibroblast/decidual divider only when the five fibroblast subtypes are shown separately
            ax.axhline(n_fib_rows - 0.5, color=fg, lw=0.8)
        for i, k in enumerate(n.values):
            ax.text(nc - 0.5 + 0.15, i, f"n = {int(k):,}", ha="left", va="center", fontsize=7, color=grey)
        for j, v in enumerate(cmax.values):
            ax.text(j, nr - 0.5 + 0.12, f"{v:.2f}" if v >= 0.01 else f"{v:.3f}", ha="center", va="top", fontsize=6.5, color=grey, rotation=90)
        ax.text(-0.5 - 0.15, nr - 0.5 + 0.12, "99th pct\n(per cell)" if scale == "p99" else "column max\n(mean AUCell)", ha="right", va="top", fontsize=6.5, color=grey)
        fig.text(0.01, 0.015, f"mouse fibroblast/decidual cells, " + (f"{timepoint.replace('_', ' ')} only (both lanes of the section pooled)" if timepoint else "all timepoints pooled") + f", n = {int(n.sum()):,}", ha="left", va="bottom", fontsize=7, color=grey)
        x0 = (L + nc * cw + 0.95) / W
        cax = fig.add_axes([x0, (B + nr * cw * 0.25) / H, 0.12 / W, nr * cw * 0.5 / H])
        cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0, 1), cmap=CMAP), cax=cax, ticks=[0, 1])
        cb.ax.set_yticklabels(["0", "99th pct" if scale == "p99" else "max"]); cb.outline.set_visible(False); cb.ax.tick_params(color=fg, labelcolor=fg, length=2)
        cb.solids.set_rasterized(False); cb.solids.set_edgecolor("face")   # colour bar as vector too
        cax.set_title("mean AUCell\n(0 \u2013 99th pct\nof cells)" if scale == "p99" else "mean AUCell\n(0 \u2013 column max)", loc="left", fontsize=8, color=fg, pad=6)
    return fig, {"heatmap_scale": scale, ("ceiling_99th_pct_per_cell" if scale == "p99" else "column_max_mean_auc"): {m: round(float(v), 4) for m, v in zip(ORDER_, cmax.values)},
                 "n_per_subtype": {k: int(v) for k, v in n.items()},
                 "mean_auc": {s_: {m: round(float(M.loc[s_, f"auc_{m}"]), 4) for m in ORDER_} for s_ in subtypes}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="."); ap.add_argument("--force", action="store_true"); ap.add_argument("--timepoint", default="3.5_dpi")
    ap.add_argument("--plot-table", required=True); ap.add_argument("--suffix", default="")
    ap.add_argument("--dpi", type=int, default=300); ap.add_argument("--formats", default="png,svg")
    ap.add_argument("--cmap", default="magma"); ap.add_argument("--figures", default="spatial,heatmap")
    ap.add_argument("--heatmap-scale", default="colmax", choices=["colmax", "p99"])   # p99 adds _p99 to heatmap names
    ap.add_argument("--modules", default="")   # comma list subset of columns, e.g. PDZ,TDZ,SDZ,VSZ,MD,udStr; adds _<n>zones to heatmap names
    ap.add_argument("--heatmap-timepoint", default="")   # e.g. 3.5_dpi; heatmap from that timepoint only; adds _<tp> to heatmap names
    ap.add_argument("--section-rotation", default="")   # resources/figures/mouse_section_rotation.yaml: spatial panels turned by the timepoint's angle; adds _rot<deg>
    ap.add_argument("--rotate", type=float, default=0)   # by hand only: degrees CCW about the section centre; adds _rot<deg>
    ap.add_argument("--drop-apoptotic", action="store_true")   # spatial panels: drop decAp before scaling/plotting; adds _noapop
    ap.add_argument("--pool-fib", action="store_true")   # one fibroblast row; adds _fibpooled to heatmap names
    ap.add_argument("--vocab", default="resources/annotation_vocabulary.yaml")   # legacy_aliases for older plot tables
    ap.add_argument("--show-flagged", action="store_true")   # also draw embedding_outlier nuclei in spatial panels
    a = ap.parse_args()
    global CMAP; CMAP = a.cmap; ctag = "" if a.cmap == "magma" else f"_{a.cmap}"
    want = {f.strip() for f in a.figures.split(",")}
    fam = style(); os.makedirs(a.out, exist_ok=True)
    vocab = yaml.safe_load(open(a.vocab))
    a.timepoint = current_label(a.timepoint, vocab, "timepoint")   # older names (e.g. 3_dpi) are accepted
    if a.heatmap_timepoint:
        a.heatmap_timepoint = current_label(a.heatmap_timepoint, vocab, "timepoint")
    g, D, uns, srcs, scored_tps, partial = load(a.timepoint, a.plot_table, vocab)
    # Display rotation: from the shared table, or --rotate by hand. Applied to the whole placed array before the
    # display filters below, so hidden nuclei sit where they belong and do not depend on the angle.
    if a.section_rotation and a.rotate % 360:
        sys.exit("give --section-rotation or --rotate, not both")
    rot = section_rotation.load_rotation_table(a.section_rotation)[a.timepoint] if a.section_rotation else a.rotate
    g = rotate_section(g, a.timepoint, rot)
    rot_src = (f"{a.section_rotation} (sha256 {section_rotation.table_record(a.section_rotation)['sha256'][:12]})"
               if a.section_rotation else "--rotate")
    win = (f"AUCell window {uns['auc_threshold']} (top {uns['rank_cutoff']} genes, seed {uns['seed']}, {uns['normalization']}; "
           f"scored timepoints {uns['timepoints']}; table {os.path.basename(uns['scores_file'])})")
    n_hid, hid_note = 0, ""
    if not a.show_flagged:   # spatial panels only; heatmaps (D) keep these nuclei
        n_hid = int(g["embedding_outlier"].sum()); g = g[~g["embedding_outlier"]]
        hid_note = f"; {n_hid} embedding-outlier nuclei hidden (kept in heatmaps)"
    tag = a.timepoint.replace("_", "")  # 3.5_dpi -> 3.5dpi
    n_ap = 0
    if a.drop_apoptotic:
        n_ap = int((g["cell_subtype"] == "decAp").sum()); g = g[g["cell_subtype"] != "decAp"]; tag += "_noapop"
    if rot % 360:
        tag += f"_rot{rot:g}"
    sp_note = (f"AUCell module scores, {win}, timepoint {a.timepoint}, lanes pooled {sorted(g['sample'].unique())}, "
               f"{len(g)} cells with QC'd coordinates" + hid_note + (f" ({n_ap} decAp removed)" if a.drop_apoptotic else "") + (f"; display rotated {rot:g} deg CCW about the centre of the section's bounding box ({rot_src})" if rot % 360 else "") + f"; colour 0 -> per-module 99th percentile; font {fam}")
    mods = [m.strip() for m in a.modules.split(",") if m.strip()] or None
    assert not mods or set(mods) <= set(ORDER), mods
    htp = a.heatmap_timepoint or None
    if "heatmap" in want:
        assert htp is not None or not partial, ("the scores cover timepoints %s only; pass --heatmap-timepoint "
                                                          "(an 'all timepoints pooled' heatmap would be mislabelled)" % scored_tps)
        assert htp is None or htp in scored_tps, f"heatmap timepoint {htp!r} is not scored ({scored_tps})"
    ptag = (f"_{htp.replace('_', '')}" if htp else "") + (f"_{len(mods)}zones" if mods else "") + ("_fibpooled" if a.pool_fib else "") + ("_p99" if a.heatmap_scale == "p99" else "")
    Dh = D[D["timepoint"] == htp] if htp else D
    hm_note = lambda st: (f"mean AUCell per fibroblast/decidual subtype ({len(st)} subtypes), " + (f"timepoint {htp} only" if htp else "all timepoints pooled") + ", "
                          f"{int(Dh.cell_subtype.isin(st).sum())} cells; {win}; colour 0 -> " + ("99th pct of per-cell scores" if a.heatmap_scale == "p99" else "largest subtype mean") + f" per module; font {fam}")
    jobs = [("pdz_sdz_spatial_%s" % tag, lambda dark: fig_two(g, a.timepoint, dark), sp_note),
            ("zone_modules_spatial_%s" % tag, lambda dark: fig_all(g, a.timepoint, dark), sp_note),
            *([("zone_modules_spatial_%s_%dzones" % (tag, len(mods)), lambda dark: fig_six(g, a.timepoint, dark, mods), sp_note)] if mods else []),
            ("zone_modules_heatmap_subtypes" + ptag, lambda dark: fig_heatmap(D, dark, SUBTYPES, a.heatmap_scale, a.pool_fib, mods, htp), hm_note(SUBTYPES)),
            ("zone_modules_heatmap_subtypes_with_apoptotic" + ptag, lambda dark: fig_heatmap(D, dark, SUBTYPES_APOP, a.heatmap_scale, a.pool_fib, mods, htp), hm_note(SUBTYPES_APOP))]
    jobs = [j for j in jobs if ("heatmap" if "heatmap" in j[0] else "spatial") in want]
    for stem, fn, jnote in jobs:
        for dark in (True, False):
            base = stem + ctag + a.suffix + ("_dark" if dark else ""); fmts = [f.strip() for f in a.formats.split(",") if f.strip()]
            paths = {f: os.path.join(a.out, base + "." + f) for f in fmts}
            for path in paths.values():
                if os.path.exists(path) and not a.force:
                    sys.exit(f"refusing to overwrite {path} (use --force)")
            fig, info = fn(dark)
            for f, path in paths.items():
                fig.savefig(path, dpi=a.dpi)   # svg: rasterized scatters embed at this dpi, text stays vector
            plt.close(fig); path = paths[fmts[0]]
            side = {"figure": os.path.abspath(path), "files": [os.path.abspath(p) for p in paths.values()], "dpi": a.dpi, "written": datetime.datetime.now().isoformat(timespec="seconds"),
                    "note": jnote, "cmap": a.cmap, "aucell_scores": uns,
                    "embedding_outliers": {"hidden_in": "spatial panels only" if not a.show_flagged else "nowhere (--show-flagged)", "n_hidden_in_spatial_set": n_hid}, **info,
                    "script": os.path.abspath(__file__), "sources": srcs}
            json.dump(side, open(os.path.join(a.out, base + ".provenance.json"), "w"), indent=1)
            print("wrote", *paths.values())


if __name__ == "__main__":
    main()
