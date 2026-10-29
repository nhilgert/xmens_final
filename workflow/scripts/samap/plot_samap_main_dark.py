"""Main-figure SAMap panels on black: mouse pseudotime transferred onto human nuclei.
Called by rule samap_main_panels (workflow/rules/10c_samap_figures.smk).

  umap_dark_transferred_pseudotime            joint SAMap UMAP; human nuclei colored by transferred pseudotime,
                                              mouse nuclei dark grey underneath for context
  umap_dark_transferred_pseudotime_humanonly  same, without the mouse nuclei
  <date>-Human-<phase>-transferred_pseudotime-spatial   bare spatial panel per human section (no text, black)
  <date>-spatial_key_transferred_pseudotime   horizontal colorbar matching the spatial panels

Transferred pseudotime = SAMap transfer of the mouse global rank pseudotime (pt_harmonized_rank), i.e. the
`transferred_mm_dpt_pseudotime_ranknorm` column of samap_transferred_pt_harmonized.h5ad. Values are shown as stored on
a fixed 0-1 viridis scale; nothing is re-ranked (the Oct 14 version re-ranked them within the human nuclei).
Each human phase is one physical section: its three 10x libraries are drawn together in their shared x_um_qc/y_um_qc
frame. Frames and rotations are those of the Oct 14 panels. Nuclei flagged embedding_outlier are omitted from both
plot types (display rule used for all final pseudotime figures).

usage: plot_samap_main_dark.py --joint PARQUET --human TABLE --out DIR [--date 20261027]
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D

TRANSFERRED = "transferred_mm_dpt_pseudotime_ranknorm"
CMAP = mpl.colormaps["viridis"]
NORM = mpl.colors.Normalize(0, 1)

# Dark-theme colors, as in the Oct 14 panels
MOUSE_CONTEXT = "#3a3a3a"     # mouse nuclei under the human ones
NO_VALUE = "#9a9a9a"          # human nuclei without a mouse SAMap neighbor carrying pseudotime (UMAP)
NO_VALUE_SPATIAL = "#B0B0B0"  # same, spatial panels
TEXT, AXIS_ARROWS = "white", "#bbbbbb"

UMAP_DOT = 0.5                # points^2
SPATIAL_DOT = 8               # points^2 on a 6 in, 300 dpi bare panel
SPATIAL_IN, SPATIAL_DPI = 6, 300

# Per human phase: crop (x, y limits in um, or None for the full section) and rotation (degrees, counterclockwise)
SECTIONS = {
    "secretory": dict(limits=[(400, 5500), (750, 5250)], rotation=-90),
    "menstrual": dict(limits=None, rotation=90),
}


def set_font() -> None:
    helvetica = Path("/System/Library/Fonts/HelveticaNeue.ttc")
    if helvetica.exists():
        font_manager.fontManager.addfont(str(helvetica))
    mpl.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "DejaVu Sans"],
                         "svg.fonttype": "none", "pdf.fonttype": 42})


# ---- UMAP ----

def square_limits(umap: pd.DataFrame) -> tuple:
    # Square frame around the full joint embedding (both species), so every UMAP panel shares one frame.
    x, y = umap["X_umap__0"], umap["X_umap__1"]
    side = max(x.max() - x.min(), y.max() - y.min()) * 1.06
    cx, cy = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2
    return (cx - side / 2, cx + side / 2), (cy - side / 2, cy + side / 2), side


def umap_figure(title: str, xlim, ylim, side):
    fig = plt.figure(figsize=(90 / 25.4, 90 / 25.4), facecolor="black")
    ax = fig.add_axes([0.06, 0.06, 0.88, 0.84])
    ax.set_facecolor("black")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(0.06, 0.965, title, ha="left", va="top", fontsize=9, color=TEXT)

    # Short axis arrows in place of a box
    length = 0.14 * side
    arrow = dict(arrowprops=dict(arrowstyle="-|>", lw=0.8, color=AXIS_ARROWS, mutation_scale=6, shrinkA=0, shrinkB=0),
                 annotation_clip=False)
    ax.annotate("", xy=(xlim[0] + length, ylim[0]), xytext=(xlim[0], ylim[0]), **arrow)
    ax.annotate("", xy=(xlim[0], ylim[0] + length), xytext=(xlim[0], ylim[0]), **arrow)
    ax.text(xlim[0] + 0.015 * side, ylim[0] - 0.012 * side, "UMAP 1", ha="left", va="top", fontsize=7, color=AXIS_ARROWS)
    ax.text(xlim[0] - 0.015 * side, ylim[0] + 0.015 * side, "UMAP 2", ha="right", va="bottom", rotation=90,
            fontsize=7, color=AXIS_ARROWS)
    return fig, ax


def scatter_umap(ax, cells: pd.DataFrame, colors, zorder: int, size: float = UMAP_DOT) -> None:
    ax.scatter(cells["X_umap__0"], cells["X_umap__1"], s=size, c=colors, linewidths=0, rasterized=True, zorder=zorder)


def colorbar_vertical(fig) -> None:
    cax = fig.add_axes([0.80, 0.66, 0.025, 0.20])
    bar = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, ticks=[0, 0.5, 1])
    bar.ax.set_yticklabels(["0 early", "0.5", "1 late"], fontsize=7, color=TEXT)
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=1.5, pad=1.5, color=TEXT)
    cax.set_title("pseudotime\n(rank)", fontsize=7, loc="left", pad=4, color=TEXT)


def legend_bottom_right(fig, entries) -> None:
    handles = [Line2D([], [], ls="", marker="o", ms=4, mfc=color, mec="none", label=label) for color, label in entries]
    kax = fig.add_axes([0.56, 0.0, 0.44, 0.16])
    kax.axis("off")
    kax.legend(handles=handles, loc="lower left", fontsize=7, handletextpad=0.3, labelspacing=0.4, borderaxespad=0,
               labelcolor=TEXT, frameon=False)


# ---- spatial ----

def rotate(x, y, degrees, cx, cy):
    if degrees % 360 == 0:
        return x, y
    t = np.deg2rad(degrees)
    return (x - cx) * np.cos(t) - (y - cy) * np.sin(t) + cx, (x - cx) * np.sin(t) + (y - cy) * np.cos(t) + cy


def section_frame(cells: pd.DataFrame, limits, degrees) -> pd.DataFrame:
    x, y = cells["x_um_qc"].to_numpy(float), cells["y_um_qc"].to_numpy(float)
    if limits is not None:
        (x0, x1), (y0, y1) = limits
        keep = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        cells, x, y = cells[keep], x[keep], y[keep]
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    else:
        cx, cy = (x.min() + x.max()) / 2, (y.min() + y.max()) / 2
    xr, yr = rotate(x, y, degrees, cx, cy)
    return cells.assign(xr=xr, yr=yr)


def marker_padding(span_x, span_y, dot_size) -> tuple:
    # Pad the data range by one marker radius (+2 px) so edge dots are not clipped on the bare, axis-free panel.
    radius_px = np.sqrt(dot_size / np.pi) / 72 * SPATIAL_DPI + 2
    side_px = SPATIAL_IN * SPATIAL_DPI
    pad = lambda span: span * (radius_px / side_px) / (1 - 2 * radius_px / side_px) + 0.005 * span
    return pad(span_x), pad(span_y)


def spatial_panel(frame: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(SPATIAL_IN, SPATIAL_IN), facecolor="black")
    ax.set_position([0, 0, 1, 1])
    ax.set_facecolor("black")
    pad_x, pad_y = marker_padding(np.ptp(frame.xr), np.ptp(frame.yr), SPATIAL_DOT)
    ax.set_xlim(frame.xr.min() - pad_x, frame.xr.max() + pad_x)
    ax.set_ylim(frame.yr.min() - pad_y, frame.yr.max() + pad_y)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    values = frame[TRANSFERRED].to_numpy()
    missing = ~np.isfinite(values)
    ax.scatter(frame.xr[missing], frame.yr[missing], s=SPATIAL_DOT, color=NO_VALUE_SPATIAL, linewidths=0, zorder=1)
    ax.scatter(frame.xr[~missing], frame.yr[~missing], s=SPATIAL_DOT, c=CMAP(NORM(values[~missing])), linewidths=0, zorder=2)
    fig.savefig(path, dpi=SPATIAL_DPI, facecolor="black")
    fig.savefig(path.with_suffix(".pdf"), dpi=SPATIAL_DPI, facecolor="black")
    plt.close(fig)


def spatial_key(path: Path) -> None:
    fig = plt.figure(figsize=(3.0, 1.0), facecolor="black")
    cax = fig.add_axes([0.12, 0.45, 0.76, 0.2])
    bar = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
    bar.ax.set_xticklabels(["0 early", "0.5", "1 late"], color=TEXT, fontsize=8)
    bar.outline.set_visible(False)
    bar.ax.tick_params(color=TEXT, length=2)
    cax.set_title("transferred pseudotime (rank)", color=TEXT, fontsize=8, pad=4)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=300, facecolor="black")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--joint", required=True, help="joint SAMap export: species, X_umap__0/1, transferred columns")
    parser.add_argument("--human", required=True, help="human release plot table (x_um_qc, y_um_qc, timepoint, embedding_outlier)")
    parser.add_argument("--out", required=True)
    parser.add_argument("--date", default="20261027")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    set_font()

    joint = pd.read_parquet(args.joint)
    human_table = pd.read_parquet(args.human, columns=["x_um_qc", "y_um_qc", "timepoint", "embedding_outlier"])
    # timepoint and sample come from the release table; the joint export carries its own copies.
    human = joint[joint.species == "hs"].drop(columns=["timepoint"]).join(human_table, how="inner")
    if len(human) != (joint.species == "hs").sum():
        raise ValueError("human nuclei in the joint export are missing from the release table")
    flagged = set(human.index[human.embedding_outlier.astype(bool)])
    mouse = joint[joint.species == "mm"]
    human = human[~human.index.isin(flagged)]
    has_value = human[TRANSFERRED].notna()

    # UMAP panels. Seeded random order so neither early nor late nuclei systematically sit on top.
    xlim, ylim, side = square_limits(joint)
    shuffled = human[has_value].sample(frac=1, random_state=3)
    no_value_label = f"human, no mouse\nneighbor ({int((~has_value).sum())})"
    variants = {
        "umap_dark_transferred_pseudotime": [(MOUSE_CONTEXT, "mouse cells"), (NO_VALUE, no_value_label)],
        "umap_dark_transferred_pseudotime_humanonly": [(NO_VALUE, no_value_label)],
    }
    for stem, legend in variants.items():
        fig, ax = umap_figure("Pseudotime transferred onto human cells", xlim, ylim, side)
        if "humanonly" not in stem:
            scatter_umap(ax, mouse, MOUSE_CONTEXT, zorder=1)
        # Human nuclei without a value sit behind the colored ones (above the mouse context layer).
        scatter_umap(ax, human[~has_value], NO_VALUE, zorder=2, size=UMAP_DOT * 3)
        scatter_umap(ax, shuffled, CMAP(NORM(shuffled[TRANSFERRED].to_numpy())), zorder=3)
        colorbar_vertical(fig)
        legend_bottom_right(fig, legend)
        for ext in ("png", "svg", "pdf"):
            fig.savefig(out / f"{stem}.{ext}", dpi=600, facecolor="black")
        plt.close(fig)

    # Spatial panels, one per physical human section
    counts = {}
    for phase, spec in SECTIONS.items():
        cells = human[(human.timepoint.astype(str) == phase) & human.x_um_qc.notna()]
        frame = section_frame(cells, spec["limits"], spec["rotation"])
        spatial_panel(frame, out / f"{args.date}-Human-{phase}-transferred_pseudotime-spatial.png")
        counts[phase] = dict(libraries=sorted(frame["sample"].astype(str).unique()), drawn=int(len(frame)),
                             no_value=int(frame[TRANSFERRED].isna().sum()),
                             median_transferred=round(float(frame[TRANSFERRED].median()), 4),
                             limits_um=spec["limits"], rotation_deg=spec["rotation"])
    spatial_key(out / f"{args.date}-spatial_key_transferred_pseudotime.png")

    summary = dict(transferred_column=TRANSFERRED, rescaled_or_reranked=False, colormap="viridis, fixed 0-1",
                   embedding_outliers_omitted=len(flagged), umap_human_with_value=int(has_value.sum()),
                   umap_human_no_value=int((~has_value).sum()), umap_mouse_context=int(len(mouse)), spatial=counts)
    (out / "samap_main_dark_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
