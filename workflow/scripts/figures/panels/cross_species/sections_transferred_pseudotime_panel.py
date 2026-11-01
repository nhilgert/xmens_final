"""Fig. 6E: X-Mens pseudotime transferred onto the human secretory (top) and menstrual (bottom) sections, on black,
in the orientation of the main text.

Called by rule panels_cross_species_transferred_pseudotime (workflow/rules/91_panels_cross_species.smk).
Reads:  the joint SAMap table (transferred pseudotime of every nucleus; rule samap_figure_tables) and the human
        fibroblast/decidual table that rule samap_fig6e reads (QC'd positions, phase, embedding_outlier; a plot-table
        parquet or a release .h5ad).
Writes: the panel as PDF, SVG and 300-dpi PNG on black, and a provenance JSON with the counts and median per section;
        each section and the key also on their own (--parts-prefix P, --parts-tag T: P_secretory_T, P_menstrual_T,
        P_key_T, each as PDF, SVG and PNG) at the scale of the panel, for assembling the figure by hand.

Nuclei, values, color map, crop and key scale are those of rule samap_fig6e
(notebooks/samap/fig6e_sections_transferred_pseudotime.ipynb). Orientation is that of the main text and of Fig. 6D:
crop, mirror and rotation of each section come from spatial_efib_dec_panel.SECTIONS, so each section appears
mirrored left-right relative to the samap_fig6e panels. The key is labeled without samap_fig6e's "(rank)" and
"early"/"late".

usage: sections_transferred_pseudotime_panel.py --joint T --human-table T --pdf F --svg F --png F --provenance F
                                                --parts-prefix P --parts-tag T
"""
import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import cross_species_panels as common
from spatial_efib_dec_panel import CELL_MM, GAP_MM, SECTIONS, place


# Constants ----

# The transfer of record and its display, as in samap_fig6e: one fixed 0-1 viridis scale with values as stored
# (never re-ranked); nuclei without a valued mouse SAMap neighbor are drawn in grey underneath, not imputed.
TRANSFERRED = "transferred_mm_dpt_pseudotime_ranknorm"
CMAP = mpl.colormaps["viridis"]
NORM = mpl.colors.Normalize(0, 1)
NO_VALUE = "#B0B0B0"

# Sections from top to bottom. Each reuses the human row of Fig. 6D, so 6D and 6E share crop and orientation.
PHASES = {
    "secretory": SECTIONS[("human", "decidualized")],
    "menstrual": SECTIONS[("human", "menstrual")],
}

# Layout in mm: the square cells of 6D stacked with 6D's gap, then a band for the key.
KEY_BAND_MM = 12
KEY_TITLE_GAP_MM = 5.0  # cell bottom to the top of the key bar; the title sits in this gap
KEY_BAR_HEIGHT_MM = 2.2
KEY_BAR_SHARE = 0.6  # key bar width as a share of the panel width

# Dot area of the samap_fig6e panels (8 pt^2 on a 6 in square) scaled to a 6D cell, so dots cover the same share
# of the section as there.
NOTEBOOK_DOT_PT2 = 8
NOTEBOOK_PANEL_IN = 6
DOT_PT2 = NOTEBOOK_DOT_PT2 * (CELL_MM / 25.4 / NOTEBOOK_PANEL_IN) ** 2

# Key: the 0-1 scale of samap_fig6e, at the 7 pt of the other keys in this figure; title and tick labels without
# samap_fig6e's "(rank)" and "early"/"late".
KEY_TITLE = "transferred pseudotime"
KEY_TICKS = {0: "0", 0.5: "0.5", 1: "1"}
KEY_SIZE = 7

# Formats of the files with one section or the key alone.
PART_FORMATS = ("pdf", "svg", "png")

TABLE_COLUMNS = ["x_um_qc", "y_um_qc", "spatial_qc_reason", "timepoint", "embedding_outlier"]

LOG = logging.getLogger("sections_transferred_pseudotime_panel")


# Load ----

def load_human_nuclei(joint_path: Path, table_path: Path) -> pd.DataFrame:
    """Human nuclei of the joint table with QC'd position, phase and embedding_outlier from the human table."""
    joint = pd.read_parquet(joint_path, columns=["species", "sample", TRANSFERRED])
    if table_path.suffix == ".h5ad":
        table = common.load_release_obs(table_path, TABLE_COLUMNS)
    else:
        table = pd.read_parquet(table_path, columns=TABLE_COLUMNS)

    is_human = joint["species"].astype(str) == "hs"
    human = joint[is_human].join(table, how="inner")
    n_missing = int(is_human.sum()) - len(human)
    if n_missing:
        raise ValueError(f"{n_missing} human nuclei of the joint table are not in {table_path}")
    return human


# Select and place ----

def section_nuclei(human: pd.DataFrame, phase: str, spec: dict) -> tuple:
    """Drawn nuclei of one phase's section with display coordinates (xr, yr), and the counts behind each filter."""
    of_phase = human[human["timepoint"].astype(str) == phase]

    # Positions of record: nuclei the release marks 'kept'. samap_fig6e selects nuclei with QC'd coordinates;
    # the two selections must coincide, so the panels draw the same nuclei.
    is_kept = of_phase["spatial_qc_reason"].astype(str) == "kept"
    has_position = of_phase["x_um_qc"].notna() & of_phase["y_um_qc"].notna()
    if not (is_kept == has_position).all():
        n_disagree = int((is_kept != has_position).sum())
        raise ValueError(f"{phase}: 'kept' and QC'd coordinates disagree for {n_disagree} nuclei")
    kept = of_phase[is_kept]

    # Embedding outliers are a release display filter: kept in statistics, never drawn.
    is_outlier = kept["embedding_outlier"].astype(bool)
    candidates = kept[~is_outlier]

    if spec["crop_um"] is not None:
        (x0, x1), (y0, y1) = spec["crop_um"]
        in_crop = candidates["x_um_qc"].between(x0, x1) & candidates["y_um_qc"].between(y0, y1)
    else:
        in_crop = pd.Series(True, index=candidates.index)
    drawn = candidates[in_crop].copy()

    drawn["xr"], drawn["yr"] = place(
        drawn["x_um_qc"].to_numpy(float),
        drawn["y_um_qc"].to_numpy(float),
        spec["mirror_x"],
        spec["rotation_deg_ccw"],
    )
    counts = {
        "libraries": sorted(drawn["sample"].astype(str).unique()),
        "nuclei_of_phase": len(of_phase),
        "kept": len(kept),
        "embedding_outliers_omitted": int(is_outlier.sum()),
        "outside_crop": int((~in_crop).sum()),
        "drawn": len(drawn),
        "no_value": int(drawn[TRANSFERRED].isna().sum()),
        "median_transferred": round(float(drawn[TRANSFERRED].median()), 4),
    }
    return drawn, counts


# Figure ----

def draw_section(ax, drawn: pd.DataFrame) -> float:
    """Dots of one section fitted to its square cell; returns the micrometres spanned by the cell's side."""
    values = drawn[TRANSFERRED].to_numpy(float)
    has_value = np.isfinite(values)

    # Nuclei without a value in grey underneath, the others colored by transferred pseudotime on top. Table order
    # within each layer, as in samap_fig6e.
    ax.scatter(
        drawn["xr"][~has_value],
        drawn["yr"][~has_value],
        s=DOT_PT2,
        color=NO_VALUE,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    ax.scatter(
        drawn["xr"][has_value],
        drawn["yr"][has_value],
        s=DOT_PT2,
        c=CMAP(NORM(values[has_value])),
        linewidths=0,
        rasterized=True,
        zorder=2,
    )

    # Square limits around the section, widened by one dot radius (plus 1%) so edge dots are not clipped (as 6D).
    half = max(np.ptp(drawn["xr"]), np.ptp(drawn["yr"])) / 2
    radius_in = np.sqrt(DOT_PT2 / np.pi) / 72
    half = half / (1 - 2 * radius_in / (CELL_MM / 25.4)) * 1.01
    cx = (drawn["xr"].min() + drawn["xr"].max()) / 2
    cy = (drawn["yr"].min() + drawn["yr"].max()) / 2
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
    return 2 * half


def draw_key(fig, fig_w_mm: float, fig_h_mm: float, band_top_mm: float) -> None:
    """Horizontal color key centered under the sections, on the 0-1 scale of samap_fig6e."""
    bar_w_mm = KEY_BAR_SHARE * fig_w_mm
    bar_top_mm = band_top_mm + KEY_TITLE_GAP_MM
    cax = fig.add_axes([
        (fig_w_mm - bar_w_mm) / 2 / fig_w_mm,
        1 - (bar_top_mm + KEY_BAR_HEIGHT_MM) / fig_h_mm,
        bar_w_mm / fig_w_mm,
        KEY_BAR_HEIGHT_MM / fig_h_mm,
    ])
    bar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP),
        cax=cax,
        orientation="horizontal",
        ticks=list(KEY_TICKS),
    )
    bar.ax.set_xticklabels(list(KEY_TICKS.values()), color=common.TEXT, fontsize=KEY_SIZE)
    bar.outline.set_visible(False)
    bar.ax.tick_params(color=common.TEXT, length=2, width=0.6, pad=1.5)
    cax.set_title(KEY_TITLE, color=common.TEXT, fontsize=KEY_SIZE, pad=3)


def draw_panel(sections: dict, fig_w_mm: float, fig_h_mm: float):
    fig = plt.figure(figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4), facecolor=common.BACKGROUND)
    scale = {}
    for row, (phase, (drawn, _)) in enumerate(sections.items()):
        top_mm = row * (CELL_MM + GAP_MM)
        ax = fig.add_axes([0, 1 - (top_mm + CELL_MM) / fig_h_mm, 1, CELL_MM / fig_h_mm])
        ax.set_facecolor(common.BACKGROUND)
        ax.axis("off")
        scale[phase] = draw_section(ax, drawn) / CELL_MM

    draw_key(fig, fig_w_mm, fig_h_mm, band_top_mm=len(sections) * CELL_MM + (len(sections) - 1) * GAP_MM)
    return fig, scale


def draw_section_alone(drawn: pd.DataFrame):
    """One section in a figure of one square cell: the cell, crop and scale it has in the panel."""
    fig = plt.figure(figsize=(CELL_MM / 25.4, CELL_MM / 25.4), facecolor=common.BACKGROUND)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(common.BACKGROUND)
    ax.axis("off")
    draw_section(ax, drawn)
    return fig


def draw_key_alone():
    """The key in a figure of the key band, placed as under the panel."""
    fig = plt.figure(figsize=(CELL_MM / 25.4, KEY_BAND_MM / 25.4), facecolor=common.BACKGROUND)
    draw_key(fig, CELL_MM, KEY_BAND_MM, band_top_mm=0)
    return fig


# Run ----

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--joint", required=True, type=Path)
    parser.add_argument("--human-table", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--svg", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--parts-prefix", required=True, type=Path)
    parser.add_argument("--parts-tag", required=True, help="release version tag closing each part's file name")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    # fontTools logs every subset glyph at INFO when the PDF embeds the font.
    logging.getLogger("fontTools").setLevel(logging.WARNING)

    human = load_human_nuclei(args.joint, args.human_table)
    sections = {}
    for phase, spec in PHASES.items():
        sections[phase] = section_nuclei(human, phase, spec)
        LOG.info("%s: %s", phase, sections[phase][1])

    fig_w_mm = CELL_MM
    fig_h_mm = len(PHASES) * CELL_MM + (len(PHASES) - 1) * GAP_MM + KEY_BAND_MM
    font = common.set_font()
    fig, scale = draw_panel(sections, fig_w_mm, fig_h_mm)
    paths = {"pdf": args.pdf, "svg": args.svg, "png": args.png}
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    common.save_panel(fig, paths)
    plt.close(fig)
    LOG.info("Wrote %s (font %s)", [str(path) for path in paths.values()], font)

    # Each section and the key on their own, for assembling the figure by hand.
    part_paths = {}
    for part in [*sections, "key"]:
        part_fig = draw_key_alone() if part == "key" else draw_section_alone(sections[part][0])
        part_paths[part] = {ext: Path(f"{args.parts_prefix}_{part}_{args.parts_tag}.{ext}") for ext in PART_FORMATS}
        common.save_panel(part_fig, part_paths[part])
        plt.close(part_fig)
    LOG.info("Wrote the sections and the key alone: %s", list(part_paths))

    inputs = {"joint": args.joint, "human_table": args.human_table}
    record = common.base_record(Path(__file__), inputs, paths, font)
    record["panel"] = "6E"
    record["transferred_column"] = TRANSFERRED
    record["rescaled_or_reranked"] = False
    record["colormap"] = "viridis, fixed 0-1"
    record["no_value_color"] = NO_VALUE
    record["filters"] = [
        "human nuclei of the joint SAMap table, phase from the human table's timepoint",
        "spatial_qc_reason == 'kept' (identical to having x_um_qc/y_um_qc, the samap_fig6e selection; checked)",
        "embedding_outlier == False (flags of the human table given)",
        "inside crop_um, where a crop is given (raw coordinates, before mirror and rotation)",
    ]
    record["orientation"] = (
        "crop, mirror and rotation of the human row of Fig. 6D (spatial_efib_dec_panel.SECTIONS); "
        "mirror then rotation equals the samap_fig6e rotation followed by a left-right mirror"
    )
    record["layout_mm"] = {
        "cell": CELL_MM,
        "gap": GAP_MM,
        "key_band": KEY_BAND_MM,
        "figure": [fig_w_mm, fig_h_mm],
        "section_alone": [CELL_MM, CELL_MM],
        "key_alone": [CELL_MM, KEY_BAND_MM],
        "dot_pt2": round(DOT_PT2, 4),
    }

    # Per section: its definition, the counts behind each filter, and the scale of its cell.
    record["sections"] = {}
    for phase, (_, counts) in sections.items():
        scale_record = {"um_per_mm_of_panel": round(scale[phase], 2)}
        record["sections"][phase] = PHASES[phase] | counts | scale_record
    record["parts"] = {part: {ext: str(path) for ext, path in files.items()} for part, files in part_paths.items()}
    common.write_json(args.provenance, record)
    LOG.info("Wrote %s", args.provenance)


if __name__ == "__main__":
    main()
