"""Fig. 6D: eFib and dec nuclei on the tissue in X-Mens (top) and human (bottom), decidualized (left) and menstrual
(right) sections, on black.

Called by rule panels_cross_species_spatial (workflow/rules/91_panels_cross_species.smk).
Reads:  the release fibroblast/decidual objects of both species (obs only: cell_type, sample, tissue_section,
        spatial_qc_reason, x_um_qc, y_um_qc, embedding_outlier; timepoint for mouse), and the mouse section rotation
        table (--section-rotation, resources/figures/mouse_section_rotation.yaml).
Writes: the panel as PDF, SVG and 300-dpi PNG, and a provenance JSON with each section's nuclei, filters,
        transform and counts.

usage: spatial_efib_dec_panel.py --mouse-obj H5AD --human-obj H5AD --section-rotation YAML --pdf F --svg F --png F
                                 --provenance F
"""
import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import cross_species_panels as common

# The mouse rotation reader is shared by every mouse spatial panel and lives in workflow/scripts/figures/.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import section_rotation  # noqa: E402


# Constants ----

# Rows and columns as in the manuscript panel; the keys double as the printed labels.
ROWS = ("mouse", "human")
COLUMNS = ("decidualized", "menstrual")

# One physical Slide-tags section per cell of the grid. Libraries of one section share a coordinate frame and are
# drawn together (as in workflow/scripts/pseudotime/plot_pt_spatial_sections.py). Sections, crop, mirror and
# rotation reproduce the manuscript's Fig. 6D: each section was matched to its panel in the figure image by
# correlating binned eFib and dec densities over all rotations (5-degree grid, then 1 degree) with and without a
# mirror. Correlation of the chosen placement: mouse 1.5 dpi late r = 0.97 (the 3.5 dpi section reaches 0.26 at
# its best rotation, so the manuscript's decidualized mouse panel is 1.5 dpi late), mouse 4.5 dpi r = 0.98, human
# secretory r = 0.94, human menstrual r = 0.91 (the human quarter turns are within 0.003 of the best angle).
# The mouse sections are turned by the table every mouse spatial panel reads (--section-rotation,
# resources/figures/mouse_section_rotation.yaml), keyed by obs['timepoint']: 303 and 81 degrees, within 1 degree of
# the fitted 302 and 81, and not mirrored. The human crop is that of Fig. 6E. Human mirror (x -> -x) comes before
# the counter-clockwise rotation; both are display-only.
SECTIONS = {
    ("mouse", "decidualized"): {
        "label": "1.5 dpi late",
        "tissue_sections": ["1_dpi_late_mck12a__A", "1_dpi_late_mck12b__A"],
        "crop_um": None,
        "timepoint": "1.5_dpi_late",
    },
    ("mouse", "menstrual"): {
        "label": "4.5 dpi bleeding",
        "tissue_sections": ["4_dpi_bleeding_mck11a__A"],
        "crop_um": None,
        "timepoint": "4.5_dpi_bleeding",
    },
    ("human", "decidualized"): {
        "label": "secretory",
        "tissue_sections": ["secretory_mck6a__A", "secretory_mck6b__A", "secretory_mck6c__A"],
        "crop_um": [(400, 5500), (750, 5250)],
        "mirror_x": True,
        "rotation_deg_ccw": 90,
    },
    ("human", "menstrual"): {
        "label": "menstrual",
        "tissue_sections": ["menstrual_mck5a__A", "menstrual_mck5b__A", "menstrual_mck5c__A"],
        "crop_um": None,
        "mirror_x": True,
        "rotation_deg_ccw": 270,
    },
}

# Layout in mm: a band for the row labels on the left, one for the column headers on top, then a 2 x 2 grid of
# square cells. Each section fills its own cell, as in the manuscript, so the spatial scale differs between cells;
# the provenance record gives it per cell.
CELL_MM = 42
GAP_MM = 2
LABEL_BAND_MM = 6
HEADER_BAND_MM = 7
HEADER_SIZE = 9

ORDER_SEED = 3  # seeded random draw order, so neither eFib nor dec always sits on top

OBS_COLUMNS = [
    "cell_type",
    "sample",
    "tissue_section",
    "spatial_qc_reason",
    "x_um_qc",
    "y_um_qc",
    "embedding_outlier",
]

LOG = logging.getLogger("spatial_efib_dec_panel")


# Select and place ----

def section_nuclei(obs: pd.DataFrame, spec: dict, angles: dict) -> tuple:
    """Drawn nuclei of one section with display coordinates (xr, yr), and the counts behind each filter.

    angles: the mouse rotation table (section_rotation.load_rotation_table), used for sections keyed by timepoint.
    """
    libraries = sorted({name.split("__")[0] for name in spec["tissue_sections"]})
    of_libraries = obs[obs["sample"].astype(str).isin(libraries)]
    in_section = of_libraries[of_libraries["tissue_section"].astype(str).isin(spec["tissue_sections"])]

    # Positions of record: QC'd Slide-tags coordinates of nuclei the release marks 'kept'.
    is_kept = in_section["spatial_qc_reason"].astype(str) == "kept"
    kept = in_section[is_kept]
    if kept[["x_um_qc", "y_um_qc"]].isna().any().any():
        raise ValueError(f"{spec['label']}: 'kept' nuclei without QC'd coordinates")

    # Embedding outliers are a release display filter: kept in statistics, never drawn.
    is_outlier = kept["embedding_outlier"].astype(bool)
    candidates = kept[~is_outlier]

    if spec["crop_um"] is not None:
        (x0, x1), (y0, y1) = spec["crop_um"]
        in_crop = candidates["x_um_qc"].between(x0, x1) & candidates["y_um_qc"].between(y0, y1)
    else:
        in_crop = pd.Series(True, index=candidates.index)
    drawn = candidates[in_crop].copy()

    if "timepoint" in spec:
        # Mouse: the shared table, keyed by obs['timepoint']. The section's whole array (every library, outliers
        # included) turns as one body; the drawn nuclei then take their places in it.
        timepoints = set(kept["timepoint"].astype(str))
        if timepoints != {spec["timepoint"]}:
            raise ValueError(f"{spec['label']}: expected timepoint {spec['timepoint']}, the section holds {timepoints}")
        x_turned, y_turned = section_rotation.rotate_by_timepoint(
            kept["x_um_qc"], kept["y_um_qc"], kept["timepoint"], angles
        )
        turned = pd.DataFrame({"xr": x_turned, "yr": y_turned}, index=kept.index)
        drawn["xr"], drawn["yr"] = turned.loc[drawn.index, "xr"], turned.loc[drawn.index, "yr"]
    else:
        drawn["xr"], drawn["yr"] = place(
            drawn["x_um_qc"].to_numpy(float),
            drawn["y_um_qc"].to_numpy(float),
            spec["mirror_x"],
            spec["rotation_deg_ccw"],
        )
    counts = {
        "libraries": libraries,
        "nuclei_of_libraries": len(of_libraries),
        "spatial_qc_reason_of_libraries": count_values(of_libraries["spatial_qc_reason"]),
        "nuclei_in_section_kept": len(kept),
        "not_kept_in_section": int((~is_kept).sum()),
        "embedding_outliers_omitted": int(is_outlier.sum()),
        "outside_crop": int((~in_crop).sum()),
        "drawn": len(drawn),
        "drawn_by_cell_type": count_values(drawn["cell_type"].map(common.CELL_TYPE_LABELS)),
    }
    return drawn, counts


def count_values(values: pd.Series) -> dict:
    return {str(key): int(n) for key, n in values.astype(str).value_counts().sort_index().items()}


def place(x: np.ndarray, y: np.ndarray, mirror_x: bool, degrees_ccw: float) -> tuple:
    """Display coordinates: mirror x if asked, then rotate counter-clockwise about the bounding-box center."""
    dx = x - (x.min() + x.max()) / 2
    dy = y - (y.min() + y.max()) / 2
    if mirror_x:
        dx = -dx
    theta = np.deg2rad(degrees_ccw)
    return dx * np.cos(theta) - dy * np.sin(theta), dx * np.sin(theta) + dy * np.cos(theta)


def marker_area(n_nuclei: int, cell_in: float) -> float:
    """Dot area (pt^2) so dots cover a similar share of each cell whatever the section's size (as in Fig. 4D)."""
    return float(np.clip(0.45 * (cell_in * 72) ** 2 / n_nuclei, 0.4, 6.0))


# Figure ----

def cell_axes(fig, row: int, column: int, fig_w_mm: float, fig_h_mm: float):
    left = LABEL_BAND_MM + column * (CELL_MM + GAP_MM)
    top = HEADER_BAND_MM + row * (CELL_MM + GAP_MM)
    ax = fig.add_axes([left / fig_w_mm, 1 - (top + CELL_MM) / fig_h_mm, CELL_MM / fig_w_mm, CELL_MM / fig_h_mm])
    ax.set_facecolor(common.BACKGROUND)
    ax.axis("off")
    return ax


def draw_section(ax, drawn: pd.DataFrame, dot_area: float) -> float:
    """Dots of one section fitted to its square cell; returns the micrometres spanned by the cell's side."""
    # Color = cell type. One shuffled layer, so neither type systematically covers the other.
    order = np.random.default_rng(ORDER_SEED).permutation(len(drawn))
    shuffled = drawn.iloc[order]
    ax.scatter(
        shuffled["xr"],
        shuffled["yr"],
        s=dot_area,
        c=shuffled["cell_type"].astype(str).map(common.CELL_TYPE_COLORS).to_list(),
        linewidths=0,
        rasterized=True,
    )

    # Square limits around the section, widened by one dot radius (plus 1%) so edge dots are not clipped.
    half = max(np.ptp(drawn["xr"]), np.ptp(drawn["yr"])) / 2
    radius_in = np.sqrt(dot_area / np.pi) / 72
    half = half / (1 - 2 * radius_in / (CELL_MM / 25.4)) * 1.01
    cx = (drawn["xr"].min() + drawn["xr"].max()) / 2
    cy = (drawn["yr"].min() + drawn["yr"].max()) / 2
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
    return 2 * half


def draw_panel(sections: dict, fig_w_mm: float, fig_h_mm: float):
    fig = plt.figure(figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4), facecolor=common.BACKGROUND)
    scale = {}
    for (species, phase), (drawn, counts) in sections.items():
        ax = cell_axes(fig, ROWS.index(species), COLUMNS.index(phase), fig_w_mm, fig_h_mm)
        span_um = draw_section(ax, drawn, counts["dot_pt2"])
        scale[(species, phase)] = span_um / CELL_MM

    # Column headers above the grid and row labels left of it, bold white as in the manuscript.
    header = dict(ha="center", va="center", fontsize=HEADER_SIZE, fontweight="bold", color=common.TEXT)
    for column, label in enumerate(COLUMNS):
        x_mm = LABEL_BAND_MM + column * (CELL_MM + GAP_MM) + CELL_MM / 2
        fig.text(x_mm / fig_w_mm, 1 - (HEADER_BAND_MM / 2) / fig_h_mm, label, **header)
    for row, label in enumerate(ROWS):
        y_mm = HEADER_BAND_MM + row * (CELL_MM + GAP_MM) + CELL_MM / 2
        fig.text((LABEL_BAND_MM / 2) / fig_w_mm, 1 - y_mm / fig_h_mm, label, rotation=90, **header)
    return fig, scale


# Run ----

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mouse-obj", required=True, type=Path)
    parser.add_argument("--human-obj", required=True, type=Path)
    parser.add_argument("--section-rotation", required=True, type=Path,
                        help="resources/figures/mouse_section_rotation.yaml: display rotation of the mouse sections")
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--svg", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    # fontTools logs every subset glyph at INFO when the PDF embeds the font.
    logging.getLogger("fontTools").setLevel(logging.WARNING)

    obs = {
        "mouse": common.load_release_obs(args.mouse_obj, OBS_COLUMNS + ["timepoint"]),
        "human": common.load_release_obs(args.human_obj, OBS_COLUMNS),
    }
    angles = section_rotation.load_rotation_table(args.section_rotation)
    cell_in = CELL_MM / 25.4
    sections = {}
    for (species, phase), spec in SECTIONS.items():
        drawn, counts = section_nuclei(obs[species], spec, angles)
        if "timepoint" in spec:
            counts["mirror_x"] = False   # the rotation table has no mirror
            counts["rotation_deg_ccw"] = angles[spec["timepoint"]]
        counts["dot_pt2"] = round(marker_area(len(drawn), cell_in), 3)
        sections[(species, phase)] = (drawn, counts)
        LOG.info("%s %s (%s): %s", species, phase, spec["label"], counts)

    fig_w_mm = LABEL_BAND_MM + 2 * CELL_MM + GAP_MM
    fig_h_mm = HEADER_BAND_MM + 2 * CELL_MM + GAP_MM
    font = common.set_font()
    fig, scale = draw_panel(sections, fig_w_mm, fig_h_mm)
    paths = {"pdf": args.pdf, "svg": args.svg, "png": args.png}
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    common.save_panel(fig, paths)
    plt.close(fig)
    LOG.info("Wrote %s (font %s)", [str(path) for path in paths.values()], font)

    inputs = {"mouse_object": args.mouse_obj, "human_object": args.human_obj, "section_rotation": args.section_rotation}
    record = common.base_record(Path(__file__), inputs, paths, font)
    record["panel"] = "6D"
    record["colors"] = {common.CELL_TYPE_LABELS[key]: color for key, color in common.CELL_TYPE_COLORS.items()}
    record["filters"] = [
        "nuclei of the section's tissue_section values",
        "spatial_qc_reason == 'kept' (positions x_um_qc, y_um_qc)",
        "embedding_outlier == False",
        "inside crop_um, where a crop is given (raw coordinates, before mirror and rotation)",
    ]
    record["mouse_rotation"] = section_rotation.table_record(args.section_rotation)
    record["layout_mm"] = {
        "cell": CELL_MM,
        "gap": GAP_MM,
        "label_band": LABEL_BAND_MM,
        "header_band": HEADER_BAND_MM,
        "figure": [fig_w_mm, fig_h_mm],
    }

    # Per section: its definition, the counts behind each filter, and the scale of its cell.
    record["sections"] = {}
    for key, (_, counts) in sections.items():
        scale_record = {"um_per_mm_of_panel": round(scale[key], 2)}
        record["sections"]["/".join(key)] = SECTIONS[key] | counts | scale_record
    record["order_seed"] = ORDER_SEED
    common.write_json(args.provenance, record)
    LOG.info("Wrote %s", args.provenance)


if __name__ == "__main__":
    main()
