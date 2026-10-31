"""Draw the spatial grid of Fig. S13: one column per major cell class (magenta) over all other nuclei (grey).

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for S13.
Reads cells.csv.gz and extract_provenance.json from --tables-dir, and the display rotation of each stage from
--section-rotation (resources/figures/mouse_section_rotation.yaml, read through figures/section_rotation.py).
Writes <out-prefix>.pdf, .svg, .png and .provenance.json.
"""

import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mouse_tissue_style as style

# The rotation reader is shared by every mouse spatial panel and lives in workflow/scripts/figures/.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import section_rotation  # noqa: E402


# Square panel per stage x class, as in the manuscript's grid; label bands on the left and top.
PANEL_IN = 0.95
GAP_IN = 0.04
ROW_LABEL_IN = 0.3
COLUMN_TITLE_IN = 0.42

# Point area (pt^2): about one point across, the manuscript's dot size at this panel size.
POINT_AREA_PT2 = 0.9

# Room around each section, as a fraction of its larger extent.
SECTION_PAD = 0.03


def placed_nuclei(cells: pd.DataFrame, angles: dict) -> pd.DataFrame:
    """Spatial QC 'kept' nuclei with display coordinates x_plot/y_plot.

    Each stage's whole array (every library, both uninduced pieces, embedding outliers included) turns as one rigid
    body, so what is drawn afterwards sits exactly where it does in the array.
    """
    placed = cells[cells["spatial_qc_reason"] == "kept"].copy()
    placed["x_plot"], placed["y_plot"] = section_rotation.rotate_by_timepoint(
        placed["x_um_qc"], placed["y_um_qc"], placed["timepoint"], angles
    )
    return placed


def stage_nuclei(placed: pd.DataFrame, stage: str, uninduced_section: str) -> pd.DataFrame:
    """The nuclei drawn for one stage: placed, not an embedding outlier, one section."""
    nuclei = placed[placed["timepoint"] == stage]
    if stage == "uninduced":
        # Two sections share the uninduced array; the manuscript draws one of them.
        nuclei = nuclei[nuclei["tissue_section"].str.endswith(f"__{uninduced_section}")]
    return nuclei[~nuclei["embedding_outlier"].astype(bool)]


def square_limits(nuclei: pd.DataFrame):
    """Square x and y limits around a rotated section, so every column of a row shows it at one scale."""
    x_mid = (nuclei["x_plot"].min() + nuclei["x_plot"].max()) / 2
    y_mid = (nuclei["y_plot"].min() + nuclei["y_plot"].max()) / 2
    extent = max(np.ptp(nuclei["x_plot"]), np.ptp(nuclei["y_plot"])) * (1 + 2 * SECTION_PAD)
    return (x_mid - extent / 2, x_mid + extent / 2), (y_mid - extent / 2, y_mid + extent / 2), extent


def main() -> None:
    parser = style.panel_arguments(__doc__)
    parser.add_argument("--uninduced-section", required=True, choices=["A", "B"])
    parser.add_argument("--section-rotation", required=True, type=Path,
                        help="resources/figures/mouse_section_rotation.yaml: display rotation per stage")
    args = parser.parse_args()
    style.apply_panel_style()

    # Load ----

    cells = style.read_cells(args.tables_dir)
    classes = [c for c in style.SUPP_CLASS_ORDER if c in set(cells["cell_type"])]
    angles = section_rotation.load_rotation_table(args.section_rotation)
    placed = placed_nuclei(cells, angles)
    n_hidden = int(placed["embedding_outlier"].astype(bool).sum())

    # Figure ----

    # Fixed layout in inches: stages are rows, classes columns, as in the manuscript.
    n_rows, n_cols = len(style.STAGE_ORDER), len(classes)
    fig_w = ROW_LABEL_IN + n_cols * PANEL_IN + (n_cols - 1) * GAP_IN + 0.05
    fig_h = COLUMN_TITLE_IN + n_rows * PANEL_IN + (n_rows - 1) * GAP_IN + 0.05
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="black")

    per_stage = {}
    for row, stage in enumerate(style.STAGE_ORDER):
        nuclei = stage_nuclei(placed, stage, args.uninduced_section)
        xlim, ylim, extent_um = square_limits(nuclei)
        per_stage[stage] = {
            "n_drawn": int(len(nuclei)),
            "rotation_deg_ccw": angles[stage],
            "panel_extent_um": round(float(extent_um)),
            "n_per_class": nuclei["cell_type"].value_counts().reindex(classes, fill_value=0).to_dict(),
        }

        for col, cell_type in enumerate(classes):
            left = (ROW_LABEL_IN + col * (PANEL_IN + GAP_IN)) / fig_w
            bottom = 1 - (COLUMN_TITLE_IN + (row + 1) * PANEL_IN + row * GAP_IN) / fig_h
            ax = fig.add_axes([left, bottom, PANEL_IN / fig_w, PANEL_IN / fig_h])
            ax.set_facecolor("black")
            ax.set_axis_off()

            # All other nuclei in grey first, the highlighted class in magenta on top.
            is_class = nuclei["cell_type"] == cell_type
            for subset, color in ((nuclei[~is_class], style.GRID_OTHERS), (nuclei[is_class], style.GRID_HIGHLIGHT)):
                ax.scatter(
                    subset["x_plot"],
                    subset["y_plot"],
                    s=POINT_AREA_PT2,
                    color=color,
                    linewidths=0,
                    rasterized=True,
                )

            # y up, the convention of the rotation table: with the table's angles each section matches its
            # manuscript S13 map, which the unrotated image orientation (y down) did not; one scale per row.
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.set_aspect("equal")

            if row == 0:
                ax.text(0.5, 1.04, style.GRID_CLASS_TITLES[cell_type], transform=ax.transAxes,
                        ha="center", va="bottom", color="white", fontweight="bold")
            if col == 0:
                ax.text(-0.06, 0.5, style.GRID_STAGE_LABELS[stage], transform=ax.transAxes,
                        ha="right", va="center", rotation=90, color="white", fontweight="bold")

    logging.info("Drew %s nuclei per stage", {s: v["n_drawn"] for s, v in per_stage.items()})

    # Save ----

    provenance = style.base_provenance(
        "S13",
        args,
        inputs=[args.tables_dir / "cells.csv.gz", args.section_rotation],
        script=__file__,
    )
    provenance.update({
        "columns": classes,
        "column_titles": [style.GRID_CLASS_TITLES[c].replace("\n", " ") for c in classes],
        "coordinates": ("x_um_qc / y_um_qc rotated counter-clockwise by the stage's angle about the centre of the "
                        "stage's whole array (display only), drawn with y increasing upwards"),
        "rotation": section_rotation.table_record(args.section_rotation),
        "filters": {
            "spatial_qc_reason": "kept",
            "embedding_outlier_hidden_among_kept": n_hidden,
            "uninduced_section": f"tissue_section *__{args.uninduced_section} only (both lanes of the array)",
        },
        "scale": "each row (stage) scaled to fill its panels; columns of a row share one scale",
        "colors": {"highlight": style.GRID_HIGHLIGHT, "others": style.GRID_OTHERS},
        "per_stage": per_stage,
    })
    style.save_panel(fig, args.out_prefix, provenance, facecolor="black")
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
