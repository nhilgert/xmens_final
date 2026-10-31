"""Figs. S15B and S22: one gene's expression on the mouse Slide-tags sections, black background.

Reads the whole-tissue plot table of one release (extract_panel_table.py) and writes <out>.{pdf,svg,png}
and <out>.provenance.json. One map per stage, every placed nucleus of the section (spatial QC 'kept',
embedding outliers omitted), coloured by log-normalized expression divided by its maximum over the drawn
nuclei, so all maps share one 0-1 colour scale as in the manuscript panels. Expressing nuclei are drawn on
top. Each section is turned by its stage's angle in --section-rotation (the table of every mouse spatial panel)
and scaled to fill its map (no common spatial scale, no scale bar), as in the manuscript.
Rules: panels_mouse_stroma_s15b, panels_mouse_stroma_s22 (workflow/rules/91_panels_mouse_stroma.smk).
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

import panel_style as ps

MAP_IN = 2.0  # side of one section map (in); each section fills its map, as in the manuscript panels
TITLE_IN = 0.3  # band above each map for its stage title
COLORBAR_IN = 0.75  # band right of the maps for the colour key
MARGIN_IN = 0.08
CMAP = "viridis"  # perceptually uniform; 0 = not detected, 1 = the highest drawn value
MAP_PAD = 1.03  # maps clear the section by 3% on every side
MARKER_COVERAGE = 0.35  # summed marker area / map area; lower than the subtype grid so single nuclei stay distinct


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--table-prov", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--stages", required=True, help="comma-separated release timepoints, in reading order")
    parser.add_argument("--ncols", required=True, type=int)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--section-rotation", required=True, type=Path,
                        help="resources/figures/mouse_section_rotation.yaml: display rotation per stage")
    args = parser.parse_args()
    stages = args.stages.split(",")

    foreground = ps.apply_style(dark=True)
    table, input_record = ps.load_table(args.table, args.table_prov)
    angles = ps.section_rotation.load_rotation_table(args.section_rotation)
    cells, counts = ps.placed_nuclei(table, stages, angles)

    # Colour scale: expression over its maximum among the drawn nuclei (shared by every map).
    expression = cells[args.gene].to_numpy(dtype=float)
    max_drawn = float(expression.max())
    cells["colour_value"] = expression / max_drawn
    for stage in stages:
        in_stage = cells["timepoint"] == stage
        counts[stage]["nuclei_detected"] = int((cells.loc[in_stage, args.gene] > 0).sum())

    # Each section fills its map (the manuscript convention); physical extents go to the provenance.
    extent_mm = {
        stage: [round(float(np.ptp(g["x_plot"])) / 1000, 2), round(float(np.ptp(g["y_plot"])) / 1000, 2)]
        for stage, g in cells.groupby("timepoint")
    }

    # One marker size for every map, set by the median section, so densities compare across stages.
    marker_area = ps.marker_area(int(cells.groupby("timepoint").size().median()), MAP_IN, coverage=MARKER_COVERAGE)

    n_rows = int(np.ceil(len(stages) / args.ncols))
    fig_w = 2 * MARGIN_IN + args.ncols * MAP_IN + COLORBAR_IN
    fig_h = 2 * MARGIN_IN + n_rows * (MAP_IN + TITLE_IN)
    fig = plt.figure(figsize=(fig_w, fig_h))

    for i, stage in enumerate(stages):
        row, column = divmod(i, args.ncols)
        left = MARGIN_IN + column * MAP_IN
        bottom = fig_h - MARGIN_IN - (row + 1) * (MAP_IN + TITLE_IN)
        ax = fig.add_axes([left / fig_w, bottom / fig_h, MAP_IN / fig_w, MAP_IN / fig_h])

        # Position = tissue coordinates (display-rotated by the stage's angle); colour = scaled expression;
        # ascending draw order keeps the few expressing nuclei visible above the non-expressing majority.
        section = cells[cells["timepoint"] == stage].sort_values("colour_value", kind="stable")
        ax.scatter(
            section["x_plot"],
            section["y_plot"],
            c=section["colour_value"],
            cmap=CMAP,
            vmin=0,
            vmax=1,
            s=marker_area,
            linewidths=0,
            rasterized=True,
        )
        half_width_um = MAP_PAD * max(np.ptp(section["x_plot"]), np.ptp(section["y_plot"])) / 2
        centre_x = (section["x_plot"].min() + section["x_plot"].max()) / 2
        centre_y = (section["y_plot"].min() + section["y_plot"].max()) / 2
        ax.set_xlim(centre_x - half_width_um, centre_x + half_width_um)
        ax.set_ylim(centre_y - half_width_um, centre_y + half_width_um)
        ax.set_aspect("equal")
        ax.axis("off")
        fig.text((left + MAP_IN / 2) / fig_w, (bottom + MAP_IN + 0.06) / fig_h, ps.STAGE_LABELS[stage],
                 ha="center", va="bottom", color=foreground, fontsize=ps.SIZE_TITLE, fontweight="bold")

    # Colour key right of the maps, vertically centred, label reading upwards with the gene in italics.
    key_h_in = min(1.0, fig_h * 0.4)
    key_left_in = MARGIN_IN + args.ncols * MAP_IN + 0.1
    key_bottom_in = (fig_h - key_h_in) / 2
    cax = fig.add_axes([key_left_in / fig_w, key_bottom_in / fig_h, 0.08 / fig_w, key_h_in / fig_h])
    colorbar = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=CMAP), cax=cax, ticks=[])
    colorbar.outline.set_visible(False)
    key_centre_x = (key_left_in + 0.04) / fig_w
    fig.text(key_centre_x, (key_bottom_in + key_h_in + 0.04) / fig_h, "1", ha="center", va="bottom",
             color=foreground, fontsize=ps.SIZE_TITLE)
    fig.text(key_centre_x, (key_bottom_in - 0.04) / fig_h, "0", ha="center", va="top",
             color=foreground, fontsize=ps.SIZE_TITLE)
    label_parts = [(args.gene, "italic"), (" expression (log normalized)", "normal")]
    ps.mixed_style_text(fig, (key_left_in + 0.3) / fig_w, 0.5 - 0.95 / fig_h, label_parts,
                        rotation=90, color=foreground, fontsize=ps.SIZE_TITLE)

    provenance = {
        "panel": args.panel,
        "script": Path(__file__).name,
        "inputs": input_record,
        "gene": args.gene,
        "stages": stages,
        "nuclei": "spatial_qc_reason == 'kept' in the stage's section (lanes pooled; uninduced = piece B); embedding outliers omitted (release rule)",
        "coordinates": ("x_um_qc, y_um_qc rotated counter-clockwise by the stage's angle about the centre of the "
                        "stage's whole array (display only), y up"),
        "rotation": ps.section_rotation.table_record(args.section_rotation),
        "colour": f"X (log1p counts per 10k) / {max_drawn:.4f} (maximum over the drawn nuclei); viridis 0-1",
        "colour_max_log_norm": round(max_drawn, 5),
        "draw_order": "ascending colour value",
        "marker_area_pt2": round(marker_area, 3),
        "per_stage": counts,
        "drawn_nuclei_by_class": cells["coarse_refined"].value_counts().to_dict(),
        "drawn_nuclei_by_subtype": cells["cell_subtype"].value_counts().to_dict(),
        "pending_labels_drawn": ps.pending_label_notes(cells["cell_subtype"]),
        "spatial_scale": "each section scaled to fill its map (no common scale, no scale bar), as in the manuscript",
        "section_extent_mm_xy": extent_mm,
    }
    ps.save_panel(fig, args.out, provenance)


if __name__ == "__main__":
    main()
