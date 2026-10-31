"""Fig. S16: each fibroblast/decidual subtype (columns) on the mouse section of each stage (rows).

Reads the fib/dec plot table (extract_panel_table.py) and writes <out>.{pdf,svg,png} and
<out>.provenance.json. In every map the highlighted subtype is magenta, drawn over the stage's other
fibroblast and decidual nuclei in grey; black background, as in the manuscript panel. Nuclei: spatial QC
'kept' in the stage's section, embedding outliers omitted (release rule). Each section is turned by its stage's
angle in --section-rotation (the table of every mouse spatial panel) and scaled to fill its maps (no common
spatial scale, no scale bar), as in the manuscript.
Rule: panels_mouse_stroma_s16 (workflow/rules/91_panels_mouse_stroma.smk).
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import panel_style as ps

MAP_IN = 0.72  # side of one map (in); 9 columns fill the 180 mm figure width with the row labels
ROW_LABEL_IN = 0.3
HEADER_IN = 0.25
MARGIN_IN = 0.04
MAP_PAD = 1.03
HIGHLIGHT = "#FF00FF"  # magenta, as in the manuscript S13/S16 panels
OTHERS = "#A0A0A0"  # light grey context nuclei


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--table-prov", required=True, type=Path)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--section-rotation", required=True, type=Path,
                        help="resources/figures/mouse_section_rotation.yaml: display rotation per stage")
    args = parser.parse_args()

    foreground = ps.apply_style(dark=True)
    table, input_record = ps.load_table(args.table, args.table_prov)
    angles = ps.section_rotation.load_rotation_table(args.section_rotation)
    cells, counts = ps.placed_nuclei(table, ps.STAGES, angles)

    unknown = set(cells["cell_subtype"]) - set(ps.SUBTYPES)
    if unknown:
        raise ValueError(f"subtypes outside the S16 columns: {sorted(unknown)}")

    extent_mm = {
        stage: [round(float(np.ptp(g["x_plot"])) / 1000, 2), round(float(np.ptp(g["y_plot"])) / 1000, 2)]
        for stage, g in cells.groupby("timepoint")
    }
    # One marker size for every map, set by the median section, so densities compare across stages.
    area = ps.marker_area(int(cells.groupby("timepoint").size().median()), MAP_IN, coverage=0.55)

    fig_w = 2 * MARGIN_IN + ROW_LABEL_IN + len(ps.SUBTYPES) * MAP_IN
    fig_h = 2 * MARGIN_IN + HEADER_IN + len(ps.STAGES) * MAP_IN
    fig = plt.figure(figsize=(fig_w, fig_h))

    per_map = {}
    for row, stage in enumerate(ps.STAGES):
        section = cells[cells["timepoint"] == stage]
        half_width_um = MAP_PAD * max(np.ptp(section["x_plot"]), np.ptp(section["y_plot"])) / 2
        centre_x = (section["x_plot"].min() + section["x_plot"].max()) / 2
        centre_y = (section["y_plot"].min() + section["y_plot"].max()) / 2
        bottom = fig_h - MARGIN_IN - HEADER_IN - (row + 1) * MAP_IN

        for column, subtype in enumerate(ps.SUBTYPES):
            left = MARGIN_IN + ROW_LABEL_IN + column * MAP_IN
            ax = fig.add_axes([left / fig_w, bottom / fig_h, MAP_IN / fig_w, MAP_IN / fig_h])
            is_subtype = section["cell_subtype"] == subtype

            # Grey context first, highlighted subtype on top so sparse subtypes stay visible.
            ax.scatter(section.loc[~is_subtype, "x_plot"], section.loc[~is_subtype, "y_plot"],
                       s=area, c=OTHERS, linewidths=0, rasterized=True)
            ax.scatter(section.loc[is_subtype, "x_plot"], section.loc[is_subtype, "y_plot"],
                       s=area, c=HIGHLIGHT, linewidths=0, rasterized=True)
            ax.set_xlim(centre_x - half_width_um, centre_x + half_width_um)
            ax.set_ylim(centre_y - half_width_um, centre_y + half_width_um)
            ax.set_aspect("equal")
            ax.axis("off")
            per_map.setdefault(stage, {})[subtype] = int(is_subtype.sum())

            if row == 0:
                fig.text((left + MAP_IN / 2) / fig_w, (fig_h - MARGIN_IN - HEADER_IN / 2) / fig_h, subtype,
                         ha="center", va="center", color=foreground, fontsize=ps.SIZE_TITLE, fontweight="bold")

        label = ps.STAGE_LABELS[stage].replace("4.5 dpi menstrual", "4.5 dpi\nmenstrual")
        fig.text((MARGIN_IN + ROW_LABEL_IN / 2) / fig_w, (bottom + MAP_IN / 2) / fig_h, label, rotation=90,
                 ha="center", va="center", color=foreground, fontsize=ps.SIZE_TITLE, fontweight="bold",
                 linespacing=1.0)

    provenance = {
        "panel": args.panel,
        "script": Path(__file__).name,
        "inputs": input_record,
        "rows": ps.STAGES,
        "columns": ps.SUBTYPES,
        "nuclei": "fibroblast/decidual nuclei with spatial_qc_reason == 'kept' in the stage's section (lanes pooled; uninduced = piece B); embedding outliers omitted (release rule)",
        "grey": "the stage's other drawn fibroblast/decidual nuclei",
        "coordinates": ("x_um_qc, y_um_qc rotated counter-clockwise by the stage's angle about the centre of the "
                        "stage's whole array (display only), y up"),
        "rotation": ps.section_rotation.table_record(args.section_rotation),
        "per_stage": counts,
        "highlighted_nuclei": per_map,
        "marker_area_pt2": round(area, 3),
        "spatial_scale": "each section scaled to fill its maps (no common scale, no scale bar), as in the manuscript",
        "section_extent_mm_xy": extent_mm,
    }
    ps.save_panel(fig, args.out, provenance)


if __name__ == "__main__":
    main()
