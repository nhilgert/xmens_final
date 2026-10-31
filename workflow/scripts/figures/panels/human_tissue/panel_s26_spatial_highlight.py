"""Fig. S26: one spatial map per cell class (rows) and sample (columns), the class in magenta over gray.

Called by rule `panels_human_tissue_s26`, once per release.
Reads:  the plot table and extraction provenance of one release (extract_human_tissue.py).
Writes: <out-stem>.pdf/.svg/.png and <out-stem>.provenance.json.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

import panel_common as pc


# Row order of the published S26.
ROW_ORDER = [
    "epithelial_cell",
    "endometrial_fibroblast",
    "decidual_cell",
    "immune_cell",
    "endothelial_cell",
    "perivascular_cell",
    "smooth_muscle_cell",
]

# Geometry (in): seven rows of two square maps, row labels on the left and sample names on top,
# the same proportions as the published page-height supplement.
MAP_IN = 1.3
GAP_IN = 0.04
ROW_LABEL_IN = 0.34
HEADER_IN = 0.22
MARGIN_IN = 0.04
FIG_W_IN = ROW_LABEL_IN + 2 * MAP_IN + GAP_IN + MARGIN_IN
FIG_H_IN = HEADER_IN + len(ROW_ORDER) * MAP_IN + (len(ROW_ORDER) - 1) * GAP_IN + MARGIN_IN

# Highlight colours of the published S26: the class in magenta, every other nucleus in light gray.
HIGHLIGHT = "#FF00FF"
BACKGROUND_CELLS = "#AFAFAF"
FOREGROUND = "white"

# Marker area (pt^2) for about 24k nuclei per 1.3 in map; the highlighted class uses the same size,
# so its density is read against the gray tissue rather than inflated by larger marks.
MARKER_AREA = 0.5


def draw_highlight_map(ax, section, cell_class: str, half_width_um: float) -> None:
    # Gray first, then the highlighted class on top, so sparse classes (pv: about 110 nuclei per
    # section) remain visible.
    is_class = section["cell_type"] == cell_class
    for subset, color in ((section[~is_class], BACKGROUND_CELLS), (section[is_class], HIGHLIGHT)):
        ax.scatter(subset["x_plot"], subset["y_plot"], c=color, s=MARKER_AREA, linewidths=0, rasterized=True)
    ax.set_xlim(-half_width_um, half_width_um)
    ax.set_ylim(-half_width_um, half_width_um)
    ax.set_aspect("equal")
    ax.set_facecolor("black")
    ax.axis("off")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--extract-provenance", required=True, type=Path)
    parser.add_argument("--release", required=True)
    parser.add_argument("--out-stem", required=True, type=Path)
    args = parser.parse_args()

    pc.set_style()

    # Load: the nuclei of the 5B maps, in the same orientation and spatial scale ----
    cells = pc.load_cells(args.cells)
    drawn, counts = pc.spatial_cells(cells)
    half_width_um = pc.shared_half_width_um(drawn)

    # Figure: rows = classes, columns = samples, one shared spatial scale ----
    fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN), facecolor="black")
    for row, cell_class in enumerate(ROW_ORDER):
        top_in = HEADER_IN + row * (MAP_IN + GAP_IN)
        bottom = 1 - (top_in + MAP_IN) / FIG_H_IN
        for column, sample in enumerate(pc.SAMPLES):
            left = (ROW_LABEL_IN + column * (MAP_IN + GAP_IN)) / FIG_W_IN
            ax = fig.add_axes([left, bottom, MAP_IN / FIG_W_IN, MAP_IN / FIG_H_IN])
            draw_highlight_map(ax, drawn[drawn["timepoint"] == sample], cell_class, half_width_um)
            if row == len(ROW_ORDER) - 1 and column == 0:
                pc.draw_scale_bar(ax, half_width_um, FOREGROUND)

        # Row label, rotated, centred on the row.
        fig.text(
            (ROW_LABEL_IN / 2) / FIG_W_IN,
            bottom + (MAP_IN / 2) / FIG_H_IN,
            pc.FULL_NAMES[cell_class],
            rotation=90,
            ha="center",
            va="center",
            color=FOREGROUND,
            fontweight="bold",
            linespacing=1.1,
        )

    for column, sample in enumerate(pc.SAMPLES):
        fig.text(
            (ROW_LABEL_IN + column * (MAP_IN + GAP_IN) + MAP_IN / 2) / FIG_W_IN,
            1 - (HEADER_IN / 2) / FIG_H_IN,
            sample,
            ha="center",
            va="center",
            color=FOREGROUND,
            fontweight="bold",
        )

    outputs = pc.save_panel(fig, args.out_stem, facecolor="black")

    # Provenance ----
    scripts = [Path(__file__), Path(pc.__file__)]
    provenance = pc.base_provenance("S26", args.release, args.cells, args.extract_provenance, scripts)
    provenance.update({
        "filters": [
            "obs spatial_qc_reason == 'kept' (x_um_qc/y_um_qc present)",
            "obs embedding_outlier hidden (release rule: omitted from embedding and spatial plots)",
        ],
        "counts_by_sample": counts,
        "highlighted_by_sample_and_class": drawn.groupby(["timepoint", "cell_type"]).size().unstack(fill_value=0).to_dict(orient="index"),
        "row_order": ROW_ORDER,
        "display_transform": {s: {"mirror_x": m, "rotate_ccw_deg": r} for s, (m, r) in pc.DISPLAY_TRANSFORM.items()},
        "map_side_um": round(2 * half_width_um, 1),
        "scale_bar_um": pc.SCALE_BAR_UM,
        "marker_area_pt2": MARKER_AREA,
        "outputs": outputs,
    })
    provenance["colors"].update({"highlight": HIGHLIGHT, "other_nuclei": BACKGROUND_CELLS})
    pc.write_json(Path(f"{args.out_stem}.provenance.json"), provenance)


if __name__ == "__main__":
    main()
