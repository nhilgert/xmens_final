"""Fig. 5B: spatial mapping of all cell classes in the secretory and menstrual sections, on black.

Called by rule `panels_human_tissue_5b`, once per release.
Reads:  the plot table and extraction provenance of one release (extract_human_tissue.py).
Writes: <out-stem>.pdf/.svg/.png and <out-stem>.provenance.json.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import panel_common as pc


# Geometry (in): two square maps side by side, legend band below; 3.4 in matches the width of
# panel B in the published Fig. 5 (about half of the 7.25 in page).
MAP_IN = 1.62
GAP_IN = 0.08
MARGIN_IN = 0.04
LEGEND_IN = 0.26
FIG_W_IN = 2 * MAP_IN + GAP_IN + 2 * MARGIN_IN
FIG_H_IN = MAP_IN + LEGEND_IN + MARGIN_IN

# Marker area (pt^2): about 24k nuclei per section in a 1.6 in map; at 0.9 pt^2 neighbouring nuclei
# touch, so the tissue reads as continuous as in the published panel.
MARKER_AREA = 0.9

FOREGROUND = "white"


def draw_map(ax, section, half_width_um: float) -> None:
    # Position = tissue coordinates in the registered display orientation; colour = cell class.
    # Rasterized so the PDF/SVG stay light, while labels and legend remain vector text.
    ax.scatter(
        section["x_plot"],
        section["y_plot"],
        c=section["cell_type"].map(pc.CLASS_COLORS),
        s=MARKER_AREA,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xlim(-half_width_um, half_width_um)
    ax.set_ylim(-half_width_um, half_width_um)
    ax.set_aspect("equal")
    ax.set_facecolor("black")
    ax.axis("off")


def add_class_legend(fig) -> None:
    # One row of swatches in the Fig. 5B order, text in the manuscript abbreviations.
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=4.2, markeredgewidth=0, color=pc.CLASS_COLORS[c])
        for c in pc.CLASS_ORDER
    ]
    labels = [pc.ABBREVIATIONS[c] for c in pc.CLASS_ORDER]
    fig.legend(
        handles,
        labels,
        loc="center",
        bbox_to_anchor=(0.5, (LEGEND_IN / 2) / FIG_H_IN),
        ncol=len(handles),
        handletextpad=0.05,
        columnspacing=0.9,
        labelcolor=FOREGROUND,
        fontsize=7,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--extract-provenance", required=True, type=Path)
    parser.add_argument("--release", required=True)
    parser.add_argument("--out-stem", required=True, type=Path)
    args = parser.parse_args()

    pc.set_style()

    # Load: one point per nucleus that passed spatial QC, embedding outliers hidden ----
    cells = pc.load_cells(args.cells)
    drawn, counts = pc.spatial_cells(cells)
    half_width_um = pc.shared_half_width_um(drawn)
    for sample, sample_counts in counts.items():
        pc.log.info(f"{args.release} {sample}: {sample_counts}")

    # Figure: both sections at one spatial scale, so their sizes compare directly ----
    fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN), facecolor="black")
    for i, sample in enumerate(pc.SAMPLES):
        left_in = MARGIN_IN + i * (MAP_IN + GAP_IN)
        ax = fig.add_axes([left_in / FIG_W_IN, LEGEND_IN / FIG_H_IN, MAP_IN / FIG_W_IN, MAP_IN / FIG_H_IN])
        draw_map(ax, drawn[drawn["timepoint"] == sample], half_width_um)
        if i == 0:
            pc.draw_scale_bar(ax, half_width_um, FOREGROUND)

    add_class_legend(fig)
    outputs = pc.save_panel(fig, args.out_stem, facecolor="black")

    # Provenance ----
    scripts = [Path(__file__), Path(pc.__file__)]
    provenance = pc.base_provenance("5B", args.release, args.cells, args.extract_provenance, scripts)
    provenance.update({
        "filters": [
            "obs spatial_qc_reason == 'kept' (x_um_qc/y_um_qc present)",
            "obs embedding_outlier hidden (release rule: omitted from embedding and spatial plots)",
        ],
        "counts_by_sample": counts,
        "drawn_by_sample_and_class": drawn.groupby(["timepoint", "cell_type"]).size().unstack(fill_value=0).to_dict(orient="index"),
        "display_transform": {s: {"mirror_x": m, "rotate_ccw_deg": r} for s, (m, r) in pc.DISPLAY_TRANSFORM.items()},
        "map_side_um": round(2 * half_width_um, 1),
        "scale_bar_um": pc.SCALE_BAR_UM,
        "marker_area_pt2": MARKER_AREA,
        "draw_order_seed": pc.SEED,
        "outputs": outputs,
    })
    provenance["colors"]["classes"] = {pc.ABBREVIATIONS[c]: pc.CLASS_COLORS[c] for c in pc.CLASS_ORDER}
    pc.write_json(Path(f"{args.out_stem}.provenance.json"), provenance)


if __name__ == "__main__":
    main()
