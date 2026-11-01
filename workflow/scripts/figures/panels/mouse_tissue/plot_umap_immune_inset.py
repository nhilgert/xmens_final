"""Draw the mouse whole-tissue UMAP (Fig. 3A) with the immune-subtype UMAP (Fig. 3C) as an inset.

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for panel 3Ainset.
Reads cells.csv.gz and extract_provenance.json from --tables-dir; the immune subtypes' order, names,
colors and hand-set label directions from --vocabulary and --immune-display.
Writes <out-prefix>.pdf, .svg, .png and .provenance.json.

The atlas is drawn as in 3A and the inset as in 3C: the same embeddings, colors, label names and
embedding-outlier rule (taken from plot_umap.py), with the inset's own hand-set label directions. A box around the
atlas's immune nuclei is joined to the inset's frame by two connector lines. The inset sits in a
corner of the panel; each corner is tried, the atlas is drawn as large as it can be while leaving
that corner empty, and the corner allowing the largest atlas is kept.
"""

import logging

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.transforms import Bbox

import mouse_tissue_style as style
from plot_umap import LABEL_MARGIN, with_display


# Sizes in inches. The panel takes the width of panels 3A and 3C together in the manuscript layout
# (two quarters of a 180 mm page). The inset has the proportions of panel 3C (2.3 x 1.5 in, 8 pt
# labels) at 7 pt labels, so 3C's label margin and hand-set label directions carry over.
FIG_SIZE_IN = (3.6, 2.4)
INSET_SIZE_IN = (1.9, 1.25)
EDGE_GAP_IN = 0.02        # between the inset frame and the panel edge
CLEARANCE_IN = 0.05       # kept free of atlas points around the inset frame
LABEL_ROOM_IN = 0.2       # left around the atlas for its outside labels

# Atlas scales tried, as fractions of the largest scale that fits the panel; for each corner the
# first scale that leaves the corner empty is used.
SCALE_FRACTIONS = np.linspace(1.0, 0.5, 26)
CORNERS = ("lower right", "lower left", "upper right", "upper left")

# Inset points smaller than in 3C (2.0) in proportion to the smaller frame; inset labels one font
# step below the atlas labels. Frame, island box and connectors share one thin grey line.
INSET_POINT_SIZE = 1.0

# Label directions set by hand for the inset: display level "inset" of the display table, chosen by eye
# on the rendered v9 panel; without that table the directions of 3C (level "immune") are used.
LINE_COLOR = "#808080"
LINE_WIDTH_PT = 0.5
IMMUNE_BOX_PAD_PX = 3.0
# The immune box spans this central share of the atlas's immune nuclei in each direction, so a few
# immune nuclei sitting among other lineages do not stretch it.
IMMUNE_BOX_PERCENTILES = (0.5, 99.5)


def drawn_nuclei(cells, level, seed):
    """Coordinates and groups of the nuclei an embedding shows, in a seeded random order.

    Flagged embedding outliers sit among another lineage's nuclei; the release rule (uns
    embedding_outlier_rule) keeps them in every count and leaves them off embeddings. The shuffle
    keeps any one group from being painted over the others wholesale where islands touch.
    """
    shown = cells[~cells["embedding_outlier"].astype(bool)]
    shown = shown.iloc[np.random.default_rng(seed).permutation(len(shown))]
    return shown[level["coords"]].to_numpy(), shown[level["group"]].to_numpy(), len(cells) - len(shown)


def corner_layout(atlas_xy, corner):
    """Largest atlas scale (inches per UMAP unit) that leaves `corner` free for the inset.

    The atlas is pushed into the opposite corner, LABEL_ROOM_IN from the panel edges. Returns
    {scale, offset_in, inset_in} or None when no tried scale clears the corner.
    """
    width, height = FIG_SIZE_IN
    inset_w, inset_h = INSET_SIZE_IN
    low = atlas_xy.min(axis=0)
    span = atlas_xy.max(axis=0) - low
    largest = min((width - 2 * LABEL_ROOM_IN) / span[0], (height - 2 * LABEL_ROOM_IN) / span[1])

    right, upper = "right" in corner, "upper" in corner
    inset_in = Bbox.from_bounds(
        width - EDGE_GAP_IN - inset_w if right else EDGE_GAP_IN,
        height - EDGE_GAP_IN - inset_h if upper else EDGE_GAP_IN,
        inset_w,
        inset_h,
    )
    for fraction in SCALE_FRACTIONS:
        scale = largest * fraction
        size_in = span * scale
        offset_in = np.array([
            LABEL_ROOM_IN if right else width - LABEL_ROOM_IN - size_in[0],
            LABEL_ROOM_IN if upper else height - LABEL_ROOM_IN - size_in[1],
        ])
        points_in = (atlas_xy - low) * scale + offset_in
        if inset_in.padded(CLEARANCE_IN).count_contains(points_in) == 0:
            return {"scale": scale, "fraction_of_largest": float(fraction), "low": low,
                    "offset_in": offset_in, "inset_in": inset_in}
    return None


def facing_corners(box, other):
    """The two corners of `box` that face `other` (both Bboxes in one coordinate system)."""
    dx = (other.x0 + other.x1 - box.x0 - box.x1) / 2
    dy = (other.y0 + other.y1 - box.y0 - box.y1) / 2
    if abs(dx) >= abs(dy):
        x = box.x1 if dx > 0 else box.x0
        return [(x, box.y0), (x, box.y1)]
    y = box.y1 if dy > 0 else box.y0
    return [(box.x0, y), (box.x1, y)]


def main() -> None:
    parser = style.panel_arguments(__doc__)
    args = parser.parse_args()
    style.apply_panel_style()

    # Load ----

    cells = style.read_cells(args.tables_dir)
    immune_cells = cells[cells["cell_type"] == "immune_cell"]
    atlas = with_display("class", cells["cell_type"], args)
    immune = with_display("immune", immune_cells["cell_subtype"], args)
    inset_angles = immune["display"]["angles"].get("inset", immune["angles"])
    atlas_xy, atlas_groups, atlas_hidden = drawn_nuclei(cells, atlas, args.seed)
    immune_xy, immune_groups, immune_hidden = drawn_nuclei(immune_cells, immune, args.seed)
    logging.info("Atlas: %d nuclei drawn, %d outliers hidden; inset: %d drawn, %d hidden",
                 len(atlas_xy), atlas_hidden, len(immune_xy), immune_hidden)

    # Layout: the corner that allows the largest atlas; ties go to the corner nearest the immune island.

    island_centre = np.median(atlas_xy[atlas_groups == "immune_cell"], axis=0)
    candidates = []
    for corner in CORNERS:
        layout = corner_layout(atlas_xy, corner)
        if layout is None:
            continue
        island_in = (island_centre - layout["low"]) * layout["scale"] + layout["offset_in"]
        distance_in = float(np.hypot(*(island_in - layout["inset_in"].get_points().mean(axis=0))))
        candidates.append((-layout["scale"], distance_in, corner, layout))
    if not candidates:
        raise RuntimeError("no corner can be cleared for the inset; enlarge FIG_SIZE_IN")
    _, distance_in, corner, layout = min(candidates, key=lambda c: c[:2])
    logging.info("Inset corner: %s (atlas at %.2f of its largest scale)", corner, layout["fraction_of_largest"])

    # Atlas ----

    width, height = FIG_SIZE_IN
    fig = plt.figure(figsize=FIG_SIZE_IN)
    ax = fig.add_axes([0, 0, 1, 1])
    data_low = layout["low"] - layout["offset_in"] / layout["scale"]
    ax.set_xlim(data_low[0], data_low[0] + width / layout["scale"])
    ax.set_ylim(data_low[1], data_low[1] + height / layout["scale"])
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    ax.scatter(atlas_xy[:, 0], atlas_xy[:, 1], c=[atlas["colors"][g] for g in atlas_groups],
               s=atlas["point_size"], linewidths=0, rasterized=True)
    fig.canvas.draw()

    # Box around the atlas's immune nuclei (all their islands), in display pixels and in data.
    immune_px = ax.transData.transform(atlas_xy[atlas_groups == "immune_cell"])
    low_px, high_px = np.percentile(immune_px, IMMUNE_BOX_PERCENTILES, axis=0)
    immune_box_px = Bbox([low_px, high_px]).padded(IMMUNE_BOX_PAD_PX)
    immune_box = immune_box_px.transformed(ax.transData.inverted())
    ax.add_patch(Rectangle((immune_box.x0, immune_box.y0), immune_box.width, immune_box.height,
                           fill=False, edgecolor=LINE_COLOR, linewidth=LINE_WIDTH_PT))

    # Atlas labels as in 3A, kept off the inset (with its clearance) and off the immune box.
    inset_px = layout["inset_in"].transformed(fig.dpi_scale_trans)
    atlas_labels = style.place_direct_labels(
        ax, atlas_xy, atlas_groups, atlas["order"], atlas["labels"],
        atlas["angles"],
        n_orders=atlas["label_orders"], seed=args.seed,
        keep_out=[inset_px.padded(CLEARANCE_IN * fig.dpi), immune_box_px],
    )

    # Inset ----

    inset = layout["inset_in"]
    ax_inset = fig.add_axes([inset.x0 / width, inset.y0 / height, inset.width / width, inset.height / height])
    ax_inset.scatter(immune_xy[:, 0], immune_xy[:, 1], c=[immune["colors"][g] for g in immune_groups],
                     s=INSET_POINT_SIZE, linewidths=0, rasterized=True)
    ax_inset.set_aspect("equal", adjustable="datalim")
    ax_inset.margins(LABEL_MARGIN["immune"])
    ax_inset.set_xticks([])
    ax_inset.set_yticks([])
    for spine in ax_inset.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor(LINE_COLOR)
        spine.set_linewidth(LINE_WIDTH_PT)

    # Inset labels as in 3C, one font step smaller and inside the frame.
    with matplotlib.rc_context({"font.size": style.FONT_SIZE_SMALL}):
        fig.canvas.draw()
        inset_labels = style.place_direct_labels(
            ax_inset, immune_xy, immune_groups, immune["order"], immune["labels"],
            inset_angles,
            n_orders=immune["label_orders"], seed=args.seed,
            bounds=ax_inset.bbox.padded(-2.0),
        )

    # Connectors from the immune box's corners to the facing corners of the inset frame, in data
    # and axes-fraction coordinates so they stay put when the files are written at other dpi. They
    # belong to the atlas axes, beneath its points, so they pass behind any island they cross.
    to_atlas_data = ax.transData.inverted().transform
    to_inset_fraction = ax_inset.transAxes.inverted().transform
    for start, end in zip(facing_corners(immune_box_px, inset_px), facing_corners(inset_px, immune_box_px)):
        ax.add_artist(ConnectionPatch(
            xyA=to_atlas_data(start), coordsA=ax.transData,
            xyB=to_inset_fraction(end).round(6), coordsB=ax_inset.transAxes,
            color=LINE_COLOR, linewidth=LINE_WIDTH_PT, zorder=0.5, clip_on=False,
        ))

    # Save ----

    provenance = style.base_provenance("3Ainset", args, inputs=[args.tables_dir / "cells.csv.gz"], script=__file__)
    provenance.update({
        "embedding": {"atlas": atlas["embedding"], "inset": immune["embedding"]},
        "nuclei": {"atlas": "all released nuclei", "inset": "released immune nuclei only"},
        "filters": {"embedding_outlier_hidden": {"atlas": atlas_hidden, "inset": immune_hidden}},
        "n_drawn": {"atlas": int(len(atlas_xy)), "inset": int(len(immune_xy))},
        "layout": {
            "inset_corner": corner,
            "atlas_scale_fraction_of_largest": layout["fraction_of_largest"],
            "inset_in": [round(float(v), 3) for v in inset.bounds],
            "island_to_inset_in": round(distance_in, 3),
        },
        "colors_source": "as panels 3A and 3C (style module; immune subtypes: the --immune-display table)",
        "immune_display": str(args.immune_display),
        "introduced": style.introduced_labels(set(atlas_groups) | set(immune_groups), args.immune_display),
        **immune["record"],
        "labels": {"atlas": atlas_labels, "inset": inset_labels},
    })
    style.save_panel(fig, args.out_prefix, provenance)
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
