"""Fig. S26B-C: the human epithelium split into inner and outer layers, on the tissue (B) and on a UMAP (C).

Called by rules `panels_human_epithelium_s26b`, `panels_human_epithelium_s26c` and
`panels_human_epithelium_s26bc` (workflow/rules/91_panels_human_epithelium.smk), once per release.
Reads:  the release plot table, its uns export and extraction provenance (extract_human_tissue.py); for the
        lineage-fit UMAP also the table and provenance of extract_epithelial_lineage_umap.py.
Writes: <out-stem>.pdf/.svg/.png and <out-stem>.provenance.json of one layout: the spatial maps (S26B),
        the UMAP (S26C) or both side by side (S26BC).

S26B reuses the Fig. S26 grid unchanged: sections, display transform, shared spatial scale, point size,
map geometry and fonts come from panel_common.py and panel_s26_spatial_highlight.py, so S26B can sit
under the S26 grid as one more row.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_hex, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import panel_common as pc
import panel_s26_spatial_highlight as s26


# Vocabulary ----

# Epithelial subtypes (obs cell_subtype) in legend order: the two layers, then the ciliated cells.
SUBTYPE_ORDER = ["epithelial_inner", "epithelial_outer", "epithelial_ciliated"]

# The layer of an epithelial nucleus was decided by counting these two classes among its 20 nearest
# placed nuclei (decidual cells: inner, endometrial fibroblasts: outer); the decfib background shows them.
LAYER_REFERENCE_CLASSES = ["decidual_cell", "endometrial_fibroblast"]

ROW_LABEL = "epithelial\nsubtypes"


# Colours ----

# The manuscript shows no epithelial subtypes, so these colours are new: Okabe-Ito orange and blue for
# the inner/outer pair (distinct under the common colour-vision deficiencies, and legible on the black
# maps and on the white UMAP), Okabe-Ito bluish green for the ciliated cells. None is a Fig. 5B class
# colour.
SUBTYPE_COLORS = {
    "epithelial_inner": "#E69F00",
    "epithelial_outer": "#0072B2",
    "epithelial_ciliated": "#009E73",
}

# Backgrounds of the maps. grey: every other nucleus in the light grey of the S26 grid. decfib: decidual
# cells and endometrial fibroblasts in pale tints of their Fig. 5B colours, the rest in the S26 grey, so
# the stroma behind the layer rule is visible without competing with the epithelium. Fibroblasts are
# tinted further toward white than decidual cells: tinted equally, the two collapse under deuteranopia
# (simulated CIELAB distance 11, against 19 for these tints and 16 for the full Fig. 5B colours).
OTHER_NUCLEI = s26.BACKGROUND_CELLS
TINT_TOWARD_WHITE = {"decidual_cell": 0.3, "endometrial_fibroblast": 0.6}
BACKGROUNDS = ["grey", "decfib"]

# The UMAP sits on white like the manuscript's UMAPs (Fig. 5C); its other nuclei use the map colours.
UMAP_BACKGROUND = "white"
UMAP_TEXT = "black"


# Geometry (in) ----

# Maps, gaps, row-label band and header band are those of the S26 grid (s26.*), so S26B lines up under
# it column for column. The legend band below the maps holds three 6.5 pt rows.
LEGEND_IN = 0.42
LEGEND_GAP_IN = 0.06  # between the bottom of the maps (or the UMAP) and the first legend row
LEGEND_ROW_IN = 0.12  # one 6.5 pt legend row
LEGEND_MARKER_PT = 4.2  # swatch diameter of the Fig. 5B legend

# UMAP: as tall as a map row; its width follows the embedding's aspect (equal axis scales). The padding
# around the points is a fraction of the larger data span and clears the largest dot.
UMAP_H_IN = s26.MAP_IN
UMAP_GAP_IN = 0.16  # between the black map block and the UMAP in the side-by-side layout
UMAP_PAD = 0.04
UMAP_TITLE_GAP_IN = 0.05  # between the UMAP title's baseline and the top of the UMAP
UMAP_MIN_W_IN = 1.9  # standalone S26C: room for the legend's longest label

# Marker area (pt^2) on the UMAP: about 11,700 epithelial nuclei in a 1.3 in square; at 0.6 pt^2 the
# dots overlap inside the islands and stay separate at their edges.
UMAP_MARKER_AREA = 0.6

# Locator of the atlas option: the whole embedding in miniature with the zoom window boxed. The
# epithelial island of the whole-tissue embedding is isolated, so no other nucleus falls inside the zoom
# window and the zoomed axes alone cannot show where the island sits among the other nuclei. It goes into
# the first corner of the UMAP (in LOCATOR_CORNERS order) whose footprint covers no epithelial nucleus.
LOCATOR_IN = 0.36
LOCATOR_INSET_IN = 0.03  # gap between the locator and the UMAP's edge
LOCATOR_LABEL_IN = 0.11  # band above the locator for its label
LOCATOR_LABEL_GAP_IN = 0.02
LOCATOR_LABEL = "all nuclei"
LOCATOR_LABEL_PT = 6  # the scale-bar label size of the maps, the smallest step of the type ladder
LOCATOR_MARKER_AREA = 0.02
LOCATOR_CORNERS = ["lower right", "upper left", "lower left", "upper right"]
LOCATOR_FRAME = "0.6"

# The two UMAP options of Fig. 3C: the lineage fit on which the epithelial clusters were decided, and
# the release's whole-tissue embedding zoomed to its epithelial island.
EMBEDDINGS = {
    "lineage_fit": {
        "coords": ["lineage_umap_1", "lineage_umap_2"],
        "title": "epithelial lineage UMAP",
        "description": (
            "human epithelial lineage fit (latent_human_epithelial.h5ad) obsm X_umap_scVI as stored; "
            "released epithelial nuclei only"
        ),
    },
    "atlas": {
        "coords": ["umap_1", "umap_2"],
        "title": "whole-tissue UMAP",
        "description": (
            "release obsm X_umap_scVI (whole-tissue fit) as stored, axes zoomed to the epithelial nuclei; "
            "other nuclei in the window in their map colours, and a locator of the whole embedding"
        ),
    },
}


# Data ----

def lower_first(text: str) -> str:
    return text[:1].lower() + text[1:]


def pale_tint(color: str, toward_white: float) -> str:
    rgb = np.array(to_rgb(color))
    return to_hex(rgb + (1 - rgb) * toward_white)


def background_colors(background: str) -> dict:
    """Map colour of every non-epithelial class under one background option."""
    colors = {cell_class: OTHER_NUCLEI for cell_class in pc.CLASS_ORDER if cell_class != "epithelial_cell"}
    if background == "decfib":
        for cell_class in LAYER_REFERENCE_CLASSES:
            colors[cell_class] = pale_tint(pc.CLASS_COLORS[cell_class], TINT_TOWARD_WHITE[cell_class])
    return colors


def legend_entries(uns: dict, background: str) -> tuple[list, list]:
    """(label, colour) pairs: the epithelial subtypes, and the other nuclei under one background.

    Labels are the release's display names (uns display_names), lower-case as the S26 row labels.
    """
    names = uns["display_names"]
    subtypes = [(lower_first(names["cell_subtype"][s]), SUBTYPE_COLORS[s]) for s in SUBTYPE_ORDER]
    others = []
    if background == "decfib":
        colors = background_colors(background)
        others = [(lower_first(names["cell_type"][c]), colors[c]) for c in LAYER_REFERENCE_CLASSES]
    others.append(("other nuclei", OTHER_NUCLEI))
    return subtypes, others


def load_cells(path: Path) -> pd.DataFrame:
    """Plot table of one release; fails on an epithelial subtype the panels have no colour for."""
    cells = pc.load_cells(path)
    cells["cell_subtype"] = cells["cell_subtype"].astype(str)
    is_epithelial = cells["cell_type"] == "epithelial_cell"
    unknown = set(cells.loc[is_epithelial, "cell_subtype"]) - set(SUBTYPE_ORDER)
    if unknown:
        raise ValueError(f"epithelial cell_subtype without a colour: {sorted(unknown)}")
    return cells


def subtype_counts(epithelial: pd.DataFrame) -> dict:
    """Nuclei per epithelial subtype, zeros included."""
    return epithelial["cell_subtype"].value_counts().reindex(SUBTYPE_ORDER, fill_value=0).to_dict()


def count_table(frame: pd.DataFrame, rows: str, columns: str) -> dict:
    return frame.groupby([rows, columns]).size().unstack(fill_value=0).to_dict(orient="index")


# Spatial maps (S26B) ----

def draw_map(ax, section: pd.DataFrame, colors_by_class: dict, half_width_um: float) -> None:
    # Position = tissue coordinates in the S26 display orientation; colour = epithelial subtype, or the
    # background colour of the nucleus's class. Other nuclei first and the epithelium on top, so the thin
    # glands are not covered by the denser stroma; within each layer the seeded shuffle of
    # pc.spatial_cells decides the order, so neither layer systematically covers the other.
    is_epithelial = section["cell_type"] == "epithelial_cell"
    stroma = section[~is_epithelial]
    epithelium = section[is_epithelial]
    layers = (
        (stroma, stroma["cell_type"].map(colors_by_class)),
        (epithelium, epithelium["cell_subtype"].map(SUBTYPE_COLORS)),
    )
    for nuclei, colors in layers:
        ax.scatter(
            nuclei["x_plot"],
            nuclei["y_plot"],
            c=colors.to_list(),
            s=s26.MARKER_AREA,
            linewidths=0,
            rasterized=True,
        )
    ax.set_xlim(-half_width_um, half_width_um)
    ax.set_ylim(-half_width_um, half_width_um)
    ax.set_aspect("equal")
    ax.axis("off")


def add_legend(fig, entries: list, x: float, y_top: float, text_color: str) -> None:
    """One column of round swatches and labels, top-centred on (x, y_top) in figure coordinates."""
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=LEGEND_MARKER_PT, markeredgewidth=0, color=color)
        for _, color in entries
    ]
    fig.legend(
        handles,
        [label for label, _ in entries],
        loc="upper center",
        bbox_to_anchor=(x, y_top),
        ncol=1,
        handletextpad=0.05,
        labelspacing=0.3,
        borderpad=0,
        borderaxespad=0,
        labelcolor=text_color,
    )


def spatial_block_size_in() -> tuple[float, float]:
    width = s26.ROW_LABEL_IN + 2 * s26.MAP_IN + s26.GAP_IN + s26.MARGIN_IN
    height = s26.HEADER_IN + s26.MAP_IN + LEGEND_IN + s26.MARGIN_IN
    return width, height


def draw_spatial_block(fig, fig_size_in, drawn, half_width_um, background, uns) -> None:
    """S26B in the left part of the figure: sample headers, row label, both maps, scale bar, legend."""
    fig_w, fig_h = fig_size_in
    block_w, block_h = spatial_block_size_in()

    # Black ground of the S26 grid, also when the block shares the figure with the white UMAP.
    fig.add_artist(Rectangle(
        (0, 0),
        block_w / fig_w,
        1,
        transform=fig.transFigure,
        facecolor="black",
        linewidth=0,
        zorder=-1,
    ))

    colors_by_class = background_colors(background)
    map_bottom = 1 - (s26.HEADER_IN + s26.MAP_IN) / fig_h
    map_centers_x = []
    for column, sample in enumerate(pc.SAMPLES):
        left_in = s26.ROW_LABEL_IN + column * (s26.MAP_IN + s26.GAP_IN)
        ax = fig.add_axes([left_in / fig_w, map_bottom, s26.MAP_IN / fig_w, s26.MAP_IN / fig_h])
        draw_map(ax, drawn[drawn["timepoint"] == sample], colors_by_class, half_width_um)
        if column == 0:
            pc.draw_scale_bar(ax, half_width_um, s26.FOREGROUND)
        map_centers_x.append((left_in + s26.MAP_IN / 2) / fig_w)

        # Sample name centred over its map, as the S26 column headers.
        fig.text(
            map_centers_x[-1],
            1 - (s26.HEADER_IN / 2) / fig_h,
            sample,
            ha="center",
            va="center",
            color=s26.FOREGROUND,
            fontweight="bold",
        )

    # Row label, rotated and centred on the maps, as the S26 row labels.
    fig.text(
        (s26.ROW_LABEL_IN / 2) / fig_w,
        map_bottom + (s26.MAP_IN / 2) / fig_h,
        ROW_LABEL,
        rotation=90,
        ha="center",
        va="center",
        color=s26.FOREGROUND,
        fontweight="bold",
        linespacing=1.1,
    )

    # Legend under the maps: the epithelial subtypes under the first map, the other nuclei under the second.
    subtypes, others = legend_entries(uns, background)
    legend_top = map_bottom - LEGEND_GAP_IN / fig_h
    add_legend(fig, subtypes, map_centers_x[0], legend_top, s26.FOREGROUND)
    add_legend(fig, others, map_centers_x[1], legend_top, s26.FOREGROUND)


def spatial_provenance(cells, drawn, counts, half_width_um, background) -> dict:
    epithelial = cells[cells["cell_type"] == "epithelial_cell"]
    drawn_epithelial = drawn[drawn["cell_type"] == "epithelial_cell"]
    has_position = epithelial["spatial_qc_reason"].astype(str) == pc.SPATIAL_KEPT
    is_outlier = epithelial["embedding_outlier"].astype(bool)
    colors = background_colors(background)
    return {
        "background": background,
        "filters": [
            "obs spatial_qc_reason == 'kept' (x_um_qc/y_um_qc present)",
            "obs embedding_outlier hidden (release rule: omitted from embedding and spatial plots)",
        ],
        "sections": {
            sample: sorted(drawn.loc[drawn["timepoint"] == sample, "tissue_section"].astype(str).unique())
            for sample in pc.SAMPLES
        },
        "counts_by_sample": counts,
        "epithelial_drawn_by_sample_and_subtype": count_table(drawn_epithelial, "timepoint", "cell_subtype"),
        "epithelial_drawn_by_tissue_section_and_subtype": count_table(
            drawn_epithelial.assign(tissue_section=drawn_epithelial["tissue_section"].astype(str)),
            "tissue_section",
            "cell_subtype",
        ),
        "epithelial_not_drawn_no_position_by_subtype": subtype_counts(epithelial[~has_position]),
        "epithelial_not_drawn_embedding_outlier_by_subtype": subtype_counts(epithelial[has_position & is_outlier]),
        "other_nuclei_drawn_by_sample_and_class": count_table(
            drawn[drawn["cell_type"] != "epithelial_cell"], "timepoint", "cell_type"
        ),
        "display_transform": {s: {"mirror_x": m, "rotate_ccw_deg": r} for s, (m, r) in pc.DISPLAY_TRANSFORM.items()},
        "crop": "none: whole sections, as the S26 grid",
        "map_side_um": round(2 * half_width_um, 1),
        "scale_bar_um": pc.SCALE_BAR_UM,
        "marker_area_pt2": s26.MARKER_AREA,
        "draw_order_seed": pc.SEED,
        "background_colors": {c: colors[c] for c in sorted(colors)},
        "tint_toward_white": TINT_TOWARD_WHITE if background == "decfib" else None,
    }


# UMAP (S26C) ----

def umap_nuclei(cells: pd.DataFrame, embedding: str, lineage_table: Path | None) -> tuple[pd.DataFrame, dict]:
    """Nuclei drawn on one UMAP option, with columns umap_x/umap_y, in a seeded random order."""
    is_outlier = cells["embedding_outlier"].astype(bool)
    is_epithelial = cells["cell_type"] == "epithelial_cell"
    if embedding == "lineage_fit":
        lineage = pd.read_parquet(lineage_table)
        nuclei = cells[is_epithelial & ~is_outlier].join(lineage, how="left")
        if nuclei["lineage_umap_1"].isna().any():
            raise ValueError("released epithelial nuclei without lineage-fit coordinates")
    else:
        nuclei = cells[~is_outlier].copy()
    x_col, y_col = EMBEDDINGS[embedding]["coords"]
    nuclei["umap_x"] = nuclei[x_col].astype(float)
    nuclei["umap_y"] = nuclei[y_col].astype(float)

    counts = {
        "embedding_outliers_hidden_epithelial": int((is_epithelial & is_outlier).sum()),
        "embedding_outliers_hidden_all": int(is_outlier.sum()) if embedding == "atlas" else None,
    }
    order = np.random.default_rng(pc.SEED).permutation(len(nuclei))
    return nuclei.iloc[order], counts


def umap_limits(nuclei: pd.DataFrame) -> tuple[tuple, tuple]:
    """Axis limits around the epithelial nuclei: the whole lineage fit, or the atlas zoomed to the epithelium."""
    epithelial = nuclei[nuclei["cell_type"] == "epithelial_cell"]
    x_lo, x_hi = epithelial["umap_x"].min(), epithelial["umap_x"].max()
    y_lo, y_hi = epithelial["umap_y"].min(), epithelial["umap_y"].max()
    pad = UMAP_PAD * max(x_hi - x_lo, y_hi - y_lo)
    return (x_lo - pad, x_hi + pad), (y_lo - pad, y_hi + pad)


def umap_width_in(x_lim: tuple, y_lim: tuple) -> float:
    return UMAP_H_IN * (x_lim[1] - x_lim[0]) / (y_lim[1] - y_lim[0])


def draw_umap(ax, nuclei: pd.DataFrame, x_lim: tuple, y_lim: tuple, colors_by_class: dict) -> dict:
    # Colour = epithelial subtype, the colours of the maps. Any other nucleus inside the window is drawn
    # first in its map colour, so the epithelial island reads against its few neighbours from other
    # classes. UMAP axes carry no units: no ticks or frame, as the manuscript's UMAPs.
    is_epithelial = nuclei["cell_type"] == "epithelial_cell"
    in_window = nuclei["umap_x"].between(*x_lim) & nuclei["umap_y"].between(*y_lim)
    others = nuclei[~is_epithelial & in_window]
    epithelium = nuclei[is_epithelial]
    layers = (
        (others, others["cell_type"].map(colors_by_class)),
        (epithelium, epithelium["cell_subtype"].map(SUBTYPE_COLORS)),
    )
    for subset, colors in layers:
        ax.scatter(
            subset["umap_x"],
            subset["umap_y"],
            c=colors.to_list(),
            s=UMAP_MARKER_AREA,
            linewidths=0,
            rasterized=True,
        )
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_aspect("equal")
    ax.axis("off")
    return {
        "epithelial_drawn_by_subtype": subtype_counts(epithelium),
        "other_nuclei_in_window_by_class": others["cell_type"].value_counts().to_dict(),
        "other_nuclei_outside_window": int((~is_epithelial & ~in_window).sum()),
    }


def draw_locator(fig, ax, nuclei: pd.DataFrame, x_lim: tuple, y_lim: tuple, colors_by_class: dict) -> dict:
    """Miniature of the whole atlas with the zoom window boxed, in an empty corner of the zoomed UMAP."""
    fig_w, fig_h = fig.get_size_inches()
    box = ax.get_position()
    ax_w_in, ax_h_in = box.width * fig_w, box.height * fig_h
    units_per_in = (x_lim[1] - x_lim[0]) / ax_w_in
    footprint_w = LOCATOR_IN + LOCATOR_INSET_IN
    footprint_h = LOCATOR_IN + LOCATOR_LABEL_IN + LOCATOR_INSET_IN

    # Footprint of each candidate corner in data units; the first one holding no epithelial nucleus wins.
    epithelial_xy = nuclei.loc[nuclei["cell_type"] == "epithelial_cell", ["umap_x", "umap_y"]].to_numpy()
    corner = None
    for candidate in LOCATOR_CORNERS:
        vertical, horizontal = candidate.split()
        x0 = x_lim[1] - footprint_w * units_per_in if horizontal == "right" else x_lim[0]
        y0 = y_lim[0] if vertical == "lower" else y_lim[1] - footprint_h * units_per_in
        covered = (
            (epithelial_xy[:, 0] >= x0)
            & (epithelial_xy[:, 0] <= x0 + footprint_w * units_per_in)
            & (epithelial_xy[:, 1] >= y0)
            & (epithelial_xy[:, 1] <= y0 + footprint_h * units_per_in)
        )
        if not covered.any():
            corner = (candidate, horizontal, vertical)
            break
    if corner is None:
        raise ValueError("no corner of the zoomed UMAP is free of epithelial nuclei for the locator")
    candidate, horizontal, vertical = corner

    left_in = ax_w_in - LOCATOR_IN - LOCATOR_INSET_IN if horizontal == "right" else LOCATOR_INSET_IN
    bottom_in = LOCATOR_INSET_IN if vertical == "lower" else ax_h_in - LOCATOR_IN - LOCATOR_LABEL_IN - LOCATOR_INSET_IN
    locator = fig.add_axes([
        box.x0 + left_in / fig_w,
        box.y0 + bottom_in / fig_h,
        LOCATOR_IN / fig_w,
        LOCATOR_IN / fig_h,
    ])

    # Every drawn nucleus in its map colour, other nuclei under the epithelium as on the maps; the box is
    # the zoom window of the large panel.
    is_epithelial = nuclei["cell_type"] == "epithelial_cell"
    layers = (
        (nuclei[~is_epithelial], nuclei.loc[~is_epithelial, "cell_type"].map(colors_by_class)),
        (nuclei[is_epithelial], nuclei.loc[is_epithelial, "cell_subtype"].map(SUBTYPE_COLORS)),
    )
    for subset, colors in layers:
        locator.scatter(
            subset["umap_x"],
            subset["umap_y"],
            c=colors.to_list(),
            s=LOCATOR_MARKER_AREA,
            linewidths=0,
            rasterized=True,
        )
    locator.add_patch(Rectangle(
        (x_lim[0], y_lim[0]),
        x_lim[1] - x_lim[0],
        y_lim[1] - y_lim[0],
        fill=False,
        edgecolor=UMAP_TEXT,
        linewidth=0.5,
    ))
    locator.set_aspect("equal", adjustable="datalim")
    locator.set_xticks([])
    locator.set_yticks([])
    for spine in locator.spines.values():
        spine.set_color(LOCATOR_FRAME)
        spine.set_linewidth(0.4)
    fig.text(
        box.x0 + (left_in + LOCATOR_IN / 2) / fig_w,
        box.y0 + (bottom_in + LOCATOR_IN + LOCATOR_LABEL_GAP_IN) / fig_h,
        LOCATOR_LABEL,
        ha="center",
        va="bottom",
        fontsize=LOCATOR_LABEL_PT,
        color=UMAP_TEXT,
    )
    return {"corner": candidate, "side_in": LOCATOR_IN, "nuclei": int(len(nuclei)), "label": LOCATOR_LABEL}


def add_umap_title(fig, x_left: float, umap_top: float, embedding: str) -> None:
    """Name of the embedding above the UMAP's left edge, bold as the S26 column headers."""
    fig_h = fig.get_size_inches()[1]
    fig.text(
        x_left,
        umap_top + UMAP_TITLE_GAP_IN / fig_h,
        EMBEDDINGS[embedding]["title"],
        ha="left",
        va="bottom",
        color=UMAP_TEXT,
        fontweight="bold",
    )


# Layouts ----

def figure_spatial(maps: dict, background: str, uns: dict):
    size_in = spatial_block_size_in()
    fig = plt.figure(figsize=size_in, facecolor="black")
    draw_spatial_block(fig, size_in, maps["drawn"], maps["half_width_um"], background, uns)
    return fig, {}


def figure_umap(umap: dict, embedding: str, uns: dict):
    """Standalone S26C: title, UMAP, legend of the subtypes (and of the other nuclei on the atlas)."""
    nuclei = umap["nuclei"]
    x_lim, y_lim = umap_limits(nuclei)
    umap_w = umap_width_in(x_lim, y_lim)

    # Other nuclei appear only in the atlas option's locator; the lineage fit holds epithelium alone.
    subtypes, others = legend_entries(uns, "grey")
    entries = subtypes + (others if embedding == "atlas" else [])
    fig_w = max(umap_w + 2 * s26.MARGIN_IN, UMAP_MIN_W_IN)
    fig_h = s26.HEADER_IN + UMAP_H_IN + LEGEND_GAP_IN + len(entries) * LEGEND_ROW_IN + s26.MARGIN_IN
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=UMAP_BACKGROUND)

    left = (fig_w - umap_w) / 2 / fig_w
    bottom = 1 - (s26.HEADER_IN + UMAP_H_IN) / fig_h
    ax = fig.add_axes([left, bottom, umap_w / fig_w, UMAP_H_IN / fig_h])
    counts = draw_umap(ax, nuclei, x_lim, y_lim, background_colors("grey"))
    if embedding == "atlas":
        counts["locator"] = draw_locator(fig, ax, nuclei, x_lim, y_lim, background_colors("grey"))
    add_umap_title(fig, left, bottom + UMAP_H_IN / fig_h, embedding)
    add_legend(fig, entries, 0.5, bottom - LEGEND_GAP_IN / fig_h, UMAP_TEXT)
    counts["limits"] = {"x": [round(v, 3) for v in x_lim], "y": [round(v, 3) for v in y_lim]}
    return fig, counts


def figure_both(maps: dict, background: str, umap: dict, embedding: str, uns: dict):
    """S26B and S26C side by side: the maps on black, the UMAP on white, one legend under the maps."""
    block_w, block_h = spatial_block_size_in()
    nuclei = umap["nuclei"]
    x_lim, y_lim = umap_limits(nuclei)
    umap_w = umap_width_in(x_lim, y_lim)
    fig_w = block_w + UMAP_GAP_IN + umap_w + s26.MARGIN_IN
    fig = plt.figure(figsize=(fig_w, block_h), facecolor=UMAP_BACKGROUND)
    draw_spatial_block(fig, (fig_w, block_h), maps["drawn"], maps["half_width_um"], background, uns)

    # The UMAP shares the maps' top and height, so a colour is read across from map to UMAP at one glance.
    left = (block_w + UMAP_GAP_IN) / fig_w
    bottom = 1 - (s26.HEADER_IN + UMAP_H_IN) / block_h
    ax = fig.add_axes([left, bottom, umap_w / fig_w, UMAP_H_IN / block_h])
    counts = draw_umap(ax, nuclei, x_lim, y_lim, background_colors(background))
    if embedding == "atlas":
        counts["locator"] = draw_locator(fig, ax, nuclei, x_lim, y_lim, background_colors(background))
    add_umap_title(fig, left, bottom + UMAP_H_IN / block_h, embedding)
    counts["limits"] = {"x": [round(v, 3) for v in x_lim], "y": [round(v, 3) for v in y_lim]}
    return fig, counts


# Run ----

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--uns", required=True, type=Path)
    parser.add_argument("--extract-provenance", required=True, type=Path)
    parser.add_argument("--release", required=True)
    parser.add_argument("--layout", required=True, choices=["spatial", "umap", "both"])
    parser.add_argument("--background", choices=BACKGROUNDS, help="map background (spatial and both)")
    parser.add_argument("--embedding", choices=list(EMBEDDINGS), help="UMAP option (umap and both)")
    parser.add_argument("--lineage-table", type=Path, help="extract_epithelial_lineage_umap.py table")
    parser.add_argument("--lineage-provenance", type=Path)
    parser.add_argument("--out-stem", required=True, type=Path)
    args = parser.parse_args()

    draws_maps = args.layout in ("spatial", "both")
    draws_umap = args.layout in ("umap", "both")
    needs_lineage = draws_umap and args.embedding == "lineage_fit"
    if draws_maps and args.background is None:
        parser.error("--background is required for the spatial and both layouts")
    if draws_umap and args.embedding is None:
        parser.error("--embedding is required for the umap and both layouts")
    if needs_lineage and (args.lineage_table is None or args.lineage_provenance is None):
        parser.error("--lineage-table and --lineage-provenance are required for the lineage_fit UMAP")

    panel = {
        "spatial": f"S26B_{args.background}",
        "umap": f"S26C_{args.embedding}",
        "both": f"S26BC_{args.background}_{args.embedding}",
    }[args.layout]

    pc.set_style()

    # Load: one point per nucleus. The maps keep the nuclei that passed spatial QC, in the S26
    # orientation and spatial scale; embedding outliers are hidden on maps and UMAPs (release rule) ----
    cells = load_cells(args.cells)
    uns = json.loads(args.uns.read_text())
    maps = None
    if draws_maps:
        drawn, spatial_counts = pc.spatial_cells(cells)
        maps = {"drawn": drawn, "counts": spatial_counts, "half_width_um": pc.shared_half_width_um(drawn)}
    umap = None
    if draws_umap:
        nuclei, umap_counts = umap_nuclei(cells, args.embedding, args.lineage_table)
        umap = {"nuclei": nuclei, "counts": umap_counts}

    # Figure ----
    if args.layout == "spatial":
        fig, drawn_counts = figure_spatial(maps, args.background, uns)
        facecolor = "black"
    elif args.layout == "umap":
        fig, drawn_counts = figure_umap(umap, args.embedding, uns)
        facecolor = UMAP_BACKGROUND
    else:
        fig, drawn_counts = figure_both(maps, args.background, umap, args.embedding, uns)
        facecolor = UMAP_BACKGROUND
    outputs = pc.save_panel(fig, args.out_stem, facecolor=facecolor)

    # Provenance ----
    scripts = [Path(__file__), Path(pc.__file__), Path(s26.__file__)]
    provenance = pc.base_provenance(panel, args.release, args.cells, args.extract_provenance, scripts)
    provenance["layout"] = args.layout
    provenance["uns_export"] = {"path": str(args.uns.resolve()), "sha256": pc.sha256(args.uns)}
    provenance["layer_rule"] = (
        "annotation vocabulary v7 (human_epithelial_split): a placed epithelial nucleus votes inner when decidual "
        "cells outnumber endometrial fibroblasts among its 20 nearest placed nuclei of the same section, outer for "
        "the reverse; each leiden_scVI_3.0 cluster of the epithelial lineage fit takes its majority"
    )
    provenance["colors"]["subtypes"] = SUBTYPE_COLORS
    provenance["colors"]["subtype_source"] = (
        "new (the manuscript shows no epithelial subtypes): Okabe-Ito orange, blue, bluish green"
    )
    if draws_maps:
        provenance["spatial"] = spatial_provenance(
            cells,
            maps["drawn"],
            maps["counts"],
            maps["half_width_um"],
            args.background,
        )
    if draws_umap:
        provenance["umap"] = {
            "embedding": args.embedding,
            "description": EMBEDDINGS[args.embedding]["description"],
            "marker_area_pt2": UMAP_MARKER_AREA,
            "draw_order_seed": pc.SEED,
            **umap["counts"],
            **drawn_counts,
        }
    if needs_lineage:
        lineage = json.loads(args.lineage_provenance.read_text())
        provenance["umap"]["lineage_table"] = {
            "path": str(args.lineage_table.resolve()),
            "sha256": pc.sha256(args.lineage_table),
        }
        provenance["umap"]["lineage_fit"] = {
            key: lineage["fit"][key] for key in ("path", "bytes", "sha256", "sha256_matches_release_record")
        }
        provenance["umap"]["fit_nuclei_not_released_by_cluster"] = lineage["fit_nuclei_not_released_by_cluster"]
        provenance["umap"]["clusters_by_subtype"] = lineage["clusters_by_subtype"]
    provenance["outputs"] = outputs
    pc.write_json(Path(f"{args.out_stem}.provenance.json"), provenance)


if __name__ == "__main__":
    main()
