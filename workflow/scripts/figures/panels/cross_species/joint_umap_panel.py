"""Fig. 6A and 6B: human and X-Mens eFib/dec nuclei on the joint SAMap UMAP, by species (6A) or cell type (6B).

Called by rule panels_cross_species_joint_umap (workflow/rules/91_panels_cross_species.smk).
Reads:  the joint SAMap table of Fig. 6C (species and joint UMAP of every nucleus; rule samap_figure_tables) and the
        release fibroblast/decidual objects of both species (obs only: cell_type, embedding_outlier).
Writes: the panel as PDF, SVG and 300-dpi PNG on black, and a provenance JSON with inputs, filters and counts.

Panel size, square frame, dot size, draw-order seed, axis arrows and key position copy the Fig. 6C notebook
(notebooks/samap/fig6c_umap_transferred_pseudotime.ipynb), so 6A, 6B and 6C overlay exactly.

usage: joint_umap_panel.py --panel {6A,6B} --joint T --mouse-obj H5AD --human-obj H5AD
                           --pdf F --svg F --png F --provenance F
"""
import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

import cross_species_panels as common


# Constants ----

# Geometry of the Fig. 6C notebook: a 90 mm black square, the UMAP axes inset 6% at the sides and bottom.
PANEL_MM = 90
AXES_RECT = [0.06, 0.06, 0.88, 0.84]
FRAME_MARGIN = 1.06  # square frame 6% wider than the larger extent of the joint embedding
UMAP_DOT = 0.5  # pt^2
ORDER_SEED = 3  # seeded random draw order, so neither species nor cell type always sits on top
AXIS_ARROWS = "#bbbbbb"

# Font sizes (pt) of 6C: the title, then key and axis-arrow labels one step down.
TITLE_SIZE = 9
KEY_SIZE = 7

# One entry per panel: the column that sets the color, its colors and key labels in key order, and the title.
PANELS = {
    "6A": {
        "column": "species",
        "colors": common.SPECIES_COLORS,
        "labels": common.SPECIES_LABELS,
        "title": "Human and mouse cells co-embedded by SAMap",
    },
    "6B": {
        "column": "cell_type",
        "colors": common.CELL_TYPE_COLORS,
        "labels": common.CELL_TYPE_LABELS,
        "title": "Fibroblasts and decidual cells of both species",
    },
}

LOG = logging.getLogger("joint_umap_panel")


# Load ----

def load_cells(joint_path: Path, mouse_path: Path, human_path: Path) -> pd.DataFrame:
    """One row per nucleus of the joint table: species, joint UMAP, release cell_type and embedding_outlier."""
    joint = pd.read_parquet(joint_path, columns=["species", "X_umap__0", "X_umap__1"])
    labels = pd.concat([
        common.load_release_obs(mouse_path, ["cell_type", "embedding_outlier"]),
        common.load_release_obs(human_path, ["cell_type", "embedding_outlier"]),
    ])

    # Every nucleus of the joint embedding must carry a release label, and no label may be ambiguous.
    if not labels.index.is_unique:
        raise ValueError("cell names repeat across the mouse and human release objects")
    missing = joint.index.difference(labels.index)
    if len(missing):
        raise ValueError(f"{len(missing)} joint nuclei are not in the release objects, e.g. {list(missing[:3])}")

    cells = joint.join(labels)
    cells["species"] = cells["species"].astype(str)
    cells["cell_type"] = cells["cell_type"].astype(str)
    cells["embedding_outlier"] = cells["embedding_outlier"].astype(bool)
    unknown = (set(cells["cell_type"]) - set(common.CELL_TYPE_COLORS)) | (
        set(cells["species"]) - set(common.SPECIES_COLORS)
    )
    if unknown:
        raise ValueError(f"cell types or species without a color: {sorted(unknown)}")
    return cells


def count_by(cells: pd.DataFrame, column: str) -> dict:
    return {str(key): int(value) for key, value in cells[column].value_counts().sort_index().items()}


# Figure ----

def square_limits(cells: pd.DataFrame) -> tuple:
    """Square frame around the whole joint embedding (both species, outliers included), as in Fig. 6C."""
    x, y = cells["X_umap__0"], cells["X_umap__1"]
    side = max(x.max() - x.min(), y.max() - y.min()) * FRAME_MARGIN
    cx, cy = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2
    return (cx - side / 2, cx + side / 2), (cy - side / 2, cy + side / 2), side


def umap_figure(title: str, xlim, ylim, side):
    """Black square panel with a boxless equal-aspect UMAP axes and short axis arrows."""
    fig = plt.figure(figsize=(PANEL_MM / 25.4, PANEL_MM / 25.4), facecolor=common.BACKGROUND)
    ax = fig.add_axes(AXES_RECT)
    ax.set_facecolor(common.BACKGROUND)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(0.06, 0.965, title, ha="left", va="top", fontsize=TITLE_SIZE, color=common.TEXT)

    # UMAP coordinates have no units, so two short arrows replace a box and ticks.
    length = 0.14 * side
    arrow = dict(
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color=AXIS_ARROWS, mutation_scale=6, shrinkA=0, shrinkB=0),
        annotation_clip=False,
    )
    ax.annotate("", xy=(xlim[0] + length, ylim[0]), xytext=(xlim[0], ylim[0]), **arrow)
    ax.annotate("", xy=(xlim[0], ylim[0] + length), xytext=(xlim[0], ylim[0]), **arrow)
    ax.text(
        xlim[0] + 0.015 * side,
        ylim[0] - 0.012 * side,
        "UMAP 1",
        ha="left",
        va="top",
        fontsize=KEY_SIZE,
        color=AXIS_ARROWS,
    )
    ax.text(
        xlim[0] - 0.015 * side,
        ylim[0] + 0.015 * side,
        "UMAP 2",
        ha="right",
        va="bottom",
        rotation=90,
        fontsize=KEY_SIZE,
        color=AXIS_ARROWS,
    )
    return fig, ax


def key_bottom_right(fig, entries) -> None:
    """Color key in the empty bottom-right corner of the embedding, where Fig. 6C keeps its key."""
    handles = [
        Line2D([], [], ls="", marker="o", ms=4, mfc=color, mec="none", label=label)
        for color, label in entries
    ]
    kax = fig.add_axes([0.56, 0.0, 0.44, 0.16])
    kax.axis("off")
    kax.legend(
        handles=handles,
        loc="lower left",
        fontsize=KEY_SIZE,
        handletextpad=0.3,
        labelspacing=0.4,
        borderaxespad=0,
        labelcolor=common.TEXT,
        frameon=False,
    )


def draw_panel(cells: pd.DataFrame, spec: dict, xlim, ylim, side):
    fig, ax = umap_figure(spec["title"], xlim, ylim, side)

    # One dot per nucleus, color = species or cell type. All nuclei go down in one shuffled layer: drawing one
    # group after the other would let the last group hide the first wherever they overlap.
    shuffled = cells.sample(frac=1, random_state=ORDER_SEED)
    ax.scatter(
        shuffled["X_umap__0"],
        shuffled["X_umap__1"],
        s=UMAP_DOT,
        c=shuffled[spec["column"]].map(spec["colors"]).to_list(),
        linewidths=0,
        rasterized=True,
    )

    key_bottom_right(fig, [(spec["colors"][key], label) for key, label in spec["labels"].items()])
    return fig


# Run ----

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--panel", required=True, choices=sorted(PANELS))
    parser.add_argument("--joint", required=True, type=Path)
    parser.add_argument("--mouse-obj", required=True, type=Path)
    parser.add_argument("--human-obj", required=True, type=Path)
    parser.add_argument("--pdf", required=True, type=Path)
    parser.add_argument("--svg", required=True, type=Path)
    parser.add_argument("--png", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    # fontTools logs every subset glyph at INFO when the PDF embeds the font.
    logging.getLogger("fontTools").setLevel(logging.WARNING)

    spec = PANELS[args.panel]
    cells = load_cells(args.joint, args.mouse_obj, args.human_obj)

    # Every frame shares one square fitted to the whole joint embedding, outliers included (as in 6C).
    xlim, ylim, side = square_limits(cells)

    # Nuclei flagged embedding_outlier by the release stay in every statistic but are not drawn.
    is_outlier = cells["embedding_outlier"]
    drawn = cells[~is_outlier]
    LOG.info(
        "%s: %d joint nuclei; %d embedding outliers left out (%s); %d drawn",
        args.panel,
        len(cells),
        int(is_outlier.sum()),
        count_by(cells[is_outlier], "species"),
        len(drawn),
    )

    font = common.set_font()
    fig = draw_panel(drawn, spec, xlim, ylim, side)
    paths = {"pdf": args.pdf, "svg": args.svg, "png": args.png}
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    common.save_panel(fig, paths)
    plt.close(fig)
    LOG.info("Wrote %s (font %s)", [str(path) for path in paths.values()], font)

    inputs = {"joint": args.joint, "mouse_object": args.mouse_obj, "human_object": args.human_obj}
    record = common.base_record(Path(__file__), inputs, paths, font)
    record["panel"] = args.panel
    record["color_by"] = spec["column"]
    record["colors"] = {spec["labels"][key]: spec["colors"][key] for key in spec["labels"]}
    record["filters"] = [
        "all nuclei of the joint SAMap table, labelled from the release objects",
        "nuclei with obs embedding_outlier == True are not drawn (frame still fitted to them, as in 6C)",
    ]
    record["counts"] = {
        "joint_nuclei": len(cells),
        "joint_by_species": count_by(cells, "species"),
        "embedding_outliers_omitted": int(is_outlier.sum()),
        "embedding_outliers_by_species": count_by(cells[is_outlier], "species"),
        "drawn": len(drawn),
        "drawn_by_species": count_by(drawn, "species"),
        "drawn_by_cell_type": count_by(drawn, "cell_type"),
        "drawn_by_species_and_cell_type": {
            f"{species}/{cell_type}": int(n)
            for (species, cell_type), n in drawn.groupby(["species", "cell_type"]).size().items()
        },
    }
    record["frame"] = {
        "xlim": [float(value) for value in xlim],
        "ylim": [float(value) for value in ylim],
        "side": float(side),
        "margin": FRAME_MARGIN,
    }
    record["drawing"] = {
        "panel_mm": PANEL_MM,
        "dot_pt2": UMAP_DOT,
        "order_seed": ORDER_SEED,
        "raster_dpi": common.RASTER_DPI,
        "png_dpi": common.PNG_DPI,
    }
    common.write_json(args.provenance, record)
    LOG.info("Wrote %s", args.provenance)


if __name__ == "__main__":
    main()
