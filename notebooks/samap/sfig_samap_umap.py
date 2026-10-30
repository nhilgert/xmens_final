"""Joint SAMap UMAP of panels a-c of the SAMap supplementary figure: the nuclei drawn, their order, and the frame.

Imported by the notebooks of panels a, b and c, which add their own color layer and key. load_umap_cells() reads
the three tables it is given; nothing is written.
"""
import matplotlib as mpl
import pandas as pd
from matplotlib.lines import Line2D

from sfig_samap_style import PSEUDOTIME_CMAP, PSEUDOTIME_NORM, new_panel, panel_title

DOT = 0.22                 # pt^2: about 90,000 nuclei share each 57 mm panel
OTHER_SPECIES = "#d4d4d4"  # nuclei of the species whose values a panel does not show
AXIS_GREY = "#555555"      # UMAP axis arrows and their labels


def load_umap_cells(joint_path, mouse_cells_path, human_cells_path, draw_outliers: bool) -> dict:
    """Nuclei drawn on the UMAPs, each species in a seeded random order, and the frame around them.

    Nuclei flagged embedding_outlier in the release are drawn when draw_outliers is True (S30 shows every nucleus;
    Fig. 6C and 6E leave them out). The random order keeps early or late nuclei from always sitting on top; the seeds
    and the joint object's nucleus order fix it. Returns {"mouse", "human"} DataFrames (umap1, umap2 and the
    values), "xlim", "ylim", "aspect", and the flagged nuclei and those left out, per species.
    """
    umap = (
        pd.read_parquet(joint_path, columns=["species", "X_umap__0", "X_umap__1"])
        .rename(columns={"X_umap__0": "umap1", "X_umap__1": "umap2"})
    )
    umap["species"] = umap["species"].astype(str)
    if not umap.index.is_unique:
        raise ValueError(f"{joint_path} has duplicated nucleus names")

    mouse = pd.read_csv(mouse_cells_path, index_col="cell")
    human = pd.read_csv(human_cells_path, index_col="cell")

    # Stored global rank (average rank / N over the non-apoptotic mouse nuclei; NaN for apoptotic ones).
    mouse_values = mouse[["cell_type", "cell_subtype"]].assign(pt_rank=mouse["pt_harmonized_rank"])

    # Human pseudotime (pt_human) and the transfer of record (pt_ranknorm), one row per human nucleus.
    human_values = human[["pt_transferred", "pt_harmonized_rank", "cell_type"]].rename(
        columns={"pt_transferred": "pt_ranknorm", "pt_harmonized_rank": "pt_human"})

    flagged_human = set(human.index[human["embedding_outlier"].astype(bool)])
    flagged_mouse = set(mouse.index[mouse["embedding_outlier"].astype(bool)])
    is_outlier = (
        ((umap.species == "hs") & umap.index.isin(flagged_human))
        | ((umap.species == "mm") & umap.index.isin(flagged_mouse))
    )
    drawn = umap if draw_outliers else umap[~is_outlier]
    flagged = {
        "mouse": int((umap.species[is_outlier] == "mm").sum()),
        "human": int((umap.species[is_outlier] == "hs").sum()),
    }

    mouse_drawn = drawn[drawn.species == "mm"].join(mouse_values[["pt_rank", "cell_type", "cell_subtype"]], how="left")
    human_drawn = drawn[drawn.species == "hs"].join(human_values[["pt_ranknorm", "pt_human", "cell_type"]], how="left")

    # Frame: the drawn nuclei plus 2% of their range on each side.
    xy = drawn[["umap1", "umap2"]].values
    pad = 0.02 * (xy.max(0) - xy.min(0))
    xlim = (xy[:, 0].min() - pad[0], xy[:, 0].max() + pad[0])
    ylim = (xy[:, 1].min() - pad[1], xy[:, 1].max() + pad[1])

    return {
        "mouse": mouse_drawn.sample(frac=1, random_state=1),
        "human": human_drawn.sample(frac=1, random_state=3),
        "xlim": xlim,
        "ylim": ylim,
        "aspect": (xlim[1] - xlim[0]) / (ylim[1] - ylim[0]),
        "outliers_flagged": flagged,
        "outliers_left_out": {"mouse": 0, "human": 0} if draw_outliers else flagged,
    }


def umap_panel(letter: str, title: str, cells: dict, draw, key=None, colorbar_title=None):
    """One UMAP panel: equal-aspect axes without a box, short axis arrows, and an optional colorbar and key.

    draw(ax) adds the panel's own scatter layers. key is a list of (color, label) for the grey nuclei; the colorbar
    is the fixed 0-1 pseudotime scale.
    """
    fig, width_mm, height_mm = new_panel(letter)
    top_mm, bottom_mm, left_mm = 10.0, 4.0, 6.0
    axes_h = (height_mm - top_mm - bottom_mm) / height_mm
    axes_w = axes_h * height_mm * cells["aspect"] / width_mm
    ax = fig.add_axes([left_mm / width_mm, bottom_mm / height_mm, axes_w, axes_h])
    xlim, ylim = cells["xlim"], cells["ylim"]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_title(fig, width_mm, height_mm, title)

    draw(ax)

    # UMAP coordinates have no units, so two short arrows replace a box and ticks.
    dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]
    length = 0.16 * dx
    arrow = dict(
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color=AXIS_GREY, mutation_scale=6,
                        shrinkA=0, shrinkB=0, capstyle="projecting", joinstyle="miter"),
        annotation_clip=False,
    )
    ax.annotate("", xy=(xlim[0] + length, ylim[0]), xytext=(xlim[0], ylim[0]), **arrow)
    ax.annotate("", xy=(xlim[0], ylim[0] + length), xytext=(xlim[0], ylim[0]), **arrow)
    ax.text(xlim[0] + 0.02 * dx, ylim[0] - 0.012 * dy, "UMAP 1", ha="left", va="top", fontsize=7, color=AXIS_GREY)
    ax.text(xlim[0] - 0.03 * dx, ylim[0] + 0.02 * dy, "UMAP 2", ha="right", va="bottom", rotation=90,
            fontsize=7, color=AXIS_GREY)

    right = (left_mm + axes_w * width_mm + 3.0) / width_mm
    if colorbar_title:
        cax = fig.add_axes([right, 0.42, 2.2 / width_mm, 0.30])
        bar = fig.colorbar(mpl.cm.ScalarMappable(norm=PSEUDOTIME_NORM, cmap=PSEUDOTIME_CMAP), cax=cax,
                           ticks=[0, 0.5, 1])
        bar.ax.set_yticklabels(["0 early", "0.5", "1 late"])
        bar.outline.set_visible(False)
        bar.ax.tick_params(length=1.5, pad=1.5)
        cax.set_title(colorbar_title, fontsize=8, loc="left", pad=8)
    if key:
        top = 0.33 if colorbar_title else 1 - top_mm / height_mm
        kax = fig.add_axes([right - 1.2 / width_mm, bottom_mm / height_mm, 1 - right, top - bottom_mm / height_mm])
        kax.axis("off")
        kax.legend(
            handles=[Line2D([], [], ls="none", marker="o", ms=3.2, color=color, label=label) for color, label in key],
            loc="upper left",
            frameon=False,
            fontsize=8,
            handletextpad=0.2,
            labelspacing=0.6,
            borderaxespad=0,
            handlelength=0.9,
        )
    return fig
