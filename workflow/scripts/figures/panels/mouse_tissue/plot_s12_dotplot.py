"""Draw the class marker dot plot of Fig. S12B: detection fraction and scaled mean expression.

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for S12B.
Reads marker_stats.csv and extract_provenance.json from --tables-dir.
Writes <out-prefix>.pdf, .svg, .png and .provenance.json.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize

import mouse_tissue_style as style


# Half of a 180 mm page, as in the manuscript's S12 layout.
FIG_SIZE_IN = (3.7, 1.95)

# Dot area (pt^2) of a gene detected in every nucleus of a class; area scales with the fraction.
# The manuscript's largest dots are about 7 pt across.
MAX_DOT_AREA_PT2 = 48.0
LEGEND_FRACTIONS = [100, 50]

# White-to-dark-red ramp of the manuscript for the per-gene scaled mean (0-1).
CMAP = colormaps["Reds"]
DOT_EDGE = "#7f7f7f"


def dotplot_genes(classes: list) -> list:
    """The manuscript's markers, plus the perivascular marker when that class is present (v9)."""
    genes = list(style.DOTPLOT_GENES)
    if "perivascular_cell" in classes:
        at = genes.index(style.PERIVASCULAR_GENES_AFTER) + 1
        genes[at:at] = style.PERIVASCULAR_GENES
    return genes


def main() -> None:
    args = style.panel_arguments(__doc__).parse_args()
    style.apply_panel_style()

    # Load and reshape ----

    stats = pd.read_csv(args.tables_dir / "marker_stats.csv")
    classes = [c for c in style.SUPP_CLASS_ORDER if c in set(stats["cell_type"])]
    unknown = set(stats["cell_type"]) - set(style.SUPP_CLASS_ORDER)
    if unknown:
        raise ValueError(f"cell types without a dot-plot row: {sorted(unknown)}")
    genes = dotplot_genes(classes)
    frac = stats.pivot(index="cell_type", columns="gene", values="frac_detected").loc[classes, genes]
    mean = stats.pivot(index="cell_type", columns="gene", values="mean_lognorm").loc[classes, genes]

    # "Normalised mean expression in cluster": each gene's class means rescaled to 0-1 across
    # the classes shown, so every column spans the full colour ramp.
    scaled = (mean - mean.min()) / (mean.max() - mean.min())
    logging.info("Dot plot: %d classes x %d genes", len(classes), len(genes))

    # Figure ----

    fig = plt.figure(figsize=FIG_SIZE_IN)
    width_share = len(genes) / (len(genes) + 6.5)
    ax = fig.add_axes([0.1, 0.27, 0.86 * width_share, 0.7])

    # One dot per class x gene: area = fraction of nuclei detecting the gene, colour = scaled mean.
    x, y = np.meshgrid(np.arange(len(genes)), np.arange(len(classes)))
    ax.scatter(
        x.ravel(),
        y.ravel(),
        s=MAX_DOT_AREA_PT2 * frac.to_numpy().ravel(),
        c=scaled.to_numpy().ravel(),
        cmap=CMAP,
        vmin=0,
        vmax=1,
        edgecolors=DOT_EDGE,
        linewidths=0.3,
    )

    # Boxed frame as in the manuscript; classes top to bottom in legend order.
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(classes) - 0.5, -0.5)
    ax.set_xticks(np.arange(len(genes)), labels=genes, rotation=90)
    ax.set_yticks(np.arange(len(classes)), labels=[style.DOTPLOT_CLASS_LABELS[c] for c in classes])

    # Legends ----

    # Size key: dots for 100% and 50%, drawn at the same scale as the plot.
    legend_x = 0.1 + 0.86 * width_share + 0.04
    key = fig.add_axes([legend_x, 0.62, 0.95 - legend_x, 0.35])
    key.set_axis_off()
    key.set_xlim(0, 1)
    key.set_ylim(0, 1)
    key.text(0.5, 1.0, "fraction of cells\nin cluster (%)", ha="center", va="top", fontsize=style.FONT_SIZE_SMALL)
    for i, pct in enumerate(LEGEND_FRACTIONS):
        y_dot = 0.42 - 0.3 * i
        key.scatter([0.4], [y_dot], s=MAX_DOT_AREA_PT2 * pct / 100, color=CMAP(1.0), edgecolors=DOT_EDGE, linewidths=0.3)
        key.text(0.52, y_dot, f"{pct}", ha="left", va="center", fontsize=style.FONT_SIZE_SMALL)

    # Colour key: the 0-1 ramp of the scaled mean.
    fig.text(
        (legend_x + 0.95) / 2,
        0.5,
        "normalized mean\nexpression in cluster",
        ha="center",
        va="top",
        fontsize=style.FONT_SIZE_SMALL,
    )
    bar = fig.add_axes([legend_x + 0.03, 0.31, 0.95 - legend_x - 0.06, 0.05])
    colorbar = ColorbarBase(bar, cmap=CMAP, norm=Normalize(0, 1), orientation="horizontal", ticks=[])
    colorbar.outline.set_visible(False)
    bar.text(-0.06, 0.5, "0", transform=bar.transAxes, ha="right", va="center", fontsize=style.FONT_SIZE_SMALL)
    bar.text(1.06, 0.5, "1", transform=bar.transAxes, ha="left", va="center", fontsize=style.FONT_SIZE_SMALL)

    # Save ----

    provenance = style.base_provenance(
        "S12B",
        args,
        inputs=[args.tables_dir / "marker_stats.csv"],
        script=__file__,
    )
    provenance.update({
        "rows": classes,
        "row_labels": [style.DOTPLOT_CLASS_LABELS[c] for c in classes],
        "genes": genes,
        "genes_added_for_v9": [g for g in genes if g in style.PERIVASCULAR_GENES],
        "size": "fraction of the class's nuclei with the gene detected (X > 0)",
        "color": "mean of X (log1p normalised) per class, min-max scaled per gene across the rows",
        "n_nuclei_per_row": stats.drop_duplicates("cell_type").set_index("cell_type")["n_nuclei"].reindex(classes).to_dict(),
        "frac_detected": frac.round(4).to_dict(orient="index"),
        "scaled_mean": scaled.round(4).to_dict(orient="index"),
    })
    style.save_panel(fig, args.out_prefix, provenance)
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
