"""Fig. S15A: dot plot of one gene per major cell class (rows) and stage (columns), mouse whole tissue.

Reads the whole-tissue plot table of one release (extract_panel_table.py) and writes <out>.{pdf,svg,png}
and <out>.provenance.json. Dot size = % of nuclei with the gene detected (X > 0); colour = mean
log-normalized expression over all nuclei of the group, min-max scaled across the plotted groups
(scanpy dotplot, standard_scale='var'), as in the manuscript panel.
Rule: panels_mouse_stroma_s15a (workflow/rules/91_panels_mouse_stroma.smk).
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

import panel_style as ps

FIG_SIZE_IN = (3.3, 2.15)  # S15A at the manuscript's two-panel row (180 mm figure width)
DOTS_RECT = [0.17, 0.33, 0.43, 0.56]  # dot-plot axes [left, bottom, width, height] in figure fraction
LARGEST_DOT_AREA = 80.0  # pt^2 at dot_max
SIZE_KEYS_PCT = [5, 10, 15, 20]  # legend keys of the manuscript panel
CMAP = "Reds"  # white = lowest group mean, dark red = highest, as in the manuscript
EDGE_GREY = "#808080"


def class_stage_stats(table, gene: str, classes: list[str]):
    """Fraction detected (%), mean expression and n per class x stage, over every nucleus of the object."""
    table = table.assign(cls=table["coarse_refined"].map(ps.CLASS_ABBREV))
    if table["cls"].isna().any():
        raise ValueError(f"coarse_refined values without an abbreviation: {sorted(table.loc[table['cls'].isna(), 'coarse_refined'].unique())}")
    grouped = table.groupby(["cls", "timepoint"])[gene]
    fraction_pct = grouped.apply(lambda v: (v > 0).mean() * 100).unstack().loc[classes, ps.STAGES]
    mean = grouped.mean().unstack().loc[classes, ps.STAGES]
    n = grouped.size().unstack().loc[classes, ps.STAGES]
    return fraction_pct, mean, n


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--table-prov", required=True, type=Path)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--classes", required=True, help="comma-separated abbreviations, top to bottom")
    parser.add_argument("--panel", required=True)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    classes = args.classes.split(",")

    ps.apply_style(dark=False)
    table, input_record = ps.load_table(args.table, args.table_prov)

    # Statistics: every nucleus of the class and stage, placed on the tissue or not; embedding outliers are
    # retained in analyses, so they count here.
    fraction_pct, mean, n = class_stage_stats(table, args.gene, classes)
    scaled = (mean - mean.values.min()) / (mean.values.max() - mean.values.min())
    dot_max = ps.scanpy_dot_max(fraction_pct.values.max() / 100)

    fig = plt.figure(figsize=FIG_SIZE_IN)
    ax = fig.add_axes(DOTS_RECT)

    # One dot per class x stage: area = detection fraction, colour = scaled group mean.
    xs, ys = np.meshgrid(np.arange(len(ps.STAGES)), np.arange(len(classes)))
    ax.scatter(
        xs.ravel(),
        ys.ravel(),
        s=ps.dot_area(fraction_pct.values.ravel() / 100, dot_max, LARGEST_DOT_AREA),
        c=scaled.values.ravel(),
        cmap=CMAP,
        vmin=0,
        vmax=1,
        edgecolors=EDGE_GREY,
        linewidths=0.4,
    )

    # Rows top to bottom in the manuscript order; the box frame stays closed as in the manuscript panel.
    ax.set_xlim(-0.6, len(ps.STAGES) - 0.4)
    ax.set_ylim(len(classes) - 0.5, -0.5)
    ax.set_xticks(range(len(ps.STAGES)))
    ax.set_xticklabels([ps.STAGE_LABELS[s].replace("4.5 dpi menstrual", "4.5 dpi\nmenstrual") for s in ps.STAGES],
                       rotation=45, ha="right", rotation_mode="anchor", fontsize=ps.SIZE_TITLE, fontweight="bold")
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes, fontsize=ps.SIZE_TITLE, fontweight="bold")
    for side in ("top", "right"):
        ax.spines[side].set_visible(True)

    # Title with the gene name in italics.
    ps.mixed_style_text(fig, DOTS_RECT[0] + 0.07, DOTS_RECT[1] + DOTS_RECT[3] + 0.04,
                        [(args.gene, "italic"), (" expression", "normal")], fontsize=ps.SIZE_TITLE, fontweight="bold")

    # Size key: the manuscript's legend values on the same area scale.
    key_ax = fig.add_axes([0.66, 0.50, 0.12, 0.33])
    key_ax.scatter(np.zeros(len(SIZE_KEYS_PCT)), np.arange(len(SIZE_KEYS_PCT)),
                   s=ps.dot_area(np.array(SIZE_KEYS_PCT) / 100, dot_max, LARGEST_DOT_AREA), c="black")
    for i, pct in enumerate(SIZE_KEYS_PCT):
        key_ax.text(0.9, i, str(pct), va="center", fontsize=ps.SIZE_LEGEND)
    key_ax.set_xlim(-0.8, 2.5)
    key_ax.set_ylim(len(SIZE_KEYS_PCT) - 0.5, -0.6)
    key_ax.axis("off")
    fig.text(0.665 + 0.06, 0.87, "fraction of\ncells (%)", ha="center", va="bottom", fontsize=ps.SIZE_LEGEND)

    # Colour key: 0 = lowest and 1 = highest group mean of this panel.
    cbar_ax = fig.add_axes([0.66, 0.24, 0.27, 0.035])
    colorbar = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=CMAP), cax=cbar_ax, orientation="horizontal", ticks=[])
    colorbar.outline.set_linewidth(0.4)
    fig.text(0.655, 0.255, "0", ha="right", va="center", fontsize=ps.SIZE_LEGEND)
    fig.text(0.935, 0.255, "1", ha="left", va="center", fontsize=ps.SIZE_LEGEND)
    fig.text(0.795, 0.30, "mean expression\n(log normalized)", ha="center", va="bottom", fontsize=ps.SIZE_LEGEND)

    provenance = {
        "panel": args.panel,
        "script": Path(__file__).name,
        "inputs": input_record,
        "gene": args.gene,
        "rows": classes,
        "row_source": "obs['coarse_refined'] mapped to manuscript abbreviations",
        "columns": ps.STAGES,
        "nuclei": "all nuclei of each class and stage (placed or not; embedding outliers included, as the release keeps them in analyses)",
        "size": f"% nuclei with X > 0; area = (fraction / {dot_max})^1.5 * {LARGEST_DOT_AREA} pt^2 (scanpy dotplot rule)",
        "colour": "group mean of X (log1p counts per 10k), min-max scaled across the plotted groups (cmap Reds)",
        "fraction_pct": fraction_pct.round(3).to_dict(orient="index"),
        "mean_log_norm": mean.round(5).to_dict(orient="index"),
        "colour_scaled": scaled.round(4).to_dict(orient="index"),
        "n_nuclei": n.to_dict(orient="index"),
        "pending_labels_in_rows": ps.pending_label_notes(
            table.loc[table["coarse_refined"].map(ps.CLASS_ABBREV).isin(classes), "cell_subtype"]
        ),
        "classes_not_shown": sorted(set(table["coarse_refined"].map(ps.CLASS_ABBREV)) - set(classes)),
        "class_composition": {
            cls: table.loc[table["coarse_refined"].map(ps.CLASS_ABBREV) == cls, "cell_subtype"].value_counts().to_dict()
            for cls in classes
        },
    }
    ps.save_panel(fig, args.out, provenance)


if __name__ == "__main__":
    main()
