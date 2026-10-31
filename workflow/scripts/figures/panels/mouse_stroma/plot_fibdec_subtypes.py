"""Fig. S17A-C: mouse fibroblast and decidual subtypes across the X-Mens stages.

Reads the fib/dec plot table (extract_panel_table.py) and writes, each with a .provenance.json:
  --out-a  A: the fib/dec object's own UMAP (obsm X_umap) coloured by subtype, direct labels
  --out-b  B: % of each subtype among the fibroblast/decidual nuclei of each stage (stacked bars)
  --out-c  C: marker dot plot, subtypes x genes (size = % detected, colour = mean scaled 0-1 per gene)
White background; subtype colours from the manuscript (panel_style.SUBTYPE_COLORS).
Rule: panels_mouse_stroma_s17 (workflow/rules/91_panels_mouse_stroma.smk).
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

import panel_style as ps

SEED = 0  # UMAP draw order is shuffled so no subtype systematically covers another

# Panel A: the manuscript shows this embedding turned 90 degrees clockwise (eFib3 top left, dec3 bottom
# left); the same display rotation is applied here. UMAP axes carry no meaning, so rotation changes nothing.
UMAP_FIG_IN = (2.7, 2.6)
UMAP_MARKER_AREA = 0.6  # pt^2
# Label offset from the subtype's median (rotated UMAP units) and whether the label sits outside the cloud
# with a leader line; a leader label is placed dx left of the cloud's left edge, at the median's height + dy.
UMAP_LABELS = {
    "eFib3": ((-0.6, 3.2), False),
    "eFib2": ((1.6, 2.6), False),
    "eFibP": ((-2.6, 0.6), False),
    "eFib1": ((2.3, 0.4), False),
    "eFib2i": ((3.2, -0.6), False),
    "dec1": ((-0.3, 0.0), True),
    "dec2": ((4.2, 0.4), False),
    "dec3": ((-2.4, -1.6), False),
    "decAp": ((2.6, 0.0), False),
}

# Panel B: segment order and bar order of the manuscript panel.
BAR_FIG_IN = (3.7, 2.4)
BAR_SEGMENT_ORDER = ["eFib3", "eFib2", "eFib2i", "eFib1", "eFibP", "dec2", "dec1", "dec3", "decAp"]

# Panel C: 33 marker genes in the manuscript order; the manuscript labels 'Col1a3' and 'Kcdn2' are not mouse
# genes and correspond to Col1a2 and Kcnd2 (Col1a2 reproduces the manuscript dots; see provenance).
DOT_FIG_IN = (7.1, 2.75)
DOT_RECT = [0.075, 0.28, 0.74, 0.69]
DOT_LARGEST_AREA = 95.0  # pt^2 at dot_max
DOT_SIZE_KEYS_PCT = [20, 60, 100]
DOT_CMAP = "Reds"
EDGE_GREY = "#808080"


def rotated_umap(table: pd.DataFrame) -> pd.DataFrame:
    """Display coordinates: 90 degrees clockwise, x = UMAP2, y = -UMAP1."""
    return pd.DataFrame({"x": table["X_umap_2"], "y": -table["X_umap_1"]}, index=table.index)


def draw_umap(table: pd.DataFrame, input_record: dict, panel: str, out: Path) -> None:
    ps.apply_style(dark=False)
    shown = table[~table["embedding_outlier"]]
    coords = rotated_umap(shown)
    order = np.random.default_rng(SEED).permutation(len(shown))

    fig = plt.figure(figsize=UMAP_FIG_IN)
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])

    # One point per nucleus, colour = subtype; shuffled draw order.
    ax.scatter(
        coords["x"].to_numpy()[order],
        coords["y"].to_numpy()[order],
        c=shown["cell_subtype"].map(ps.SUBTYPE_COLORS).to_numpy()[order],
        s=UMAP_MARKER_AREA,
        linewidths=0,
        rasterized=True,
    )

    # Direct labels beside each subtype's median. dec1 lies inside the decidual cloud, so its label sits outside
    # with a leader to the cloud's left edge (15th percentile of x at the median y), as in the manuscript.
    medians = coords.groupby(shown["cell_subtype"]).median()
    left_edges = coords.groupby(shown["cell_subtype"])["x"].quantile(0.15)
    for subtype, ((dx, dy), has_leader) in UMAP_LABELS.items():
        anchor = medians.loc[subtype]
        if has_leader:
            position = (coords["x"].min() + dx, anchor["y"] + dy)
            target = (left_edges.loc[subtype], anchor["y"])
            ax.annotate(subtype, xy=target, xytext=position, ha="right", va="center", fontsize=ps.SIZE_TITLE,
                        arrowprops=dict(arrowstyle="-", lw=0.8, color="black", shrinkA=1, shrinkB=0))
        else:
            position = (anchor["x"] + dx, anchor["y"] + dy)
            ax.text(*position, subtype, ha="center", va="center", fontsize=ps.SIZE_TITLE)

    ax.set_aspect("equal")
    ax.margins(0.04)
    ax.axis("off")

    provenance = {
        "panel": panel,
        "script": Path(__file__).name,
        "inputs": input_record,
        "embedding": "obsm['X_umap'] of the fib/dec object (its own UMAP; identical in the 2026_10_22, v8 and v9 objects)",
        "display": "rotated 90 deg clockwise (x = UMAP2, y = -UMAP1) to match the manuscript layout",
        "nuclei_drawn": int(len(shown)),
        "embedding_outliers_hidden": int(table["embedding_outlier"].sum()),
        "nuclei_by_subtype_drawn": shown["cell_subtype"].value_counts().reindex(ps.SUBTYPES).to_dict(),
        "colours": ps.SUBTYPE_COLORS,
        "colour_source": "sampled from the manuscript Fig. S17B bars (release objects carry no mouse label_colors)",
        "draw_order": f"numpy default_rng({SEED}) permutation",
    }
    ps.save_panel(fig, out, provenance)


def draw_bars(table: pd.DataFrame, input_record: dict, panel: str, out: Path) -> None:
    ps.apply_style(dark=False)

    # One bar per stage: every fibroblast/decidual nucleus of the stage (placed or not, outliers included).
    counts = pd.crosstab(table["timepoint"], table["cell_subtype"]).reindex(index=ps.STAGES, columns=BAR_SEGMENT_ORDER, fill_value=0)
    percent = counts.div(counts.sum(axis=1), axis=0) * 100

    fig = plt.figure(figsize=BAR_FIG_IN)
    ax = fig.add_axes([0.2, 0.12, 0.76, 0.78])

    # Segments stacked left to right in the manuscript order; colour = subtype as in panel A.
    y = np.arange(len(ps.STAGES))
    left = np.zeros(len(ps.STAGES))
    for subtype in BAR_SEGMENT_ORDER:
        ax.barh(y, percent[subtype], left=left, height=0.8, color=ps.SUBTYPE_COLORS[subtype], linewidth=0)
        left += percent[subtype].to_numpy()

    ax.set_ylim(len(ps.STAGES) - 0.5, -0.5)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.set_yticks(y)
    labels = [ps.STAGE_LABELS[s].replace("4.5 dpi menstrual", "4.5 dpi\nmenstrual") for s in ps.STAGES]
    ax.set_yticklabels(labels, rotation=30, ha="right", rotation_mode="anchor", fontsize=ps.SIZE_LEGEND)
    ax.tick_params(axis="x", labelsize=ps.SIZE_LEGEND)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("relative cell subtype abundance", loc="center", fontweight="bold", fontsize=ps.SIZE_TITLE, pad=4)

    provenance = {
        "panel": panel,
        "script": Path(__file__).name,
        "inputs": input_record,
        "unit": "% of all fibroblast/decidual nuclei of the stage (all nuclei of the fib/dec object; embedding outliers included)",
        "segment_order": BAR_SEGMENT_ORDER,
        "n_per_stage": counts.sum(axis=1).to_dict(),
        "counts": counts.to_dict(orient="index"),
        "percent": percent.round(3).to_dict(orient="index"),
        "colours": ps.SUBTYPE_COLORS,
    }
    ps.save_panel(fig, out, provenance)


def draw_dotplot(table: pd.DataFrame, genes: list[str], input_record: dict, panel: str, out: Path) -> None:
    ps.apply_style(dark=False)

    # Per subtype: % nuclei with the gene detected and mean expression over all its nuclei (all stages).
    grouped = table.groupby("cell_subtype")
    fraction = (table[genes] > 0).groupby(table["cell_subtype"]).mean().loc[ps.SUBTYPES]
    mean = grouped[genes].mean().loc[ps.SUBTYPES]
    scaled = (mean - mean.min()) / (mean.max() - mean.min())
    dot_max = ps.scanpy_dot_max(float(fraction.values.max()))

    fig = plt.figure(figsize=DOT_FIG_IN)
    ax = fig.add_axes(DOT_RECT)
    xs, ys = np.meshgrid(np.arange(len(genes)), np.arange(len(ps.SUBTYPES)))
    ax.scatter(
        xs.ravel(),
        ys.ravel(),
        s=ps.dot_area(fraction.values.ravel(), dot_max, DOT_LARGEST_AREA),
        c=scaled.values.ravel(),
        cmap=DOT_CMAP,
        vmin=0,
        vmax=1,
        edgecolors=EDGE_GREY,
        linewidths=0.4,
    )
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(ps.SUBTYPES) - 0.5, -0.6)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=ps.SIZE_TITLE, fontweight="bold")
    ax.set_yticks(range(len(ps.SUBTYPES)))
    ax.set_yticklabels(ps.SUBTYPES, fontsize=ps.SIZE_TITLE, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(True)

    # Size key with the manuscript's legend values.
    key_ax = fig.add_axes([0.85, 0.55, 0.06, 0.28])
    key_ax.scatter(np.zeros(len(DOT_SIZE_KEYS_PCT)), np.arange(len(DOT_SIZE_KEYS_PCT)),
                   s=ps.dot_area(np.array(DOT_SIZE_KEYS_PCT) / 100, dot_max, DOT_LARGEST_AREA), c="black")
    for i, pct in enumerate(DOT_SIZE_KEYS_PCT):
        key_ax.text(1.2, i, str(pct), va="center", fontsize=ps.SIZE_LEGEND)
    key_ax.set_xlim(-1, 3)
    key_ax.set_ylim(len(DOT_SIZE_KEYS_PCT) - 0.5, -0.6)
    key_ax.axis("off")
    fig.text(0.895, 0.86, "fraction of\ncells (%)", ha="center", va="bottom", fontsize=ps.SIZE_LEGEND)

    # Colour key: per gene, 0 = lowest and 1 = highest subtype mean.
    cbar_ax = fig.add_axes([0.85, 0.33, 0.1, 0.035])
    colorbar = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap=DOT_CMAP), cax=cbar_ax, orientation="horizontal", ticks=[])
    colorbar.outline.set_linewidth(0.4)
    fig.text(0.845, 0.347, "0", ha="right", va="center", fontsize=ps.SIZE_LEGEND)
    fig.text(0.955, 0.347, "1", ha="left", va="center", fontsize=ps.SIZE_LEGEND)
    fig.text(0.9, 0.39, "mean expression\n(log normalized)", ha="center", va="bottom", fontsize=ps.SIZE_LEGEND)

    provenance = {
        "panel": panel,
        "script": Path(__file__).name,
        "inputs": input_record,
        "genes": genes,
        "gene_label_corrections": {"Col1a3": "Col1a2", "Kcdn2": "Kcnd2"},
        "nuclei": "all fibroblast/decidual nuclei of the object, all stages (embedding outliers included)",
        "size": f"% nuclei with X > 0; area = (fraction / {dot_max})^1.5 * {DOT_LARGEST_AREA} pt^2 (scanpy dotplot rule)",
        "colour": "subtype mean of X (log1p counts per 10k), min-max scaled per gene across subtypes (cmap Reds)",
        "n_nuclei": table["cell_subtype"].value_counts().reindex(ps.SUBTYPES).to_dict(),
        "fraction": fraction.round(4).to_dict(orient="index"),
        "mean_log_norm": mean.round(5).to_dict(orient="index"),
        "colour_scaled": scaled.round(4).to_dict(orient="index"),
    }
    ps.save_panel(fig, out, provenance)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--table-prov", required=True, type=Path)
    parser.add_argument("--genes", required=True, help="comma-separated, in display order")
    parser.add_argument("--out-a", required=True, type=Path)
    parser.add_argument("--out-b", required=True, type=Path)
    parser.add_argument("--out-c", required=True, type=Path)
    args = parser.parse_args()

    table, input_record = ps.load_table(args.table, args.table_prov)
    unknown = set(table["cell_subtype"]) - set(ps.SUBTYPES)
    if unknown:
        raise ValueError(f"subtypes outside the S17 vocabulary: {sorted(unknown)}")
    genes = args.genes.split(",")
    missing = [g for g in genes if g not in table.columns]
    if missing:
        raise KeyError(f"genes not in {args.table.name}: {missing}")

    draw_umap(table, input_record, "S17A", args.out_a)
    draw_bars(table, input_record, "S17B", args.out_b)
    draw_dotplot(table, genes, input_record, "S17C", args.out_c)


if __name__ == "__main__":
    main()
