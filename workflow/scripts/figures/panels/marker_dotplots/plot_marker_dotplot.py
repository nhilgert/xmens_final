"""Draw one final-label marker dot plot of release v9 on the Table M1 v9 genes.

Called by rule panels_marker_dotplots_panel (workflow/rules/91_panels_marker_dotplots.smk), once per species
and plot (marker_panels.PLOTS). Reads <tables-dir>/<species>_marker_stats.csv, _table_m1_rule_check.csv and
_extract_provenance.json, Table M1 and the annotation vocabulary. Writes <out-prefix>.svg, .png (300 dpi) and
.provenance.json.

The plot is scanpy's DotPlot with the layout, legends and fonts of the earlier immune dot plots (scanpy
defaults: 0.37 in per gene, 0.35 in per row, 1.5 in of legends, labels 'small' on a 10 pt base). Dot size is
the % of the row's nuclei detecting the gene, colour the row's mean log1p CP10k scaled 0-1 per gene across
the rows shown. Columns are each row's Table M1 genes in Table M1 order (marker_panels.plot_genes): a gene the
Methods list under two rows of the plot is shown under both.
"""

import argparse
import json
import logging
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

import dotplot_style as style
import marker_panels as panels


# Encoding of the earlier dot plots: white-to-dark-red ramp, thin black dot edges, dot area growing with the
# detected fraction to the power 1.5 (scanpy's default, which keeps sparse genes visible but small).
CMAP = "Reds"
DOT_EDGE_COLOR = "black"
DOT_EDGE_WIDTH = 0.2
LARGEST_DOT_PT2 = 200.0
SIZE_EXPONENT = 1.5

# Size key fixed at 0-100%, so every plot's dots share one scale and one legend (20, 40, 60, 80, 100).
DOT_MIN, DOT_MAX = 0.0, 1.0

# Figure size of the earlier plots (inches): scanpy's per-gene and per-row sizes plus the legend column,
# passed explicitly so the main plot fills the height (no spacer row above it). The 2.5 in floor is scanpy's
# minimum height, below which the size key and colour bar overlap.
INCHES_PER_GENE = 0.37
INCHES_PER_ROW = 0.35
LABEL_HEIGHT_IN = 1.0
LEGENDS_WIDTH_IN = 1.5
MIN_HEIGHT_IN = 2.5
SIZE_TITLE = "Fraction of nuclei\nin group (%)"
COLOR_TITLE = "Mean expression\n(scaled per gene)"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, choices=list(panels.PLOTS))
    parser.add_argument("--plot", required=True)
    parser.add_argument("--tables-dir", required=True, type=Path)
    parser.add_argument("--table-m1", required=True, type=Path)
    parser.add_argument("--vocabulary", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    return parser.parse_args()


def plot_matrices(stats: pd.DataFrame, level: str, rows: list, genes: list):
    """Rows x genes (each gene once): fraction detecting (0-1), mean log1p CP10k, and the mean scaled 0-1 per gene."""
    genes = list(dict.fromkeys(genes))
    at_level = stats[stats.level == level]
    fraction = at_level.pivot(index="label", columns="gene", values="pct_detected").loc[rows, genes] / 100.0
    mean = at_level.pivot(index="label", columns="gene", values="mean_log1p_cp10k").loc[rows, genes]

    # scanpy's standard_scale='var': each gene's row means rescaled to span 0-1 over the rows shown.
    scaled = mean - mean.min()
    scaled = (scaled / scaled.max()).fillna(0.0)
    return fraction, mean, scaled


def column_keys(columns: pd.DataFrame) -> list:
    """Unique column names for scanpy, which rejects repeated var names: a repeated gene gets its row appended."""
    keys, seen = [], set()
    for gene, row in zip(columns.gene, columns.m1_cell_type):
        keys.append(f"{gene} ({row})" if gene in seen else gene)
        seen.add(gene)
    return keys


def draw(fraction: pd.DataFrame, scaled: pd.DataFrame, row_names: list, gene_labels: list):
    """scanpy DotPlot of precomputed values; returns the figure."""
    # A one-nucleus-per-row AnnData carries the row and gene names; the values come from the tables.
    holder = ad.AnnData(
        X=np.zeros(fraction.shape, dtype=np.float32),
        obs=pd.DataFrame({"row": pd.Categorical(row_names, categories=row_names)}, index=row_names),
        var=pd.DataFrame(index=list(fraction.columns)),
    )
    size_table = fraction.set_axis(row_names, axis=0)
    color_table = scaled.set_axis(row_names, axis=0)
    figsize = (
        INCHES_PER_GENE * fraction.shape[1] + LEGENDS_WIDTH_IN,
        max(MIN_HEIGHT_IN, INCHES_PER_ROW * fraction.shape[0] + LABEL_HEIGHT_IN),
    )
    dotplot = sc.pl.DotPlot(
        holder,
        list(fraction.columns),
        "row",
        categories_order=row_names,
        figsize=figsize,
        dot_size_df=size_table,
        dot_color_df=color_table,
        standard_scale=None,
        vmin=0.0,
        vmax=1.0,
    )
    dotplot.style(
        cmap=CMAP,
        dot_edge_color=DOT_EDGE_COLOR,
        dot_edge_lw=DOT_EDGE_WIDTH,
        smallest_dot=0.0,
        largest_dot=LARGEST_DOT_PT2,
        size_exponent=SIZE_EXPONENT,
        dot_min=DOT_MIN,
        dot_max=DOT_MAX,
    )
    dotplot.legend(size_title=SIZE_TITLE, colorbar_title=COLOR_TITLE)
    dotplot.make_figure()

    # Tick labels are the gene names (a repeated gene's column key carries its row); scanpy's rotation and size
    # stay on the tick labels.
    main_ax = dotplot.get_axes()["mainplot_ax"]
    main_ax.set_xticklabels(gene_labels)
    return dotplot.fig


def main() -> None:
    args = parse_arguments()
    definition = panels.PLOTS[args.species][args.plot]
    matplotlib.rcParams.update({**style.FONT_RC, "font.size": style.BASE_FONT_SIZE})

    # Load ----

    stats = pd.read_csv(args.tables_dir / f"{args.species}_marker_stats.csv")
    checks = pd.read_csv(args.tables_dir / f"{args.species}_table_m1_rule_check.csv")
    extract = json.loads((args.tables_dir / f"{args.species}_extract_provenance.json").read_text())
    table_m1 = panels.read_table_m1(args.table_m1)
    vocabulary = yaml.safe_load(args.vocabulary.read_text())

    # Rows are release labels; columns the Table M1 genes of those rows that the object has.
    rows = definition["rows"]
    names = panels.display_names(vocabulary)
    columns = panels.plot_genes(args.species, args.plot, table_m1)
    missing = list(dict.fromkeys(g for g in columns.gene if g in extract["table_m1_genes_missing_from_var"]))
    columns = columns[~columns.gene.isin(missing)].reset_index(drop=True)
    genes = list(columns.gene)
    fraction, mean, scaled = plot_matrices(stats, definition["level"], rows, genes)
    n_nuclei = (
        stats[stats.level == definition["level"]]
        .drop_duplicates("label")
        .set_index("label")["n_nuclei"]
        .loc[rows]
    )
    logging.info("%s %s: %d rows x %d gene columns; missing from var: %s",
                 args.species, args.plot, len(rows), len(genes), missing)

    # Figure ----

    keys = column_keys(columns)
    fig = draw(
        fraction[genes].set_axis(keys, axis=1),
        scaled[genes].set_axis(keys, axis=1),
        [names[label] for label in rows],
        genes,
    )
    n_anchored = style.anchor_rotated_text(fig)

    # Save ----

    outputs = style.save_svg_png(fig, args.out_prefix)
    plt.close(fig)

    # Rule status of each plotted column for its own Table M1 row (the row it is listed under in this plot).
    own_row = checks.set_index(["m1_cell_type", "gene"])
    gene_records = []
    for column in columns.itertuples():
        check = own_row.loc[(column.m1_cell_type, column.gene)]
        gene_records.append({
            "gene": column.gene,
            "table_m1_row": column.m1_cell_type,
            "table_m1_class": column.m1_class,
            "criterion": check["criterion"],
            "pct_detected_type": round(float(check["pct_detected_type"]), 2),
            "pct_detected_comparison": round(float(check["pct_detected_comparison"]), 2),
            "comparison_nuclei": check["comparison_nuclei"],
            "meets_criterion_v9": bool(check["meets_criterion"]),
        })
    provenance = {
        "plot": args.plot,
        "title": definition["title"],
        "species": args.species,
        "version": "v9",
        "release": extract["release"],
        "object": extract["object"],
        "object_sha256": extract["object_sha256"],
        "vocabulary_schema": extract["vocabulary_schema"],
        "rows": [
            {"label": label, "shown_as": names[label], "level": definition["level"], "n_nuclei": int(n_nuclei[label]),
             **({"pending": panels.PENDING_SIGN_OFF[label]} if label in panels.PENDING_SIGN_OFF else {})}
            for label in rows
        ],
        "genes": gene_records,
        "repeated_genes": sorted(set(columns.gene[columns.gene.duplicated()])),
        "genes_missing_from_var": missing,
        "dot_size": ("% of the row's nuclei with >= 1 count of the gene "
                     "(scanpy size scale: area ~ fraction^1.5, 0-100%)"),
        "dot_color": "mean log1p CP10k of the row's nuclei, min-max scaled 0-1 per gene across the rows shown",
        "embedding_outliers": (f"included ({extract['n_embedding_outliers_included']} in the object; "
                               "hidden on UMAPs only)"),
        "fraction_detected": fraction.round(4).to_dict(orient="index"),
        "mean_log1p_cp10k": mean.round(4).to_dict(orient="index"),
        "scaled_mean": scaled.round(4).to_dict(orient="index"),
        "rotated_labels_anchored": n_anchored,
        "outputs": [str(p) for p in outputs],
        "inputs": {
            str(path): style.sha256_of(path)
            for path in [
                args.tables_dir / f"{args.species}_marker_stats.csv",
                args.tables_dir / f"{args.species}_table_m1_rule_check.csv",
                args.table_m1,
                args.vocabulary,
            ]
        },
        "script_sha256": style.sha256_of(Path(__file__)),
        "marker_panels_sha256": style.sha256_of(Path(panels.__file__)),
        "style_sha256": style.sha256_of(Path(style.__file__)),
        "scanpy": version("scanpy"),
    }
    style.write_provenance(Path(f"{args.out_prefix}.provenance.json"), provenance)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
