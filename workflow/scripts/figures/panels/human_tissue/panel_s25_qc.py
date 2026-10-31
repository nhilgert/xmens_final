"""Fig. S25A-D: sample table, cell-class marker dot plot and per-sample QC violins (human Slide-tags).

Called by rule `panels_human_tissue_s25`, once per release.
Reads:  the plot table and extraction provenance of one release (extract_human_tissue.py) and the
        sample sheet (config/samples.tsv, 10x lanes per sample).
Writes: <out-dir>/panel_S25{A,B,C,D}_<release>.pdf/.svg/.png and .provenance.json, one set per sub-panel.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

import panel_common as pc


# Panel A: sample table ----

# Column headers as in the published table; the first reads "cycle phase" as in the S25 legend
# (the published image says "timepoint").
TABLE_HEADERS = ["cycle phase", "section\nthickness", "#10X\nlanes", "number of\ncells", "spatially\nmapped (%)"]

# Section thickness is not recorded in config/samples.tsv; these are the values of the published table.
SECTION_THICKNESS_UM = {"secretory": 20, "menstrual": 20}

TABLE_W_IN, TABLE_H_IN = 3.15, 0.78
HEADER_LINESPACING = 1.15  # two-line headers; at 1.0 the descender of "spatially" touches "(%)" in Nimbus Sans
TABLE_COLUMN_X_IN = [0.04, 1.05, 1.63, 2.17, 2.8]  # first column left-aligned, the others centred


# Panel B: dot plot ----

# Row order and abbreviations of the published S25B; its legend writes perivascular cells as "pvc".
DOTPLOT_ROWS = [
    "epithelial_cell",
    "endometrial_fibroblast",
    "decidual_cell",
    "immune_cell",
    "endothelial_cell",
    "perivascular_cell",
    "smooth_muscle_cell",
]
ROW_LABELS = {**pc.ABBREVIATIONS, "perivascular_cell": "pvc"}

# Encodings follow the scanpy dot plot of the published panel: colour = mean log1p(CP10k) per class,
# min-max scaled across classes for each gene (standard_scale="var"); dot area = (fraction of nuclei
# with at least one count) ** 1.5, the scanpy size exponent, at most MAX_DOT_AREA.
CMAP = "Reds"
SIZE_EXPONENT = 1.5
MAX_DOT_AREA = 75  # pt^2; the 100% dot is 8.7 pt wide, just under the 9.7 pt column pitch
DOT_EDGE = {"edgecolors": "black", "linewidths": 0.2}
LEGEND_FRACTIONS = [1.0, 0.5]

COL_PITCH_IN, ROW_PITCH_IN = 0.135, 0.17
DOT_LEFT_IN, DOT_BOTTOM_IN, DOT_TOP_IN = 0.34, 0.5, 0.06
KEY_GAP_IN, KEY_W_IN = 0.25, 1.05


# Panels C and D: QC violins ----

# Natural log of (1 + value), as the published axes (secretory counts span about 6.6 to 10.2).
QC_METRICS = {
    "C": ("total_counts", "counts per cell (log)"),
    "D": ("n_genes", "genes per cell (log)"),
}
VIOLIN_W_IN, VIOLIN_H_IN = 1.5, 1.5


def save_subpanel(fig, out_dir: Path, letter: str, release: str, facecolor: str = "white") -> list[str]:
    return pc.save_panel(fig, out_dir / f"panel_S25{letter}_{release}", facecolor=facecolor)


def count_lanes(samples_path: Path) -> dict:
    """10x lanes per human sample (rows of the sample sheet)."""
    samples = pd.read_csv(samples_path, sep="\t")
    human = samples[samples["species"] == "human"]
    return human.groupby("timepoint").size().reindex(pc.SAMPLES).astype(int).to_dict()


def table_rows(cells: pd.DataFrame, lanes: dict) -> pd.DataFrame:
    # Number of cells = every nucleus of the release object (embedding outliers included, as in all
    # statistics); spatially mapped = nuclei with QC'd coordinates (spatial_qc_reason 'kept').
    rows = []
    for sample in pc.SAMPLES:
        in_sample = cells["timepoint"] == sample
        reason = cells.loc[in_sample, "spatial_qc_reason"].astype(str)
        n_cells = int(in_sample.sum())
        n_kept = int((reason == pc.SPATIAL_KEPT).sum())
        rows.append({
            "sample": sample,
            "section_thickness_um": SECTION_THICKNESS_UM[sample],
            "lanes": lanes[sample],
            "n_cells": n_cells,
            "n_spatial_kept": n_kept,
            "n_positioned_before_qc": int((reason != "unmapped").sum()),
            "pct_spatially_mapped": 100 * n_kept / n_cells,
        })
    return pd.DataFrame(rows)


def draw_table(rows: pd.DataFrame):
    fig = plt.figure(figsize=(TABLE_W_IN, TABLE_H_IN), facecolor="white")

    def place(x_in: float, y_in_from_top: float, text: str, column: int, **kwargs) -> None:
        ha = "left" if column == 0 else "center"
        fig.text(x_in / TABLE_W_IN, 1 - y_in_from_top / TABLE_H_IN, text, ha=ha, **kwargs)

    # A single rule above the header, as in the published table.
    fig.add_artist(plt.Line2D([0.0, 1.0], [1 - 0.02 / TABLE_H_IN] * 2, color="black", lw=0.8))

    # Headers share their bottom line, so the one-line header sits level with the two-line ones.
    for column, (x_in, header) in enumerate(zip(TABLE_COLUMN_X_IN, TABLE_HEADERS)):
        place(x_in, 0.32, header, column, va="bottom", fontweight="bold", linespacing=HEADER_LINESPACING)

    for i, row in enumerate(rows.itertuples()):
        values = [
            row.sample,
            f"{row.section_thickness_um} µm",
            f"{row.lanes}",
            f"{row.n_cells:,}",
            f"{row.pct_spatially_mapped:.2f}",
        ]
        for column, (x_in, value) in enumerate(zip(TABLE_COLUMN_X_IN, values)):
            place(x_in, 0.47 + i * 0.2, value, column, va="center")
    return fig


def dotplot_stats(cells: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Per class and gene: n, fraction detected, mean log1p(CP10k) and its per-gene 0-1 scaling."""
    expr = cells[[f"expr_{g}" for g in genes]].set_axis(genes, axis=1)
    groups = expr.groupby(cells["cell_type"])
    mean = groups.mean().reindex(DOTPLOT_ROWS)
    fraction = (expr > 0).groupby(cells["cell_type"]).mean().reindex(DOTPLOT_ROWS)
    scaled = (mean - mean.min()) / (mean.max() - mean.min())

    stats = pd.concat(
        {"mean_log1p_cp10k": mean.stack(), "fraction_detected": fraction.stack(), "scaled_mean": scaled.stack()},
        axis=1,
    )
    stats.index.names = ["cell_type", "gene"]
    stats = stats.reset_index()
    stats["n_nuclei"] = stats["cell_type"].map(cells["cell_type"].value_counts())
    return stats


def draw_dotplot_keys(fig, left_in: float, norm: Normalize) -> None:
    """Size and colour keys stacked in one column, positioned in inches from the figure's top."""
    fig_w_in, fig_h_in = fig.get_size_inches()
    key_ax = fig.add_axes([left_in / fig_w_in, 0, KEY_W_IN / fig_w_in, 1])
    key_ax.set_xlim(0, KEY_W_IN)
    key_ax.set_ylim(fig_h_in, 0)
    key_ax.axis("off")

    # Size key: the 100% and 50% dots, drawn with the same area rule as the data.
    key_ax.text(0, 0.06, "fraction of cells\nin cluster (%)", va="top", fontsize=6.5, linespacing=HEADER_LINESPACING)
    for y_in, fraction in zip((0.43, 0.67), LEGEND_FRACTIONS):
        key_ax.scatter(
            0.1,
            y_in,
            s=MAX_DOT_AREA * fraction**SIZE_EXPONENT,
            color=plt.get_cmap(CMAP)(0.95),
            clip_on=False,
            **DOT_EDGE,
        )
        key_ax.text(0.25, y_in, f"{100 * fraction:.0f}", va="center", fontsize=6.5)

    # Colour key: the per-gene scaled mean, 0 = lowest class, 1 = highest.
    key_ax.text(0, 0.93, "normalized mean\nexpression in cluster", va="top", fontsize=6.5, linespacing=HEADER_LINESPACING)
    cbar_ax = fig.add_axes([left_in / fig_w_in, 1 - 1.3 / fig_h_in, 0.9 / fig_w_in, 0.07 / fig_h_in])
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbar_ax, orientation="horizontal", ticks=[0, 1])
    colorbar.outline.set_linewidth(0.4)
    colorbar.ax.tick_params(labelsize=6.5, length=2, width=0.5, pad=1.5)


def draw_dotplot(stats: pd.DataFrame, genes: list[str]):
    dot_w_in = len(genes) * COL_PITCH_IN
    dot_h_in = len(DOTPLOT_ROWS) * ROW_PITCH_IN
    fig_w_in = DOT_LEFT_IN + dot_w_in + KEY_GAP_IN + KEY_W_IN
    fig_h_in = DOT_TOP_IN + dot_h_in + DOT_BOTTOM_IN
    fig = plt.figure(figsize=(fig_w_in, fig_h_in), facecolor="white")
    ax = fig.add_axes([DOT_LEFT_IN / fig_w_in, DOT_BOTTOM_IN / fig_h_in, dot_w_in / fig_w_in, dot_h_in / fig_h_in])

    # Dots: x = gene, y = class (first row on top); colour and area as defined in the constants.
    x = stats["gene"].map({g: i for i, g in enumerate(genes)})
    y = stats["cell_type"].map({c: len(DOTPLOT_ROWS) - 1 - i for i, c in enumerate(DOTPLOT_ROWS)})
    norm = Normalize(0, 1)
    ax.scatter(
        x,
        y,
        s=MAX_DOT_AREA * stats["fraction_detected"] ** SIZE_EXPONENT,
        c=stats["scaled_mean"],
        cmap=CMAP,
        norm=norm,
        zorder=2,
        **DOT_EDGE,
    )

    # Boxed axes with one tick per gene and class, as in the published dot plot.
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(-0.6, len(DOTPLOT_ROWS) - 0.4)
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(DOTPLOT_ROWS)), [ROW_LABELS[c] for c in reversed(DOTPLOT_ROWS)], fontsize=7)
    ax.tick_params(length=2.5, width=0.6, pad=1.5)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.6)

    draw_dotplot_keys(fig, DOT_LEFT_IN + dot_w_in + KEY_GAP_IN, norm)
    return fig


def violin_summary(values: pd.Series) -> dict:
    quantiles = values.quantile([0, 0.25, 0.5, 0.75, 1]).round(4).tolist()
    return {"n": int(values.size), "min_q25_median_q75_max": quantiles}


def draw_violins(cells: pd.DataFrame, column: str, label: str) -> tuple:
    fig = plt.figure(figsize=(VIOLIN_W_IN, VIOLIN_H_IN), facecolor="white")
    ax = fig.add_axes([0.42 / VIOLIN_W_IN, 0.24 / VIOLIN_H_IN, 1.03 / VIOLIN_W_IN, 1.2 / VIOLIN_H_IN])

    # One violin per sample over all its nuclei: Gaussian KDE (Scott bandwidth) cut at the observed
    # range, filled with the sample colour of the published panel.
    log_values = [np.log1p(cells.loc[cells["timepoint"] == s, column].to_numpy()) for s in pc.SAMPLES]
    bodies = ax.violinplot(log_values, positions=range(len(pc.SAMPLES)), widths=0.8, showextrema=False)["bodies"]
    for body, sample in zip(bodies, pc.SAMPLES):
        body.set_facecolor(pc.SAMPLE_COLORS[sample])
        body.set_edgecolor("black")
        body.set_linewidth(0.5)
        body.set_alpha(1)

    # The y range is the data range plus 4% padding; ticks every 0.5 log units inside it.
    low = min(v.min() for v in log_values)
    high = max(v.max() for v in log_values)
    padding = 0.04 * (high - low)
    ax.set_ylim(low - padding, high + padding)
    ax.set_yticks(np.arange(np.ceil(2 * (low - padding)) / 2, high + padding, 0.5))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.1f"))

    ax.set_xticks(range(len(pc.SAMPLES)), pc.SAMPLES, fontsize=7)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_ylabel(label, fontweight="bold")
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)

    summary = {s: violin_summary(pd.Series(v)) for s, v in zip(pc.SAMPLES, log_values)}
    return fig, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--extract-provenance", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path, help="config/samples.tsv")
    parser.add_argument("--genes", required=True, help="comma-separated dot plot genes, in column order")
    parser.add_argument("--release", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    genes = args.genes.split(",")

    pc.set_style()
    cells = pc.load_cells(args.cells)
    scripts = [Path(__file__), Path(pc.__file__)]

    def provenance_for(letter: str) -> dict:
        provenance = pc.base_provenance(f"S25{letter}", args.release, args.cells, args.extract_provenance, scripts)
        provenance["nuclei"] = "all nuclei of the release object (embedding outliers included)"
        provenance["n_nuclei"] = int(len(cells))
        return provenance

    def write_provenance(letter: str, provenance: dict, outputs: list[str]) -> None:
        provenance["outputs"] = outputs
        pc.write_json(args.out_dir / f"panel_S25{letter}_{args.release}.provenance.json", provenance)

    # Panel A ----
    lanes = count_lanes(args.samples)
    rows = table_rows(cells, lanes)
    pc.log.info(f"{args.release} table:\n{rows.to_string(index=False)}")
    outputs = save_subpanel(draw_table(rows), args.out_dir, "A", args.release)
    provenance = provenance_for("A")
    provenance.update({
        "samples_sheet": {"path": str(args.samples.resolve()), "sha256": pc.sha256(args.samples)},
        "sources": {
            "lanes": "rows of config/samples.tsv per human timepoint",
            "section_thickness": "published S25A table (not in the sample sheet)",
            "number_of_cells": "nuclei per obs timepoint in the release object",
            "spatially_mapped": "100 x nuclei with spatial_qc_reason 'kept' / number of cells",
        },
        "header": TABLE_HEADERS,
        "rows": rows.round(4).to_dict(orient="records"),
    })
    write_provenance("A", provenance, outputs)

    # Panel B ----
    stats = dotplot_stats(cells, genes)
    outputs = save_subpanel(draw_dotplot(stats, genes), args.out_dir, "B", args.release)
    provenance = provenance_for("B")
    provenance.update({
        "expression": "mean of log1p(counts / total_counts * 1e4) from layers['counts'], per obs cell_type",
        "color": "per gene, (mean - min over classes) / (max - min over classes)",
        "size": f"area = {MAX_DOT_AREA} pt^2 x fraction_detected ** {SIZE_EXPONENT}; detected = count > 0",
        "genes": genes,
        "row_labels": {ROW_LABELS[c]: c for c in DOTPLOT_ROWS},
        "stats": stats.round(5).to_dict(orient="records"),
    })
    provenance["colors"]["colormap"] = CMAP
    write_provenance("B", provenance, outputs)

    # Panels C and D ----
    for letter, (column, label) in QC_METRICS.items():
        fig, summary = draw_violins(cells, column, label)
        outputs = save_subpanel(fig, args.out_dir, letter, args.release)
        provenance = provenance_for(letter)
        provenance.update({
            "metric": f"natural log of (1 + {column}); {column} from layers['counts'] of the release object",
            "summary_log_scale": summary,
        })
        provenance["colors"]["samples"] = pc.SAMPLE_COLORS
        write_provenance(letter, provenance, outputs)


if __name__ == "__main__":
    main()
