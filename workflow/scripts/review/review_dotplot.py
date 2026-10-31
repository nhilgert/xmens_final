"""Review dot plot of one annotation step: the layout of the 29 Oct all-steps dot plots, with more QC.

Rows are the step's clusters grouped by label (removed clusters last, in a red EXCLUDE block); columns are
marker genes grouped by panel. Dot size = % of the cluster's nuclei detecting the gene; colour = mean
log1p CP10k, min-max scaled per gene over the clusters shown. Side columns (bars with values): nuclei,
median genes, median UMIs, mean % mitochondrial reads, mean doublet score, % of nuclei expressing two
lineage panels and the share of the largest array. Text columns on the right carry flags, verdicts and exclusion reasons.
"""
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROW_H = 0.155          # inches per cluster row
COL_W = 0.155          # inches per gene column
SIDE_W = 0.56          # each QC side column
SIDE_GAP = 0.07
GAP = 0.10
CHAR_W6 = 0.047        # approximate width of a 6 pt character
CHAR_W7 = 0.055
TITLE_SIZE = 8.5
CHAR_W_TITLE = 0.066
MAX_DOT_AREA = 80.0
TEXT_MAX_CHARS = 200   # text columns are one line per row; longer texts end in "..."
EXCLUDE_COLOUR = "#C0392B"   # red is reserved for EXCLUDE (and for flag text)
DOT_CMAP = "viridis_r"
GROUP_PALETTE = ["#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
                 "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#C5B0D5", "#C49C94"]
SIDE_COLUMNS = [("n_nuclei", "nuclei", "{:,.0f}"),
                ("median_genes", "median\ngenes", "{:,.0f}"),
                ("median_umis", "median\nUMIs", "{:,.0f}"),
                ("mean_pct_mito", "mean\n% mito", "{:.2f}"),
                ("mean_doublet", "mean\ndoublet score", "{:.1f}"),
                ("pct_two_lineages", "% two\nlineages", "{:.0f}"),
                ("top_sample_share", "largest\narray share", "{:.0%}")]
FIGURE_RC = {
    "axes.labelsize": 8.0, "axes.linewidth": 0.6, "axes.spines.right": False, "axes.spines.top": False,
    "figure.dpi": 200.0, "font.size": 8.0, "legend.fontsize": 7.0, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "tight", "savefig.dpi": 300.0, "svg.fonttype": "none",
    "xtick.labelsize": 6.0, "ytick.labelsize": 6.0, "xtick.major.size": 3.0, "ytick.major.size": 3.0,
}


def wrap_to_lines(text, max_lines, min_chars=24):
    """Wrap at >= min_chars per line, widening until the text needs at most max_lines lines."""
    for width in range(min(min_chars, len(text)), len(text) + 1):
        lines = textwrap.wrap(text, width=width, break_long_words=False)
        if len(lines) <= max_lines:
            return lines
    return [text]


def panel_columns(panels, available):
    """(category, group, gene) for every panel gene present, each gene once, in panel order."""
    columns, seen = [], set()
    for panel in panels:
        for group in panel["groups"]:
            for gene in group["genes"]:
                if gene in available and gene not in seen:
                    columns.append((panel["category"], group["label"], gene))
                    seen.add(gene)
    return columns


def render(rows, gene_stats, panels, label_order, display, title, subtitle, notes=(), text_columns=()):
    """Draw one review dot plot.

    rows: one row per cluster with columns cluster, label, excluded and the SIDE_COLUMNS values.
    gene_stats: long table cluster, gene, pct_detected, mean_log1p_cp10k (annotation_review.marker_stats).
    panels: list of {"category", "groups": [{"label", "genes"}]}.
    label_order: kept labels in figure order; display: label -> text shown for its row group.
    text_columns: (header, {cluster: text}, colour) shown right of the side columns.
    Returns the figure and the long table of plotted values.
    """
    rank = {l: i for i, l in enumerate(list(label_order) + ["EXCLUDE"])}
    rows = rows.assign(_rank=rows.label.map(rank).fillna(len(rank)),
                       _num=pd.to_numeric(rows.cluster, errors="coerce"))
    rows = rows.sort_values(["_rank", "_num", "cluster"]).reset_index(drop=True)
    n_rows = len(rows)

    columns = panel_columns(panels, set(gene_stats.gene))
    genes = [g for _, _, g in columns]
    pct = gene_stats.pivot(index="cluster", columns="gene", values="pct_detected").loc[rows.cluster, genes]
    mean = gene_stats.pivot(index="cluster", columns="gene", values="mean_log1p_cp10k").loc[rows.cluster, genes]
    span = (mean.max() - mean.min()).replace(0, np.nan)
    scaled = ((mean - mean.min()) / span).fillna(0.0)

    group_text = {l: display.get(l, l) for l in label_order}
    group_text["EXCLUDE"] = "EXCLUDE"
    colours = {l: GROUP_PALETTE[i % len(GROUP_PALETTE)] for i, l in enumerate(label_order)}
    group_lines = {}
    for label, group_rows in rows.groupby("label", sort=False):
        group_lines[label] = wrap_to_lines(group_text.get(label, label), max(1, int(len(group_rows) * ROW_H / 0.105)))
    label_w = max(max(len(line) for line in lines) for lines in group_lines.values()) * CHAR_W7 + 0.15

    text_widths = []
    for header, texts, _ in text_columns:
        longest = max([len(t) for t in texts.values()] + [len(header)]) if texts else len(header)
        text_widths.append(min(longest, TEXT_MAX_CHARS) * CHAR_W6 + 0.2)

    # Category brackets: a label that would run into the previous one at the same height is raised a level.
    category_spans, start = [], 0
    for j in range(1, len(columns) + 1):
        if j == len(columns) or columns[j][0] != columns[start][0]:
            category_spans.append((start, j - 1, columns[start][0]))
            start = j
    level_end, category_levels = [], []
    for first, last, category in category_spans:
        x_start, level = first - 0.4, 0
        while level < len(level_end) and x_start < level_end[level] + 0.5:
            level += 1
        if level == len(level_end):
            level_end.append(-np.inf)
        level_end[level] = x_start + len(category) * CHAR_W_TITLE / COL_W
        category_levels.append(level)
    n_levels = max(category_levels) + 1

    longest_group = max(len(g) for _, g, _ in columns)
    longest_gene = max(len(g) for g in genes)
    header_h = 0.45 + 0.19 * n_levels + longest_group * CHAR_W_TITLE + 0.12 + longest_gene * CHAR_W6
    title_h, footer_h = 0.72, 1.45
    body_h = n_rows * ROW_H
    left_w = label_w + 0.14 + 0.30
    dots_w = len(genes) * COL_W
    side_w = len(SIDE_COLUMNS) * (SIDE_W + SIDE_GAP)
    fig_w = 0.2 + left_w + dots_w + GAP + side_w + sum(GAP + w for w in text_widths) + 0.2
    fig_w = max(fig_w, 7.2)
    fig_h = title_h + header_h + body_h + footer_h
    fig = plt.figure(figsize=(fig_w, fig_h))

    def axes_at(x, y_top, w, h):
        return fig.add_axes([x / fig_w, 1 - (y_top + h) / fig_h, w / fig_w, h / fig_h])

    body_top = title_h + header_h
    x_dots = 0.2 + left_w

    ax = axes_at(x_dots, body_top, dots_w, body_h)
    grid_x, grid_y = np.meshgrid(np.arange(len(genes)), np.arange(n_rows))
    ax.scatter(grid_x.ravel(), grid_y.ravel(), s=MAX_DOT_AREA * pct.values.ravel() / 100.0,
               c=scaled.values.ravel(), cmap=DOT_CMAP, vmin=0, vmax=1, edgecolors="#555555", linewidths=0.2)
    ax.set_xlim(-0.5, len(genes) - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    for side in ["left", "bottom", "top", "right"]:
        ax.spines[side].set_visible(False)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontstyle="italic")
    ax.xaxis.tick_top()
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(rows.cluster)
    ax.tick_params(length=0, pad=2)
    for i in range(1, n_rows):
        if rows.label[i] != rows.label[i - 1]:
            heavy = rows.label[i] == "EXCLUDE"
            ax.axhline(i - 0.5, color=EXCLUDE_COLOUR if heavy else "#BBBBBB", lw=0.8 if heavy else 0.4)
    group_keys = [(c, g) for c, g, _ in columns]
    for j in range(1, len(columns)):
        if group_keys[j] != group_keys[j - 1]:
            heavy = columns[j][0] != columns[j - 1][0]
            ax.axvline(j - 0.5, color="#999999" if heavy else "#DDDDDD", lw=0.6 if heavy else 0.4, zorder=0)

    # Row-group labels: colour block and display name, EXCLUDE in red.
    ax_lab = axes_at(0.2, body_top, label_w + 0.14, body_h)
    ax_lab.set_xlim(0, label_w + 0.14)
    ax_lab.set_ylim(n_rows - 0.5, -0.5)
    ax_lab.axis("off")
    for label, group_rows in rows.groupby("label", sort=False):
        top, bottom = group_rows.index.min() - 0.5, group_rows.index.max() + 0.5
        excluded = label == "EXCLUDE"
        colour = EXCLUDE_COLOUR if excluded else colours.get(label, "#999999")
        ax_lab.add_patch(mpl.patches.Rectangle((label_w + 0.02, top + 0.08), 0.1, bottom - top - 0.16, color=colour, lw=0))
        ax_lab.text(label_w - 0.04, (top + bottom) / 2, "\n".join(group_lines[label]), ha="right", va="center",
                    fontsize=7, color=EXCLUDE_COLOUR if excluded else "black",
                    fontweight="bold" if excluded else "normal", linespacing=1.05)

    # Header: category brackets (horizontal labels) and group brackets (vertical labels).
    head_h = header_h - longest_gene * CHAR_W6 - 0.08
    ax_head = axes_at(x_dots, title_h, dots_w, head_h)
    ax_head.set_xlim(-0.5, len(genes) - 0.5)
    ax_head.set_ylim(0, 1)
    ax_head.axis("off")
    y_group_bracket = 0.04 / head_h
    y_category_bracket = (0.08 + longest_group * CHAR_W_TITLE + 0.06) / head_h
    start = 0
    for j in range(1, len(columns) + 1):
        if j == len(columns) or group_keys[j] != group_keys[start]:
            ax_head.plot([start - 0.35, j - 1 + 0.35], [y_group_bracket] * 2, color="black", lw=0.6)
            ax_head.text((start + j - 1) / 2, y_group_bracket + 0.03 / head_h, columns[start][1], rotation=90,
                         ha="center", va="bottom", fontsize=TITLE_SIZE)
            start = j
    for (first, last, category), level in zip(category_spans, category_levels):
        y_line = y_category_bracket + 0.19 * level / head_h
        ax_head.plot([first - 0.4, last + 0.4], [y_line] * 2, color="black", lw=0.8)
        ax_head.text(first - 0.4, y_line + 0.03 / head_h, category, ha="left", va="bottom", fontsize=TITLE_SIZE)

    # QC side columns: bars with values, removed clusters in red.
    x_side = x_dots + dots_w + GAP
    for k, (column, header, fmt) in enumerate(SIDE_COLUMNS):
        ax_side = axes_at(x_side + k * (SIDE_W + SIDE_GAP), body_top, SIDE_W, body_h)
        values = rows[column].to_numpy(dtype=float)
        bar_colours = [EXCLUDE_COLOUR if e else "#B0B0B0" for e in rows.excluded]
        ax_side.barh(np.arange(n_rows), values, height=0.62, color=bar_colours, alpha=0.55, lw=0)
        top_value = np.nanmax(values) if np.nanmax(values) > 0 else 1.0
        ax_side.set_xlim(0, top_value * 1.9)
        ax_side.set_ylim(n_rows - 0.5, -0.5)
        for i, v in enumerate(values):
            ax_side.text(v + top_value * 0.03, i, fmt.format(v), va="center", ha="left", fontsize=6)
        ax_side.axis("off")
        ax_side.text(0, -0.9, header, ha="left", va="bottom", fontsize=6)

    # Text columns (flags, verdicts, reasons), one line per row.
    x_text = x_side + side_w
    for (header, texts, colour), width in zip(text_columns, text_widths):
        x_text += GAP
        ax_text = axes_at(x_text, body_top, width, body_h)
        ax_text.set_xlim(0, 1)
        ax_text.set_ylim(n_rows - 0.5, -0.5)
        ax_text.axis("off")
        ax_text.text(0, -0.9, header, ha="left", va="bottom", fontsize=6, color=colour)
        for i, cluster in enumerate(rows.cluster):
            if texts.get(cluster):
                ax_text.text(0, i, textwrap.shorten(texts[cluster], width=TEXT_MAX_CHARS, placeholder="..."),
                             ha="left", va="center", fontsize=6, color=colour)
        x_text += width

    fig.text(0.2 / fig_w, 1 - 0.24 / fig_h, title, ha="left", va="center", fontsize=TITLE_SIZE)
    fig.text(0.2 / fig_w, 1 - 0.50 / fig_h, subtitle, ha="left", va="center", fontsize=6, color="#444444")

    # Footer: dot-size key, colour bar and how to read the plot.
    footer_top = body_top + body_h + 0.25
    ax_size = axes_at(0.2, footer_top + 0.15, 1.9, 0.35)
    ax_size.axis("off")
    for k, p in enumerate([10, 25, 50, 75, 100]):
        ax_size.scatter(k, 0, s=MAX_DOT_AREA * p / 100, color="#777777", edgecolors="#555555", linewidths=0.2)
        ax_size.text(k, -0.9, f"{p}%", ha="center", va="top", fontsize=6)
    ax_size.set_xlim(-0.6, 4.6)
    ax_size.set_ylim(-1.6, 0.8)
    ax_size.text(-0.6, 1.0, "nuclei with >= 1 count", ha="left", va="bottom", fontsize=6)
    ax_bar = axes_at(2.35, footer_top + 0.22, 1.5, 0.09)
    mpl.colorbar.ColorbarBase(ax_bar, cmap=mpl.colormaps[DOT_CMAP], orientation="horizontal",
                              norm=mpl.colors.Normalize(0, 1), ticks=[0, 0.5, 1])
    ax_bar.tick_params(labelsize=6, length=2)
    ax_bar.set_title("mean log1p(CP10k), scaled 0-1 per gene\nacross the clusters of this plot", fontsize=6, loc="left")
    footnote = ("Dot size: % of the cluster's nuclei with counts > 0 for the gene. Colour: mean of log(1 + 10,000 x counts / "
                f"total counts of the nucleus), min-max scaled per gene over the {n_rows} clusters shown. Side columns: nuclei; "
                "median detected genes and UMIs per nucleus; mean % mitochondrial reads; mean DoubletDetection score; "
                "% of nuclei expressing two or more lineage panels (>= 2 genes each); share of the cluster's nuclei from its "
                "largest array. Each gene is shown once.")
    width_chars = int((fig_w - 4.3) / CHAR_W6)
    lines = textwrap.wrap(footnote, width=width_chars) + [l for n in notes for l in textwrap.wrap(n, width=width_chars)]
    fig.text(4.15 / fig_w, 1 - (footer_top + 0.05) / fig_h, "\n".join(lines), ha="left", va="top", fontsize=6, color="#333333")

    plotted = (pd.DataFrame({"cluster": np.repeat(rows.cluster.to_numpy(), len(genes)), "gene": genes * n_rows})
               .assign(pct_detected=pct.to_numpy().ravel(), mean_log1p_cp10k=mean.to_numpy().ravel(),
                       scaled_mean_0_1=scaled.to_numpy().ravel())
               .merge(rows[["cluster", "label", "excluded"]], on="cluster"))
    plotted["panel"] = plotted.gene.map({g: c for c, _, g in columns})
    plotted["panel_group"] = plotted.gene.map({g: grp for _, grp, g in columns})
    return fig, plotted
