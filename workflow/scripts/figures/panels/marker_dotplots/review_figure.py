"""Leiden review dot plot of one annotation step with explicit gene columns.

The review kit's figure (workflow/scripts/review/review_dotplot.py, render(): clusters grouped by label with
removed clusters in a red EXCLUDE block, genes under group and block headers, QC side columns, text columns,
footer) with two changes for the marker dot plots of release v9:
  - the caller lists every column as (block, group, gene), so a gene the Methods list under two Table M1 rows
    is shown under both, where the kit shows each gene once;
  - a gap of one column with a dark rule can separate the first block from the rest (a lineage's own Table M1
    rows from the other lineages' class genes).
Without a gap and with the kit's footer sentence, the figure is the kit's. Sizes, colours and side columns are
imported from the kit.
"""

import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from review_dotplot import (
    CHAR_W6,
    CHAR_W7,
    CHAR_W_TITLE,
    COL_W,
    DOT_CMAP,
    EXCLUDE_COLOUR,
    GAP,
    GROUP_PALETTE,
    MAX_DOT_AREA,
    ROW_H,
    SIDE_COLUMNS,
    SIDE_GAP,
    SIDE_W,
    TEXT_MAX_CHARS,
    TITLE_SIZE,
    wrap_to_lines,
)


GAP_COLUMNS = 1.0              # width of the gap between the first block and the rest, in gene columns
BLOCK_RULE_COLOUR = "#555555"  # the rule in that gap, darker than the kit's block separators (#999999)
KIT_GENE_NOTE = "Each gene is shown once."


def column_positions(n_columns: int, gap_before) -> np.ndarray:
    """x of each gene column: 0, 1, 2, ..., shifted right by GAP_COLUMNS from column gap_before on."""
    x = np.arange(n_columns)
    if gap_before is not None:
        x = x + np.where(x >= gap_before, GAP_COLUMNS, 0.0)
    return x


def render(rows, gene_stats, columns, label_order, display, title, subtitle, notes=(), text_columns=(),
           gap_before=None, gene_note=KIT_GENE_NOTE):
    """Draw one review dot plot.

    rows: one row per cluster with columns cluster, label, excluded and the kit's SIDE_COLUMNS values.
    gene_stats: long table cluster, gene, pct_detected, mean_log1p_cp10k (annotation_review.marker_stats).
    columns: (block, group, gene) per gene column, left to right; a gene may appear more than once.
    label_order: kept labels in figure order; display: label -> text shown for its row group.
    text_columns: (header, {cluster: text}, colour) shown right of the side columns.
    gap_before: index of the first column after the gap, or None; gene_note: last sentence of the footer.
    Returns the figure and the long table of plotted values (one row per cluster and column).
    """
    rank = {l: i for i, l in enumerate(list(label_order) + ["EXCLUDE"])}
    rows = rows.assign(_rank=rows.label.map(rank).fillna(len(rank)),
                       _num=pd.to_numeric(rows.cluster, errors="coerce"))
    rows = rows.sort_values(["_rank", "_num", "cluster"]).reset_index(drop=True)
    n_rows = len(rows)

    # Values per cluster and gene, scaled per gene over the clusters shown; a repeated column repeats them.
    genes = [g for _, _, g in columns]
    unique_genes = list(dict.fromkeys(genes))
    pct = gene_stats.pivot(index="cluster", columns="gene", values="pct_detected").loc[rows.cluster, unique_genes]
    mean = gene_stats.pivot(index="cluster", columns="gene", values="mean_log1p_cp10k").loc[rows.cluster, unique_genes]
    span = (mean.max() - mean.min()).replace(0, np.nan)
    scaled = ((mean - mean.min()) / span).fillna(0.0)
    pct, mean, scaled = pct[genes], mean[genes], scaled[genes]
    x = column_positions(len(columns), gap_before)
    x_right = x[-1] + 0.5

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

    # Block brackets: a label that would run into the previous one at the same height is raised a level.
    block_spans, start = [], 0
    for j in range(1, len(columns) + 1):
        if j == len(columns) or columns[j][0] != columns[start][0]:
            block_spans.append((start, j - 1, columns[start][0]))
            start = j
    level_end, block_levels = [], []
    for first, last, block in block_spans:
        x_start, level = x[first] - 0.4, 0
        while level < len(level_end) and x_start < level_end[level] + 0.5:
            level += 1
        if level == len(level_end):
            level_end.append(-np.inf)
        level_end[level] = x_start + len(block) * CHAR_W_TITLE / COL_W
        block_levels.append(level)
    n_levels = max(block_levels) + 1

    longest_group = max(len(g) for _, g, _ in columns)
    longest_gene = max(len(g) for g in genes)
    header_h = 0.45 + 0.19 * n_levels + longest_group * CHAR_W_TITLE + 0.12 + longest_gene * CHAR_W6
    title_h, footer_h = 0.72, 1.45
    body_h = n_rows * ROW_H
    left_w = label_w + 0.14 + 0.30
    dots_w = (x[-1] + 1) * COL_W
    side_w = len(SIDE_COLUMNS) * (SIDE_W + SIDE_GAP)
    fig_w = 0.2 + left_w + dots_w + GAP + side_w + sum(GAP + w for w in text_widths) + 0.2
    fig_w = max(fig_w, 7.2)
    fig_h = title_h + header_h + body_h + footer_h
    fig = plt.figure(figsize=(fig_w, fig_h))

    def axes_at(x0, y_top, w, h):
        return fig.add_axes([x0 / fig_w, 1 - (y_top + h) / fig_h, w / fig_w, h / fig_h])

    body_top = title_h + header_h
    x_dots = 0.2 + left_w

    ax = axes_at(x_dots, body_top, dots_w, body_h)
    grid_x, grid_y = np.meshgrid(x, np.arange(n_rows))
    ax.scatter(grid_x.ravel(), grid_y.ravel(), s=MAX_DOT_AREA * pct.values.ravel() / 100.0,
               c=scaled.values.ravel(), cmap=DOT_CMAP, vmin=0, vmax=1, edgecolors="#555555", linewidths=0.2)
    ax.set_xlim(-0.5, x_right)
    ax.set_ylim(n_rows - 0.5, -0.5)
    for side in ["left", "bottom", "top", "right"]:
        ax.spines[side].set_visible(False)
    ax.set_xticks(x)
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
        if j == gap_before:
            ax.axvline((x[j - 1] + x[j]) / 2, color=BLOCK_RULE_COLOUR, lw=0.8, zorder=0)
        elif group_keys[j] != group_keys[j - 1]:
            heavy = columns[j][0] != columns[j - 1][0]
            ax.axvline((x[j - 1] + x[j]) / 2, color="#999999" if heavy else "#DDDDDD", lw=0.6 if heavy else 0.4,
                       zorder=0)

    # Row-group labels: colour block and display name, EXCLUDE in red.
    ax_lab = axes_at(0.2, body_top, label_w + 0.14, body_h)
    ax_lab.set_xlim(0, label_w + 0.14)
    ax_lab.set_ylim(n_rows - 0.5, -0.5)
    ax_lab.axis("off")
    for label, group_rows in rows.groupby("label", sort=False):
        top, bottom = group_rows.index.min() - 0.5, group_rows.index.max() + 0.5
        excluded = label == "EXCLUDE"
        colour = EXCLUDE_COLOUR if excluded else colours.get(label, "#999999")
        ax_lab.add_patch(mpl.patches.Rectangle((label_w + 0.02, top + 0.08), 0.1, bottom - top - 0.16,
                                               color=colour, lw=0))
        ax_lab.text(label_w - 0.04, (top + bottom) / 2, "\n".join(group_lines[label]), ha="right", va="center",
                    fontsize=7, color=EXCLUDE_COLOUR if excluded else "black",
                    fontweight="bold" if excluded else "normal", linespacing=1.05)

    # Header: block brackets (horizontal labels) and group brackets (vertical labels).
    head_h = header_h - longest_gene * CHAR_W6 - 0.08
    ax_head = axes_at(x_dots, title_h, dots_w, head_h)
    ax_head.set_xlim(-0.5, x_right)
    ax_head.set_ylim(0, 1)
    ax_head.axis("off")
    y_group_bracket = 0.04 / head_h
    y_block_bracket = (0.08 + longest_group * CHAR_W_TITLE + 0.06) / head_h
    start = 0
    for j in range(1, len(columns) + 1):
        if j == len(columns) or group_keys[j] != group_keys[start]:
            ax_head.plot([x[start] - 0.35, x[j - 1] + 0.35], [y_group_bracket] * 2, color="black", lw=0.6)
            ax_head.text((x[start] + x[j - 1]) / 2, y_group_bracket + 0.03 / head_h, columns[start][1], rotation=90,
                         ha="center", va="bottom", fontsize=TITLE_SIZE)
            start = j
    for (first, last, block), level in zip(block_spans, block_levels):
        y_line = y_block_bracket + 0.19 * level / head_h
        ax_head.plot([x[first] - 0.4, x[last] + 0.4], [y_line] * 2, color="black", lw=0.8)
        ax_head.text(x[first] - 0.4, y_line + 0.03 / head_h, block, ha="left", va="bottom", fontsize=TITLE_SIZE)

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

    # Text columns (decisions, reasons), one line per row.
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
    footnote = ("Dot size: % of the cluster's nuclei with counts > 0 for the gene. "
                "Colour: mean of log(1 + 10,000 x counts / "
                f"total counts of the nucleus), min-max scaled per gene over the {n_rows} clusters shown. "
                "Side columns: nuclei; "
                "median detected genes and UMIs per nucleus; mean % mitochondrial reads; mean DoubletDetection score; "
                "% of nuclei expressing two or more lineage panels (>= 2 genes each); "
                "share of the cluster's nuclei from its "
                "largest array." + (f" {gene_note}" if gene_note else ""))
    width_chars = int((fig_w - 4.3) / CHAR_W6)
    lines = textwrap.wrap(footnote, width=width_chars) + [l for n in notes for l in textwrap.wrap(n, width=width_chars)]
    fig.text(4.15 / fig_w, 1 - (footer_top + 0.05) / fig_h, "\n".join(lines), ha="left", va="top", fontsize=6,
             color="#333333")

    # One row per cluster and column, in plot order (a left merge keeps it), with the column's block and group.
    plotted = (pd.DataFrame({"cluster": np.repeat(rows.cluster.to_numpy(), len(genes)), "gene": genes * n_rows})
               .assign(pct_detected=pct.to_numpy().ravel(), mean_log1p_cp10k=mean.to_numpy().ravel(),
                       scaled_mean_0_1=scaled.to_numpy().ravel())
               .merge(rows[["cluster", "label", "excluded"]], on="cluster", how="left", validate="many_to_one"))
    plotted["panel"] = [b for b, _, _ in columns] * n_rows
    plotted["panel_group"] = [g for _, g, _ in columns] * n_rows
    return fig, plotted
