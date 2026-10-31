"""Draw per-stage violins of log counts (Fig. S12C) or log detected genes (Fig. S12D) per nucleus.

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for S12C, S12D.
Reads cells.csv.gz and extract_provenance.json from --tables-dir.
Writes <out-prefix>.pdf, .svg, .png and .provenance.json.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np

import mouse_tissue_style as style


# Half of a 180 mm page, as in the manuscript's S12 layout.
FIG_SIZE_IN = (3.6, 2.0)

# Natural log of (1 + value), the scale of the manuscript's axes (counts 6-11, genes 5.5-9).
METRICS = {
    "total_counts": {"panel": "S12C", "label": "counts per cell (log)"},
    "n_genes": {"panel": "S12D", "label": "genes per cell (log)"},
}
VIOLIN_EDGE = "#333333"


def main() -> None:
    parser = style.panel_arguments(__doc__)
    parser.add_argument("--metric", required=True, choices=list(METRICS))
    args = parser.parse_args()
    metric = METRICS[args.metric]
    style.apply_panel_style()

    # Load ----

    # One value per released nucleus (embedding outliers included: this is library QC).
    cells = style.read_cells(args.tables_dir)
    values = {stage: np.log1p(cells.loc[cells["timepoint"] == stage, args.metric].to_numpy()) for stage in style.STAGE_ORDER}
    empty = [stage for stage, v in values.items() if len(v) == 0]
    if empty:
        raise ValueError(f"stages without nuclei: {empty}")

    # Figure ----

    fig, ax = plt.subplots(figsize=FIG_SIZE_IN)
    positions = np.arange(len(style.STAGE_ORDER))

    # Kernel density violins spanning each stage's range, filled with the stage colour.
    parts = ax.violinplot(
        [values[stage] for stage in style.STAGE_ORDER],
        positions=positions,
        widths=0.8,
        showextrema=False,
    )
    for body, stage in zip(parts["bodies"], style.STAGE_ORDER):
        body.set_facecolor(style.STAGE_COLORS[stage])
        body.set_edgecolor(VIOLIN_EDGE)
        body.set_linewidth(0.6)
        body.set_alpha(1.0)

    # Labels ----

    ax.set_xticks(positions, labels=[style.VIOLIN_STAGE_LABELS[s] for s in style.STAGE_ORDER])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(positions) - 0.4)
    ax.set_ylabel(metric["label"], fontweight="bold")

    # Save ----

    provenance = style.base_provenance(
        metric["panel"],
        args,
        inputs=[args.tables_dir / "cells.csv.gz"],
        script=__file__,
    )
    raw = {stage: cells.loc[cells["timepoint"] == stage, args.metric] for stage in style.STAGE_ORDER}
    provenance.update({
        "metric": f"log1p({args.metric}), {args.metric} summed from the release's layers/counts",
        "nuclei": "all released nuclei of each stage",
        "kde": "matplotlib violinplot defaults (Gaussian KDE, Scott bandwidth, cut at the data range)",
        "per_stage": {
            stage: {
                "n": int(len(raw[stage])),
                "median": float(raw[stage].median()),
                "min": float(raw[stage].min()),
                "max": float(raw[stage].max()),
                "median_log1p": round(float(np.median(values[stage])), 4),
            }
            for stage in style.STAGE_ORDER
        },
        "colors": {s: style.STAGE_COLORS[s] for s in style.STAGE_ORDER},
        "colors_source": "manuscript Fig. S12C-D; the release's uns['timepoint_colors'] are matplotlib tab10 defaults",
    })
    style.save_panel(fig, args.out_prefix, provenance)
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
