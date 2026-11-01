"""Draw relative abundance per stage as stacked bars: cell classes (Fig. 3B) or immune subtypes (Fig. 3D).

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for 3B and 3D.
Reads cells.csv.gz and extract_provenance.json from --tables-dir; the immune subtypes' order and colors
from --vocabulary and --immune-display.
Writes <out-prefix>.pdf, .svg, .png and .provenance.json.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import mouse_tissue_style as style


# Same footprint as the neighbouring UMAP panels of Fig. 3.
FIG_SIZE_IN = (2.1, 1.8)

# Bar thickness as a fraction of the row pitch; the white gap separates stages as in the manuscript.
BAR_HEIGHT = 0.78

# The class palette is the style module's; the immune subtypes' order and colors come from the display
# table (style.immune_display).
LEVELS = {
    "class": {
        "panel": "3B",
        "group": "cell_type",
        "title": "relative cell type abundance",
    },
    "immune": {
        "panel": "3D",
        "group": "cell_subtype",
        "title": "relative immune cell abundance",
    },
}


def main() -> None:
    parser = style.panel_arguments(__doc__)
    parser.add_argument("--level", required=True, choices=list(LEVELS))
    args = parser.parse_args()
    level = dict(LEVELS[args.level])
    style.apply_panel_style()

    # Load and reshape ----

    cells = style.read_cells(args.tables_dir)
    record = {}
    if args.level == "immune":
        cells = cells[cells["cell_type"] == "immune_cell"]
        display = style.immune_display(cells[level["group"]], args.vocabulary, args.immune_display)
        level.update(order=display["order"], colors=display["colors"])
        record = {"pending_sign_off": display["pending_sign_off"], "display_defaults": display["display_defaults"]}
    else:
        level.update(order=style.CLASS_ORDER, colors=style.CLASS_COLORS)
        unknown = set(cells[level["group"]]) - set(level["order"])
        if unknown:
            raise ValueError(f"{level['group']} values without a color: {sorted(unknown)}")

    # One bar per stage: the share of all released nuclei of that stage (embedding outliers
    # included; the flag only governs embedding and spatial plots).
    counts = pd.crosstab(cells["timepoint"], cells[level["group"]])
    groups = [g for g in level["order"] if g in counts.columns]
    counts = counts.reindex(index=style.STAGE_ORDER, columns=groups, fill_value=0)
    percent = counts.div(counts.sum(axis=1), axis=0) * 100
    logging.info("Stages %s; nuclei per stage %s", list(counts.index), counts.sum(axis=1).tolist())

    # Figure ----

    fig, ax = plt.subplots(figsize=FIG_SIZE_IN)
    y = np.arange(len(style.STAGE_ORDER))

    # Segments stacked left to right in the manuscript's group order, so a colour sits at the
    # same place in every bar.
    left = np.zeros(len(y))
    for group in groups:
        ax.barh(
            y,
            percent[group].to_numpy(),
            left=left,
            height=BAR_HEIGHT,
            color=level["colors"][group],
            linewidth=0,
        )
        left += percent[group].to_numpy()

    # Scales and labels ----

    # Shares of 100%; the earliest stage on top, reading down the time course.
    ax.set_xlim(0, 100)
    ax.set_xticks(np.arange(0, 101, 20))
    ax.set_ylim(len(y) - 0.5, -0.5)
    ax.set_yticks(y, labels=[style.BAR_STAGE_LABELS[s] for s in style.STAGE_ORDER])

    # Stage names slanted as in the manuscript, so they fit beside bars this thin.
    for tick in ax.get_yticklabels():
        tick.set_rotation(30)
        tick.set_ha("right")
        tick.set_rotation_mode("anchor")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_title(level["title"], loc="center", fontweight="bold")

    # Save ----

    provenance = style.base_provenance(
        level["panel"],
        args,
        inputs=[args.tables_dir / "cells.csv.gz"],
        script=__file__,
    )
    provenance.update({
        "nuclei": "released immune nuclei" if args.level == "immune" else "all released nuclei",
        "filters": {"embedding_outlier_hidden": 0, "note": "embedding outliers counted"},
        "stack_order": groups,
        "counts": {stage: counts.loc[stage].to_dict() for stage in style.STAGE_ORDER},
        "percent": {stage: percent.loc[stage].round(2).to_dict() for stage in style.STAGE_ORDER},
        "n_per_stage": counts.sum(axis=1).to_dict(),
        "colors": {g: level["colors"][g] for g in groups},
        "immune_display": str(args.immune_display),
        "introduced": style.introduced_labels(groups, args.immune_display),
        **record,
    })
    style.save_panel(fig, args.out_prefix, provenance)
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
