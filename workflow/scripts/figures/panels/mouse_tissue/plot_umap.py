"""Draw a mouse UMAP with direct labels: cell classes (Fig. 3A) or immune subtypes (Fig. 3C).

The immune subtypes are drawn either on the immune lineage fit's own embedding (--level immune,
panel 3C) or where the immune nuclei sit in the whole-tissue embedding, zoomed to its immune
islands (--level immune_atlas, panel 3Catlas).

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for 3A, 3C and
3Catlas.
Reads cells.csv.gz and extract_provenance.json from --tables-dir; the immune subtypes' order, names,
colors and hand-set label directions from --vocabulary and --immune-display.
Writes <out-prefix>.pdf, .svg, .png and .provenance.json.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np

import mouse_tissue_style as style


# About a quarter of a 180 mm page wide, as panels 3A and 3C are in the manuscript layout; 3C a
# little wider for its longer subtype names. The height follows the embedding's aspect (equal axis scales), within these bounds, so a wide
# embedding gives a short panel instead of empty bands.
FIG_WIDTH_IN = {"class": 2.0, "immune": 2.3, "immune_atlas": 2.3}
FIG_HEIGHT_BOUNDS_IN = (1.5, 2.2)

# Room around the points for outside labels, as a fraction of the data range: wider for the
# immune panel, whose subtype names are longer and whose islands are smaller.
LABEL_MARGIN = {"class": 0.12, "immune": 0.2, "immune_atlas": 0.2}

# 3A draws the whole-tissue embedding of the release (~58,000 nuclei, so small points); 3C draws
# the immune lineage fit's embedding of the released immune nuclei (~2,100, so larger points);
# 3Catlas draws the same immune nuclei at their whole-tissue coordinates, so the axes close in on
# the atlas's immune islands.
# label_orders: placing orders the label planner tries; the 12 crowded immune labels need more.
# The immune levels take order, names, colors and label directions from the display table
# (with_display); label directions set by hand exist for 3C (display level "immune") only.
LEVELS = {
    "class": {
        "panel": "3A",
        "coords": ["umap_1", "umap_2"],
        "group": "cell_type",
        "point_size": 0.25,
        "label_orders": 2,
        "immune_only": False,
        "embedding": "release obsm X_umap_scVI (whole-tissue fit), as stored",
    },
    "immune": {
        "panel": "3C",
        "coords": ["immune_umap_1", "immune_umap_2"],
        "group": "cell_subtype",
        "point_size": 2.0,
        "label_orders": 24,
        "immune_only": True,
        "embedding": "immune lineage fit obsm X_umap_scVI, as stored",
    },
    "immune_atlas": {
        "panel": "3Catlas",
        "coords": ["umap_1", "umap_2"],
        "group": "cell_subtype",
        "point_size": 2.0,
        "label_orders": 24,
        "immune_only": True,
        "embedding": "release obsm X_umap_scVI (whole-tissue fit), as stored; immune nuclei only",
    },
}


def with_display(level_name: str, groups, args) -> dict:
    """The level's settings plus the order, names, colors and label directions of the groups present.

    Classes use the manuscript's palette of the style module (no hand-set directions); immune subtypes the
    display table (style.immune_display). Also returns what the provenance records about the display.
    """
    level = dict(LEVELS[level_name])
    if level["group"] == "cell_type":
        unknown = set(groups) - set(style.CLASS_ORDER)
        if unknown:
            raise ValueError(f"cell_type values without a color: {sorted(unknown)}")
        level.update(order=style.CLASS_ORDER, labels=style.CLASS_LABELS, colors=style.CLASS_COLORS, angles={},
                     record={})
        return level

    display = style.immune_display(groups, args.vocabulary, args.immune_display)
    level.update(
        order=display["order"],
        labels=display["labels"],
        colors=display["colors"],
        angles=display["angles"].get(level_name, {}),
        display=display,
        record={
            "pending_sign_off": display["pending_sign_off"],
            "display_defaults": display["display_defaults"],
        },
    )
    return level


def main() -> None:
    parser = style.panel_arguments(__doc__)
    parser.add_argument("--level", required=True, choices=list(LEVELS))
    args = parser.parse_args()
    style.apply_panel_style()

    # Load ----

    cells = style.read_cells(args.tables_dir)
    if LEVELS[args.level]["immune_only"]:
        cells = cells[cells["cell_type"] == "immune_cell"]
    level = with_display(args.level, cells[LEVELS[args.level]["group"]], args)

    # Flagged embedding outliers sit among another lineage's nuclei; the release rule
    # (uns embedding_outlier_rule) keeps them in every count and leaves them off embeddings.
    is_outlier = cells["embedding_outlier"].astype(bool)
    drawn = cells[~is_outlier]
    logging.info("Drawing %d nuclei; %d embedding outliers hidden", len(drawn), int(is_outlier.sum()))

    # Seeded shuffle so no group is painted over the others wholesale where islands touch.
    rng = np.random.default_rng(args.seed)
    drawn = drawn.iloc[rng.permutation(len(drawn))]
    xy = drawn[level["coords"]].to_numpy()
    groups = drawn[level["group"]].to_numpy()

    # Figure ----

    span = xy.max(axis=0) - xy.min(axis=0)
    width_in = FIG_WIDTH_IN[args.level]
    height_in = float(np.clip(width_in * span[1] / span[0], *FIG_HEIGHT_BOUNDS_IN))
    fig, ax = plt.subplots(figsize=(width_in, height_in))

    # One point per nucleus, coloured by group; rasterised so the vector files stay small.
    ax.scatter(
        xy[:, 0],
        xy[:, 1],
        c=[level["colors"][g] for g in groups],
        s=level["point_size"],
        linewidths=0,
        rasterized=True,
    )

    # UMAP axes carry no units, so no ticks or frame, as in the manuscript. Equal aspect keeps
    # distances comparable in both directions; the margins leave room for outside labels.
    ax.set_axis_off()
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(LABEL_MARGIN[args.level])
    fig.canvas.draw()

    # Direct labels instead of a legend, as in the manuscript.
    placements = style.place_direct_labels(
        ax,
        xy,
        groups,
        level["order"],
        level["labels"],
        level["angles"],
        n_orders=level["label_orders"],
        seed=args.seed,
    )

    # Save ----

    provenance = style.base_provenance(
        level["panel"],
        args,
        inputs=[args.tables_dir / "cells.csv.gz"],
        script=__file__,
    )
    provenance.update({
        "embedding": level["embedding"],
        "nuclei": "released immune nuclei only" if level["immune_only"] else "all released nuclei",
        "filters": {"embedding_outlier_hidden": int(is_outlier.sum())},
        "n_drawn": int(len(drawn)),
        "n_drawn_per_group": drawn[level["group"]].value_counts().to_dict(),
        "colors": {g: level["colors"][g] for g in level["order"] if g in set(groups)},
        "colors_source": "manuscript Fig. 3 (style module); immune subtypes: the --immune-display table",
        "immune_display": str(args.immune_display),
        "introduced": style.introduced_labels(groups, args.immune_display),
        **level["record"],
        "labels": placements,
    })
    style.save_panel(fig, args.out_prefix, provenance)
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
