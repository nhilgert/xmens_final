"""Draw a mouse UMAP with direct labels: cell classes (Fig. 3A) or immune subtypes (Fig. 3C).

The immune subtypes are drawn either on the immune lineage fit's own embedding (--level immune,
panel 3C) or where the immune nuclei sit in the whole-tissue embedding, zoomed to its immune
islands (--level immune_atlas, panel 3Catlas).

Called by rule panels_mouse_tissue_panel (workflow/rules/91_panels_mouse_tissue.smk) for 3A, 3C and
3Catlas.
Reads cells.csv.gz and extract_provenance.json from --tables-dir.
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

# Label directions set by hand (degrees counter-clockwise from east, multiples of 22.5) where the
# automatic layout reads ambiguously; chosen by eye on the rendered panel. The v9 immune embedding
# packs five uNK and T islands against each other, so its whole layout is given: lymphoid labels
# around the left cluster, myeloid labels around the right one.
LABEL_ANGLE_OVERRIDES = {
    ("v9", "immune"): {
        "cytotoxic_unk_cell": 112.5,
        "unk_cell": 67.5,
        "t_cell": 0,
        "gamma_delta_t_cell": -90,
        "cycling_unk_cell": -112.5,
        "plasma_cell": -90,
        "b_cell": 67.5,
        "neutrophil": 112.5,
        "cycling_macrophage": 90,
        "macrophage": 0,
        "monocyte": -45,
        "dendritic_cell": -90,
    },
}

# 3A draws the whole-tissue embedding of the release (~58,000 nuclei, so small points); 3C draws
# the immune lineage fit's embedding of the released immune nuclei (~2,100, so larger points);
# 3Catlas draws the same immune nuclei at their whole-tissue coordinates, so the axes close in on
# the atlas's immune islands.
# label_orders: placing orders the label planner tries; the 12 crowded immune labels need more.
LEVELS = {
    "class": {
        "panel": "3A",
        "coords": ["umap_1", "umap_2"],
        "group": "cell_type",
        "order": style.CLASS_ORDER,
        "labels": style.CLASS_LABELS,
        "colors": style.CLASS_COLORS,
        "point_size": 0.25,
        "label_orders": 2,
        "immune_only": False,
        "embedding": "release obsm X_umap_scVI (whole-tissue fit), as stored",
    },
    "immune": {
        "panel": "3C",
        "coords": ["immune_umap_1", "immune_umap_2"],
        "group": "cell_subtype",
        "order": style.IMMUNE_ORDER,
        "labels": style.IMMUNE_LABELS,
        "colors": style.IMMUNE_COLORS,
        "point_size": 2.0,
        "label_orders": 24,
        "immune_only": True,
        "embedding": "immune lineage fit obsm X_umap_scVI, as stored",
    },
    "immune_atlas": {
        "panel": "3Catlas",
        "coords": ["umap_1", "umap_2"],
        "group": "cell_subtype",
        "order": style.IMMUNE_ORDER,
        "labels": style.IMMUNE_LABELS,
        "colors": style.IMMUNE_COLORS,
        "point_size": 2.0,
        "label_orders": 24,
        "immune_only": True,
        "embedding": "release obsm X_umap_scVI (whole-tissue fit), as stored; immune nuclei only",
    },
}


def main() -> None:
    parser = style.panel_arguments(__doc__)
    parser.add_argument("--level", required=True, choices=list(LEVELS))
    args = parser.parse_args()
    level = LEVELS[args.level]
    style.apply_panel_style()

    # Load ----

    cells = style.read_cells(args.tables_dir)
    if level["immune_only"]:
        cells = cells[cells["cell_type"] == "immune_cell"]
    unknown = set(cells[level["group"]]) - set(level["order"])
    if unknown:
        raise ValueError(f"{level['group']} values without a colour: {sorted(unknown)}")

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
    angles_deg = LABEL_ANGLE_OVERRIDES.get((args.version, args.level), {})
    placements = style.place_direct_labels(
        ax,
        xy,
        groups,
        level["order"],
        level["labels"],
        angles_deg,
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
        "colors_source": "manuscript Fig. 3 (classes new in v9: Fig. 5 or new, see style module)",
        "new_in_v9": [g for g in style.NEW_IN_V9 if g in set(groups)],
        "pending_sign_off": {g: v for g, v in style.PENDING_SIGN_OFF.items() if g in set(groups)},
        "labels": placements,
    })
    style.save_panel(fig, args.out_prefix, provenance)
    plt.close(fig)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
