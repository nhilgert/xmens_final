"""Table M1 v9 marker genes and the definitions of the marker dot plots of release v9.

Imported by the scripts of rules panels_marker_dotplots_* (workflow/rules/91_panels_marker_dotplots.smk).
Table M1 (resources/figures/table_m1_markers_v9.csv; columns Class, Cell type, Mouse, Human) lists the marker
genes of each cell type per species exactly as the Methods table does; "—" marks a row not used in that
species, and the "*" after some cell-type names is the Methods table's footnote mark. Final-label dot plots
("A" plots) show release labels as rows and, as columns, each row's Table M1 genes in Table M1 order, so a
gene the Methods list under two rows of a plot is shown under both; the row Stromal cells (all) covers two
classes (CLASS_GROUPS), so the major-class plots show its genes before the decidual genes. Each Leiden review
dot plot of a lineage shows the Table M1 rows of its own classes first and then the class-level genes of every
other class (review_columns); the atlas review dot plot shows the class-level genes of every class.
"""

from pathlib import Path

import pandas as pd


# Table M1 ----

NOT_USED = {"—", "-", ""}
FOOTNOTE_MARKS = "†‡§¶*"
SPECIES_COLUMN = {"mouse": "Mouse", "human": "Human"}

# Release label of each Table M1 row: obs cell_type for class-level rows, cell_subtype otherwise (the
# label_in_objects column of the v9 rule check); a row of several classes has a CLASS_GROUPS key instead.
TABLE_M1_LABELS = {
    "Stromal cells (all)": "stromal_cell",
    "Decidual cells": "decidual_cell",
    "Decidual boundary (dec1)": "dec1",
    "Proliferative decidual (dec2)": "dec2",
    "Mature decidual (dec3)": "dec3",
    "Apoptotic decidual (decAp)": "decAp",
    "Endometrial fibroblasts": "endometrial_fibroblast",
    "Outer eFib (eFib1)": "eFib1",
    "Middle eFib (eFib2)": "eFib2",
    "Induced middle eFib (eFib2i)": "eFib2i",
    "Subluminal eFib (eFib3)": "eFib3",
    "Proliferative eFib (eFibP)": "eFibP",
    "Perivascular cells": "perivascular_cell",
    "Smooth muscle cells": "smooth_muscle_cell",
    "Perimetrium mesothelial cells": "perimetrial_cell",
    "Epithelial cells (all)": "epithelial_cell",
    "Ciliated epithelial cells": "epithelial_ciliated",
    "Inner epithelial cells": "epithelial_inner",
    "Outer epithelial cells": "epithelial_outer",
    "Luminal epithelial cells": "epithelial_luminal",
    "Glandular epithelial cells": "epithelial_glandular",
    "Proliferative epithelial cells": "epithelial_proliferative",
    "Endothelial cells (all)": "endothelial_cell",
    "Blood endothelial cells": "endothelial_blood",
    "Lymphatic endothelial cells": "endothelial_lymph",
    "Immune cells (all)": "immune_cell",
    "Macrophages": "macrophage",
    "Cycling macrophages": "cycling_macrophage",
    "Monocytes": "monocyte",
    "Neutrophils": "neutrophil",
    "Dendritic cells": "dendritic_cell",
    "T cells": "t_cell",
    "CD4 T cells": "cd4_t_cell",
    "CD8 T cells": "cd8_t_cell",
    "Regulatory T cells": "treg_cell",
    "γδ T cells": "gamma_delta_t_cell",
    "RORC-positive lymphocytes": "rorc_positive_lymphocyte",
    "B cells": "b_cell",
    "Plasma cells": "plasma_cell",
    "Mast cells": "mast_cell",
    "uNK cells": "unk_cell",
    "Cytotoxic uNK cells": "cytotoxic_unk_cell",
    "Cycling uNK cells": "cycling_unk_cell",
}

# Table M1 rows of a group of classes: "Stromal cells (all)" is every decidual cell and endometrial fibroblast,
# compared with the other classes like the class-level rows Epithelial cells (all) and Immune cells (all). The
# group is no release label: its genes go before the first of its classes in a plot that shows all of them.
CLASS_GROUPS = {
    "stromal_cell": ["decidual_cell", "endometrial_fibroblast"],
}

# Marker rules of Table M1 (v9 rule check): a gene marks its row when it is detected (>= 1 count) in at least
# MIN_PCT % of the row's nuclei and in at least FOLD times the share of the comparison nuclei. Genes that
# are informative but detected in few nuclei of their type follow the sparse rule; the Table M1 caption names them.
STANDARD_RULE = {"name": ">=20% and >=2x", "min_pct": 20.0, "fold": 2.0}
SPARSE_RULE = {"name": ">=5% and >=5x", "min_pct": 5.0, "fold": 5.0}
SPARSE_RULE_GENES = {
    "human": ["IGFBP1", "FOXP3", "CTLA4", "CD40LG"],
    "mouse": ["Xcr1", "Cd3e"],
}

# Perivascular cells and smooth muscle cells share mural genes, so the rule check also tests the perivascular genes
# against smooth muscle cells alone.
SIBLING_CHECKS = {"perivascular_cell": "smooth_muscle_cell"}

# Label shown pending the PI's sign-off (vocabulary v9, approved 30 Oct 2026).
PENDING_SIGN_OFF = {"neutrophil": "label approved 30 Oct 2026 for vocabulary v9; PI sign-off pending"}


def clean_cell_type(name: str) -> str:
    """Table M1 cell-type name without its footnote marks."""
    return name.rstrip(" " + FOOTNOTE_MARKS).strip()


def read_table_m1(path: Path) -> pd.DataFrame:
    """One row per Table M1 row and species with genes: class, cell type, label, species, genes (list)."""
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    rows = []
    for order, record in enumerate(table.itertuples(index=False)):
        cell_type = clean_cell_type(record[1])
        if cell_type not in TABLE_M1_LABELS:
            raise KeyError(f"Table M1 row without a release label: {record[1]!r}")
        for species, column in SPECIES_COLUMN.items():
            cell = getattr(record, column).strip()
            if cell in NOT_USED:
                continue
            rows.append({
                "m1_order": order,
                "m1_class": record.Class,
                "m1_cell_type": cell_type,
                "label": TABLE_M1_LABELS[cell_type],
                "species": species,
                "genes": [gene.strip() for gene in cell.split(",")],
            })
    return pd.DataFrame(rows)


def criterion_of(species: str, gene: str) -> dict:
    return SPARSE_RULE if gene in SPARSE_RULE_GENES[species] else STANDARD_RULE


# Final-label dot plots ----

# Every plot: the obs column of its rows (level), the release labels shown as rows from top to bottom (Table M1
# order, so the dots form a diagonal), and the cell types whose subtypes must all be rows (covers: a v9 label
# missing from a plot is an error, not a silent drop). Columns are the Table M1 genes of the rows (plot_genes).
PLOTS = {
    "mouse": {
        "major_classes": {
            "title": "major cell classes",
            "level": "cell_type",
            "rows": [
                "decidual_cell",
                "endometrial_fibroblast",
                "perivascular_cell",
                "smooth_muscle_cell",
                "perimetrial_cell",
                "epithelial_cell",
                "endothelial_cell",
                "immune_cell",
            ],
            "covers": "all",
        },
        "fib_dec_subtypes": {
            "title": "decidual and endometrial fibroblast subtypes",
            "level": "cell_subtype",
            "rows": ["dec1", "dec2", "dec3", "decAp", "eFib1", "eFib2", "eFib2i", "eFib3", "eFibP"],
            "covers": ["decidual_cell", "endometrial_fibroblast"],
        },
        "epithelial": {
            "title": "epithelial subtypes",
            "level": "cell_subtype",
            "rows": ["epithelial_luminal", "epithelial_glandular", "epithelial_proliferative"],
            "covers": ["epithelial_cell"],
        },
        "endothelial": {
            "title": "endothelial subtypes",
            "level": "cell_subtype",
            "rows": ["endothelial_blood", "endothelial_lymph"],
            "covers": ["endothelial_cell"],
        },
        "smooth_muscle_perivascular": {
            "title": "smooth muscle and perivascular cells",
            "level": "cell_type",
            "rows": ["perivascular_cell", "smooth_muscle_cell"],
            "covers": ["perivascular_cell", "smooth_muscle_cell"],
        },
        "immune_v9": {
            "title": "immune subtypes of release v9",
            "level": "cell_subtype",
            "rows": [
                "macrophage",
                "cycling_macrophage",
                "monocyte",
                "neutrophil",
                "dendritic_cell",
                "t_cell",
                "gamma_delta_t_cell",
                "b_cell",
                "plasma_cell",
                "unk_cell",
                "cytotoxic_unk_cell",
                "cycling_unk_cell",
            ],
            "covers": ["immune_cell"],
        },
        # The cell types of the earlier mouse immune dot plot, without the subtypes new in v9 (neutrophils,
        # plasma cells), with the Methods' genes in Table M1 order like every other plot.
        "immune_as_before": {
            "title": "immune subtypes of the earlier plot (no v9-only subtypes)",
            "level": "cell_subtype",
            "rows": [
                "macrophage",
                "cycling_macrophage",
                "monocyte",
                "dendritic_cell",
                "t_cell",
                "gamma_delta_t_cell",
                "b_cell",
                "unk_cell",
                "cytotoxic_unk_cell",
                "cycling_unk_cell",
            ],
            "covers": None,
        },
    },
    "human": {
        "major_classes": {
            "title": "major cell classes",
            "level": "cell_type",
            "rows": [
                "decidual_cell",
                "endometrial_fibroblast",
                "perivascular_cell",
                "smooth_muscle_cell",
                "epithelial_cell",
                "endothelial_cell",
                "immune_cell",
            ],
            "covers": "all",
        },
        "epithelial": {
            "title": "epithelial subtypes",
            "level": "cell_subtype",
            "rows": ["epithelial_ciliated", "epithelial_inner", "epithelial_outer"],
            "covers": ["epithelial_cell"],
        },
        "endothelial": {
            "title": "endothelial subtypes",
            "level": "cell_subtype",
            "rows": ["endothelial_blood", "endothelial_lymph"],
            "covers": ["endothelial_cell"],
        },
        "smooth_muscle_perivascular": {
            "title": "smooth muscle and perivascular cells",
            "level": "cell_type",
            "rows": ["perivascular_cell", "smooth_muscle_cell"],
            "covers": ["perivascular_cell", "smooth_muscle_cell"],
        },
        "immune_v9": {
            "title": "immune subtypes of release v9",
            "level": "cell_subtype",
            "rows": [
                "macrophage",
                "cd4_t_cell",
                "cd8_t_cell",
                "treg_cell",
                "rorc_positive_lymphocyte",
                "b_cell",
                "mast_cell",
                "unk_cell",
                "cytotoxic_unk_cell",
            ],
            "covers": ["immune_cell"],
        },
        # The cell types of the earlier human immune dot plot, without the subtype new in v9 (B cells), with the
        # Methods' genes in Table M1 order like every other plot. Its gamma-delta T row is the v9 label
        # rorc_positive_lymphocyte.
        "immune_as_before": {
            "title": "immune subtypes of the earlier plot (no v9-only subtypes)",
            "level": "cell_subtype",
            "rows": [
                "macrophage",
                "cd4_t_cell",
                "cd8_t_cell",
                "treg_cell",
                "rorc_positive_lymphocyte",
                "mast_cell",
                "unk_cell",
                "cytotoxic_unk_cell",
            ],
            "covers": None,
        },
    },
}

# Leiden review dot plot of each annotation step -> the A plot whose rows order its cluster groups (and whose
# genes the atlas shows).
REVIEW_PLOT = {
    "mouse": {
        "atlas": "major_classes",
        "stromal": "fib_dec_subtypes",
        "epithelial": "epithelial",
        "vascular": "endothelial",
        "perimetrium": "major_classes",
        "smc": "smooth_muscle_perivascular",
        "immune": "immune_v9",
    },
    "human": {
        "atlas": "major_classes",
        "stromal": "major_classes",
        "epithelial": "epithelial",
        "perivascular": "endothelial",
        "smc": "smooth_muscle_perivascular",
        "immune": "immune_v9",
    },
}

# Atlas lineages in the order of their Table M1 classes, so the review rows follow the class plot's diagonal.
ATLAS_LINEAGE_ORDER = {
    "mouse": ["stromal", "smc", "perimetrium", "epithelial", "vascular", "immune"],
    "human": ["stromal", "smc", "epithelial", "perivascular", "immune"],
}

# Classes (obs cell_type) that each lineage step annotated. Its review dot plot shows the Table M1 rows of these
# classes first, then the class-level genes of every other class, so a cluster carrying another lineage's RNA
# shows it. Human perivascular cells were annotated in the stromal lineage, mouse ones with the smooth muscle
# cells; the human step named "perivascular" annotated the endothelial cells.
REVIEW_OWN_CLASSES = {
    "mouse": {
        "stromal": ["decidual_cell", "endometrial_fibroblast"],
        "epithelial": ["epithelial_cell"],
        "vascular": ["endothelial_cell"],
        "perimetrium": ["perimetrial_cell"],
        "smc": ["perivascular_cell", "smooth_muscle_cell"],
        "immune": ["immune_cell"],
    },
    "human": {
        "stromal": ["decidual_cell", "endometrial_fibroblast", "perivascular_cell"],
        "epithelial": ["epithelial_cell"],
        "perivascular": ["endothelial_cell"],
        "smc": ["smooth_muscle_cell"],
        "immune": ["immune_cell"],
    },
}
OWN_BLOCK = "own lineage: Table M1 rows"
OTHER_BLOCK = "other lineages: Table M1 class genes"


def display_names(vocabulary: dict) -> dict:
    """label -> row name: the vocabulary's display names of cell types and subtypes, the same in both species."""
    return {**vocabulary["display_names"]["cell_type"], **vocabulary["display_names"]["cell_subtype"]}


def table_m1_rows_of(label: str, species: str, table_m1: pd.DataFrame) -> pd.DataFrame:
    """The Table M1 rows that give a label its genes, in Table M1 order."""
    species_rows = table_m1[table_m1.species == species]
    rows = species_rows[species_rows.label == label]
    if rows.empty:
        raise KeyError(f"no Table M1 row for {species} label {label!r}")
    return rows.sort_values("m1_order")


def groups_before(label: str, definition: dict) -> list:
    """Class groups (CLASS_GROUPS) whose genes a class-level plot shows just before the genes of this row.

    A group comes before the first of its classes, and only in a plot that has every one of its classes as a row:
    the major-class plot shows the Stromal cells (all) genes before the decidual and fibroblast genes.
    """
    rows = definition["rows"]
    if definition["level"] != "cell_type":
        return []
    return [
        group for group, classes in CLASS_GROUPS.items()
        if set(classes) <= set(rows) and label == next(r for r in rows if r in classes)
    ]


def plot_genes(species: str, plot: str, table_m1: pd.DataFrame) -> pd.DataFrame:
    """Gene columns of an A plot: gene, Table M1 row and class, row by row as the Methods list them.

    Each row contributes all of its Table M1 genes in Table M1 order, so a gene listed under two rows of the plot
    (mouse Top2a under dec2 and eFibP) is a column under each. A class group's row (Stromal cells (all)) comes
    before its first class (groups_before).
    """
    definition = PLOTS[species][plot]
    columns = []
    for label in definition["rows"]:
        for source in groups_before(label, definition) + [label]:
            for row in table_m1_rows_of(source, species, table_m1).itertuples():
                columns += [{"gene": g, "m1_cell_type": row.m1_cell_type, "m1_class": row.m1_class}
                            for g in row.genes]
    return pd.DataFrame(columns)


def class_of_labels(vocabulary: dict) -> dict:
    """label -> its class (obs cell_type): a class is its own class, a subtype that of its cell type."""
    allowed = vocabulary["allowed"]
    return {**{c: c for c in allowed}, **{s: c for c, subtypes in allowed.items() for s in subtypes}}


def classes_of_row(label: str, class_of: dict) -> set:
    """Classes (obs cell_type) a Table M1 row covers: a group's classes, else the class of its label."""
    return set(CLASS_GROUPS[label]) if label in CLASS_GROUPS else {class_of[label]}


def review_columns(species: str, step: str, table_m1: pd.DataFrame, vocabulary: dict) -> pd.DataFrame:
    """Gene columns of a step's review dot plot: block, group, gene, m1_cell_type, m1_class.

    The atlas shows the major-class plot's genes, blocked by Table M1 class and grouped by row. A lineage shows
    first the Table M1 rows of its own classes (REVIEW_OWN_CLASSES) in Table M1 order, each row's genes as the
    Methods list them (Stromal cells (all) with the decidual and fibroblast classes, Endothelial cells (all) with
    the endothelial class), then the major-class plot's genes of every other class, one group per class-level
    Table M1 row (Stromal cells (all) included), each gene once and none repeating an own-lineage gene.
    """
    if step == "atlas":
        columns = plot_genes(species, "major_classes", table_m1)
        return columns.assign(block=columns.m1_class, group=columns.m1_cell_type)

    own_classes = set(REVIEW_OWN_CLASSES[species][step])
    class_of = class_of_labels(vocabulary)
    species_rows = table_m1[table_m1.species == species].sort_values("m1_order")
    row_classes = {label: classes_of_row(label, class_of) for label in species_rows.label}
    for label, classes in row_classes.items():
        if classes & own_classes and not classes <= own_classes:
            raise ValueError(f"{species} {step}: Table M1 row {label} spans own and other classes")
    own_rows = species_rows[[bool(row_classes[label] & own_classes) for label in species_rows.label]]
    columns = [
        {"block": OWN_BLOCK, "group": row.m1_cell_type, "gene": g, "m1_cell_type": row.m1_cell_type,
         "m1_class": row.m1_class}
        for row in own_rows.itertuples()
        for g in row.genes
    ]
    shown = {column["gene"] for column in columns}

    # Other lineages: a class row is named as its class (vocabulary), a group row as its Table M1 row.
    class_names = vocabulary["display_names"]["cell_type"]
    m1_label = dict(zip(species_rows.m1_cell_type, species_rows.label))
    for column in plot_genes(species, "major_classes", table_m1).itertuples():
        label = m1_label[column.m1_cell_type]
        if row_classes[label] & own_classes or column.gene in shown:
            continue
        columns.append({"block": OTHER_BLOCK, "group": class_names.get(label, column.m1_cell_type),
                        "gene": column.gene, "m1_cell_type": column.m1_cell_type, "m1_class": column.m1_class})
        shown.add(column.gene)
    return pd.DataFrame(columns)


def check_review_classes(species: str, step: str, kept_labels: set, vocabulary: dict) -> None:
    """Fail when a lineage step kept a label outside its own classes, or none of one of them."""
    if step == "atlas":
        return
    kept_classes = {class_of_labels(vocabulary)[label] for label in kept_labels}
    own_classes = set(REVIEW_OWN_CLASSES[species][step])
    if kept_classes != own_classes:
        raise ValueError(f"{species} {step}: kept labels' classes {sorted(kept_classes)} "
                         f"differ from {sorted(own_classes)}")


def check_rows_cover(species: str, plot: str, cell_types: pd.Series, cell_subtypes: pd.Series) -> None:
    """Fail when a plot misses a release label it should show, or shows one the release lacks."""
    definition = PLOTS[species][plot]
    level_values = cell_types if definition["level"] == "cell_type" else cell_subtypes
    absent = sorted(set(definition["rows"]) - set(level_values))
    if absent:
        raise ValueError(f"{species} {plot}: rows not in the release: {absent}")
    if definition["covers"] is None:
        return
    if definition["covers"] == "all":
        expected = set(level_values)
    else:
        expected = set(level_values[cell_types.isin(definition["covers"]).to_numpy()])
    missing = sorted(expected - set(definition["rows"]))
    if missing:
        raise ValueError(f"{species} {plot}: release labels without a row: {missing}")
