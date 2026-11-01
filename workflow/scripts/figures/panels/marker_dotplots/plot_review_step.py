"""Draw the Leiden review dot plot of one annotation step on the Table M1 v9 genes of the Methods.

Called by rule panels_marker_dotplots_review_panel (workflow/rules/91_panels_marker_dotplots.smk), once per
species and step. Reads <tables-dir>/<species>_<step>_review_clusters.csv, _review_marker_stats.csv and
_review_extract_provenance.json, Table M1 and the annotation vocabulary. Writes <out-prefix>.svg, .png
(300 dpi) and .provenance.json. plot_review_pdf.py draws the same pages into each species' reading copy.

The layout is the review kit's (workflow/scripts/review/review_dotplot.py, as in the all-steps dot plots, drawn
by review_figure.render): clusters grouped by label with removed clusters in a red EXCLUDE block, QC side
columns, then each cluster's decision and, for removed clusters, the recorded reason. Columns
(marker_panels.review_columns): the atlas shows the class-level genes of every class; a lineage shows the Table
M1 rows of its own classes, each row's genes as the Methods list them, and after a gap the class-level genes of
every other class, so a cluster removed for another lineage's RNA shows that lineage's genes.
"""

import argparse
import json
import logging
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import yaml

import dotplot_style as style
import marker_panels as panels
from review_dotplot import EXCLUDE_COLOUR, FIGURE_RC
from review_figure import KIT_GENE_NOTE, render
from review_step import ATLAS_DISPLAY


DECISION_COLOUR = "#1F4E79"   # dark blue, apart from the red of EXCLUDE and the black of the data labels

# Resolution named in the title: the atlas is clustered at Leiden 2.5, every lineage at 3.0 (lineage_spec).
LEIDEN_LEVEL = {"atlas": "Leiden 2.5"}

# The review kit's sizes and layout, with the plots' font family and pinned SVG ids. Active while a page is
# drawn and saved: savefig reads the SVG and PDF settings.
REVIEW_RC = {**FIGURE_RC, **style.FONT_RC}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, choices=list(panels.REVIEW_PLOT))
    parser.add_argument("--step", required=True)
    parser.add_argument("--tables-dir", required=True, type=Path)
    parser.add_argument("--table-m1", required=True, type=Path)
    parser.add_argument("--vocabulary", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    return parser.parse_args()


def render_columns(columns: pd.DataFrame) -> list:
    """render()'s columns: (block, group, gene) per gene column, left to right."""
    return list(zip(columns.block, columns.group, columns.gene))


def row_order_and_names(species: str, step: str, labels: set, vocabulary: dict, table_m1: pd.DataFrame):
    """Kept labels in the order of the matching marker plot's rows, and the name shown for each."""
    if step == "atlas":
        order = panels.ATLAS_LINEAGE_ORDER[species]
        names = ATLAS_DISPLAY
    else:
        kept = labels - {"EXCLUDE"}
        order = panels.plot_rows(species, panels.REVIEW_PLOT[species][step], kept, table_m1, vocabulary)
        names = panels.display_names(vocabulary)
    unknown = sorted(labels - set(order) - {"EXCLUDE"})
    if unknown:
        raise ValueError(f"{species} {step}: labels outside the rows of its marker plot: {unknown}")
    return [label for label in order if label in labels], names


def decision_text(row) -> str:
    """'assigned', 'relabelled (v8 label)' or 'EXCLUDE', with a pending sign-off where flagged."""
    if row.decision == "relabelled":
        text = f"relabelled ({row.v8_reference_text})"
    else:
        text = row.decision
    if row.pending_sign_off:
        text += "; PI sign-off pending"
    return text


def footnotes(step: str, extract: dict, is_round2: bool, sparse_genes: list, own_names: list, missing: list) -> list:
    """Lines under the plot: the decision provenance and how to read the decision column and the genes."""
    reference = (
        "the label most of the cluster's nuclei carried in release v8 (round-2 re-fit: new clusters)"
        if is_round2
        else "the cluster's label in the v8 table of this clustering"
    )
    table_m1 = f"Table M1 {panels.table_m1_version(extract['table_m1'])}"
    if step == "atlas":
        genes = (f"Genes: {table_m1} rows of the {extract['marker_plot'].replace('_', ' ')} marker dot plot, "
                 "each gene once; ")
    else:
        genes = (
            f"Genes: {table_m1} (Methods). Own lineage ({', '.join(own_names)}): every Table M1 row of "
            f"{'these classes' if len(own_names) > 1 else 'this class'} in Table M1 order, with each row's genes as "
            "listed (a gene listed under two rows is shown under both). Other lineages: the genes of the class-level "
            "Table M1 rows of every other class, each gene once and none repeating an own-lineage gene. "
        )
    rules = (
        f"sparse-rule markers {', '.join(sparse_genes)} ({panels.SPARSE_RULE['name']}); other genes follow "
        if sparse_genes else "all follow "
    )
    # The atlas keeps its sentence; a lineage's gene note ends a sentence, so its rule text starts one.
    if step != "atlas":
        rules = rules[0].upper() + rules[1:]
    notes = [
        f"Decisions: {extract['decided_at']} ({extract['approval_date']}, {extract['approver']}); "
        f"{extract['decision_record']}.",
        f"Decision column: assigned = same label as {reference}; relabelled = a different label (v8 label shown); "
        "EXCLUDE = removed, with the recorded reason [source].",
        genes + rules + f"{panels.STANDARD_RULE['name']}. Missing from this object: {', '.join(missing) or 'none'}.",
    ]
    if extract["pending_sign_off"]:
        notes.append("PI sign-off pending: " + "; ".join(extract["pending_sign_off"]) + ".")
    return notes


def load_step(tables_dir: Path, species: str, step: str, table_m1: pd.DataFrame, vocabulary: dict) -> dict:
    """One step's cluster table, marker statistics, extract provenance and gene columns."""
    stem = tables_dir / f"{species}_{step}"
    clusters = pd.read_csv(f"{stem}_review_clusters.csv", dtype={"cluster": str}, keep_default_na=False)
    clusters["excluded"] = clusters.excluded.astype(str).eq("True")
    clusters["pending_sign_off"] = clusters.pending_sign_off.astype(str).eq("True")
    extract = json.loads(Path(f"{stem}_review_extract_provenance.json").read_text())
    columns = panels.review_columns(species, step, table_m1, vocabulary)
    missing = list(dict.fromkeys(g for g in columns.gene if g in extract["genes_missing_from_var"]))
    return {
        "clusters": clusters,
        "stats": pd.read_csv(f"{stem}_review_marker_stats.csv", dtype={"cluster": str}),
        "extract": extract,
        "columns": columns[~columns.gene.isin(missing)].reset_index(drop=True),
        "missing": missing,
        "is_round2": extract["fit"].startswith("round 2"),
        "inputs": [Path(f"{stem}_review_clusters.csv"), Path(f"{stem}_review_marker_stats.csv")],
        "table_m1": table_m1,
    }


def draw_step(species: str, step: str, data: dict, vocabulary: dict):
    """The review figure of one step; call inside mpl.rc_context(REVIEW_RC).

    Returns the figure, the long table of plotted values, the row groups shown and the number of rotated
    labels anchored.
    """
    clusters, extract, columns = data["clusters"], data["extract"], data["columns"]
    order, names = row_order_and_names(species, step, set(clusters.label), vocabulary, data["table_m1"])
    kept, removed = clusters[~clusters.excluded], clusters[clusters.excluded]
    panels.check_review_classes(species, step, set(kept.label), vocabulary)

    # Annotations ----

    level = LEIDEN_LEVEL.get(step, "Leiden 3.0")
    step_name = "whole-tissue atlas" if step == "atlas" else step + (", round 2" if data["is_round2"] else "")
    title = (
        f"{species.capitalize()} {step_name} ({level}): {len(kept)} clusters kept ({kept.n_nuclei.sum():,} nuclei) "
        f"as {kept.label.nunique()} label{'s' if kept.label.nunique() != 1 else ''}; {len(removed)} removed "
        f"({removed.n_nuclei.sum():,} nuclei); Table M1 {panels.table_m1_version(extract['table_m1'])} markers"
    )
    own_names = []
    if step == "atlas":
        genes_text = f"{extract['table_m1']} ({extract['marker_plot']})"
    else:
        own_names = [vocabulary["display_names"]["cell_type"][c] for c in panels.REVIEW_OWN_CLASSES[species][step]]
        genes_text = f"{extract['table_m1']} (own lineage: {', '.join(own_names)}; other lineages: class genes)"
    subtitle = (
        f"Clusters: {extract['object']} ({extract['cluster_column']}); "
        f"decisions ({extract['version']}): {extract['decisions']}; genes: {genes_text}"
    )
    decisions = {row.cluster: decision_text(row) for row in clusters.itertuples()}
    reasons = {
        row.cluster: f"{row.exclusion_reason} [{row.reason_source}]"
        for row in clusters.itertuples() if row.excluded
    }
    text_columns = [("decision", decisions, DECISION_COLOUR)]

    # A step that removed nothing gets no empty reason column.
    if reasons:
        text_columns.append(("exclusion reason [source]", reasons, EXCLUDE_COLOUR))

    # Sparse-rule genes are named in the footnote rather than marked on the plot (as in the Table M1 caption).
    sparse_genes = list(dict.fromkeys(g for g in columns.gene if g in panels.SPARSE_RULE_GENES[species]))

    # Figure ----

    # A lineage's own genes and the other lineages' genes are separated by a gap; the atlas is the kit's figure.
    other = columns.index[columns.block == panels.OTHER_BLOCK]
    fig, plotted = render(
        clusters,
        data["stats"],
        render_columns(columns),
        order,
        names,
        title,
        subtitle,
        footnotes(step, extract, data["is_round2"], sparse_genes, own_names, data["missing"]),
        text_columns,
        gap_before=int(other[0]) if len(other) else None,
        gene_note=KIT_GENE_NOTE if step == "atlas" else None,
    )

    n_anchored = style.anchor_rotated_text(fig)
    row_groups = [{"label": label, "shown_as": names.get(label, label)} for label in order] + [{"label": "EXCLUDE"}]
    return fig, plotted, row_groups, n_anchored


def main() -> None:
    args = parse_arguments()
    table_m1 = panels.read_table_m1(args.table_m1)
    vocabulary = yaml.safe_load(args.vocabulary.read_text())
    data = load_step(args.tables_dir, args.species, args.step, table_m1, vocabulary)
    clusters, extract, columns = data["clusters"], data["extract"], data["columns"]

    with mpl.rc_context(REVIEW_RC):
        fig, plotted, row_groups, n_anchored = draw_step(args.species, args.step, data, vocabulary)
        outputs = style.save_svg_png(fig, args.out_prefix)
    plt.close(fig)
    logging.info("%s %s: %d clusters x %d gene columns", args.species, args.step, len(clusters), len(columns))

    # Provenance ----

    provenance = {
        "plot": f"review_{args.step}",
        "species": args.species,
        "version": extract["version"],
        # The release whose labels these decisions set, as the step's extract records it.
        "release": extract["release"],
        "release_dir": extract["release_dir"],
        "vocabulary_schema": extract["vocabulary_schema"],
        "step": args.step,
        "fit": extract["fit"],
        "object": extract["object"],
        "object_sha256": extract["object_sha256"],
        "cluster_column": extract["cluster_column"],
        "decisions": extract["decisions"],
        "decisions_sha256": extract["decisions_sha256"],
        "v8_reference": extract["v8_reference"],
        "marker_plot": extract["marker_plot"],
        "genes": [
            {"gene": g, "block": b, "group": grp, "rule": panels.criterion_of(args.species, g)["name"],
             "table_m1_row": r, "table_m1_class": c}
            for g, b, grp, r, c in zip(columns.gene, columns.block, columns.group, columns.m1_cell_type,
                                       columns.m1_class)
        ],
        "own_classes": panels.REVIEW_OWN_CLASSES[args.species].get(args.step, "all (atlas)"),
        "genes_missing_from_var": data["missing"],
        "row_groups": row_groups,
        "clusters": [
            {"cluster": row.cluster, "label": row.label, "decision": row.decision, "n_nuclei": int(row.n_nuclei),
             "v8_reference": row.v8_reference_text, "pending_sign_off": bool(row.pending_sign_off)}
            for row in clusters.itertuples()
        ],
        "clusters_by_decision": extract["clusters_by_decision"],
        "pending_sign_off": extract["pending_sign_off"],
        "dot_size": "% of the cluster's nuclei with >= 1 count of the gene (review kit scale)",
        "dot_color": "mean log1p CP10k, min-max scaled 0-1 per gene across the clusters shown (review kit, viridis_r)",
        "rotated_labels_anchored": n_anchored,
        "outputs": [str(p) for p in outputs],
        "inputs": {str(p): style.sha256_of(p) for p in data["inputs"] + [args.table_m1, args.vocabulary]},
        "plotted_values": plotted.round(4).to_dict(orient="records"),
        "script_sha256": style.sha256_of(Path(__file__)),
        "marker_panels_sha256": style.sha256_of(Path(panels.__file__)),
        "style_sha256": style.sha256_of(Path(style.__file__)),
    }
    style.write_provenance(Path(f"{args.out_prefix}.provenance.json"), provenance)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
