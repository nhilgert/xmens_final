"""Marker statistics of every release label on the Table M1 v9 genes, and the Table M1 rule check on the release.

Called by rule panels_marker_dotplots_extract (workflow/rules/91_panels_marker_dotplots.smk), once per species.
Reads the release's All_Cells object (obs, var and the counts layer; sha256 checked against the release
MANIFEST.json) and Table M1. Writes <out-dir>/<species>_marker_stats.csv (one row per label and gene),
<species>_table_m1_rule_check.csv (one row per Table M1 row and gene) and <species>_extract_provenance.json.
"""

import argparse
import json
import logging
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import marker_panels as panels
from dotplot_style import sha256_of


TARGET_SUM = 1e4      # counts per nucleus after normalisation (log1p CP10k, as the release X)
N_X_CHECK = 2000      # nuclei whose X is compared with log1p CP10k of their counts
LEVELS = ["cell_type", "cell_subtype"]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, choices=list(panels.SPECIES_COLUMN))
    parser.add_argument("--release-dir", required=True, type=Path)
    parser.add_argument("--version", required=True, help="version tag of the release (config releases)")
    parser.add_argument("--object-name", required=True, help="file name of the All_Cells object in the release")
    parser.add_argument("--table-m1", required=True, type=Path)
    parser.add_argument("--plots", required=True, nargs="+", help="the plot names the rule file renders")
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def read_release(path: Path):
    """obs, gene names and the raw counts (CSR) of a release object; X and embeddings are not read."""
    with h5py.File(path, "r") as handle:
        obs = ad.io.read_elem(handle["obs"])
        genes = pd.Index(ad.io.read_elem(handle["var"]).index.astype(str))
        counts = sp.csr_matrix(ad.io.read_elem(handle["layers/counts"]))
        x_head = read_csr_rows(handle["X"], N_X_CHECK, len(genes))
    return obs, genes, counts, x_head


def read_csr_rows(group, n_rows: int, n_cols: int) -> sp.csr_matrix:
    """The first n_rows rows of a CSR matrix stored in an h5ad group, without reading the rest."""
    indptr = group["indptr"][: n_rows + 1]
    data = group["data"][indptr[0]:indptr[-1]]
    indices = group["indices"][indptr[0]:indptr[-1]]
    return sp.csr_matrix((data, indices, indptr - indptr[0]), shape=(n_rows, n_cols))


def log1p_cp10k(counts: sp.csr_matrix, library: np.ndarray) -> sp.csr_matrix:
    scaled = sp.diags(TARGET_SUM / np.maximum(library, 1)) @ counts.astype(np.float64)
    scaled.data = np.log1p(scaled.data)
    return scaled.tocsr()


def label_statistics(obs: pd.DataFrame, detected: np.ndarray, expression: np.ndarray, genes: list) -> pd.DataFrame:
    """Per label of each level and gene: nuclei, nuclei detecting the gene, % detecting and mean log1p CP10k."""
    rows = []
    for level in LEVELS:
        labels = obs[level].astype(str).to_numpy()
        for label in sorted(set(labels)):
            member = labels == label
            n_detecting = detected[member].sum(axis=0)
            mean_expression = expression[member].mean(axis=0)
            for j, gene in enumerate(genes):
                rows.append({
                    "level": level,
                    "label": label,
                    "n_nuclei": int(member.sum()),
                    "gene": gene,
                    "n_detecting": int(n_detecting[j]),
                    "pct_detected": 100.0 * n_detecting[j] / member.sum(),
                    "mean_log1p_cp10k": float(mean_expression[j]),
                })
    return pd.DataFrame(rows)


def label_level(label: str, obs: pd.DataFrame) -> str:
    """cell_type when the label is a class of the release or a group of classes, else cell_subtype."""
    classes = set(obs["cell_type"].astype(str))
    return "cell_type" if label in classes or label in panels.CLASS_GROUPS else "cell_subtype"


def row_nuclei(label: str, obs: pd.DataFrame) -> np.ndarray:
    """Nuclei of a Table M1 row: those of its label, or of every class of a group of classes."""
    if label in panels.CLASS_GROUPS:
        return np.isin(obs["cell_type"].astype(str).to_numpy(), panels.CLASS_GROUPS[label])
    return obs[label_level(label, obs)].astype(str).to_numpy() == label


def row_label_text(label: str) -> str:
    """The release label(s) of a Table M1 row as the rule check names them."""
    return " + ".join(panels.CLASS_GROUPS[label]) if label in panels.CLASS_GROUPS else label


def comparison_nuclei(label: str, obs: pd.DataFrame, class_types: dict):
    """Nuclei a Table M1 row is compared with, and how the rule check names them.

    A subtype is compared with the other subtypes of its cell type; a class-level type with the other types of
    its Table M1 class (decidual cells vs endometrial fibroblasts), or with every other class when it is alone;
    a group of classes (Stromal cells (all)) with every class outside the group.
    """
    cell_type = obs["cell_type"].astype(str).to_numpy()
    if label in panels.CLASS_GROUPS:
        return ~np.isin(cell_type, panels.CLASS_GROUPS[label]), "other classes"
    if label_level(label, obs) == "cell_subtype":
        subtype = obs["cell_subtype"].astype(str).to_numpy()
        own_type = pd.Series(cell_type[subtype == label]).unique()
        if len(own_type) != 1:
            raise ValueError(f"subtype {label} spans cell types {list(own_type)}")
        return (cell_type == own_type[0]) & (subtype != label), "rest of class"
    same_class = [t for t in class_types[label] if t != label]
    if same_class:
        return np.isin(cell_type, same_class), "rest of class"
    return cell_type != label, "other classes"


def meets(pct_type: float, pct_other: float, criterion: dict) -> bool:
    return pct_type >= criterion["min_pct"] and pct_type >= criterion["fold"] * pct_other


def rule_check(species: str, obs: pd.DataFrame, detected: np.ndarray, present: list,
               table_m1: pd.DataFrame) -> pd.DataFrame:
    """Each Table M1 gene of the species against its rule, on every nucleus of the release.

    detected: nuclei x present genes, True where the nucleus has at least one count of the gene.
    """
    column_of = {gene: j for j, gene in enumerate(present)}
    species_rows = table_m1[table_m1.species == species]
    class_level = species_rows[species_rows.label.isin(set(obs["cell_type"].astype(str)))]
    class_types = {
        label: sorted(set(class_level.label[class_level.m1_class == m1_class]))
        for label, m1_class in zip(class_level.label, class_level.m1_class)
    }

    rows = []
    for row in species_rows.itertuples():
        level = label_level(row.label, obs)
        member = row_nuclei(row.label, obs)
        other, comparison = comparison_nuclei(row.label, obs, class_types)
        sibling = panels.SIBLING_CHECKS.get(row.label)
        sibling_member = obs["cell_type"].astype(str).to_numpy() == sibling if sibling else None
        for gene in row.genes:
            criterion = panels.criterion_of(species, gene)
            record = {
                "m1_class": row.m1_class,
                "m1_cell_type": row.m1_cell_type,
                "label": row_label_text(row.label),
                "level": level,
                "gene": gene,
                "in_var": gene in column_of,
                "criterion": criterion["name"],
                "sparse_rule_marker": gene in panels.SPARSE_RULE_GENES[species],
                "comparison_nuclei": comparison,
                "n_type": int(member.sum()),
                "n_comparison": int(other.sum()),
            }
            if not record["in_var"]:
                rows.append(record)
                continue
            in_nucleus = detected[:, column_of[gene]]
            record["n_type_detecting"] = int(in_nucleus[member].sum())
            record["pct_detected_type"] = 100.0 * in_nucleus[member].mean()
            record["n_comparison_detecting"] = int(in_nucleus[other].sum())
            record["pct_detected_comparison"] = 100.0 * in_nucleus[other].mean()
            record["meets_criterion"] = meets(record["pct_detected_type"], record["pct_detected_comparison"], criterion)
            if sibling:
                record["sibling_label"] = sibling
                record["n_sibling"] = int(sibling_member.sum())
                record["pct_detected_sibling"] = 100.0 * in_nucleus[sibling_member].mean()
                record["meets_criterion_vs_sibling"] = meets(
                    record["pct_detected_type"], record["pct_detected_sibling"], criterion
                )
            rows.append(record)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_arguments()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    object_path = args.release_dir / args.object_name

    # Inputs ----

    # The release is immutable: the object must be the one its MANIFEST names.
    manifest = json.loads((args.release_dir / "MANIFEST.json").read_text())
    entry = next(e for e in manifest["objects"].values() if e["file"] == args.object_name)
    object_sha256 = sha256_of(object_path)
    if object_sha256 != entry["sha256"]:
        raise ValueError(f"{object_path}: sha256 {object_sha256} differs from MANIFEST {entry['sha256']}")

    table_m1 = panels.read_table_m1(args.table_m1)
    obs, gene_index, counts, x_head = read_release(object_path)
    logging.info("%s: %d nuclei x %d genes", object_path, *counts.shape)

    # The rule file and marker_panels.PLOTS must name the same plots, and every plot must show exactly the
    # labels the release has.
    if set(args.plots) != set(panels.PLOTS[args.species]):
        raise ValueError(f"plots of the rule file {sorted(args.plots)} differ from marker_panels "
                         f"{sorted(panels.PLOTS[args.species])}")
    for plot in panels.PLOTS[args.species]:
        panels.check_rows_cover(args.species, plot, obs["cell_type"].astype(str), obs["cell_subtype"].astype(str))

    # Statistics ----

    # Library size over all genes, as in the release normalisation; detection = at least one count.
    library = np.asarray(counts.sum(axis=1)).ravel()
    m1_genes = list(dict.fromkeys(g for genes in table_m1.genes[table_m1.species == args.species] for g in genes))
    present = [g for g in m1_genes if g in gene_index]
    missing = [g for g in m1_genes if g not in gene_index]
    selected = counts[:, gene_index.get_indexer(present)]
    expression = log1p_cp10k(selected, library).toarray()
    detected = selected.toarray() > 0
    stats = label_statistics(obs, detected, expression, present)

    # The plotted mean equals the mean of the release X (what scanpy averages) when X is log1p CP10k of the
    # counts; checked on the first nuclei.
    x_difference = np.abs(log1p_cp10k(counts[:N_X_CHECK], library[:N_X_CHECK]) - x_head)
    max_x_difference = float(x_difference.max()) if x_difference.nnz else 0.0

    checks = rule_check(args.species, obs, detected, present, table_m1)
    failing = checks[(checks.in_var) & (~checks.meets_criterion.astype(bool))]
    logging.info("Table M1 genes: %d; missing from var: %s; failing their rule: %d",
                 len(m1_genes), missing, len(failing))

    # Outputs ----

    stats.to_csv(args.out_dir / f"{args.species}_marker_stats.csv", index=False)
    checks.to_csv(args.out_dir / f"{args.species}_table_m1_rule_check.csv", index=False)
    provenance = {
        "species": args.species,
        "version": args.version,
        "release": manifest["release"],
        "vocabulary_schema": manifest["vocabulary_schema"],
        "object": str(object_path),
        "object_sha256": object_sha256,
        "object_sha256_matches_manifest": True,
        "n_nuclei": int(counts.shape[0]),
        "n_embedding_outliers_included": int(obs["embedding_outlier"].astype(bool).sum()),
        "table_m1": str(args.table_m1),
        "table_m1_sha256": sha256_of(args.table_m1),
        "table_m1_genes": m1_genes,
        "table_m1_genes_missing_from_var": missing,
        "table_m1_genes_failing_rule": failing[["m1_cell_type", "gene", "criterion"]].to_dict(orient="records"),
        "detection": "counts >= 1 in layers['counts']",
        "expression": "mean over a label's nuclei of log(1 + 10,000 x counts / total counts of the nucleus)",
        "x_vs_log1p_cp10k_max_abs_difference": max_x_difference,
        "x_check_nuclei": N_X_CHECK,
        # Relative to the working directory (the repository), so the record is the same wherever it is checked out.
        "script": os.path.relpath(__file__),
        "script_sha256": sha256_of(Path(__file__)),
        "marker_panels_sha256": sha256_of(Path(panels.__file__)),
    }
    (args.out_dir / f"{args.species}_extract_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
