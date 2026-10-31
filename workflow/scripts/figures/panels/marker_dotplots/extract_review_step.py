"""Cluster table and marker statistics of one annotation step, on every Table M1 gene of the species.

Called by rule panels_marker_dotplots_review_extract (workflow/rules/91_panels_marker_dotplots.smk), once per
species and step of the review bundle (99b_review_notebooks.smk: the step's clustering, Leiden column and v9
decision table). Reads the clustering (obs, var and the counts layer), the decision table, the v8 table of the
same clustering (atlas and July lineage fits) or the v8 release (round-2 re-fits), the recorded exclusion
reasons and the decision provenance. Writes <out-dir>/<species>_<step>_review_clusters.csv (one row per
cluster: decision, QC side columns), _review_marker_stats.csv (one row per cluster and gene) and
_review_extract_provenance.json, which also names the release whose labels the decisions set (folder, and the
release and vocabulary schema of its MANIFEST.json).
"""

import argparse
import json
import logging
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import marker_panels as panels
from annotation_review import cluster_order, marker_stats, nucleus_qc
from dotplot_style import sha256_of
from review_ledger import NOT_RELEASED, read_obs, v8_status


# How a decision table flags a label awaiting the PI (as in the review notebooks).
PENDING_MARK = "PI sign-off pending"
EXTENDED_ONLY = " (extended only)"

# v9 renamed a v8 label without changing its clusters; a renamed label is not a relabel.
V8_RENAMED = {"human": {"gamma_delta_t_cell": "rorc_positive_lymphocyte"}, "mouse": {}}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, choices=list(panels.REVIEW_PLOT))
    parser.add_argument("--step", required=True)
    parser.add_argument("--object", required=True, type=Path, help="the step's clustering (h5ad)")
    parser.add_argument("--cluster-column", required=True)
    parser.add_argument("--decisions", required=True, type=Path, help="v9 decision table (cluster, label[, note])")
    parser.add_argument("--v8-decisions", required=True, help="v8 table of the same clustering, or None (round 2)")
    parser.add_argument("--v8-release", required=True, type=Path, help="v8 extended all-cells object")
    parser.add_argument("--reasons", required=True, type=Path)
    parser.add_argument("--steps", required=True, type=Path)
    parser.add_argument("--table-m1", required=True, type=Path)
    parser.add_argument("--release-dir", required=True, type=Path, help="release folder of these decisions (v9)")
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def read_clustering(path: Path) -> ad.AnnData:
    """obs, genes and raw counts (as X) of a step's clustering; the normalised X and embeddings are not read."""
    with h5py.File(path, "r") as handle:
        obs = ad.io.read_elem(handle["obs"])
        var = pd.DataFrame(index=ad.io.read_elem(handle["var"]).index.astype(str))
        counts = sp.csr_matrix(ad.io.read_elem(handle["layers/counts"]))
    return ad.AnnData(X=counts, obs=obs, var=var)


def read_decisions(path: Path, clusters: set) -> pd.DataFrame:
    """The decision table, checked: one row per cluster of the clustering, nothing else."""
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    table["cluster"] = table["cluster"].str.strip()
    table["label"] = table["label"].str.strip()
    repeated = table.cluster[table.cluster.duplicated()].tolist()
    if set(table.cluster) != clusters or repeated:
        raise ValueError(
            f"{path}: clusters without a decision {sorted(clusters - set(table.cluster))}; rows for absent "
            f"clusters {sorted(set(table.cluster) - clusters)}; repeated rows {repeated}"
        )
    if "note" not in table:
        table["note"] = ""
    return table.set_index("cluster")


def cluster_qc(qc: pd.DataFrame, cluster: np.ndarray) -> pd.DataFrame:
    """The review dot plot's side columns per cluster, defined as in annotation_review.cluster_table."""
    by_cluster = qc.assign(cluster=cluster, two_lineages=qc.n_lineages >= 2).groupby("cluster")
    table = pd.DataFrame({
        "n_nuclei": by_cluster.size(),
        "median_genes": by_cluster.n_genes.median(),
        "median_umis": by_cluster.umis.median(),
        "mean_pct_mito": by_cluster.pct_mito.mean(),
        "mean_doublet": by_cluster.doublet_score.mean(),
        "pct_two_lineages": 100.0 * by_cluster.two_lineages.mean(),
        "top_sample": by_cluster["sample"].agg(lambda s: s.value_counts().index[0]),
        "top_sample_share": by_cluster["sample"].agg(lambda s: s.value_counts(normalize=True).iloc[0]),
    })
    return table.loc[cluster_order(table.index)]


def v8_reference_from_release(species: str, obs: pd.DataFrame, cluster: np.ndarray,
                              release_obs: pd.DataFrame) -> pd.DataFrame:
    """Per cluster of a round-2 re-fit: the label most of its nuclei carried in v8, with its share.

    Supplementary v8 populations count under their label ("(extended only)" is kept in the text); v8 names
    that v9 renamed count under the v9 name.
    """
    status = pd.Series(v8_status(obs, release_obs).to_numpy(), index=obs.index)
    rows = []
    for name in cluster_order(set(cluster)):
        shares = status[cluster == name].value_counts(normalize=True)
        top = sorted(shares.index, key=lambda label: (-shares[label], label))[0]
        label = top.removesuffix(EXTENDED_ONLY)
        rows.append({
            "cluster": name,
            "v8_reference": V8_RENAMED[species].get(label, label),
            "v8_reference_text": f"{top} {shares[top]:.0%}",
            "v8_released_share": float((status[cluster == name] != NOT_RELEASED).mean()),
        })
    return pd.DataFrame(rows).set_index("cluster")


def v8_reference_from_table(path: Path, clusters: list) -> pd.DataFrame:
    """Per cluster: its label in the v8 table of the same clustering."""
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    label_of = dict(zip(table.cluster.str.strip(), table.label.str.strip()))
    return pd.DataFrame(
        {"v8_reference": [label_of.get(c, "absent") for c in clusters],
         "v8_reference_text": [f"v8 table: {label_of.get(c, 'absent')}" for c in clusters]},
        index=pd.Index(clusters, name="cluster"),
    )


def decision_of(label: str, v8_label: str) -> str:
    if label == "EXCLUDE":
        return "EXCLUDE"
    return "assigned" if label == v8_label else "relabelled"


def main() -> None:
    args = parse_arguments()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.out_dir / f"{args.species}_{args.step}"

    # Load ----

    adata = read_clustering(args.object)
    if args.cluster_column not in adata.obs:
        raise KeyError(f"{args.object} has no column {args.cluster_column}")
    cluster = adata.obs[args.cluster_column].astype(str).to_numpy()
    clusters = cluster_order(set(cluster))
    decisions = read_decisions(args.decisions, set(clusters))
    logging.info("%s: %d nuclei in %d clusters of %s", args.object, adata.n_obs, len(clusters), args.cluster_column)

    # One recorded reason per EXCLUDE cluster and one provenance row for the step, as the review notebooks check.
    reasons = pd.read_csv(args.reasons, dtype=str, keep_default_na=False)
    reasons = reasons[(reasons.species == args.species) & (reasons.step == args.step)].set_index("cluster")
    excluded = set(decisions.index[decisions.label == "EXCLUDE"])
    if set(reasons.index) != excluded:
        raise ValueError(f"EXCLUDE clusters and recorded reasons differ: {sorted(excluded ^ set(reasons.index))}")
    steps = pd.read_csv(args.steps, dtype=str, keep_default_na=False)
    step_record = steps[(steps.species == args.species) & (steps.step == args.step)]
    if len(step_record) != 1:
        raise ValueError(f"{len(step_record)} rows for {args.species} {args.step} in {args.steps}")
    step_record = step_record.iloc[0]

    # Decisions ----

    # Relabelled means a label other than the v8 one: the v8 table of this clustering where there is one, the
    # majority v8 label of the cluster's nuclei for a round-2 re-fit (its clusters are new).
    is_round2 = args.v8_decisions in ("None", "")
    if is_round2:
        reference = v8_reference_from_release(args.species, adata.obs, cluster, read_obs(args.v8_release))
    else:
        reference = v8_reference_from_table(Path(args.v8_decisions), clusters)
    table = cluster_qc(nucleus_qc(adata, args.species)[0], cluster)
    table = table.join(decisions[["label", "note"]]).join(reference)
    table["excluded"] = table.label == "EXCLUDE"
    table["decision"] = [decision_of(label, v8) for label, v8 in zip(table.label, table.v8_reference)]
    table["pending_sign_off"] = table.note.str.contains(PENDING_MARK, case=False, regex=False)
    table["exclusion_reason"] = [reasons.reason.get(c, "") for c in table.index]
    table["reason_source"] = [reasons.source_tag.get(c, "") for c in table.index]
    table.index.name = "cluster"

    # Marker statistics ----

    # Every Table M1 gene of the species: the review plot shows its own lineage's rows and the other lineages'
    # class genes (marker_panels.review_columns), and the checks cover the whole table for every clustering.
    table_m1 = panels.read_table_m1(args.table_m1)
    plot = panels.REVIEW_PLOT[args.species][args.step]
    genes = list(dict.fromkeys(g for row_genes in table_m1.genes[table_m1.species == args.species] for g in row_genes))
    stats, missing = marker_stats(adata, args.cluster_column, genes)
    logging.info("%s %s: %d Table M1 genes; missing from var: %s", args.species, args.step, len(genes), missing)

    # Outputs ----

    table.drop(columns="note").to_csv(f"{stem}_review_clusters.csv")
    stats.to_csv(f"{stem}_review_marker_stats.csv", index=False)
    counts_by_decision = table.groupby("decision").agg(clusters=("label", "size"), nuclei=("n_nuclei", "sum"))
    release = json.loads((args.release_dir / "MANIFEST.json").read_text())
    provenance = {
        "species": args.species,
        "step": args.step,
        "release": release["release"],
        "release_dir": str(args.release_dir),
        "vocabulary_schema": release["vocabulary_schema"],
        "fit": "round 2 (lineage re-fit)" if is_round2 else (
            "July whole-tissue fit" if args.step == "atlas" else "July lineage fit"),
        "object": str(args.object),
        "object_sha256": sha256_of(args.object),
        "n_nuclei": int(adata.n_obs),
        "cluster_column": args.cluster_column,
        "n_clusters": len(clusters),
        "decisions": str(args.decisions),
        "decisions_sha256": sha256_of(args.decisions),
        "v8_reference": "v8 release majority label per cluster" if is_round2 else f"v8 table {args.v8_decisions}",
        "v8_release": str(args.v8_release) if is_round2 else None,
        "decided_at": step_record.decided_at,
        "approval_date": step_record.approval_date,
        "approver": step_record.approver,
        "decision_record": step_record.decision_record,
        "clusters_by_decision": {k: {"clusters": int(v.clusters), "nuclei": int(v.nuclei)}
                                 for k, v in counts_by_decision.iterrows()},
        "pending_sign_off": [f"cluster {c} ({table.label[c]})" for c in table.index[table.pending_sign_off]],
        "marker_plot": plot,
        "genes_scope": "every Table M1 gene of the species",
        "genes": genes,
        "genes_missing_from_var": missing,
        "table_m1": str(args.table_m1),
        "table_m1_sha256": sha256_of(args.table_m1),
        "script_sha256": sha256_of(Path(__file__)),
        "marker_panels_sha256": sha256_of(Path(panels.__file__)),
    }
    Path(f"{stem}_review_extract_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
