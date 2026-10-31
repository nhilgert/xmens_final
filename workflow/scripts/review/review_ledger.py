"""Removed-cluster ledger rows of one annotation step: QC, kit verdict, recorded reason and v8 origin.

Used by the review notebooks (notebooks/review/) after review_step. review_step gives the cluster table and the
per-nucleus QC; this module adds the doublet-score summaries the cluster table lacks, where every nucleus was in
v8 (the v8 extended release and, for the round-2 lineages, the July lineage fit) and, for each EXCLUDE cluster,
the recorded reason and decision provenance (resources/review/exclusion_reasons_v9.csv, decision_steps_v9.csv).
"""
import anndata as ad
import h5py
import numpy as np
import pandas as pd

DOUBLET_SCORE_CUT = 5.0   # DoubletDetection score; the same cut as the kit's doublet-score rule
NOT_RELEASED = "not released"

LEDGER_COLUMNS = [
    "species", "step", "fit", "object", "cluster_column", "cluster",
    "n_nuclei", "median_genes", "median_umis", "mean_pct_mito",
    "mean_doublet_score", "median_doublet_score", "pct_doublet_score_above_5",
    "pct_two_lineages", "largest_array", "largest_array_share",
    "kit_verdict", "kit_best_label", "kit_best_r",
    "final_reason", "reason_source", "reason_v8",
    "evidence_summary", "evidence_source", "n_released_in_v8", "n_in_v8_primary", "origin_v8",
    "decided_at", "approval_date", "approver", "decision_record", "decisions_table",
    # References as in the kit's rules: median genes and % above the doublet cut over the kept nuclei; medians of
    # the kept clusters' two-lineage % and mean doublet score.
    "kept_median_genes", "kept_pct_doublet_score_above_5",
    "kept_clusters_median_pct_two_lineages", "kept_clusters_median_mean_doublet_score",
]


def read_obs(path):
    """obs of an h5ad file, without reading its matrices."""
    with h5py.File(path, "r") as handle:
        return ad.io.read_elem(handle["obs"])


def join_key(obs):
    """lane|barcode: stable across the pipeline's objects (sample names may differ in prefix)."""
    sample = obs["sample"].astype(str)
    lane = sample.str.extract(r"(mck\w+)$", expand=False).fillna(sample)
    return lane + "|" + obs["cell_barcode"].astype(str)


def v8_status(obs, release_obs):
    """Per nucleus: its cell_subtype in the v8 extended release, or "not released".

    Nuclei of the supplementary populations (obs in_primary False; extended objects only) are marked
    "(extended only)".
    """
    release_key = join_key(release_obs)
    if release_key.duplicated().any():
        raise ValueError(f"lane|barcode is not unique in the v8 release ({int(release_key.duplicated().sum())} repeats)")
    label = release_obs["cell_subtype"].astype(str).to_numpy()
    in_primary = release_obs["in_primary"].astype(bool).to_numpy()
    shown = pd.Series(np.where(in_primary, label, label + " (extended only)"), index=release_key.to_numpy())
    return join_key(obs).map(shown).fillna(NOT_RELEASED)


def july_origin(obs, july_obs, july_cluster_column):
    """Per nucleus of a round-2 object: "July cluster N" (in the July lineage fit) or "atlas cluster N" (reinstated)."""
    july_cluster = pd.Series(obs.index, index=obs.index).map(july_obs[july_cluster_column].astype(str))
    atlas_cluster = "atlas cluster " + obs["atlas_cluster"].astype(str)
    return ("July cluster " + july_cluster).where(july_cluster.notna(), atlas_cluster)


def counts_text(values):
    """'a: 50, b: 2': counts of each value, largest first (ties by name, so the text is reproducible)."""
    counts = values.value_counts()
    order = sorted(counts.index, key=lambda name: (-counts[name], str(name)))
    return ", ".join(f"{name}: {counts[name]:,}" for name in order)


def origin_summary(status, origin=None, v8_decision=None):
    """One line on where a cluster's nuclei were in v8.

    status: v8_status of the cluster's nuclei. origin: july_origin of the same nuclei (round-2 lineages only).
    v8_decision: the cluster's label in the v8 table of the same clustering (atlas and July-fit lineages only).
    """
    parts = []
    if v8_decision is not None:
        parts.append(f"v8 table: {v8_decision}")
    released = status != NOT_RELEASED
    text = f"{int(released.sum()):,} of {len(status):,} nuclei in the v8 release"
    if released.any():
        text += f" ({counts_text(status[released])})"
    parts.append(text)
    if origin is not None:
        if (~released).any():
            parts.append(f"{int((~released).sum()):,} not released, from {counts_text(origin[~released])}")
        reinstated = origin.str.startswith("atlas cluster")
        if reinstated.any():
            parts.append(f"{int(reinstated.sum()):,} reinstated in round 2 from {counts_text(origin[reinstated])}")
    return "; ".join(parts)


def doublet_summary(qc, cluster):
    """Median DoubletDetection score and % of nuclei above DOUBLET_SCORE_CUT, per cluster (nuclei with a score)."""
    score = qc["doublet_score"].astype(float)
    grouped = score.groupby(pd.Series(cluster, index=qc.index))
    return pd.DataFrame({
        "median_doublet_score": grouped.median(),
        "pct_doublet_score_above_5": grouped.apply(lambda s: 100.0 * float((s.dropna() > DOUBLET_SCORE_CUT).mean())),
    })


def removed_cluster_rows(table, qc, cluster, status, step_info, reasons, origin=None, v8_decisions=None):
    """Ledger rows of the step's EXCLUDE clusters.

    table, qc: review_step's cluster table and per-nucleus QC; cluster: each nucleus's cluster (aligned with qc).
    status, origin: v8_status and july_origin of every nucleus (origin None outside the round-2 lineages).
    step_info: dict with species, step, fit, object, cluster_column, decisions_table and the step's row of
    decision_steps_v9.csv (decided_at, approval_date, approver, decision_record).
    reasons: the step's rows of exclusion_reasons_v9.csv, indexed by cluster.
    v8_decisions: cluster -> v8 label of the same clustering, or None.
    """
    doublets = doublet_summary(qc, cluster)
    kept_clusters = table.loc[~table.excluded, "cluster"]
    kept_scores = qc["doublet_score"].astype(float)[np.isin(cluster, kept_clusters)].dropna()
    kept_above_cut = 100.0 * float((kept_scores > DOUBLET_SCORE_CUT).mean()) if len(kept_scores) else np.nan

    rows = []
    for row in table[table.excluded].itertuples():
        member = cluster == row.cluster
        recorded = reasons.loc[row.cluster]
        rows.append({
            **{key: step_info[key] for key in ("species", "step", "fit", "object", "cluster_column")},
            "cluster": row.cluster,
            "n_nuclei": row.n_nuclei,
            "median_genes": row.median_genes,
            "median_umis": row.median_umis,
            "mean_pct_mito": row.mean_pct_mito,
            "mean_doublet_score": row.mean_doublet,
            "median_doublet_score": doublets.loc[row.cluster, "median_doublet_score"],
            "pct_doublet_score_above_5": doublets.loc[row.cluster, "pct_doublet_score_above_5"],
            "pct_two_lineages": row.pct_two_lineages,
            "largest_array": row.top_sample,
            "largest_array_share": row.top_sample_share,
            "kit_verdict": row.verdict,
            "kit_best_label": row.best_label,
            "kit_best_r": row.best_r,
            "final_reason": recorded["reason"],
            "reason_source": recorded["reason_source"],
            "reason_v8": recorded["reason_v8"],
            "evidence_summary": recorded["evidence"],
            "evidence_source": recorded["evidence_source"],
            "n_released_in_v8": int((status[member] != NOT_RELEASED).sum()),
            "n_in_v8_primary": int(((status[member] != NOT_RELEASED) & ~status[member].str.endswith("(extended only)")).sum()),
            "origin_v8": origin_summary(
                status[member],
                origin[member] if origin is not None else None,
                v8_decisions.get(row.cluster, "absent") if v8_decisions is not None else None,
            ),
            **{key: step_info[key] for key in ("decided_at", "approval_date", "approver", "decision_record", "decisions_table")},
            "kept_median_genes": row.ref_kept_median_genes,
            "kept_pct_doublet_score_above_5": kept_above_cut,
            "kept_clusters_median_pct_two_lineages": row.ref_kept_median_two_lineages,
            "kept_clusters_median_mean_doublet_score": row.ref_kept_median_doublet,
        })
    return pd.DataFrame(rows, columns=LEDGER_COLUMNS)
