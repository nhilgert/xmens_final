"""Comparison functions of the X-Mens release of 19 October 2026, copied verbatim (not edited) from the upstream pipeline.

Source: cagri_ivd_end_to_end_20260815/run_pipeline.py, SHA-256
1c7696a267f0b605c61cad4177c38e35d83a27fc0cc3c07d008f8b29c06764c1 (the upstream file the release recorded in
run_metadata.json). Each function is preceded by its line number in that file. Only the imports and the
PipelineError class below were added, so the functions run outside the upstream repository.

Used by compare_xmens_with_references.py:
  build_consensus_ortholog_support  X-Mens calls for the consensus genes with a strict one-to-one mouse ortholog
  build_wang_support                X-Mens calls for Wang et al.'s in vitro list, routed through current symbols and
                                    source feature IDs, with the oil and pregnancy directions alongside
"""
from __future__ import annotations

import importlib.metadata
import math
import re
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


class PipelineError(RuntimeError):
    """Raised when an input fails a check (the upstream pipeline's error type)."""


# line 59
def fail(message: str) -> None:
    raise PipelineError(message)

# line 75
def scalar_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value).decode("utf-8")
    return str(value).strip()

# line 83
def bool_value(value: Any) -> bool:
    if value is None or pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}

# line 50
VALID_DIRECTIONS = {"up", "down"}

# line 264
def direction(effect: Any) -> str:
    if effect is None or pd.isna(effect) or not math.isfinite(float(effect)):
        return ""
    return "up" if float(effect) > 0 else "down" if float(effect) < 0 else "zero"

# line 270
def classify_support(
    expected_direction: str,
    eligible: Any,
    estimable: Any,
    effect: Any,
    padj: Any,
    *,
    fdr: float,
    abs_lfc: float,
) -> tuple[bool, str, Any, str]:
    eligible_flag = bool_value(eligible)
    estimable_flag = bool_value(estimable)
    finite = (
        effect is not None
        and padj is not None
        and not pd.isna(effect)
        and not pd.isna(padj)
        and math.isfinite(float(effect))
        and math.isfinite(float(padj))
    )
    called = bool(
        eligible_flag
        and estimable_flag
        and finite
        and float(padj) < fdr
        and abs(float(effect)) > abs_lfc
    )
    observed = direction(effect)
    matches: Any = "" if not observed else observed == str(expected_direction)
    if not eligible_flag or not estimable_flag or not finite:
        status = "ineligible_or_unavailable"
    elif not called:
        status = "tested_not_called"
    elif matches:
        status = "supported"
    else:
        status = "reversed"
    return called, observed, matches, status

# line 714
def read_direction_list(path: Path, expected_n: int) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype={"gene": "string", "direction": "string"})
    if list(frame.columns) != ["gene", "direction"]:
        fail(f"Unexpected direction-list schema in {path}")
    if len(frame) != expected_n or frame["gene"].isna().any() or frame["gene"].duplicated().any():
        fail(f"Direction-list cardinality/key failure in {path}")
    if frame["direction"].isna().any():
        fail(f"Direction-list contains a missing direction in {path}")
    if not set(frame["direction"].dropna()).issubset(VALID_DIRECTIONS):
        fail(f"Unexpected direction value in {path}")
    return frame

# line 805
def build_consensus_ortholog_support(
    consensus_path: Path,
    orthology_path: Path,
    xmens: pd.DataFrame,
    calls: dict[str, Any],
    expected_n: int,
) -> pd.DataFrame:
    """Support calls for the consensus genes that have a strict 1:1 mouse ortholog.

    The comparison is scoped by orthology, not by a symbol-matched projection: every
    gene tested has exactly one mouse counterpart, so the comparator is never
    ambiguous. Genes in a one-to-many or many-to-many orthology component are excluded
    rather than resolved by symbol, because symbol identity picks the wrong paralog in
    exactly the cases that matter -- human PRL sits in an Ensembl 116 component with
    all 26 mouse Prl-family genes, and `Prl` names the pituitary gene while the mouse
    decidual program is Prl8a2/Prl3c1.

    Orthology comes from resources/consensus123_mouse_orthology.csv (Ensembl 116
    components cross-checked against OMA All.May2026 Euarchontoglires HOGs); see
    derive_consensus123_orthology.py.
    """
    consensus = pd.read_csv(consensus_path).rename(columns={"gene": "human_gene"})
    orthology = pd.read_csv(orthology_path)
    if sorted(orthology["human_gene"]) != sorted(consensus["human_gene"]):
        fail("Orthology resource does not cover exactly the consensus gene list")
    one_to_one = orthology.loc[orthology["strict_one_to_one"].map(bool_value)].copy()
    if len(one_to_one) != expected_n or one_to_one["human_gene"].nunique() != expected_n:
        fail(f"Strict 1:1 orthology does not yield {expected_n} genes: {len(one_to_one)}")
    if one_to_one["mouse_ensembl_gene_id"].isna().any() or (one_to_one["mouse_ensembl_gene_id"] == "").any():
        fail("A strict 1:1 ortholog row carries no mouse Ensembl gene id")
    if one_to_one["mouse_ensembl_gene_id"].duplicated().any():
        fail("Two consensus genes resolve to the same mouse ortholog")
    joined = consensus.merge(one_to_one, on="human_gene", how="inner", validate="one_to_one", sort=False)
    joined = joined.sort_values("human_gene", kind="stable").reset_index(drop=True)
    joined = joined.merge(
        xmens,
        left_on="mouse_ensembl_gene_id",
        right_on="ensembl_gene_id",
        how="left",
        validate="one_to_one",
    )
    if joined["ensembl_gene_id"].isna().any():
        missing = joined.loc[joined["ensembl_gene_id"].isna(), "human_gene"].tolist()
        fail(f"Mouse orthologs absent from the newly fitted X-Mens table: {missing}")
    rows = []
    for row in joined.to_dict("records"):
        called, observed, matches, status = classify_support(
            row["consensus_direction"],
            row["production_eligible"],
            row["estimable"],
            row["log2FoldChange_unshrunk"],
            row["padj_for_call"],
            fdr=float(calls["padj_lt"]),
            abs_lfc=float(calls["abs_unshrunk_lfc_gt"]),
        )
        rows.append(
            {
                "human_gene": row["human_gene"],
                # Output column keeps the name `mouse_gene`; the resource calls it
                # `mouse_ortholog`. Downstream readers key on `mouse_gene`.
                "mouse_gene": row["mouse_ortholog"],
                "consensus_direction": row["consensus_direction"],
                "tier": row["tier"],
                "n_mature": row["n_mature"],
                "mature_papers": row["mature_papers"],
                "mature_biphasic": bool_value(row["mature_biphasic"]),
                "lucas_onset_dir": row["lucas_onset_dir"],
                "onset_agrees_mature": "" if pd.isna(row["onset_agrees_mature"]) else bool_value(row["onset_agrees_mature"]),
                "orthology_source": "ensembl_116_components_oma_may2026_crosschecked",
                "current_human_symbol": row["current_human_symbol"],
                "human_ensembl_gene_id": row["human_ensembl_gene_id"],
                "mouse_ensembl_gene_id": row["mouse_ensembl_gene_id"],
                "ensembl116_homology_type": row["ensembl116_homology_type"],
                "oma_copy_number_class": row["oma_copy_number_class"],
                "oma_mouse_copy_number": row["oma_mouse_copy_number"],
                "oma_agrees_with_ensembl": row["oma_agrees_with_ensembl"],
                "source_id": "mouse_all_decidual_vs_fibroblast_timepoint_blocked",
                "carrier": bool_value(row["carrier"]),
                "eligible": bool_value(row["production_eligible"]),
                "estimable": bool_value(row["estimable"]),
                "log2FoldChange_unshrunk": row["log2FoldChange_unshrunk"],
                "lfcSE": row["lfcSE"],
                "pvalue": row["pvalue"],
                "padj": row["padj"],
                "padj_for_call": row["padj_for_call"],
                "observed_direction": observed,
                "positive_lfc_meaning": "decidual_minus_fibroblast",
                "strict_lfc1_call": called,
                "direction_match": matches,
                "support_status": status,
            }
        )
    return pd.DataFrame(rows)

# line 900
def add_direction_only(rows: pd.DataFrame, comparator: pd.DataFrame, prefix: str) -> pd.DataFrame:
    direction_by_gene = comparator.set_index("gene")["direction"]
    observed = rows["wang_ivd_gene"].map(direction_by_gene)
    called = observed.notna()
    status = pd.Series("not_in_published_deg_list", index=rows.index, dtype="string")
    status.loc[called & observed.eq(rows["wang_ivd_direction"])] = "supported"
    status.loc[called & ~observed.eq(rows["wang_ivd_direction"])] = "reversed"
    rows[f"{prefix}_called"] = called
    rows[f"{prefix}_direction"] = observed
    rows[f"{prefix}_effect_log2_fold_change"] = pd.Series(pd.NA, index=rows.index, dtype="Float64")
    rows[f"{prefix}_adjusted_p_value"] = pd.Series(pd.NA, index=rows.index, dtype="Float64")
    rows[f"{prefix}_mapping_evidence"] = called.map({True: "exact_mouse_gene_symbol", False: "absent_from_direction_only_deg_list"})
    rows[f"{prefix}_status"] = status
    return rows

# line 916
def xmens_candidates(ivd_genes: set[str], xmens: pd.DataFrame) -> dict[str, list[int]]:
    output: dict[str, list[int]] = {gene: [] for gene in ivd_genes}
    for index, row in xmens.iterrows():
        symbols = {scalar_text(row["gene_symbol"])} if scalar_text(row["gene_symbol"]) else set()
        for token in re.split(r"[|,;]", scalar_text(row["source_feature_ids"])):
            if token.strip():
                symbols.add(token.strip())
        for symbol in symbols & ivd_genes:
            output[symbol].append(index)
    return output

# line 928
def select_xmens_candidate(rows: pd.DataFrame) -> pd.Series:
    ranked = rows.copy()
    ranked["_estimable"] = ranked["estimable"].map(bool_value).astype(int)
    ranked["_eligible"] = ranked["production_eligible"].map(bool_value).astype(int)
    ranked["_carrier"] = ranked["carrier"].map(bool_value).astype(int)
    ranked["_count"] = pd.to_numeric(ranked["total_count"], errors="coerce").fillna(-1)
    return ranked.sort_values(
        ["_estimable", "_eligible", "_carrier", "_count", "ensembl_gene_id"],
        ascending=[False, False, False, False, True],
        kind="stable",
    ).iloc[0]

# line 941
def build_wang_support(
    ivd_path: Path,
    oil_path: Path,
    pregnancy_path: Path,
    xmens: pd.DataFrame,
    calls: dict[str, Any],
) -> pd.DataFrame:
    ivd = read_direction_list(ivd_path, 1501)
    oil = read_direction_list(oil_path, 4294)
    pregnancy = read_direction_list(pregnancy_path, 5277)
    rows = ivd.rename(columns={"gene": "wang_ivd_gene", "direction": "wang_ivd_direction"}).copy()
    rows.insert(0, "wang_ivd_source_row", range(1, len(rows) + 1))
    rows = add_direction_only(rows, oil, "wang_oil")
    rows = add_direction_only(rows, pregnancy, "wang_pregnancy")
    candidates = xmens_candidates(set(rows["wang_ivd_gene"].astype(str)), xmens)
    records: list[dict[str, Any]] = []
    for gene, expected_direction in zip(rows["wang_ivd_gene"].astype(str), rows["wang_ivd_direction"].astype(str)):
        indices = sorted(set(candidates[gene]))
        if not indices:
            records.append(
                {
                    "xmens_match_route": "unmatched",
                    "xmens_match_candidate_n": 0,
                    "xmens_candidate_ensembl_ids": pd.NA,
                    "xmens_selected_ensembl_gene_id": pd.NA,
                    "xmens_current_gene_symbol": pd.NA,
                    "xmens_source_feature_id": pd.NA,
                    "xmens_alias_mapping_evidence": "no_current_symbol_or_source_feature_id_match",
                    "xmens_candidate_resolution_note": pd.NA,
                    "xmens_production_eligible": pd.NA,
                    "xmens_estimable": pd.NA,
                    "xmens_effect_log2_fold_change_unshrunk": pd.NA,
                    "xmens_adjusted_p_value_for_call": pd.NA,
                    "xmens_direction": pd.NA,
                    "xmens_called_strict_lfc1": False,
                    "xmens_status": "ineligible_or_unmatched",
                }
            )
            continue
        candidate_rows = xmens.loc[indices]
        selected = select_xmens_candidate(candidate_rows)
        ids = sorted(candidate_rows["ensembl_gene_id"].astype(str).unique())
        current_symbol = scalar_text(selected["gene_symbol"])
        source_features = scalar_text(selected["source_feature_ids"])
        called, observed, _, status = classify_support(
            expected_direction,
            selected["production_eligible"],
            selected["estimable"],
            selected["log2FoldChange_unshrunk"],
            selected["padj_for_call"],
            fdr=float(calls["padj_lt"]),
            abs_lfc=float(calls["abs_unshrunk_lfc_gt"]),
        )
        route = "current_gene_symbol" if current_symbol == gene else "source_feature_id_alias"
        if status == "ineligible_or_unavailable":
            status = "ineligible_or_unmatched"
        records.append(
            {
                "xmens_match_route": route,
                "xmens_match_candidate_n": len(ids),
                "xmens_candidate_ensembl_ids": ";".join(ids),
                "xmens_selected_ensembl_gene_id": selected["ensembl_gene_id"],
                "xmens_current_gene_symbol": current_symbol,
                "xmens_source_feature_id": source_features,
                "xmens_alias_mapping_evidence": f"Wang_symbol={gene};source_feature_id={source_features};current_symbol={current_symbol}",
                "xmens_candidate_resolution_note": "multiple_candidates_resolved_by_estimable_production_carrier_count_ensembl_priority" if len(ids) > 1 else pd.NA,
                "xmens_production_eligible": bool_value(selected["production_eligible"]),
                "xmens_estimable": bool_value(selected["estimable"]),
                "xmens_effect_log2_fold_change_unshrunk": selected["log2FoldChange_unshrunk"],
                "xmens_adjusted_p_value_for_call": selected["padj_for_call"],
                "xmens_direction": observed,
                "xmens_called_strict_lfc1": called,
                "xmens_status": status,
            }
        )
    return pd.concat([rows.reset_index(drop=True), pd.DataFrame(records)], axis=1)

# line 1311
def normalized_partition(series: pd.Series, categories: Iterable[str]) -> dict[str, int]:
    counts = series.astype(str).value_counts().to_dict()
    return {str(category): int(counts.get(str(category), 0)) for category in categories}

# line 1339
def comparison_summary_row(
    comparison_id: str,
    label: str,
    status: pd.Series,
    unavailable_label: str,
) -> dict[str, Any]:
    counts = status.astype(str).value_counts()
    denominator = int(len(status))
    supported = int(counts.get("supported", 0))
    reversed_n = int(counts.get("reversed", 0))
    not_called_or_listed = int(
        counts.get("tested_not_called", 0) + counts.get("not_in_published_deg_list", 0)
    )
    unavailable = int(counts.get(unavailable_label, 0))
    called = supported + reversed_n
    return {
        "comparison_id": comparison_id,
        "comparison_label": label,
        "denominator_n": denominator,
        "called_n": called,
        "supported_n": supported,
        "reversed_n": reversed_n,
        "not_called_or_not_listed_n": not_called_or_listed,
        "ineligible_or_unavailable_n": unavailable,
        "support_fraction_all": supported / denominator if denominator else np.nan,
        "direction_agreement_among_called": supported / called if called else np.nan,
        "report_sentence": f"Of {denominator:,} genes, {supported:,} are supported.",
    }

# line 1488
def package_versions() -> dict[str, str]:
    names = [
        "anndata",
        "formulaic",
        "h5py",
        "matplotlib",
        "numpy",
        "pandas",
        "pyarrow",
        "pydeseq2",
        "PyYAML",
        "scikit-learn",
        "scipy",
        "statsmodels",
    ]
    versions: dict[str, str] = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not_installed"
    return versions

