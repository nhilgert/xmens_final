#!/usr/bin/env python3
"""Minimal object-to-DE-to-handoff pipeline.

The two expression inputs are the authoritative HECA cells and X-Mens mouse
H5AD objects.  This script streams raw counts from ``layers/counts``, creates
biological-unit pseudobulks, fits the two prespecified DE models, classifies
the five requested support comparisons, and writes the tabular inputs for the
final workbook.  No comparison result or frozen DE table is used as an input.
"""

from __future__ import annotations

import argparse
import atexit
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "cagri_ivd_matplotlib"))

try:
    import anndata as ad
    import h5py
except ModuleNotFoundError:  # Keep --help usable before the locked environment is installed.
    ad = None
    h5py = None
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


class PipelineError(RuntimeError):
    pass


VALID_DIRECTIONS = {"up", "down"}
VALID_SUPPORT = {
    "supported",
    "reversed",
    "tested_not_called",
    "ineligible_or_unavailable",
}


def fail(message: str) -> None:
    raise PipelineError(message)


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def scalar_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (bytes, np.bytes_)):
        return bytes(value).decode("utf-8")
    return str(value).strip()


def bool_value(value: Any) -> bool:
    if value is None or pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def read_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or config.get("schema_version") != "1.0":
        fail("config.yaml must be a schema_version 1.0 mapping")
    return config


def resolve_resource(repo: Path, spec: dict[str, Any]) -> Path:
    path = Path(str(spec["path"]))
    path = path if path.is_absolute() else repo / path
    path = path.resolve()
    if not path.is_file():
        fail(f"Required input is absent: {path}")
    observed = sha256_file(path)
    if observed != str(spec["sha256"]):
        fail(f"Input hash mismatch for {path}: expected {spec['sha256']}, observed {observed}")
    return path


def verify_object(path: Path, spec: dict[str, Any], *, check_hash: bool) -> None:
    if not path.is_file():
        fail(f"H5AD object is absent: {path}")
    if check_hash:
        observed = sha256_file(path)
        if observed != str(spec["sha256"]):
            fail(f"Object hash mismatch for {path}: expected {spec['sha256']}, observed {observed}")
    with h5py.File(path, "r") as handle:
        layer = f"layers/{spec['counts_layer']}"
        if layer not in handle:
            fail(f"{path} lacks /{layer}")
        matrix = handle[layer]
        shape = tuple(int(value) for value in matrix.attrs.get("shape", ()))
        if shape != tuple(int(value) for value in spec["shape"]):
            fail(f"Unexpected object shape for {path}: {shape}")
        if scalar_text(matrix.attrs.get("encoding-type")) != "csr_matrix":
            fail(f"/{layer} must be CSR")


def read_object_metadata(path: Path, columns: Iterable[str]) -> tuple[pd.DataFrame, pd.Index]:
    backed = ad.read_h5ad(path, backed="r")
    try:
        missing = sorted(set(columns) - set(backed.obs.columns))
        if missing:
            fail(f"{path} lacks required obs columns: {missing}")
        if not backed.obs_names.is_unique or not backed.var_names.is_unique:
            fail(f"{path} has duplicated observation or variable IDs")
        obs = backed.obs[list(columns)].copy()
        obs.index = pd.Index(backed.obs_names.astype(str), name="cell_id")
        var_names = pd.Index(backed.var_names.astype(str), name="source_feature_id")
    finally:
        backed.file.close()
    return obs, var_names


def read_feature_manifest(path: Path, var_names: pd.Index) -> tuple[np.ndarray, pd.DataFrame]:
    frame = pd.read_csv(path, sep="\t", dtype={"canonical_id": "string"})
    required = {
        "source_index",
        "source_feature_id",
        "source_symbol",
        "canonical_id",
        "excluded",
        "in_annotation",
        "gene_symbol",
        "biotype",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        fail(f"Feature manifest {path} lacks columns: {missing}")
    frame = frame.sort_values("source_index", kind="stable").reset_index(drop=True)
    if frame["source_index"].tolist() != list(range(len(frame))):
        fail(f"Feature manifest {path} source_index is not 0..n-1")
    if len(frame) != len(var_names):
        fail(f"Feature manifest/object width mismatch: {len(frame)} != {len(var_names)}")
    if frame["source_feature_id"].astype(str).tolist() != var_names.astype(str).tolist():
        fail(f"Feature manifest source_feature_id order differs from {path}")

    included = ~frame["excluded"].map(bool_value)
    included &= frame["canonical_id"].notna() & frame["canonical_id"].astype(str).ne("")
    canonical_order = list(dict.fromkeys(frame.loc[included, "canonical_id"].astype(str)))
    canonical_index = {gene: idx for idx, gene in enumerate(canonical_order)}
    source_to_canonical = np.full(len(frame), -1, dtype=np.int64)
    for source_idx, gene in frame.loc[included, ["source_index", "canonical_id"]].itertuples(index=False):
        source_to_canonical[int(source_idx)] = canonical_index[str(gene)]

    records: list[dict[str, Any]] = []
    included_frame = frame.loc[included].copy()
    for canonical_id, group in included_frame.groupby("canonical_id", sort=False, observed=True):
        symbols = [x for x in group["gene_symbol"].dropna().astype(str).unique() if x]
        biotypes = [x for x in group["biotype"].dropna().astype(str).unique() if x]
        records.append(
            {
                "ensembl_gene_id": str(canonical_id),
                "gene_symbol": symbols[0] if len(symbols) == 1 else pd.NA,
                "biotype": biotypes[0] if len(biotypes) == 1 else pd.NA,
                "source_feature_ids": ";".join(group["source_feature_id"].astype(str)),
                "source_feature_n": int(len(group)),
            }
        )
    canonical = pd.DataFrame(records)
    if canonical["ensembl_gene_id"].tolist() != canonical_order:
        fail("Canonical feature order changed while summarizing the feature map")
    return source_to_canonical, canonical


def stream_canonical_pseudobulk(
    h5ad_path: Path,
    counts_layer: str,
    assignments: np.ndarray,
    group_ids: list[str],
    source_to_canonical: np.ndarray,
    *,
    row_chunk_size: int = 4096,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Aggregate selected CSR rows directly into group-by-canonical counts."""

    if assignments.ndim != 1:
        fail("Assignments must be one-dimensional")
    output = np.zeros((len(group_ids), int(source_to_canonical.max()) + 1), dtype=np.int64)
    cell_counts = np.zeros(len(group_ids), dtype=np.int64)
    included_count_sum = 0
    with h5py.File(h5ad_path, "r") as handle:
        group = handle[f"layers/{counts_layer}"]
        n_rows, n_features = (int(x) for x in group.attrs["shape"])
        if len(assignments) != n_rows or len(source_to_canonical) != n_features:
            fail("Assignment/feature-map dimensions differ from the CSR matrix")
        indptr_ds = group["indptr"]
        indices_ds = group["indices"]
        data_ds = group["data"]
        for row_start in range(0, n_rows, row_chunk_size):
            row_end = min(n_rows, row_start + row_chunk_size)
            local_assignments = assignments[row_start:row_end]
            if not np.any(local_assignments >= 0):
                continue
            raw_ptr = np.asarray(indptr_ds[row_start : row_end + 1], dtype=np.int64)
            offset = int(raw_ptr[0])
            stop = int(raw_ptr[-1])
            ptr = raw_ptr - offset
            indices = np.asarray(indices_ds[offset:stop], dtype=np.int64)
            raw_data = np.asarray(data_ds[offset:stop])
            if raw_data.dtype.kind == "f":
                if not np.isfinite(raw_data).all() or not np.equal(raw_data, np.floor(raw_data)).all():
                    fail("Counts layer contains non-finite or fractional values")
            if np.any(raw_data < 0):
                fail("Counts layer contains negative values")
            data = raw_data.astype(np.int64, copy=False)
            for local_row, assignment in enumerate(local_assignments):
                target = int(assignment)
                if target < 0:
                    continue
                lo, hi = int(ptr[local_row]), int(ptr[local_row + 1])
                source_indices = indices[lo:hi]
                canonical_indices = source_to_canonical[source_indices]
                keep = canonical_indices >= 0
                if np.any(keep):
                    values = data[lo:hi][keep]
                    np.add.at(output[target], canonical_indices[keep], values)
                    included_count_sum += int(values.sum(dtype=np.int64))
                cell_counts[target] += 1
    if int(output.sum(dtype=np.int64)) != included_count_sum:
        fail("Pseudobulk canonical-count conservation failed")
    return output, cell_counts, included_count_sum


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, na_rep="", float_format="%.15g", lineterminator="\n")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.15g", lineterminator="\n")


def direction(effect: Any) -> str:
    if effect is None or pd.isna(effect) or not math.isfinite(float(effect)):
        return ""
    return "up" if float(effect) > 0 else "down" if float(effect) < 0 else "zero"


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


def prepare_heca_f(
    object_path: Path,
    object_spec: dict[str, Any],
    feature_manifest_path: Path,
    model: dict[str, Any],
    workdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    fields = [
        model["cohort_field"],
        model["state_field"],
        *model["person_fields"],
        "sample",
        "library",
        "Library_genotype",
        "Treatment",
        *model.get("constant_context_fields", []),
    ]
    obs, var_names = read_object_metadata(object_path, fields)
    source_to_canonical, features = read_feature_manifest(feature_manifest_path, var_names)
    states = obs[model["state_field"]].astype(str)
    relevant_state = states.isin([model["denominator"], model["numerator"]])
    person_parts = [obs[field].astype("string").str.strip() for field in model["person_fields"]]
    if any(part[relevant_state].isna().any() or part[relevant_state].eq("").any() for part in person_parts):
        fail("Relevant HECA cells contain missing person fields")
    person = person_parts[0].astype(str)
    for part in person_parts[1:]:
        person = person + str(model["person_separator"]) + part.astype(str)

    cohort = obs[model["cohort_field"]].astype(str).eq(str(model["cohort_value"]))
    selected = cohort & relevant_state
    if not selected.any():
        fail("No cells pass the HECA F cohort/state selection")

    # Scope the constancy guard to the cells that actually enter the model. Checking every
    # stromal cell in the atlas aborts on donors whose out-of-cohort biopsies legitimately
    # differ in fields that are never fitted. (PENDING_PIPELINE_FIXES fix 2, applied 2026-09-02.)
    relevant = pd.DataFrame(
        {"person": person[selected], "state": states[selected]},
        index=obs.index[selected],
    )
    for (person_id, state), index in relevant.groupby(["person", "state"], sort=False).groups.items():
        for field in model.get("constant_context_fields", []):
            if obs.loc[index, field].nunique(dropna=False) != 1:
                fail(f"HECA person-state {person_id}/{state} has inconsistent {field}")

    selected_frame = pd.DataFrame(
        {"person": person[selected], "state": states[selected]}, index=obs.index[selected]
    )
    keys = sorted(
        set(zip(selected_frame["person"], selected_frame["state"])),
        key=lambda item: (item[0], 0 if item[1] == model["denominator"] else 1),
    )
    group_ids = [f"person={person_id}__state={state}" for person_id, state in keys]
    key_to_group = {key: idx for idx, key in enumerate(keys)}
    assignments = np.full(len(obs), -1, dtype=np.int64)
    selected_positions = np.flatnonzero(selected.to_numpy())
    for position in selected_positions:
        assignments[position] = key_to_group[(str(person.iloc[position]), str(states.iloc[position]))]

    matrix, cell_counts, _ = stream_canonical_pseudobulk(
        object_path,
        str(object_spec["counts_layer"]),
        assignments,
        group_ids,
        source_to_canonical,
    )
    metadata_rows = []
    for idx, (person_id, state) in enumerate(keys):
        mask = selected & person.eq(person_id) & states.eq(state)
        datasets = obs.loc[mask, model["person_fields"][0]].astype(str).unique().tolist()
        if len(datasets) != 1:
            fail(f"HECA person {person_id} maps to multiple datasets")
        metadata_rows.append(
            {
                "pseudobulk_id": group_ids[idx],
                "person": person_id,
                "dataset": datasets[0],
                "state": state,
                "n_cells": int(cell_counts[idx]),
                "n_samples": int(obs.loc[mask, "sample"].nunique()),
                "n_libraries": int(obs.loc[mask, "library"].nunique()),
                "n_library_genotypes": int(obs.loc[mask, "Library_genotype"].nunique()),
                "n_treatments": int(obs.loc[mask, "Treatment"].nunique()),
            }
        )
    metadata = pd.DataFrame(metadata_rows)
    counts = pd.DataFrame(matrix, columns=features["ensembl_gene_id"], index=metadata.index)
    keep = metadata["n_cells"].ge(int(model["cell_floor_after_person_state_aggregation"]))
    metadata = metadata.loc[keep].reset_index(drop=True)
    counts = counts.loc[keep.to_numpy()].reset_index(drop=True)

    if len(metadata) != int(model["expected_pseudobulks"]):
        fail(f"HECA F pseudobulk count changed: {len(metadata)}")
    if metadata["person"].nunique() != int(model["expected_people"]):
        fail("HECA F person count changed")
    paired_people = int(metadata.groupby("person")["state"].nunique().eq(2).sum())
    if paired_people != int(model["expected_paired_people"]):
        fail(f"HECA F paired-person count changed: {paired_people}")
    observed_state_n = metadata["state"].value_counts().to_dict()
    expected_state_n = {
        str(model["numerator"]): int(model["expected_numerator_pseudobulks"]),
        str(model["denominator"]): int(model["expected_denominator_pseudobulks"]),
    }
    if observed_state_n != expected_state_n:
        fail(f"HECA F state pseudobulk counts changed: {observed_state_n}")
    if metadata.groupby("person")["dataset"].nunique().gt(1).any():
        fail("A HECA biological person maps to multiple datasets")
    if metadata.duplicated(["person", "state"]).any():
        fail("HECA technical replicates were not collapsed to person x state")

    # Eligibility (2026-09-02): a gene is testable if expressed (CPM >= threshold on raw
    # per-pseudobulk library sizes) in the majority of pseudobulks of at least one arm.
    # Arm-aware so switch-on decidualization genes can qualify through the decidual arm
    # alone; CPM >= 1 is the classic limma-workflow threshold (Law et al. 2014/2016).
    cpm_threshold = float(model["eligibility_cpm_threshold"])
    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    numerator_mask = metadata["state"].astype(str).eq(str(model["numerator"])).to_numpy()
    numerator_majority = (int(numerator_mask.sum()) + 1) // 2
    denominator_majority = (int((~numerator_mask).sum()) + 1) // 2
    eligible = (
        cpm.loc[numerator_mask].ge(cpm_threshold).sum(axis=0).ge(numerator_majority)
        | cpm.loc[~numerator_mask].ge(cpm_threshold).sum(axis=0).ge(denominator_majority)
    )
    fit_counts = counts.loc[:, eligible]
    fit_counts_out = fit_counts.copy()
    fit_counts_out.insert(0, "pseudobulk_id", metadata["pseudobulk_id"].to_numpy())
    write_tsv(fit_counts_out, workdir / "heca_f_fit_counts.tsv")
    write_tsv(metadata, workdir / "heca_f_metadata.tsv")
    write_tsv(
        pd.DataFrame(
            {
                "ensembl_gene_id": features["ensembl_gene_id"],
                "carrier": True,
                "production_eligible": eligible.to_numpy(),
                "eligibility_cpm_threshold": cpm_threshold,
                "numerator_pseudobulks_ge_cpm": cpm.loc[numerator_mask].ge(cpm_threshold).sum(axis=0).to_numpy(),
                "denominator_pseudobulks_ge_cpm": cpm.loc[~numerator_mask].ge(cpm_threshold).sum(axis=0).to_numpy(),
                "numerator_majority": numerator_majority,
                "denominator_majority": denominator_majority,
            }
        ),
        workdir / "heca_f_eligibility.tsv",
    )
    return (
        pd.concat(
            [
                features.reset_index(drop=True),
                pd.DataFrame(
                    {
                        "carrier": True,
                        "production_eligible": eligible.to_numpy(),
                        "total_count": counts.sum(axis=0).to_numpy(dtype=np.int64),
                        "nonzero_pseudobulks": counts.gt(0).sum(axis=0).to_numpy(dtype=np.int64),
                    }
                ),
            ],
            axis=1,
        ),
        metadata,
    )


def fit_heca_f(
    carrier: pd.DataFrame,
    metadata: pd.DataFrame,
    workdir: Path,
    repo: Path,
    rscript: str,
) -> pd.DataFrame:
    output = workdir / "heca_f_fit.tsv"
    command = [
        rscript,
        str(repo / "fit_heca_f.R"),
        str(workdir / "heca_f_fit_counts.tsv"),
        str(workdir / "heca_f_metadata.tsv"),
        str(output),
        str(workdir / "heca_f_sample_weights.tsv"),
    ]
    subprocess.run(command, check=True)
    fitted = pd.read_csv(output, sep="\t")
    if fitted["ensembl_gene_id"].duplicated().any():
        fail("HECA F fit contains duplicated Ensembl IDs")
    expected_ids = carrier.loc[carrier["production_eligible"], "ensembl_gene_id"].astype(str).tolist()
    if fitted["ensembl_gene_id"].astype(str).tolist() != expected_ids:
        fail("HECA F fit IDs/order differ from the eligible input")
    result = carrier.merge(fitted, on="ensembl_gene_id", how="left", validate="one_to_one")
    ineligible_with_fit = ~result["production_eligible"] & result["log2FoldChange_unshrunk"].notna()
    if ineligible_with_fit.any():
        fail("HECA F returned statistics for ineligible genes")
    fit_estimable = result["estimable"].astype("boolean").fillna(False).astype(bool)
    result["estimable"] = result["production_eligible"] & fit_estimable
    inference = ["log2FoldChange_unshrunk", "lfcSE", "stat", "pvalue", "padj"]
    result.loc[~result["estimable"], inference] = np.nan
    result["padj_for_call"] = result["padj"].where(np.isfinite(result["padj"]), 1.0)
    result["significant_fdr05_lfc1"] = (
        result["estimable"]
        & result["padj_for_call"].lt(0.05)
        & result["log2FoldChange_unshrunk"].abs().gt(1.0)
    )
    result["source_id"] = "F_population25_production"
    result["positive_lfc_meaning"] = "dStromal_late_minus_eStromal"
    result["n_pseudobulks_numerator"] = int((metadata["state"] == "dStromal_late").sum())
    result["n_pseudobulks_denominator"] = int((metadata["state"] == "eStromal").sum())
    return result


def prepare_xmens_broad(
    object_path: Path,
    object_spec: dict[str, Any],
    feature_manifest_path: Path,
    model: dict[str, Any],
    workdir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fields = [
        model["coarse_state_field"],
        "fine_refined",
        model["timepoint_field"],
        model["sample_field"],
    ]
    obs, var_names = read_object_metadata(object_path, fields)
    source_to_canonical, features = read_feature_manifest(feature_manifest_path, var_names)
    coarse = obs[model["coarse_state_field"]].astype(str)
    state = pd.Series(pd.NA, index=obs.index, dtype="string")
    state.loc[coarse.eq(str(model["numerator_value"]))] = "decidual"
    state.loc[coarse.eq(str(model["denominator_value"]))] = "fibroblast"
    selected = state.notna()
    if int(selected.sum()) != int(model["expected_selected_cells"]):
        fail(f"X-Mens selected cell count changed: {int(selected.sum())}")
    timepoint = obs[model["timepoint_field"]].astype(str)
    keys = sorted(
        set(zip(timepoint[selected], state[selected].astype(str))),
        key=lambda item: (
            list(model["expected_timepoints"]).index(item[0])
            if item[0] in model["expected_timepoints"]
            else len(model["expected_timepoints"]),
            0 if item[1] == "fibroblast" else 1,
        ),
    )
    group_ids = [f"timepoint={tp}__group={group}" for tp, group in keys]
    key_to_group = {key: idx for idx, key in enumerate(keys)}
    assignments = np.full(len(obs), -1, dtype=np.int64)
    for position in np.flatnonzero(selected.to_numpy()):
        assignments[position] = key_to_group[(str(timepoint.iloc[position]), str(state.iloc[position]))]
    matrix, cell_counts, _ = stream_canonical_pseudobulk(
        object_path,
        str(object_spec["counts_layer"]),
        assignments,
        group_ids,
        source_to_canonical,
    )
    metadata_rows = []
    for idx, (tp, group) in enumerate(keys):
        mask = selected & timepoint.eq(tp) & state.eq(group)
        metadata_rows.append(
            {
                "pseudobulk_id": group_ids[idx],
                "timepoint": tp,
                "group": group,
                "n_cells": int(cell_counts[idx]),
                "n_samples": int(obs.loc[mask, model["sample_field"]].nunique()),
                "samples": ";".join(sorted(obs.loc[mask, model["sample_field"]].astype(str).unique())),
            }
        )
    metadata = pd.DataFrame(metadata_rows)
    counts = pd.DataFrame(matrix, columns=features["ensembl_gene_id"], index=metadata.index)
    keep = metadata["n_cells"].ge(int(model["cell_floor"]))
    metadata = metadata.loc[keep].copy()
    counts = counts.loc[metadata.index].copy()
    expected = [str(x) for x in model["expected_timepoints"]]
    arms = metadata.groupby("timepoint")["group"].agg(lambda values: set(values.astype(str)))
    complete = [tp for tp in expected if arms.get(tp, set()) == {"decidual", "fibroblast"}]
    if complete != expected:
        fail(f"X-Mens matched timepoints changed: expected {expected}, observed {complete}")
    keep = metadata["timepoint"].isin(expected)
    metadata = metadata.loc[keep].copy()
    counts = counts.loc[metadata.index].copy()
    metadata["timepoint"] = pd.Categorical(metadata["timepoint"], categories=expected, ordered=True)
    metadata["_row_group_order"] = metadata["group"].map({"decidual": 0, "fibroblast": 1})
    if metadata["_row_group_order"].isna().any():
        fail("Unexpected X-Mens broad group")
    metadata = metadata.sort_values(["timepoint", "_row_group_order"], kind="stable").drop(
        columns="_row_group_order"
    )
    counts = counts.loc[metadata.index]
    metadata["group"] = pd.Categorical(
        metadata["group"], categories=["fibroblast", "decidual"], ordered=True
    )
    metadata = metadata.reset_index(drop=True)
    counts = counts.reset_index(drop=True)
    if len(metadata) != int(model["expected_pseudobulks"]):
        fail(f"X-Mens pseudobulk count changed: {len(metadata)}")
    if metadata.duplicated(["timepoint", "group"]).any():
        fail("X-Mens technical libraries were not collapsed to one timepoint x group pseudobulk")
    expected_cells = {
        (str(timepoint_key), str(group_key)): int(value)
        for timepoint_key, groups in model["expected_pseudobulk_cells"].items()
        for group_key, value in groups.items()
    }
    observed_cells = {
        (str(row.timepoint), str(row.group)): int(row.n_cells)
        for row in metadata.itertuples(index=False)
    }
    if observed_cells != expected_cells:
        fail(f"X-Mens pseudobulk cell counts changed: {observed_cells}")
    expected_sample_n = {str(key): int(value) for key, value in model["expected_pseudobulk_sample_n"].items()}
    for row in metadata.itertuples(index=False):
        if int(row.n_samples) != expected_sample_n[str(row.timepoint)]:
            fail(f"X-Mens sample provenance changed for {row.timepoint}/{row.group}")

    # Same arm-aware majority rule as HECA F (uniform across species, 2026-09-02).
    cpm_threshold = float(model["eligibility_cpm_threshold"])
    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    numerator_mask = metadata["group"].astype(str).eq("decidual").to_numpy()
    numerator_majority = (int(numerator_mask.sum()) + 1) // 2
    denominator_majority = (int((~numerator_mask).sum()) + 1) // 2
    eligible = (
        cpm.loc[numerator_mask].ge(cpm_threshold).sum(axis=0).ge(numerator_majority)
        | cpm.loc[~numerator_mask].ge(cpm_threshold).sum(axis=0).ge(denominator_majority)
    )
    carrier = pd.concat(
        [
            features.reset_index(drop=True),
            pd.DataFrame(
                {
                    "carrier": True,
                    "production_eligible": eligible.to_numpy(),
                    "total_count": counts.sum(axis=0).to_numpy(dtype=np.int64),
                    "nonzero_pseudobulks": counts.gt(0).sum(axis=0).to_numpy(dtype=np.int64),
                }
            ),
        ],
        axis=1,
    )
    counts_out = counts.copy()
    counts_out.insert(0, "pseudobulk_id", metadata["pseudobulk_id"].to_numpy())
    write_tsv(counts_out, workdir / "xmens_broad_carrier_counts.tsv")
    write_tsv(metadata.assign(timepoint=metadata["timepoint"].astype(str), group=metadata["group"].astype(str)), workdir / "xmens_broad_metadata.tsv")
    return carrier, counts, metadata, eligible


def fit_xmens_broad(
    carrier: pd.DataFrame,
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    eligible: pd.Series,
    *,
    model: dict[str, Any],
    workdir: Path,
    repo: Path,
    rscript: str,
) -> pd.DataFrame:
    """Fit the X-Mens broad contrast with the same limma engine as HECA F.

    Engine convergence (2026-09-02, user decision): both species use
    voomLmFit(sample.weights = TRUE). The mouse design has no blocking term
    because the object carries no animal identifier, so this remains
    timepoint-level inference. The former PyDESeq2 fit is in git history
    (releases up to cagri_ivd_end_to_end_20260815_ortholog114).
    """
    if str(model["model_formula"]) != "~ timepoint + group":
        fail("The frozen X-Mens formula must be ~ timepoint + group")
    genes = carrier.loc[eligible.to_numpy(), "ensembl_gene_id"].astype(str).tolist()
    fit_counts = counts.loc[:, eligible.to_numpy()].copy()
    fit_counts.columns = genes
    fit_counts_out = fit_counts.copy()
    fit_counts_out.insert(0, "pseudobulk_id", metadata["pseudobulk_id"].to_numpy())
    write_tsv(fit_counts_out, workdir / "xmens_broad_fit_counts.tsv")
    output = workdir / "xmens_broad_fit.tsv"
    command = [
        rscript,
        str(repo / "fit_xmens_broad.R"),
        str(workdir / "xmens_broad_fit_counts.tsv"),
        str(workdir / "xmens_broad_metadata.tsv"),
        str(output),
        str(workdir / "xmens_broad_sample_weights.tsv"),
    ]
    subprocess.run(command, check=True)
    fitted = pd.read_csv(output, sep="\t")
    if fitted["ensembl_gene_id"].duplicated().any():
        fail("X-Mens fit contains duplicated Ensembl IDs")
    if fitted["ensembl_gene_id"].astype(str).tolist() != genes:
        fail("X-Mens fit IDs/order differ from the eligible input")
    result = carrier.merge(fitted, on="ensembl_gene_id", how="left", validate="one_to_one")
    ineligible_with_fit = ~result["production_eligible"] & result["log2FoldChange_unshrunk"].notna()
    if ineligible_with_fit.any():
        fail("X-Mens fit returned statistics for ineligible genes")
    fit_estimable = result["estimable"].astype("boolean").fillna(False).astype(bool)
    result["estimable"] = result["production_eligible"] & fit_estimable
    inference_columns = ["log2FoldChange_unshrunk", "lfcSE", "stat", "pvalue", "padj"]
    result.loc[~result["estimable"], inference_columns] = np.nan
    result["padj_for_call"] = result["padj"].where(np.isfinite(result["padj"]), 1.0)
    result["significant_fdr05_lfc1"] = (
        result["estimable"]
        & result["padj_for_call"].lt(0.05)
        & result["log2FoldChange_unshrunk"].abs().gt(1.0)
    )
    result["source_id"] = "mouse_all_decidual_vs_fibroblast_timepoint_blocked"
    result["positive_lfc_meaning"] = "decidual_minus_fibroblast"
    result["analysis_method"] = "limma_voomLmFit_sample_weights"
    result["n_pseudobulks_numerator"] = int((metadata["group"].astype(str) == "decidual").sum())
    result["n_pseudobulks_denominator"] = int((metadata["group"].astype(str) == "fibroblast").sum())
    return result


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


def build_heca_support(
    consensus_path: Path,
    mapping_path: Path,
    heca: pd.DataFrame,
    calls: dict[str, Any],
) -> pd.DataFrame:
    consensus = pd.read_csv(consensus_path)
    mapping = pd.read_csv(mapping_path)
    if len(consensus) != 123 or consensus["gene"].nunique() != 123:
        fail("Consensus input is not the frozen 123-gene set")
    if len(mapping) != 123 or mapping["human_gene"].tolist() != consensus["gene"].tolist():
        fail("Consensus-to-HECA mapping differs from the 123-gene source/order")
    joined = consensus.merge(mapping, left_on="gene", right_on="human_gene", validate="one_to_one")
    de_columns = [
        "ensembl_gene_id",
        "gene_symbol",
        "carrier",
        "production_eligible",
        "estimable",
        "log2FoldChange_unshrunk",
        "lfcSE",
        "pvalue",
        "padj",
        "padj_for_call",
    ]
    joined = joined.merge(
        heca[de_columns],
        left_on="human_ensembl_gene_id",
        right_on="ensembl_gene_id",
        how="left",
        validate="many_to_one",
    )
    rows = []
    for row in joined.to_dict("records"):
        called, observed, matches, status = classify_support(
            row["consensus_direction"],
            row.get("production_eligible"),
            row.get("estimable"),
            row.get("log2FoldChange_unshrunk"),
            row.get("padj_for_call"),
            fdr=float(calls["padj_lt"]),
            abs_lfc=float(calls["abs_unshrunk_lfc_gt"]),
        )
        rows.append(
            {
                "human_gene": row["gene"],
                "consensus_direction": row["consensus_direction"],
                "tier": row["tier"],
                "n_mature": row["n_mature"],
                "mature_papers": row["mature_papers"],
                "mature_biphasic": bool_value(row["mature_biphasic"]),
                "lucas_onset_dir": row["lucas_onset_dir"],
                "onset_agrees_mature": "" if pd.isna(row["onset_agrees_mature"]) else bool_value(row["onset_agrees_mature"]),
                "human_ensembl_gene_id": row["human_ensembl_gene_id"],
                "current_gene_symbol": row["current_gene_symbol"],
                "current_mapping_status": row["current_mapping_status"],
                "source_id": "F_population25_production",
                "carrier": False if pd.isna(row.get("carrier")) else bool_value(row.get("carrier")),
                "production_eligible": False if pd.isna(row.get("production_eligible")) else bool_value(row.get("production_eligible")),
                "estimable": False if pd.isna(row.get("estimable")) else bool_value(row.get("estimable")),
                "log2FoldChange_unshrunk": row.get("log2FoldChange_unshrunk"),
                "lfcSE": row.get("lfcSE"),
                "pvalue": row.get("pvalue"),
                "padj": row.get("padj"),
                "padj_for_call": row.get("padj_for_call"),
                "observed_direction": observed,
                "strict_lfc1_call": called,
                "direction_match": matches,
                "support_status": status,
                "source_state_detail": "",
            }
        )
    result = pd.DataFrame(rows)
    if set(result["support_status"]) - VALID_SUPPORT:
        fail("Unexpected HECA support status")
    return result


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


def assign_zone_labels(
    subtype: pd.Series,
    coarse: pd.Series,
    zone_by_subtype: dict[str, str],
    decidual_coarse_value: str,
    nondecidual_value: str,
) -> pd.Series:
    decidual = coarse.astype(str).eq(decidual_coarse_value)
    subtype = subtype.astype(str)
    unknown = sorted(set(subtype[decidual]) - set(zone_by_subtype))
    if unknown:
        fail(f"Unknown fine_refined values among decidual cells: {unknown}")
    outside_decidual_with_zone_label = (~decidual) & subtype.isin(zone_by_subtype)
    if outside_decidual_with_zone_label.any():
        fail("A zone-mapped fine subtype occurs outside the coarse decidual population")
    zone_order = ["not_applicable", "secondary_like", "primary_like", "zone_unresolved"]
    zone = pd.Series(nondecidual_value, index=subtype.index, dtype="string")
    zone.loc[decidual] = subtype.loc[decidual].map(zone_by_subtype)
    return pd.Series(pd.Categorical(zone, categories=zone_order, ordered=True), index=subtype.index)


def build_zone_outputs(
    object_path: Path,
    object_spec: dict[str, Any],
    mapping_path: Path,
    markers_path: Path,
    zone_config: dict[str, Any],
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    obs, var_names = read_object_metadata(
        object_path,
        [zone_config["source_annotation_column"], zone_config["coarse_state_field"], "sample", "timepoint"],
    )
    mapping = pd.read_csv(mapping_path, sep="\t")
    markers = pd.read_csv(markers_path, sep="\t")
    if list(mapping.columns) != ["fine_refined", "putative_decidual_zone", "interpretation"]:
        fail("Unexpected zone-mapping schema")
    if mapping["fine_refined"].astype(str).tolist() != ["dec1", "dec2", "dec3", "decAp"]:
        fail("Zone mapping must preserve the frozen decidual-subtype order")
    if mapping["putative_decidual_zone"].astype(str).tolist() != [
        "secondary_like", "secondary_like", "primary_like", "zone_unresolved"
    ]:
        fail("Zone mapping values changed")
    if mapping["fine_refined"].duplicated().any() or len(mapping) != 4:
        fail("Zone mapping must contain four unique decidual subtypes")
    if list(markers.columns) != ["module", "rank", "gene"]:
        fail("Unexpected zone-marker schema")
    marker_contract = markers.groupby("module", sort=False)["rank"].apply(list).to_dict()
    if marker_contract != {"PDZ": list(range(1, 21)), "SDZ": list(range(1, 21))}:
        fail("Zone markers must contain ranked PDZ 1..20 followed by SDZ 1..20")
    if len(markers) != 40 or markers.duplicated(["module", "rank"]).any() or markers["gene"].duplicated().any():
        fail("Zone marker resource must contain 40 unique ranked genes")
    subtype = obs[zone_config["source_annotation_column"]].astype(str)
    coarse = obs[zone_config["coarse_state_field"]].astype(str)
    zone_by_subtype = mapping.set_index("fine_refined")["putative_decidual_zone"].to_dict()
    zone = assign_zone_labels(
        subtype,
        coarse,
        zone_by_subtype,
        str(zone_config["decidual_coarse_value"]),
        str(zone_config["nondecidual_value"]),
    )
    cell_table = pd.DataFrame(
        {
            "cell_id": obs.index,
            "sample": obs["sample"].astype(str),
            "timepoint": obs["timepoint"].astype(str),
            "coarse_refined": coarse,
            "fine_refined": subtype,
            "putative_decidual_zone": zone,
        }
    )
    observed_zone_counts = cell_table["putative_decidual_zone"].value_counts().to_dict()
    expected_zone_counts = {str(key): int(value) for key, value in zone_config["expected_counts"].items()}
    if observed_zone_counts != expected_zone_counts:
        fail(f"Zone counts changed: {observed_zone_counts}")

    subtype_order = mapping["fine_refined"].astype(str).tolist()
    subtype_to_group = {value: index for index, value in enumerate(subtype_order)}
    assignments = np.full(len(obs), -1, dtype=np.int64)
    for position, value in enumerate(subtype):
        assignments[position] = subtype_to_group.get(value, -1)
    identity = np.arange(len(var_names), dtype=np.int64)
    subtype_counts, subtype_cell_n, _ = stream_canonical_pseudobulk(
        object_path,
        str(object_spec["counts_layer"]),
        assignments,
        subtype_order,
        identity,
    )
    expected_subtype_counts = {
        str(key): int(value) for key, value in zone_config["expected_subtype_counts"].items()
    }
    observed_subtype_counts = {
        subtype_order[index]: int(value) for index, value in enumerate(subtype_cell_n)
    }
    if observed_subtype_counts != expected_subtype_counts:
        fail(f"Decidual subtype counts changed: {observed_subtype_counts}")
    library_sizes = subtype_counts.sum(axis=1).astype(float)
    if not np.isfinite(library_sizes).all() or np.any(library_sizes <= 0):
        fail("Zone subtype library sizes must be finite and positive")
    expected_library_sizes = {
        str(key): int(value)
        for key, value in zone_config["expected_subtype_library_sizes"].items()
    }
    observed_library_sizes = {
        subtype_order[index]: int(value) for index, value in enumerate(library_sizes)
    }
    if observed_library_sizes != expected_library_sizes:
        fail(f"Zone subtype library sizes changed: {observed_library_sizes}")
    marker_rows = []
    heat = []
    for marker in markers.sort_values(["module", "rank"], kind="stable").to_dict("records"):
        gene = str(marker["gene"])
        if gene not in var_names:
            fail(f"Zone marker is absent from the X-Mens object: {gene}")
        index = int(var_names.get_loc(gene))
        raw = subtype_counts[:, index].astype(float)
        log_cpm = np.log2(1.0 + 1_000_000.0 * raw / library_sizes)
        sd = float(log_cpm.std(ddof=0))
        zscore = np.zeros_like(log_cpm) if sd == 0 else (log_cpm - log_cpm.mean()) / sd
        record: dict[str, Any] = {"marker_module": marker["module"], "rank": int(marker["rank"]), "gene": gene}
        for position, value in enumerate(subtype_order):
            record[f"{value}_summed_count"] = int(raw[position])
            record[f"{value}_cpm"] = float(1_000_000.0 * raw[position] / library_sizes[position])
            record[f"{value}_log2_cpm_plus_1"] = float(log_cpm[position])
            record[f"{value}_gene_zscore"] = float(zscore[position])
        marker_rows.append(record)
        heat.append(zscore)
    heatmap_values = pd.DataFrame(marker_rows)
    matrix = np.vstack(heat)
    if not np.isfinite(matrix).all():
        fail("Zone heatmap contains a non-finite value")
    fig, ax = plt.subplots(figsize=(7.8, 14.0))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(subtype_order)))
    zone_labels = mapping.set_index("fine_refined")["putative_decidual_zone"]
    ax.set_xticklabels([f"{value}\n{zone_labels[value].replace('_', '-')}" for value in subtype_order])
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(markers.sort_values(["module", "rank"], kind="stable")["gene"], fontsize=8)
    pdz_n = int((markers["module"] == "PDZ").sum())
    ax.axhline(pdz_n - 0.5, color="#333333", linewidth=1.2)
    ax.set_title("Putative primary- and secondary-zone marker expression\nPooled-cell log2(CPM+1), standardized within gene (descriptive)")
    ax.set_xlabel("X-Mens decidual subtype and inferred zone")
    ax.set_ylabel("Frozen zone-marker gene")
    bar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    bar.set_label("Gene-wise z-score")
    fig.tight_layout()
    fig.savefig(output_dir / "decidual_zone_marker_heatmap.png", dpi=240, bbox_inches="tight")
    fig.savefig(output_dir / "decidual_zone_marker_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    write_csv(cell_table, output_dir / "decidual_zone_cell_metadata.csv")
    write_csv(heatmap_values, output_dir / "decidual_zone_marker_heatmap_values.csv")
    zone_counts = (
        cell_table.groupby(["putative_decidual_zone", "fine_refined"], observed=True)
        .size()
        .rename("cell_n")
        .reset_index()
    )
    if int(zone_counts["cell_n"].sum()) != len(obs):
        fail("Zone count table does not reconcile to the object row count")
    write_csv(zone_counts, output_dir / "decidual_zone_counts.csv")
    return cell_table, zone_counts, heatmap_values


def one_value(frame: pd.DataFrame, column: str) -> Any:
    values = frame[column].dropna().unique()
    if len(values) != 1:
        fail(f"Expected one nonmissing value for {column}; observed {values.tolist()}")
    return values[0]


def assert_close(observed: Any, expected: Any, label: str, *, tolerance: float = 1e-10) -> None:
    if not math.isclose(float(observed), float(expected), rel_tol=tolerance, abs_tol=tolerance):
        fail(f"{label} changed: expected {expected}, observed {observed}")


def validate_heca_result(
    result: pd.DataFrame,
    metadata: pd.DataFrame,
    model: dict[str, Any],
) -> list[dict[str, Any]]:
    if result["ensembl_gene_id"].duplicated().any():
        fail("HECA carrier result contains duplicated Ensembl IDs")
    exact = {
        "carrier_genes": (len(result), int(model["expected_carrier_genes"])),
        "eligible_genes": (int(result["production_eligible"].sum()), int(model["expected_eligible_genes"])),
        "estimable_genes": (int(result["estimable"].sum()), int(model["expected_eligible_genes"])),
        "carrier_count_sum": (int(result["total_count"].sum()), int(model["expected_carrier_count_sum"])),
        "fit_count_sum": (
            int(result.loc[result["production_eligible"], "total_count"].sum()),
            int(model["expected_fit_count_sum"]),
        ),
        "fdr05_calls": (int((result["estimable"] & result["padj_for_call"].lt(0.05)).sum()), int(model["expected_fdr05_calls"])),
        "strict_calls": (int(result["significant_fdr05_lfc1"].sum()), int(model["expected_strict_calls"])),
        "pseudobulks": (len(metadata), int(model["expected_pseudobulks"])),
        "people": (metadata["person"].nunique(), int(model["expected_people"])),
    }
    checks: list[dict[str, Any]] = []
    for name, (observed, expected) in exact.items():
        if int(observed) != int(expected):
            fail(f"HECA {name} changed: expected {expected}, observed {observed}")
        checks.append({"check": f"heca_{name}", "observed": observed, "expected": expected, "pass": True})
    numeric = {
        "duplicate_correlation_final": "expected_duplicate_correlation_final",
        "design_rank": "expected_design_rank",
        "residual_df": "expected_residual_df",
        "design_condition_number": "expected_design_condition_number",
        "fit_library_min": "expected_fit_library_min",
        "fit_library_median": "expected_fit_library_median",
        "fit_library_max": "expected_fit_library_max",
        "tmm_factor_min": "expected_tmm_factor_min",
        "tmm_factor_median": "expected_tmm_factor_median",
        "tmm_factor_max": "expected_tmm_factor_max",
        "voom_weight_min": "expected_voom_weight_min",
        "voom_weight_median": "expected_voom_weight_median",
        "voom_weight_max": "expected_voom_weight_max",
        "sample_weight_min": "expected_sample_weight_min",
        "sample_weight_median": "expected_sample_weight_median",
        "sample_weight_max": "expected_sample_weight_max",
    }
    for column, config_key in numeric.items():
        observed = one_value(result, column)
        expected = model[config_key]
        assert_close(observed, expected, f"HECA {column}")
        checks.append({"check": f"heca_{column}", "observed": observed, "expected": expected, "pass": True})
    if int(one_value(result, "design_column_n")) != int(model["expected_design_rank"]):
        fail("HECA design column count changed")
    if int(one_value(result, "warning_count")) != 0:
        fail("HECA limma fit emitted a warning")
    estimable_padj = result.loc[result["estimable"], "padj"]
    if not estimable_padj.dropna().between(0, 1, inclusive="both").all():
        fail("HECA estimable adjusted p-values are outside [0,1]")
    return checks


def validate_xmens_result(
    result: pd.DataFrame,
    metadata: pd.DataFrame,
    model: dict[str, Any],
) -> list[dict[str, Any]]:
    if result["ensembl_gene_id"].duplicated().any():
        fail("X-Mens carrier result contains duplicated Ensembl IDs")
    exact = {
        "carrier_genes": (len(result), int(model["expected_carrier_genes"])),
        "eligible_genes": (int(result["production_eligible"].sum()), int(model["expected_eligible_genes"])),
        "estimable_genes": (int(result["estimable"].sum()), int(model["expected_eligible_genes"])),
        "carrier_count_sum": (int(result["total_count"].sum()), int(model["expected_carrier_count_sum"])),
        "fit_count_sum": (
            int(result.loc[result["production_eligible"], "total_count"].sum()),
            int(model["expected_fit_count_sum"]),
        ),
        "strict_calls": (int(result["significant_fdr05_lfc1"].sum()), int(model["expected_strict_calls"])),
        "pseudobulks": (len(metadata), int(model["expected_pseudobulks"])),
    }
    checks: list[dict[str, Any]] = []
    for name, (observed, expected) in exact.items():
        if int(observed) != int(expected):
            fail(f"X-Mens {name} changed: expected {expected}, observed {observed}")
        checks.append({"check": f"xmens_{name}", "observed": observed, "expected": expected, "pass": True})
    assert_close(one_value(result, "design_rank"), model["expected_design_rank"], "X-Mens design rank")
    assert_close(
        one_value(result, "design_condition_number"),
        model["expected_design_condition_number"],
        "X-Mens design condition number",
    )
    numeric = {
        "fit_library_min": "expected_fit_library_min",
        "fit_library_median": "expected_fit_library_median",
        "fit_library_max": "expected_fit_library_max",
        "tmm_factor_min": "expected_tmm_factor_min",
        "tmm_factor_median": "expected_tmm_factor_median",
        "tmm_factor_max": "expected_tmm_factor_max",
        "sample_weight_min": "expected_sample_weight_min",
        "sample_weight_median": "expected_sample_weight_median",
        "sample_weight_max": "expected_sample_weight_max",
        "voom_weight_min": "expected_voom_weight_min",
        "voom_weight_median": "expected_voom_weight_median",
        "voom_weight_max": "expected_voom_weight_max",
    }
    for column, config_key in numeric.items():
        assert_close(one_value(result, column), model[config_key], f"X-Mens {column}")
        checks.append({"check": f"xmens_{column}", "observed": one_value(result, column), "expected": model[config_key], "pass": True})
    if int(one_value(result, "warning_count")) != 0:
        fail("X-Mens limma fit emitted a warning")
    estimable_padj = result.loc[result["estimable"], "padj"]
    if not estimable_padj.dropna().between(0, 1, inclusive="both").all():
        fail("X-Mens estimable adjusted p-values are outside [0,1]")
    return checks


def normalized_partition(series: pd.Series, categories: Iterable[str]) -> dict[str, int]:
    counts = series.astype(str).value_counts().to_dict()
    return {str(category): int(counts.get(str(category), 0)) for category in categories}


def validate_partition(
    label: str,
    series: pd.Series,
    expected: dict[str, Any],
) -> list[dict[str, Any]]:
    expected_int = {str(key): int(value) for key, value in expected.items()}
    unexpected = sorted(set(series.astype(str)) - set(expected_int))
    if unexpected:
        fail(f"{label} has unexpected status values: {unexpected}")
    observed = normalized_partition(series, expected_int)
    if observed != expected_int:
        fail(f"{label} status partition changed: expected {expected_int}, observed {observed}")
    return [
        {
            "check": f"{label}_{status}",
            "observed": observed[status],
            "expected": expected_int[status],
            "pass": True,
        }
        for status in expected_int
    ]


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


def validate_wang_mapping(rows: pd.DataFrame, model: dict[str, Any]) -> list[dict[str, Any]]:
    matched = int(rows["xmens_match_candidate_n"].gt(0).sum())
    unmatched = int(rows["xmens_match_candidate_n"].eq(0).sum())
    estimable = int(rows["xmens_estimable"].map(bool_value).sum())
    alias_n = int(rows["xmens_match_route"].eq("source_feature_id_alias").sum())
    called = int(rows["xmens_called_strict_lfc1"].map(bool_value).sum())
    exact = {
        "wang_matched": (matched, int(model["expected_wang_match_n"])),
        "wang_unmatched": (unmatched, int(model["expected_wang_unmatched_n"])),
        "wang_estimable": (estimable, int(model["expected_wang_estimable_n"])),
        "wang_alias_route": (alias_n, int(model["expected_wang_alias_route_n"])),
        "wang_strict_called": (called, int(model["expected_wang_strict_called_n"])),
    }
    checks = []
    for label, (observed, expected) in exact.items():
        if observed != expected:
            fail(f"{label} changed: expected {expected}, observed {observed}")
        checks.append({"check": label, "observed": observed, "expected": expected, "pass": True})
    alias_expected = {
        "2310039H08Rik": "tested_not_called",
        "Arntl": "tested_not_called",
        "Bves": "supported",
        "C77080": "supported",
        "Gyg": "tested_not_called",
        "Mettl7a1": "supported",
    }
    alias_rows = rows.set_index("wang_ivd_gene")
    observed_alias = {
        gene: str(alias_rows.loc[gene, "xmens_status"])
        for gene in alias_expected
        if str(alias_rows.loc[gene, "xmens_match_route"]) == "source_feature_id_alias"
    }
    if observed_alias != alias_expected:
        fail(f"X-Mens alias-only evidence changed: {observed_alias}")
    nrg1 = alias_rows.loc["Nrg1"]
    if int(nrg1["xmens_match_candidate_n"]) != 2 or str(nrg1["xmens_selected_ensembl_gene_id"]) != "ENSMUSG00000062991":
        fail("Nrg1 collision resolution changed")
    return checks


def write_zone_derivative(
    source_path: Path,
    cell_table: pd.DataFrame,
    mapping_path: Path,
    markers_path: Path,
    output_path: Path,
) -> None:
    source = ad.read_h5ad(source_path)
    source_shape = source.shape
    source_obs_names = source.obs_names.astype(str).tolist()
    source_var_names = source.var_names.astype(str).tolist()
    source_layers = sorted(source.layers.keys())
    source_obsm = sorted(source.obsm.keys())
    source_obsp = sorted(source.obsp.keys())
    source_varm = sorted(source.varm.keys())
    source_varp = sorted(source.varp.keys())
    source_uns = sorted(source.uns.keys())
    source_raw_shape = None if source.raw is None else source.raw.shape
    source_x_shape = None if source.X is None else source.X.shape
    source_obs_columns = source.obs.columns.astype(str).tolist()
    try:
        if "putative_decidual_zone" in source.obs.columns:
            fail("Source X-Mens object already contains putative_decidual_zone")
        if "putative_decidual_zone_annotation" in source.uns:
            fail("Source X-Mens object already contains putative_decidual_zone_annotation")
        if source_obs_names != cell_table["cell_id"].astype(str).tolist():
            fail("Zone cell table order differs from the source object")
        zone_order = ["not_applicable", "secondary_like", "primary_like", "zone_unresolved"]
        source.obs["putative_decidual_zone"] = pd.Categorical(
            cell_table["putative_decidual_zone"].astype(str).to_numpy(),
            categories=zone_order,
            ordered=True,
        )
        source.uns["putative_decidual_zone_annotation"] = {
            "schema_version": "1.0",
            "interpretation": "Literature-informed subtype analogy; not a measured embryo-distance label.",
            "zone_mapping_sha256": sha256_file(mapping_path),
            "zone_markers_sha256": sha256_file(markers_path),
            "heatmap_method": "pooled raw-count subtype totals; log2(CPM+1); gene-wise population z-score",
        }
        source.write_h5ad(output_path, compression="gzip")
    finally:
        del source
    check = ad.read_h5ad(output_path, backed="r")
    try:
        if check.shape != source_shape:
            fail("Annotated X-Mens object shape changed")
        if check.obs_names.astype(str).tolist() != source_obs_names:
            fail("Annotated X-Mens object observation order changed")
        if check.var_names.astype(str).tolist() != source_var_names:
            fail("Annotated X-Mens object variable order changed")
        if sorted(check.layers.keys()) != source_layers or sorted(check.obsm.keys()) != source_obsm:
            fail("Annotated X-Mens object layers/embeddings changed")
        if sorted(check.obsp.keys()) != source_obsp or sorted(check.varm.keys()) != source_varm:
            fail("Annotated X-Mens object graph/variable-matrix schema changed")
        if sorted(check.varp.keys()) != source_varp:
            fail("Annotated X-Mens object variable-pair schema changed")
        expected_uns = sorted(source_uns + ["putative_decidual_zone_annotation"])
        if sorted(check.uns.keys()) != expected_uns:
            fail("Annotated X-Mens object unstructured metadata schema changed unexpectedly")
        check_raw_shape = None if check.raw is None else check.raw.shape
        check_x_shape = None if check.X is None else check.X.shape
        if check_raw_shape != source_raw_shape or check_x_shape != source_x_shape:
            fail("Annotated X-Mens object X/raw shape changed")
        expected_obs_columns = source_obs_columns + ["putative_decidual_zone"]
        if check.obs.columns.astype(str).tolist() != expected_obs_columns:
            fail("Annotated X-Mens object observation schema changed unexpectedly")
        observed = check.obs["putative_decidual_zone"].astype(str).value_counts().to_dict()
        expected = cell_table["putative_decidual_zone"].astype(str).value_counts().to_dict()
        if observed != expected:
            fail("Annotated X-Mens zone column failed readback")
        if check.obs["putative_decidual_zone"].astype(str).tolist() != cell_table[
            "putative_decidual_zone"
        ].astype(str).tolist():
            fail("Annotated X-Mens zone assignments changed by cell")
    finally:
        check.file.close()


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


def artifact_inventory(root: Path) -> tuple[list[dict[str, Any]], str]:
    files = sorted(path for path in root.rglob("*") if path.is_file() and path.name != "release_manifest.json")
    entries = []
    canonical_lines = []
    for path in files:
        relative = path.relative_to(root).as_posix()
        digest = sha256_file(path)
        size = path.stat().st_size
        entries.append({"relative_path": relative, "sha256": digest, "size_bytes": size})
        canonical_lines.append(f"{relative}\t{digest}\t{size}\n")
    tree = hashlib.sha256("".join(canonical_lines).encode("utf-8")).hexdigest()
    return entries, tree


def write_release_manifest(
    root: Path,
    config_path: Path,
    object_paths: dict[str, Path],
    object_hashes: dict[str, str],
    execution_hashes: dict[str, str],
) -> None:
    entries, tree = artifact_inventory(root)
    frozen_config = read_config(config_path)
    manifest = {
        "schema_version": "1.0",
        "release_status": "complete_immutable",
        "analysis_id": "cagri_ivd_end_to_end_20260815",
        "created_utc": utc_now(),
        "artifact_tree_algorithm": "sha256(sorted(relative_path\\tsha256\\tsize_bytes\\n)); manifest excluded",
        "artifact_tree_sha256": tree,
        "config_sha256": execution_hashes["config"],
        "pipeline_sha256": execution_hashes["pipeline"],
        "heca_fitter_sha256": execution_hashes["heca_fitter"],
        "xmens_fitter_sha256": execution_hashes["xmens_fitter"],
        "object_hashes_verified": True,
        "objects": {
            key: {
                "path": str(path),
                "sha256": object_hashes[key],
                "hash_verified": True,
            }
            for key, path in object_paths.items()
        },
        "resource_inputs": frozen_config["resources"],
        "files": entries,
    }
    path = root / "release_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def verify_completed_release(root: Path) -> None:
    manifest_path = root / "release_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("release_status") != "complete_immutable":
        fail("Release manifest is not complete_immutable")
    entries, tree = artifact_inventory(root)
    if entries != manifest.get("files") or tree != manifest.get("artifact_tree_sha256"):
        fail("Completed release inventory/tree verification failed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--heca-object", type=Path)
    parser.add_argument("--xmens-object", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--rscript", default=shutil.which("Rscript") or "Rscript")
    parser.add_argument("--n-cpus", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if ad is None or h5py is None:
        fail("Install the pinned Python requirements before running the analysis")
    repo = args.repo_root.resolve()
    config_path = (args.config or repo / "config.yaml").resolve()
    config = read_config(config_path)
    if config["analysis_id"] != "cagri_ivd_end_to_end_20260815":
        fail("Unexpected analysis_id")
    if platform.python_version() != "3.11.11":
        fail(f"Expected the pinned Python 3.11.11; observed {platform.python_version()}")
    execution_inputs = {
        "config": config_path,
        "pipeline": repo / "run_pipeline.py",
        "heca_fitter": repo / "fit_heca_f.R",
        "xmens_fitter": repo / "fit_xmens_broad.R",
    }
    execution_hashes_at_start = {
        key: sha256_file(path) for key, path in execution_inputs.items()
    }
    resource_paths = {
        key: resolve_resource(repo, spec) for key, spec in config["resources"].items()
    }
    heca_object = (args.heca_object or Path(config["objects"]["heca_cells"]["path"])).resolve()
    xmens_object = (args.xmens_object or Path(config["objects"]["xmens_mouse"]["path"])).resolve()
    verify_object(heca_object, config["objects"]["heca_cells"], check_hash=True)
    verify_object(xmens_object, config["objects"]["xmens_mouse"], check_hash=True)

    output = args.output_dir or Path(config["output_directory"])
    output = output if output.is_absolute() else repo / output
    output = output.resolve()
    staging = output.with_name(output.name + ".inprogress")
    if output.exists() or staging.exists():
        fail(f"Refusing existing output/staging path: {output} or {staging}")
    staging.mkdir(parents=True)

    def _discard_staging() -> None:
        # Otherwise any failure below blocks every subsequent rerun on the
        # "Refusing existing output/staging path" gate. On success the directory has
        # already been renamed, so this is a no-op. (Does not cover SIGKILL or an
        # uncaught SIGTERM, e.g. a Slurm timeout.) PENDING_PIPELINE_FIXES fix 3.
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)

    atexit.register(_discard_staging)
    intermediate = staging / "intermediate"
    intermediate.mkdir()
    checks: list[dict[str, Any]] = []

    heca_carrier, heca_metadata = prepare_heca_f(
        heca_object,
        config["objects"]["heca_cells"],
        resource_paths["heca_feature_manifest"],
        config["heca_f"],
        intermediate,
    )
    heca_result = fit_heca_f(heca_carrier, heca_metadata, intermediate, repo, args.rscript)
    checks.extend(validate_heca_result(heca_result, heca_metadata, config["heca_f"]))

    xmens_carrier, xmens_counts, xmens_metadata, xmens_eligible = prepare_xmens_broad(
        xmens_object,
        config["objects"]["xmens_mouse"],
        resource_paths["xmens_feature_manifest"],
        config["xmens_broad"],
        intermediate,
    )
    xmens_result = fit_xmens_broad(
        xmens_carrier,
        xmens_counts,
        xmens_metadata,
        xmens_eligible,
        model=config["xmens_broad"],
        workdir=intermediate,
        repo=repo,
        rscript=args.rscript,
    )
    checks.extend(validate_xmens_result(xmens_result, xmens_metadata, config["xmens_broad"]))

    heca_result.to_parquet(staging / "heca_F_all_genes.parquet", index=False)
    xmens_result.to_parquet(staging / "xmens_broad_all_genes.parquet", index=False)
    write_csv(heca_result, staging / "heca_F_all_genes.csv")
    write_csv(xmens_result, staging / "xmens_broad_all_genes.csv")

    calls = config["support_calls"]
    heca_support = build_heca_support(
        resource_paths["consensus123"],
        resource_paths["consensus123_heca_mapping"],
        heca_result,
        calls,
    )
    consensus_ortholog_support = build_consensus_ortholog_support(
        resource_paths["consensus123"],
        resource_paths["consensus123_mouse_orthology"],
        xmens_result,
        calls,
        int(config["acceptance"]["headline_supported_counts"]["consensus_ortholog_vs_xmens_broad"][1]),
    )
    wang_support = build_wang_support(
        resource_paths["wang_ivd"],
        resource_paths["wang_oil"],
        resource_paths["wang_pregnancy"],
        xmens_result,
        calls,
    )
    expected_partitions = config["acceptance"]["status_partitions"]
    checks.extend(validate_partition("consensus123_vs_heca_F", heca_support["support_status"], expected_partitions["consensus123_vs_heca_F"]))
    checks.extend(validate_partition("consensus_ortholog_vs_xmens_broad", consensus_ortholog_support["support_status"], expected_partitions["consensus_ortholog_vs_xmens_broad"]))
    checks.extend(validate_partition("wang_ivd_vs_oil", wang_support["wang_oil_status"], expected_partitions["wang_ivd_vs_oil"]))
    checks.extend(validate_partition("wang_ivd_vs_pregnancy", wang_support["wang_pregnancy_status"], expected_partitions["wang_ivd_vs_pregnancy"]))
    checks.extend(validate_partition("wang_ivd_vs_xmens_broad", wang_support["xmens_status"], expected_partitions["wang_ivd_vs_xmens_broad"]))
    checks.extend(validate_wang_mapping(wang_support, config["xmens_broad"]))

    write_csv(heca_support, staging / "consensus123_hecaF_support_rows.csv")
    write_csv(consensus_ortholog_support, staging / "consensus_ortholog_xmens_broad_support_rows.csv")
    write_csv(wang_support, staging / "wang_ivd_support_rows.csv")
    summary = pd.DataFrame(
        [
            comparison_summary_row("consensus123_vs_heca_F", "123 human consensus IVD genes vs HECA F", heca_support["support_status"], "ineligible_or_unavailable"),
            comparison_summary_row("wang_ivd_vs_oil", "Wang IVD genes vs Wang oil/artificial in vivo", wang_support["wang_oil_status"], "ineligible_or_unmatched"),
            comparison_summary_row("wang_ivd_vs_pregnancy", "Wang IVD genes vs Wang pregnancy in vivo", wang_support["wang_pregnancy_status"], "ineligible_or_unmatched"),
            comparison_summary_row("wang_ivd_vs_xmens_broad", "Wang IVD genes vs broad X-Mens", wang_support["xmens_status"], "ineligible_or_unmatched"),
            comparison_summary_row("consensus_ortholog_vs_xmens_broad", "Consensus IVD genes with a strict 1:1 mouse ortholog vs broad X-Mens", consensus_ortholog_support["support_status"], "ineligible_or_unavailable"),
        ]
    )
    write_csv(summary, staging / "ivd_support_summary.csv")

    cell_table, zone_counts, heatmap_values = build_zone_outputs(
        xmens_object,
        config["objects"]["xmens_mouse"],
        resource_paths["zone_mapping"],
        resource_paths["zone_markers"],
        config["zones"],
        staging,
    )
    write_zone_derivative(
        xmens_object,
        cell_table,
        resource_paths["zone_mapping"],
        resource_paths["zone_markers"],
        staging / "final_mouse_with_putative_decidual_zone.h5ad",
    )

    write_csv(pd.DataFrame(checks), staging / "validation_checks.csv")
    execution_hashes_at_end = {
        key: sha256_file(path) for key, path in execution_inputs.items()
    }
    if execution_hashes_at_end != execution_hashes_at_start:
        fail(
            "Config or execution code changed during the run: "
            f"start={execution_hashes_at_start}, end={execution_hashes_at_end}"
        )
    run_metadata = {
        "analysis_id": config["analysis_id"],
        "created_utc": utc_now(),
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "analysis_n_cpus": int(args.n_cpus),
        "packages": package_versions(),
        "rscript": str(args.rscript),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "config_sha256": execution_hashes_at_start["config"],
        "pipeline_sha256": execution_hashes_at_start["pipeline"],
        "heca_fitter_sha256": execution_hashes_at_start["heca_fitter"],
        "object_hashes_verified": True,
        "annotated_object_written": True,
        "support_rule": "eligible and estimable and padj<0.05 and abs(unshrunk log2FC)>1; supported requires same direction",
        "zone_interpretation_boundary": "Literature-informed subtype analogy; no embryo-distance measurement; heatmap descriptive only.",
        "frozen_de_references_for_independent_validation": config["acceptance"].get(
            "frozen_de_references", "retired_2026-09-02_with_unified_engine_promotion"
        ),
        "xmens_fitter_sha256": execution_hashes_at_start["xmens_fitter"],
    }
    (staging / "run_metadata.json").write_text(json.dumps(run_metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if sha256_file(heca_object) != config["objects"]["heca_cells"]["sha256"]:
        fail("HECA source object changed during the run")
    if sha256_file(xmens_object) != config["objects"]["xmens_mouse"]["sha256"]:
        fail("X-Mens source object changed during the run")
    execution_hashes_before_manifest = {
        key: sha256_file(path) for key, path in execution_inputs.items()
    }
    if execution_hashes_before_manifest != execution_hashes_at_start:
        fail(
            "Config or execution code changed before release finalization: "
            f"start={execution_hashes_at_start}, final={execution_hashes_before_manifest}"
        )
    write_release_manifest(
        staging,
        config_path,
        {"heca_cells": heca_object, "xmens_mouse": xmens_object},
        {
            "heca_cells": config["objects"]["heca_cells"]["sha256"],
            "xmens_mouse": config["objects"]["xmens_mouse"]["sha256"],
        },
        execution_hashes_at_start,
    )
    staging.rename(output)
    verify_completed_release(output)
    print(json.dumps({"status": "complete_immutable", "output": str(output), "summary": summary.to_dict("records")}, indent=2))


if __name__ == "__main__":
    main()
