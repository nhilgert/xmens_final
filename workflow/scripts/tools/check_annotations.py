#!/usr/bin/env python3
"""Validate the frozen decision tables, or one annotated object, against the annotation vocabulary.

Decision-table mode (no --h5ad) checks every table's clusters and labels against its fit and the vocabulary. Object
mode (--h5ad, --species, --role) checks one object: schema (canonical columns and values, subtype parents, the
vocabulary's schema version and display names) and count gates taken from the vocabulary (expected_coarse_counts,
coarse_subtype_fill, the epithelial splits, the endothelial content tasks and the round-2 lineage label counts).
Roles all_deposit and fibdec_deposit are the release objects (data/deliverables/<release>) and their deposit copies;
there the legacy label columns the releases keep (vocabulary changelog v8 -> v9, item 5) are allowed. Rule
check_release_object (04g) runs it on every object of a cut release. Exit 0 pass, 1 schema failure, 2 count gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd
import yaml


LEGACY_OBS = {
    "coarse_refined", "fine_refined", "subtype", "annotation_level",
    "cell_type_annotations", "cell_subtype_annotations", "cell_type_subannotations",
}
# Legacy label columns the release objects keep for the earlier vocabularies; the deposit copies drop them.
RELEASE_LEGACY_OBS = {"coarse_refined", "fine_refined", "subtype", "annotation_level"}
RELEASE_ROLES = {"all_deposit", "fibdec_deposit"}
ALL_CELLS_ROLES = {"refined", "final", "all_deposit"}
FIBDEC = {"decidual_cell", "endometrial_fibroblast"}


def canonical_sets(vocab):
    types = set(vocab["allowed"])
    subtype_parent = {}
    for cell_type, subtypes in vocab["allowed"].items():
        for subtype in subtypes or []:
            if subtype in subtype_parent and subtype_parent[subtype] != cell_type:
                raise ValueError(f"duplicate subtype parent for {subtype}")
            subtype_parent[subtype] = cell_type
    return types, subtype_parent


def round2_settings(repo: Path) -> dict:
    """config annotation_round2 (04b2/04c2), or {} when the repo has none."""
    path = repo / "config/config.yaml"
    if not path.exists():
        return {}
    return (yaml.safe_load(path.read_text()) or {}).get("annotation_round2") or {}


def validate_csvs(repo: Path, species: list[str], spec: dict, vocab: dict):
    failures = []
    canonical_types, subtype_parent = canonical_sets(vocab)
    canonical = canonical_types | set(subtype_parent)
    rows_out = []
    for sp in species:
        print(f"=== {sp} decision records ===")
        pairs = [
            ("(atlas)", spec["atlas"][sp]["res"],
             repo / f"resources/annotations/{sp}_lineage.csv",
             repo / f"data/interim/latent_{sp}.h5ad", True),
        ]
        for lineage, cfg in spec["species"][sp]["lineages"].items():
            pairs.append((
                lineage, cfg.get("res", spec.get("defaults", {}).get("res")),
                repo / f"resources/annotations/{sp}_{lineage}.csv",
                repo / f"data/interim/latent_{sp}_{lineage}.h5ad", False,
            ))
        round2 = round2_settings(repo)
        if sp in round2.get("lineages", {}):
            # round 2 (v9): the round-2 atlas table and the re-annotated lineages' tables and fits
            pairs.append(("(atlas, round 2)", spec["atlas"][sp]["res"], repo / round2["atlas_decisions"][sp],
                          repo / f"data/interim/latent_{sp}.h5ad", True))
            for lineage in round2["lineages"][sp]:
                cfg = spec["species"][sp]["lineages"][lineage]
                pairs.append((f"{lineage} (round 2)", cfg.get("res", spec.get("defaults", {}).get("res")),
                              repo / round2["decisions_dir"] / f"{sp}_{lineage}.csv",
                              repo / round2["outdir"] / f"latent_{sp}_{lineage}.h5ad", False))
        for lineage, resolution, csv_path, h5ad_path, atlas in pairs:
            table = pd.read_csv(csv_path, dtype=str)
            if not {"cluster", "label"} <= set(table):
                failures.append(f"{csv_path}: missing cluster/label")
                continue
            if table.cluster.duplicated().any():
                failures.append(f"{csv_path}: duplicate cluster keys")
            csv_clusters = set(table.cluster.str.strip())
            labels = set(table.label.str.strip()) - {"EXCLUDE"}
            if atlas:
                valid = set(spec["atlas"][sp]["labels"])
            else:
                valid = canonical
            invalid_labels = sorted(labels - valid)
            if invalid_labels:
                failures.append(f"{csv_path}: invalid labels {invalid_labels}")
            if h5ad_path.exists():
                import anndata as ad
                obs = ad.read_h5ad(h5ad_path, backed="r").obs
                if resolution not in obs:
                    failures.append(f"{h5ad_path}: missing {resolution}")
                    object_clusters = set()
                else:
                    object_clusters = set(obs[resolution].astype(str).unique())
                    if object_clusters != csv_clusters:
                        failures.append(
                            f"{csv_path}: cluster set differs from {h5ad_path.name}; "
                            f"missing={sorted(object_clusters-csv_clusters)} "
                            f"extra={sorted(csv_clusters-object_clusters)}")
            rows_out.append({
                "species": sp, "lineage": lineage, "resolution": resolution,
                "csv": str(csv_path), "n_clusters": len(csv_clusters),
                "n_labels": len(labels), "status": "PASS" if not invalid_labels else "FAIL",
            })
            print(f"  {lineage:14} {resolution}: {len(csv_clusters)} clusters; "
                  f"labels={sorted(labels)}")
    return failures, rows_out


def validate_h5ad(path: Path, species: str, role: str, vocab: dict):
    import anndata as ad
    obj = ad.read_h5ad(path, backed="r")
    obs = obj.obs
    failures = []
    gates = []

    def gate(check, expected, actual, *, info=False):
        status = "INFO" if info else ("PASS" if expected == actual else "SUSPECT")
        gates.append({"check": check, "expected": expected, "actual": actual,
                      "delta": None if info else actual - expected, "status": status})

    for col in ("cell_type", "cell_subtype"):
        if col not in obs:
            failures.append(f"missing required obs column {col}")
    if failures:
        return obj, failures, gates
    allowed_legacy = RELEASE_LEGACY_OBS if role in RELEASE_ROLES else set()
    present_legacy = sorted((LEGACY_OBS - allowed_legacy) & set(obs.columns))
    if present_legacy:
        failures.append(f"legacy annotation columns remain: {present_legacy}")

    types = obs.cell_type.astype(str)
    subtypes = obs.cell_subtype.astype(str)
    if obs.cell_type.isna().any() or obs.cell_subtype.isna().any():
        failures.append("canonical annotations contain NaN")
    for label, values in (("cell_type", types), ("cell_subtype", subtypes)):
        whitespace = sorted({x for x in set(values) if x != x.strip()})
        if whitespace:
            failures.append(f"{label} has surrounding whitespace: {whitespace}")
    missing_literals = sorted(set(subtypes.str.lower()) & {"nan", "none", "na", "<na>"})
    if missing_literals:
        failures.append(f"cell_subtype contains missing-value literals: {missing_literals}")

    allowed_species = set(vocab["expected_coarse_counts"][species])
    invalid_types = sorted(set(types) - allowed_species)
    if invalid_types:
        failures.append(f"invalid {species} cell_type values: {invalid_types}")
    canonical_types, subtype_parent = canonical_sets(vocab)
    empty_count = int((subtypes == "").sum())
    if empty_count:
        failures.append(
            f"cell_subtype contains {empty_count} empty values; expected coarse-fill")
    invalid_subtypes = sorted(set(subtypes) - canonical_types - set(subtype_parent))
    if invalid_subtypes:
        failures.append(f"invalid cell_subtype values: {invalid_subtypes}")
    invalid_pairs = sorted({
        (cell_type, subtype)
        for cell_type, subtype in zip(types, subtypes)
        if subtype != cell_type and subtype_parent.get(subtype) != cell_type
    })
    if invalid_pairs:
        failures.append(f"invalid subtype-parent pairs: {invalid_pairs}")

    coarse_expected = vocab["expected_coarse_counts"][species]
    actual_counts = types.value_counts()
    if role in ALL_CELLS_ROLES:
        for label, expected in coarse_expected.items():
            gate(f"cell_type:{label}", int(expected), int(actual_counts.get(label, 0)))
        gate("total_cells", int(vocab["expected_coarse_counts"]["totals"][species]), obj.n_obs)
    if role in {"fibdec", "fibdec_deposit"}:
        if set(types) != FIBDEC:
            failures.append(f"fib/dec object has cell_type values {sorted(set(types))}")
        for label in sorted(FIBDEC):
            gate(f"cell_type:{label}", int(coarse_expected[label]), int(actual_counts.get(label, 0)))
        gate("fibdec_total", sum(int(coarse_expected[label]) for label in FIBDEC), obj.n_obs)

    if role in ALL_CELLS_ROLES:
        gate("cell_subtype_empty", 0, empty_count)

        # Coarse fills and epithelial splits come from the vocabulary.
        for entry in vocab["content_tasks"]["coarse_subtype_fill"]["applies_to"][species]:
            label = entry["cell_type"]
            gate(f"{label}_coarse_subtype", int(entry["expected_coarse_fill_count"]),
                 int(((types == label) & (subtypes == label)).sum()))
        split_task = {"human": "human_epithelial_split", "mouse": "mouse_epithelial_split"}[species]
        for label, expected in vocab["content_tasks"][split_task]["expected_counts"].items():
            gate(label, int(expected), int((subtypes == label).sum()))

        # Endothelial and perivascular content tasks of the vocabulary: endothelial_cell and perivascular_cell are
        # cell types, the other labels subtypes. The human endothelial subtypes have no expected count (INFO).
        task = {"human": "human_perivascular_split", "mouse": "mouse_endothelial_blood_only"}[species]
        for label, expected in vocab["content_tasks"][task]["expected_counts"].items():
            values = types if label in canonical_types else subtypes
            gate(f"{task}:{label}", int(expected), int((values == label).sum()))
        if species == "human":
            for label in ("endothelial_blood", "endothelial_lymph"):
                gate(label, None, int((subtypes == label).sum()), info=True)

        # Round-2 lineages: the kept label counts of each approved table (EXCLUDE is not in the object).
        if "origin_lineage" in obs:
            lineage = obs["origin_lineage"].astype(str)
            round2 = vocab["content_tasks"].get("round2_lineages", {}).get("expected_counts", {})
            for key, counts in round2.items():
                key_species, key_lineage = key.split("_", 1)
                if key_species != species:
                    continue
                found = subtypes[lineage == key_lineage].value_counts()
                for label, expected in counts.items():
                    if label != "EXCLUDE":
                        gate(f"round2 {key}:{label}", int(expected), int(found.get(label, 0)))
                gate(f"round2 {key}:total", sum(int(n) for k, n in counts.items() if k != "EXCLUDE"),
                     int((lineage == key_lineage).sum()))

    if role != "reannotated":
        if obj.uns.get("annotation_schema_version") != vocab["schema_version"]:
            failures.append("annotation_schema_version missing or wrong")
        failures += display_name_failures(obj, vocab, role, gate)
        if "" in obj.uns.get("display_names", {}).get("cell_subtype", {}):
            failures.append("display_names contains forbidden empty-string key")

    return obj, failures, gates


def display_name_failures(obj, vocab, role, gate):
    """Compare uns['display_names'] with the vocabulary.

    All-cells objects carry the whole vocabulary map. A fib/dec subset may lack entries for labels it does not contain
    (the v9 fib/dec objects lack neutrophil, plasma_cell and rorc_positive_lymphocyte): every stored entry must match
    the vocabulary and every label in the object that the vocabulary names must have one; the missing entries are
    reported as INFO.
    """
    stored = obj.uns.get("display_names", {})
    if role not in {"fibdec", "fibdec_deposit"}:
        return [] if stored == vocab["display_names"] else ["uns['display_names'] differs from locked vocabulary"]
    failures = []
    for column, names in vocab["display_names"].items():
        have = dict(stored.get(column, {}))
        wrong = sorted(label for label in have if str(have[label]) != str(names.get(label)))
        present = set(obj.obs[column].astype(str)) if column in obj.obs else set()
        # Labels the vocabulary names in this column (a subtype column also holds cell-type fill values).
        uncovered = sorted(label for label in present if label in names and label not in have)
        if wrong or uncovered:
            failures.append(f"uns['display_names'][{column}] differs from locked vocabulary: "
                            f"wrong {wrong}, labels without a name {uncovered}")
        absent = sorted(set(names) - set(have) - present)
        if absent:
            gate(f"display_names {column} entries for labels absent from the object: {', '.join(absent)}",
                 None, len(absent), info=True)
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument("--species", nargs="+", choices=["human", "mouse"],
                    default=["human", "mouse"])
    ap.add_argument("--spec", default="resources/lineage_spec.yaml")
    ap.add_argument("--vocab", default="resources/annotation_vocabulary.yaml")
    ap.add_argument("--h5ad", default=None)
    ap.add_argument("--role", choices=["reannotated", "refined", "final", "fibdec",
                                       "all_deposit", "fibdec_deposit"], default="refined")
    ap.add_argument("--report", default=None)
    ap.add_argument("--nonfatal-counts", action="store_true")
    args = ap.parse_args()
    repo = Path(args.repo).resolve()
    vocab = yaml.safe_load((repo / args.vocab).read_text())
    spec = yaml.safe_load((repo / args.spec).read_text())

    if args.h5ad:
        if len(args.species) != 1:
            raise SystemExit("--h5ad requires exactly one --species")
        path = Path(args.h5ad)
        obj, failures, gates = validate_h5ad(path, args.species[0], args.role, vocab)
        count_suspects = [row for row in gates if row["status"] == "SUSPECT"]
        payload = {
            "path": str(path), "species": args.species[0], "role": args.role,
            "shape": [obj.n_obs, obj.n_vars], "schema_failures": failures,
            "count_gates": gates,
            "status": "FAIL" if failures else ("SUSPECT" if count_suspects else "PASS"),
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        if args.report:
            report = Path(args.report)
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        if failures:
            return 1
        if count_suspects and not args.nonfatal_counts:
            return 2
        return 0

    failures, rows = validate_csvs(repo, args.species, spec, vocab)
    payload = {"schema_version": vocab["schema_version"], "failures": failures,
               "records": rows, "status": "FAIL" if failures else "PASS"}
    if args.report:
        report = Path(args.report)
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if failures:
        print("FAILURES:")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print("ALL OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
