#!/usr/bin/env python3
"""Mechanically migrate frozen annotation CSV labels to xmens.annotation.v1.

Raw collaborator exports (``*_fine.csv``) are immutable inputs and are never
edited.  Atlas decision CSVs are converted to internal lineage-routing slugs;
04c decision CSVs are converted to canonical cell types/subtypes.  The two
source-defined content tasks are explicit and count-gated downstream.

The command is idempotent and defaults to a dry run.  It refuses every label
that is not covered by the locked vocabulary or an explicit content rule.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import tempfile

import yaml


RAW_EXPORTS = {"human_immune_fine.csv", "mouse_immune_fine.csv"}
IGNORED = {"_BROKEN_fib_metadata_15972of43052_rows.csv"}
ATLAS_ROUTE_IMPORT = {
    "human": {
        "stromal": "stromal", "epithelial cells": "epithelial",
        "immune cells": "immune", "perivascular cells": "perivascular",
        "smooth muscle cells": "smc",
    },
    "mouse": {
        "stromal": "stromal", "epithelial cells": "epithelial",
        "immune cells": "immune", "vascular cells": "vascular",
        "smooth muscle cells": "smc", "perimetrium cells": "perimetrium",
    },
}


def atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    original_mode = path.stat().st_mode & 0o777
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.chmod(tmp_name, original_mode)
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def canonical_values(vocab: dict) -> set[str]:
    out = set(vocab["allowed"])
    for values in vocab["allowed"].values():
        out.update(values or [])
    return out


def mapping_for_file(path: Path, species: str, vocab: dict, atlas_routes: dict[str, str]):
    canonical = canonical_values(vocab)
    type_map = vocab["import_mapping"][species]["cell_type"]
    subtype_map = vocab["import_mapping"][species]["cell_subtype"]

    def migrate(label: str) -> str:
        label = label.strip()
        if label == "EXCLUDE" or label in canonical:
            return label
        if path.name == f"{species}_lineage.csv":
            if label in ATLAS_ROUTE_IMPORT[species]:
                return ATLAS_ROUTE_IMPORT[species][label]
            if label in set(atlas_routes.values()):
                return label
            raise KeyError(f"atlas routing label {label!r} is not in lineage_spec atlas.labels")

        # Content tasks are source-defined, never inferred from expression.
        if path.name == "human_stromal.csv" and label == "perivascular cells":
            return "perivascular_cell"
        if path.name == "mouse_vascular.csv" and label == "vascular cells":
            return "endothelial_blood"

        if label in subtype_map:
            return subtype_map[label]
        if label in type_map:
            return type_map[label]
        raise KeyError(f"unmapped stored label {label!r}")

    return migrate


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument("--vocab", default="resources/annotation_vocabulary.yaml")
    ap.add_argument("--spec", default="resources/lineage_spec.yaml")
    ap.add_argument("--apply", action="store_true", help="atomically replace covered decision CSVs")
    ap.add_argument("--report", default=None, help="optional JSON coverage/change report")
    args = ap.parse_args()

    repo = Path(args.repo).resolve()
    vocab = yaml.safe_load((repo / args.vocab).read_text())
    spec = yaml.safe_load((repo / args.spec).read_text())
    ann_dir = repo / "resources/annotations"
    report: dict[str, object] = {
        "schema_version": vocab["schema_version"],
        "mode": "apply" if args.apply else "dry-run",
        "raw_exports_untouched": sorted(RAW_EXPORTS),
        "files": {},
        "unmapped": [],
    }
    pending_writes: list[tuple[Path, list[str], list[dict[str, str]]]] = []

    for path in sorted(ann_dir.glob("*.csv")):
        if path.name in RAW_EXPORTS or path.name in IGNORED:
            continue
        species = path.name.split("_", 1)[0]
        if species not in {"human", "mouse"}:
            continue
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])
        if "label" not in fieldnames:
            continue
        migrate = mapping_for_file(path, species, vocab, spec["atlas"][species]["labels"])
        changes: dict[str, int] = {}
        bad: dict[str, str] = {}
        for row in rows:
            old = (row.get("label") or "").strip()
            try:
                new = migrate(old)
            except KeyError as exc:
                bad[old] = str(exc)
                continue
            row["label"] = new
            if new != old:
                changes[f"{old} -> {new}"] = changes.get(f"{old} -> {new}", 0) + 1
        report["files"][str(path.relative_to(repo))] = {
            "rows": len(rows), "changes": changes, "unmapped": bad,
        }
        if bad:
            report["unmapped"].append({"file": str(path.relative_to(repo)), "labels": bad})
        elif args.apply and changes:
            pending_writes.append((path, fieldnames, rows))

    if args.report:
        out = repo / args.report
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["unmapped"]:
        print("REFUSING: at least one stored label lacks an explicit mapping")
        return 2
    # No file is replaced until coverage for every decision record has passed.
    for path, fieldnames, rows in pending_writes:
        atomic_write_csv(path, fieldnames, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
