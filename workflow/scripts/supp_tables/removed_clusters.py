"""Supplementary Table S13 - every cluster removed during annotation.

Source of record: the removed-cluster ledger written by rule removed_cluster_ledger
(results/annotation_review/removed_cluster_ledger_{version}.csv), which already carries one row per removed cluster
with its quality metrics, the comparison with the clusters kept in the same step, the reason and the approval.

This script selects and renames columns for the supplement; it changes no value. The column list and the header
names live in config["supp_tables"]["s13"], so the selection can be changed without touching code.

The row count and the nucleus total are checked against config["supp_tables"]["s13"]["expected"], i.e. the numbers
the Methods quote ("38 clusters (12,274 nuclei)").
"""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    src = Path(snakemake.input.ledger)                        # noqa: F821  (Snakemake injects this)
    out_csv = Path(snakemake.output.csv)                      # noqa: F821
    out_json = Path(snakemake.output.provenance)              # noqa: F821
    spec = dict(snakemake.params.spec)                        # noqa: F821

    ledger = pd.read_csv(src)
    columns = list(spec["columns"])
    missing = [c for c in columns if c not in ledger.columns]
    if missing:
        raise SystemExit(f"[s13] the ledger has no column(s) {missing}; "
                         f"available: {sorted(ledger.columns)}")

    table = ledger[columns].copy()
    rename = dict(spec.get("header_names", {}))
    if rename:
        table = table.rename(columns=rename)
    sort_by = [c for c in spec.get("sort_by", []) if c in ledger.columns]
    if sort_by:
        table = table.assign(**{c: ledger[c] for c in sort_by if c not in table.columns})
        table = table.sort_values(sort_by, kind="stable").drop(
            columns=[c for c in sort_by if c not in [rename.get(x, x) for x in columns]])

    expected = dict(spec.get("expected", {}))
    n_rows, want_rows = int(len(table)), expected.get("n_clusters")
    if want_rows is not None and n_rows != int(want_rows):
        raise SystemExit(f"[s13] the ledger has {n_rows} removed clusters, the Methods say {want_rows}")
    count_column = spec.get("nuclei_column", "n_nuclei")
    n_nuclei, want_nuclei = int(ledger[count_column].sum()), expected.get("n_nuclei")
    if want_nuclei is not None and n_nuclei != int(want_nuclei):
        raise SystemExit(f"[s13] the removed clusters hold {n_nuclei:,} nuclei, "
                         f"the Methods say {int(want_nuclei):,}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_csv, index=False)
    provenance = {
        "table": "S13",
        "description": "Clusters removed during annotation, selected from the removed-cluster ledger.",
        "n_clusters": n_rows,
        "n_nuclei": n_nuclei,
        "by_species": {k: int(v) for k, v in ledger.groupby("species")[count_column].sum().items()},
        "columns": list(table.columns),
        "source": {"path": str(src), "sha256": hashlib.sha256(src.read_bytes()).hexdigest()},
    }
    out_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"[s13] {n_rows} clusters, {n_nuclei:,} nuclei -> {out_csv}", file=sys.stderr)


main()
