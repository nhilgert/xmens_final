"""Supplementary Table S5 - per-sample QC cutoffs and counts.

Source of record: the per-sample QC logs written by scripts/ingest_qc/qc.py (logs/01_qc_{sample}.log). qc.py prints
every cutoff and every per-stage count but writes only the filtered object, so the logs are the only record of the
numbers the Methods quote. This script parses them; it recomputes nothing, so it cannot drift from what was applied.

The species totals are checked against config["supp_tables"]["s5"]["methods_totals"], i.e. the numbers written in the
Methods. A parse that does not reproduce them stops the rule instead of shipping a wrong table.

One line of the log is deliberately NOT used: "3-MAD <metric> (log10): ... flag N" counts cells the bound flags,
which double-counts cells already removed by an earlier filter. The per-filter "removed N" lines are what the
Methods sentence adds up, and they are what this script reads.
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

PATTERNS = {
    "cellbender_called": re.compile(r"loaded (\d+) cells x \d+ genes \(species=(\w+)\)"),
    "umi_floor": re.compile(r"min (\d+) counts: removed (\d+), (\d+)/(\d+) remain"),
    "mad_total": re.compile(r"MAD total_counts \(log10\): bounds=\[([\d.]+), ([\d.]+)\]"),
    "mad_genes": re.compile(r"MAD n_genes_by_counts \(log10\): bounds=\[([\d.]+), ([\d.]+)\]"),
    "rm_mito": re.compile(r"filter pct_counts_mt: removed (\d+), (\d+) remain"),
    "rm_total": re.compile(r"filter total_counts: removed (\d+), (\d+) remain"),
    "rm_genes": re.compile(r"filter n_genes_by_counts: removed (\d+), (\d+) remain"),
    "after_thresh": re.compile(r"after threshold QC: (\d+)/(\d+) cells"),
    "doublets": re.compile(r"flagged (\d+) doublets"),
    "final": re.compile(r"final: (\d+) cells x"),
}

COLUMNS = [
    "sample", "species", "cellbender_called", "removed_umi_floor_lt500", "after_floor",
    "mad_total_counts_low", "mad_total_counts_high", "mad_n_genes_low", "mad_n_genes_high",
    "removed_mito_gt8pct", "removed_total_counts_mad", "removed_n_genes_mad",
    "after_threshold_filters", "removed_doublets", "retained",
]

TOTAL_COLUMNS = {
    "cellbender_called": "cellbender_called",
    "umi_floor_removed": "removed_umi_floor_lt500",
    "mito_removed": "removed_mito_gt8pct",
    "total_counts_mad_removed": "removed_total_counts_mad",
    "n_genes_mad_removed": "removed_n_genes_mad",
    "doublets_removed": "removed_doublets",
    "retained": "retained",
}


def need(pattern: re.Pattern, text: str, sample: str, what: str) -> re.Match:
    """First match of pattern, or a hard stop naming the sample and the missing line."""
    match = pattern.search(text)
    if match is None:
        raise SystemExit(f"[s5] {sample}: no '{what}' line in the QC log; the log format changed")
    return match


def parse_log(path: Path) -> dict:
    """One row of Table S5 from one per-sample QC log."""
    sample = path.stem.replace("01_qc_", "")
    text = path.read_text()
    called, species = need(PATTERNS["cellbender_called"], text, sample, "loaded").groups()
    floor = need(PATTERNS["umi_floor"], text, sample, "UMI floor").groups()
    tlo, thi = need(PATTERNS["mad_total"], text, sample, "total_counts MAD bounds").groups()
    glo, ghi = need(PATTERNS["mad_genes"], text, sample, "n_genes MAD bounds").groups()
    row = {
        "sample": sample,
        "species": species,
        "cellbender_called": int(called),
        "removed_umi_floor_lt500": int(floor[1]),
        "after_floor": int(floor[2]),
        "mad_total_counts_low": float(tlo),
        "mad_total_counts_high": float(thi),
        "mad_n_genes_low": float(glo),
        "mad_n_genes_high": float(ghi),
        "removed_mito_gt8pct": int(need(PATTERNS["rm_mito"], text, sample, "mito filter").group(1)),
        "removed_total_counts_mad": int(need(PATTERNS["rm_total"], text, sample, "count MAD filter").group(1)),
        "removed_n_genes_mad": int(need(PATTERNS["rm_genes"], text, sample, "gene MAD filter").group(1)),
        "after_threshold_filters": int(need(PATTERNS["after_thresh"], text, sample, "after threshold").group(1)),
        "removed_doublets": int(need(PATTERNS["doublets"], text, sample, "doublets").group(1)),
        "retained": int(need(PATTERNS["final"], text, sample, "final").group(1)),
    }
    if int(floor[3]) != row["cellbender_called"]:
        raise SystemExit(f"[s5] {sample}: floor line disagrees with the loaded count")
    stages = (row["after_floor"] - row["removed_mito_gt8pct"] - row["removed_total_counts_mad"]
              - row["removed_n_genes_mad"])
    if stages != row["after_threshold_filters"]:
        raise SystemExit(f"[s5] {sample}: filters remove {row['after_floor'] - stages} but "
                         f"after-threshold is {row['after_threshold_filters']}")
    if row["after_threshold_filters"] - row["removed_doublets"] != row["retained"]:
        raise SystemExit(f"[s5] {sample}: after-threshold minus doublets is not the retained count")
    return row


def check_totals(frame: pd.DataFrame, expected: dict) -> dict:
    """Species sums against the Methods numbers; every disagreement is reported before stopping."""
    sums, wrong = {}, []
    for name, column in TOTAL_COLUMNS.items():
        got = {sp: int(v) for sp, v in frame.groupby("species")[column].sum().items()}
        sums[name] = got
        want = expected.get(name)
        if want is not None and got != {k: int(v) for k, v in want.items()}:
            wrong.append(f"  {name}: table {got}, Methods {want}")
    if wrong:
        raise SystemExit("[s5] the table does not reproduce the Methods numbers:\n" + "\n".join(wrong))
    return sums


def main() -> None:
    logs = [Path(p) for p in snakemake.input.logs]            # noqa: F821  (Snakemake injects this)
    out_csv = Path(snakemake.output.csv)                      # noqa: F821
    out_json = Path(snakemake.output.provenance)              # noqa: F821
    expected = dict(snakemake.params.methods_totals)          # noqa: F821

    frame = pd.DataFrame([parse_log(p) for p in logs], columns=COLUMNS)
    if frame["sample"].duplicated().any():
        raise SystemExit("[s5] two logs parsed to the same sample name")
    sums = check_totals(frame, expected)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_csv, index=False)
    provenance = {
        "table": "S5",
        "description": "Per-sample QC cutoffs and counts, parsed from the QC logs of record.",
        "n_samples": int(len(frame)),
        "species_totals": sums,
        "source_logs": [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                         "mtime": p.stat().st_mtime} for p in logs],
    }
    out_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"[s5] {len(frame)} samples -> {out_csv}", file=sys.stderr)


main()
