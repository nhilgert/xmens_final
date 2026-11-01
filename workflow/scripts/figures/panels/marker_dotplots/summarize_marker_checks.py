"""Every Table M1 v9 gene missing from an object's var and every gene that fails its rule on the release drawn.

Called by rule panels_marker_dotplots_checks (workflow/rules/91_panels_marker_dotplots.smk). Reads the release
rule checks and extract provenance of both species (panels_marker_dotplots_extract) and the extract provenance
of every review step (panels_marker_dotplots_review_extract). Writes one CSV (one row per finding; header only
when there is none) and a short Markdown summary.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

import marker_panels as panels


FINDING_COLUMNS = ["finding", "species", "object", "table_m1_row", "gene", "criterion", "detail"]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rule-checks", required=True, nargs="+", type=Path)
    parser.add_argument("--release-provenance", required=True, nargs="+", type=Path)
    parser.add_argument("--review-provenance", required=True, nargs="+", type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--out-md", required=True, type=Path)
    return parser.parse_args()


def rule_findings(path: Path, species: str, object_path: str) -> list:
    """Genes of the release that fail their Table M1 rule, or the extra check against smooth muscle."""
    checks = pd.read_csv(path)
    findings = []
    for row in checks.itertuples():
        if not row.in_var:
            continue
        detail = (f"{row.pct_detected_type:.1f}% of {row.n_type:,} {row.label} nuclei vs "
                  f"{row.pct_detected_comparison:.1f}% of {row.comparison_nuclei}")
        if not row.meets_criterion:
            findings.append(["fails rule", species, object_path, row.m1_cell_type, row.gene, row.criterion, detail])
        fails_sibling = (
            "meets_criterion_vs_sibling" in checks and pd.notna(row.meets_criterion_vs_sibling)
            and not row.meets_criterion_vs_sibling
        )
        if fails_sibling:
            findings.append(["fails rule vs " + row.sibling_label, species, object_path, row.m1_cell_type, row.gene,
                             row.criterion, f"{row.pct_detected_type:.1f}% vs {row.pct_detected_sibling:.1f}%"])
    return findings


def main() -> None:
    args = parse_arguments()
    findings, objects_checked, releases, table_m1_versions = [], [], set(), set()

    # Each species' rule check is the file named after it (<species>_table_m1_rule_check.csv).
    check_of = {path.name.split("_", 1)[0]: path for path in args.rule_checks}
    for provenance_path in sorted(args.release_provenance):
        provenance = json.loads(provenance_path.read_text())
        check_path = check_of[provenance["species"]]
        objects_checked.append(f"{provenance['species']} release {provenance['release']}: {provenance['object']}")
        releases.add(f"{provenance['version']} ({provenance['release']})")
        table_m1_versions.add(panels.table_m1_version(provenance["table_m1"]))
        findings += [["missing from var", provenance["species"], provenance["object"], "", gene, "", ""]
                     for gene in provenance["table_m1_genes_missing_from_var"]]
        findings += rule_findings(check_path, provenance["species"], provenance["object"])

    # Review steps: each clustering is checked for every Table M1 gene of its species.
    for provenance_path in sorted(args.review_provenance):
        provenance = json.loads(provenance_path.read_text())
        objects_checked.append(f"{provenance['species']} {provenance['step']}: {provenance['object']}")
        findings += [["missing from var", provenance["species"], provenance["object"], provenance["marker_plot"],
                      gene, "", ""]
                     for gene in provenance["genes_missing_from_var"]]

    table = pd.DataFrame(findings, columns=FINDING_COLUMNS)
    table.to_csv(args.out_csv, index=False)

    lines = [f"# Table M1 {', '.join(sorted(table_m1_versions))} genes on release {', '.join(sorted(releases))}", "",
             f"Objects checked ({len(objects_checked)}):", ""] + [f"- {o}" for o in objects_checked] + [""]
    for finding, rows in [("missing from var", table[table.finding == "missing from var"]),
                          ("failing their rule", table[table.finding != "missing from var"])]:
        lines.append(f"Genes {finding}: " + ("none" if rows.empty else str(len(rows))))
        lines += [f"- {r.species} {r.gene} ({r.table_m1_row or r.object}): {r.finding}; {r.detail}".rstrip("; ")
                  for r in rows.itertuples()]
        lines.append("")
    args.out_md.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
