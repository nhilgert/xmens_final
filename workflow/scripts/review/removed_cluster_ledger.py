"""Assemble the removed-cluster ledger from the review notebooks' per-step rows, with a short summary.

Rule removed_cluster_ledger (workflow/rules/99b_review_notebooks.smk). Reads the {species}_{step}_removed_clusters.csv
files the review notebooks write (one row per EXCLUDE cluster) and writes the ledger CSV, in step order, and a
Markdown summary: removals per step, kit verdicts that did not support a removal, reworded reasons and how many
removed nuclei were in the v8 release. Every number in the summary is computed here.
"""
import logging

import pandas as pd

removed_paths = list(snakemake.input.removed)
csv_path = snakemake.output.csv
md_path = snakemake.output.md
version = snakemake.params.version
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format="%(message)s")


def markdown_table(frame):
    """A pandas frame as a Markdown table (no tabulate dependency)."""
    header = "| " + " | ".join(frame.columns) + " |"
    rule = "|" + "|".join("---" for _ in frame.columns) + "|"
    body = ["| " + " | ".join(str(value) for value in row) + " |" for row in frame.itertuples(index=False)]
    return "\n".join([header, rule, *body])


# Per-step rows in the rules' step order; species sorted, steps within a species as listed.
steps = []
for path in removed_paths:
    rows = pd.read_csv(path, dtype={"cluster": str, "approval_date": str}, keep_default_na=False)
    steps.append(rows)
    logging.info(f"{path}: {len(rows)} removed clusters")
ledger = pd.concat(steps, ignore_index=True)
step_rank = {(sp, st): i for i, (sp, st) in enumerate(dict.fromkeys(zip(ledger.species, ledger.step)))}
ledger = (
    ledger.assign(_species=ledger.species, _step=[step_rank[k] for k in zip(ledger.species, ledger.step)],
                  _cluster=ledger.cluster.astype(int))
    .sort_values(["_species", "_step", "_cluster"])
    .drop(columns=["_species", "_step", "_cluster"])
)
ledger.to_csv(csv_path, index=False)
logging.info(f"{csv_path}: {len(ledger)} clusters, {ledger.n_nuclei.sum():,} nuclei")

# Removals per step, including steps without any (their files are header-only).
per_step = []
for path, rows in zip(removed_paths, steps):
    species, step = path.rsplit("/", 1)[1].replace("_removed_clusters.csv", "").split("_", 1)
    per_step.append({"species": species, "step": step, "clusters removed": len(rows),
                     "nuclei removed": f"{rows.n_nuclei.sum():,}",
                     "decided at": rows.decided_at.iloc[0] if len(rows) else "",
                     "approval date": rows.approval_date.iloc[0] if len(rows) else ""})
per_step = pd.DataFrame(per_step).sort_values("species", kind="stable")

# Kit verdict classes: "exclusion supported" or not (reinstatement candidate / needs a look).
supported = ledger.kit_verdict.str.startswith("exclusion supported")
beyond_rules = ledger.loc[~supported, ["species", "step", "cluster", "n_nuclei", "kit_verdict", "final_reason"]]
beyond_rules = beyond_rules.assign(n_nuclei=beyond_rules.n_nuclei.map("{:,}".format))
reworded = ledger.loc[ledger.reason_v8 != "", ["species", "step", "cluster", "final_reason", "reason_v8"]]
released = ledger.n_released_in_v8.astype(int)
in_primary = ledger.n_in_v8_primary.astype(int)

lines = [
    f"# Removed-cluster ledger {version}",
    "",
    f"{len(ledger)} clusters removed ({ledger.n_nuclei.sum():,} nuclei) across {len(per_step)} annotation steps "
    f"({int((per_step['clusters removed'] > 0).sum())} with removals). One row per EXCLUDE cluster in "
    f"removed_cluster_ledger_{version}.csv: QC, doublet scores, kit verdict, final reason with its source, recorded "
    "evidence, v8 origin and decision provenance.",
    "",
    "## Removals per step",
    "",
    markdown_table(per_step),
    "",
    "## Kit verdicts",
    "",
    f"The fixed rules support {int(supported.sum())} of the {len(ledger)} removals. The other {len(beyond_rules)} "
    "are decisions beyond the fixed rules; the recorded final reason of each is shown:",
    "",
    markdown_table(beyond_rules) if len(beyond_rules) else "None.",
    "",
    "## Reworded reasons",
    "",
    markdown_table(reworded) if len(reworded) else "None.",
    "",
    "## v8 origin",
    "",
    f"{int(released.sum()):,} of the {ledger.n_nuclei.sum():,} removed nuclei were in the v8 release "
    f"({int((released > 0).sum())} clusters): {int(in_primary.sum()):,} in the primary objects and "
    f"{int((released - in_primary).sum()):,} only in the extended objects (supplementary populations). The rest were "
    "already removed in v8. Per cluster: columns n_released_in_v8, n_in_v8_primary and origin_v8.",
    "",
]
with open(md_path, "w") as handle:
    handle.write("\n".join(lines))
logging.info(f"{md_path} written")
