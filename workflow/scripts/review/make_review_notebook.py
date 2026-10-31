"""Write the papermill templates of the review notebooks, one per annotation step.

    python workflow/scripts/review/make_review_notebook.py --config config/config.yaml \
        --spec resources/lineage_spec.yaml --out-dir notebooks/review

Steps: the whole-tissue atlas and every lineage of config `lineages`, per species. A lineage listed in config
`annotation_round2.lineages` is reviewed on its round-2 fit and round-2 table; the others on their July fit and
table. Each template's cell tagged `parameters` holds that step's defaults (paths relative to REPO); the rule
render_review_notebook (workflow/rules/99b_review_notebooks.smk) sets every parameter when it renders.
Writes {out_dir}/{species}_{step}_review.ipynb.
"""
import argparse
import json
from pathlib import Path

import nbformat
import yaml
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

V8_RELEASE_DIR = "data/deliverables/2026_10_29b"   # the v8 release, where the ledger's origin column looks
V8_RELEASE_OBJECT = {"human": "Cevrim_Human_All_Cells_extended.h5ad", "mouse": "Cevrim_XMens_All_Cells_extended.h5ad"}

INTRO = """# {Species} {step_title}: annotation decisions, v9

The decision table labels every cluster of `{cluster_col}` ({level}) in `{object_path}` ({fit});
a removed cluster is labelled `EXCLUDE`. This notebook shows the evidence for each decision with the numbers of
the review dot plot: marker detection per cluster, QC (genes, UMIs, mitochondrial reads, doublet scores, array of
origin) and the fixed per-nucleus checks of 29 Oct. It writes the step's review dot plot, cluster table, marker
statistics and its rows of the removed-cluster ledger to `OUT_DIR`.

Fixed flag rules, against the kept clusters of this step (workflow/scripts/review/annotation_review.py):

| flag | rule |
|---|---|
| mixed | >= 15% of nuclei express two lineage panels (>= 2 genes each) and >= 3x the kept median |
| other lineage | (lineage steps) a panel kept clusters rarely express (median < 15%) is in >= 40% of nuclei and >= 4x the kept median, or one of its genes is detected in >= 50% of nuclei against < 10% (and 5x less) in kept clusters |
| doublet score | mean DoubletDetection score >= 5 and >= 5x the kept clusters' median |
| doublets | >= 20% of nuclei called doublets (rare: called doublets were removed at QC) |
| low RNA | median detected genes < half the kept median |
| small | fewer than 20 nuclei |
| fits another label | another label's mean profile correlates better (r >= 0.85, ahead by >= 0.05) |

An EXCLUDE cluster is "supported" when it is mixed, of another lineage, doublet-like or low in RNA; otherwise the
kit calls it a reinstatement candidate or asks for a look. The recorded reason next to each removed cluster comes
from `REASONS_PATH`, its decision provenance from `STEPS_PATH`.
"""

PARAMETERS = """# Paths are relative to REPO. The rule render_review_notebook (workflow/rules/99b_review_notebooks.smk) sets
# every parameter; these defaults let the template run from notebooks/review/.
REPO = "../.."
SPECIES = {species}
STEP = {step}
OBJECT_PATH = {object_path}
CLUSTER_COLUMN = {cluster_col}
DECISIONS_PATH = {decisions}
VOCABULARY_PATH = "resources/annotation_vocabulary.yaml"
PANELS_PATH = "resources/review/panels_v9.json"
REASONS_PATH = "resources/review/exclusion_reasons_v9.csv"
STEPS_PATH = "resources/review/decision_steps_v9.csv"
V8_RELEASE_PATH = {v8_release}    # v8 extended all-cells object: where each nucleus was in v8
V8_DECISIONS_PATH = {v8_decisions}    # v8 table of the same clustering (None for a round-2 re-fit)
JULY_OBJECT_PATH = {july_object}    # July fit of a round-2 lineage (None otherwise)
OUT_DIR = "results/annotation_review\""""

SETUP = '''import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from IPython.display import Image, display

repo = Path(REPO).resolve()
sys.path.insert(0, str(repo / "workflow" / "scripts" / "review"))
from review_ledger import NOT_RELEASED, counts_text, july_origin, read_obs, removed_cluster_rows, v8_status  # noqa: E402
from review_step import ATLAS_ORDER, draw_decision_umap, label_order_and_display, read_decisions, review_step  # noqa: E402

PENDING_MARK = "PI sign-off pending"   # how a decision table's note flags a label awaiting the PI
PENDING_COLOUR = "#1F4E79"             # dark blue, apart from the red of EXCLUDE and the amber of the checks

# Verdicts, reasons and evidence are sentences; the tables below show them whole.
pd.set_option("display.max_colwidth", None, "display.width", 250, "display.max_rows", 200)


def in_repo(path):
    """Absolute path of a repository-relative path; None (input not used by this step) stays None."""
    return None if path in (None, "", "None") else repo / path'''

LOAD = '''adata = ad.read_h5ad(in_repo(OBJECT_PATH))
decisions = pd.read_csv(in_repo(DECISIONS_PATH), dtype=str, keep_default_na=False)
labels = read_decisions(decisions)
vocabulary = yaml.safe_load(in_repo(VOCABULARY_PATH).read_text())

# Every cluster has exactly one decision, and every label is one this step may use.
clusters = set(adata.obs[CLUSTER_COLUMN].astype(str))
repeated = decisions.cluster[decisions.cluster.duplicated()].tolist()
if clusters != set(labels) or repeated:
    raise ValueError(f"clusters without a decision: {sorted(clusters - set(labels))}; rows for absent clusters: "
                     f"{sorted(set(labels) - clusters)}; repeated rows: {repeated}")
if STEP == "atlas":
    allowed = set(ATLAS_ORDER[SPECIES])
else:
    allowed = {label for cell_type, subtypes in vocabulary["allowed"].items() for label in [cell_type, *(subtypes or [])]}
unknown = sorted(set(labels.values()) - allowed - {"EXCLUDE"})
if unknown:
    raise ValueError(f"labels not allowed at this step: {unknown}")

# One recorded reason per EXCLUDE cluster, and one provenance row for the step.
reasons = pd.read_csv(in_repo(REASONS_PATH), dtype=str, keep_default_na=False)
reasons = reasons[(reasons.species == SPECIES) & (reasons.step == STEP)].set_index("cluster")
excluded = sorted((c for c, label in labels.items() if label == "EXCLUDE"), key=int)
if set(reasons.index) != set(excluded):
    raise ValueError(f"EXCLUDE clusters without a reason: {sorted(set(excluded) - set(reasons.index))}; "
                     f"reasons for kept clusters: {sorted(set(reasons.index) - set(excluded))}")
steps = pd.read_csv(in_repo(STEPS_PATH), dtype=str, keep_default_na=False)
step_rows = steps[(steps.species == SPECIES) & (steps.step == STEP)]
if len(step_rows) != 1:
    raise ValueError(f"{len(step_rows)} rows for {SPECIES} {STEP} in {STEPS_PATH}")
step_record = step_rows.iloc[0]

notes = decisions["note"] if "note" in decisions else pd.Series("", index=decisions.index)
pending = decisions[notes.str.contains(PENDING_MARK, case=False, regex=False)]

print(f"{adata.n_obs:,} nuclei in {len(clusters)} clusters of {CLUSTER_COLUMN}; every cluster has one decision")
print(f"{len(excluded)} EXCLUDE clusters, each with a recorded reason")
print(f"Decided at {step_record.decided_at} ({step_record.approval_date}, {step_record.approver}): {step_record.decision_record}")
print("PI sign-off pending: " + (", ".join(f"cluster {r.cluster} ({r.label})" for r in pending.itertuples()) or "none"))'''

DOTPLOT = '''order, display_names = label_order_and_display(SPECIES, STEP, labels, vocabulary)
reason_text = {cluster: f"{row.reason} [{row.source_tag}]" for cluster, row in reasons.iterrows()}
pending_text = {row.cluster: f"{PENDING_MARK}: {display_names.get(row.label, row.label)}" for row in pending.itertuples()}

# Footnotes: how the decisions were made, what changed against the v8 table of the same clustering, pending items.
footnotes = [f"Decisions: {step_record.decided_at} ({step_record.approval_date}, {step_record.approver}); "
             f"{step_record.decision_record}."]
v8_decisions = read_decisions(in_repo(V8_DECISIONS_PATH)) if in_repo(V8_DECISIONS_PATH) else None
if v8_decisions is not None:
    changed = [f"{c} {v8_decisions.get(c, 'absent')} -> {labels[c]}" for c in sorted(labels, key=int)
               if v8_decisions.get(c) != labels[c]]
    footnotes.append("Changes against the v8 table of this clustering: " + ("; ".join(changed) or "none") + ".")
if pending_text:
    footnotes.append(f"{PENDING_MARK}: " + "; ".join(f"cluster {c} ({labels[c]})" for c in pending_text) + ".")

step_name = "whole-tissue atlas" if STEP == "atlas" else (f"{STEP}, round 2" if in_repo(JULY_OBJECT_PATH) else STEP)
result = review_step(
    SPECIES,
    STEP,
    OBJECT_PATH,
    CLUSTER_COLUMN,
    decisions,
    in_repo(OUT_DIR),
    in_repo(VOCABULARY_PATH),
    in_repo(PANELS_PATH),
    subtitle=f"Clusters: {OBJECT_PATH} ({CLUSTER_COLUMN}); decisions (v9): {DECISIONS_PATH}; reasons: {REASONS_PATH}",
    notes=footnotes,
    reasons=reason_text,
    reasons_header="exclusion reason [source]",
    extra_text=("PI sign-off", pending_text, PENDING_COLOUR) if pending_text else None,
    step_name=step_name,
    adata=adata,
)
plt.close(result["figure"])
table = result["table"]
if result["missing_genes"]:
    print("panel genes not in this object:", result["missing_genes"])
display(Image(result["paths"]["png"], width=1400))'''

REMOVED = '''cluster = adata.obs[CLUSTER_COLUMN].astype(str).to_numpy()
status = v8_status(adata.obs, read_obs(in_repo(V8_RELEASE_PATH)))
origin = None
if in_repo(JULY_OBJECT_PATH):
    origin = july_origin(adata.obs, read_obs(in_repo(JULY_OBJECT_PATH)), CLUSTER_COLUMN)

if origin is not None:
    fit = "round 2 (lineage re-fit)"
else:
    fit = "July whole-tissue fit" if STEP == "atlas" else "July lineage fit"
step_info = dict(species=SPECIES, step=STEP, fit=fit, object=OBJECT_PATH, cluster_column=CLUSTER_COLUMN,
                 decisions_table=DECISIONS_PATH,
                 **step_record[["decided_at", "approval_date", "approver", "decision_record"]].to_dict())
ledger = removed_cluster_rows(table, result["qc"], cluster, status, step_info, reasons, origin, v8_decisions)
ledger_path = in_repo(OUT_DIR) / f"{SPECIES}_{STEP}_removed_clusters.csv"
ledger.to_csv(ledger_path, index=False)

kept_reference = table.iloc[0]
print(f"{len(ledger)} clusters removed ({ledger.n_nuclei.sum():,} of {adata.n_obs:,} nuclei). Kept-cluster references: "
      f"median genes {kept_reference.ref_kept_median_genes:,.0f}; two-lineage median {kept_reference.ref_kept_median_two_lineages:.1f}%; "
      f"doublet-score median {kept_reference.ref_kept_median_doublet:.2f}")
ledger[["cluster", "n_nuclei", "median_genes", "median_umis", "mean_pct_mito", "mean_doublet_score",
        "median_doublet_score", "pct_doublet_score_above_5", "pct_two_lineages", "largest_array", "largest_array_share",
        "kit_verdict", "final_reason", "reason_source", "origin_v8"]].round(2)'''

EVIDENCE = '''ledger[["cluster", "kit_best_label", "kit_best_r", "evidence_summary", "evidence_source", "reason_v8"]].round(2)'''

FLAGS = '''flagged = table[~table.excluded & (table.verdict != "fine")]
print(f"{len(flagged)} of {int((~table.excluded).sum())} kept clusters flagged by the fixed rules")
flagged[["cluster", "label", "n_nuclei", "flags", "median_genes", "mean_doublet", "pct_two_lineages",
         "best_label", "best_r", "own_label_r", "neighbour_labels", "top_genes"]].round(2)'''

V8_ORIGIN = '''rows = []
for name in sorted(set(cluster), key=int):
    member = cluster == name
    row = {"cluster": name, "decision": labels[name], "n_nuclei": int(member.sum()),
           "released_in_v8": int((status[member] != NOT_RELEASED).sum()), "v8_status": counts_text(status[member])}
    if origin is not None:
        row["origin"] = counts_text(origin[member])
    rows.append(row)
v8_by_cluster = pd.DataFrame(rows)

is_removed = v8_by_cluster.decision.eq("EXCLUDE")
for group, members in [("kept", v8_by_cluster[~is_removed]), ("removed", v8_by_cluster[is_removed])]:
    if len(members):
        print(f"{group} clusters: {members.released_in_v8.sum():,} of {members.n_nuclei.sum():,} nuclei released in v8")
v8_by_cluster'''

UMAP = '''draw_decision_umap(adata, CLUSTER_COLUMN, labels, order, display_names)
plt.show()'''

COUNTS = '''counts = table.groupby("label").agg(clusters=("cluster", "size"), nuclei=("n_nuclei", "sum"))
counts = counts.sort_values("nuclei", ascending=False)
counts["display name"] = [display_names.get(label, label) for label in counts.index]
counts.loc["total"] = [counts.clusters.sum(), counts.nuclei.sum(), ""]
counts'''

FILES = '''for name, path in {**result["paths"], "removed_clusters": ledger_path}.items():
    print(f"{name}: {Path(path).relative_to(repo)}")'''


def python_literal(value):
    """Source text of a parameter value: a double-quoted string, or None."""
    return "None" if value is None else json.dumps(value)


def review_steps(config, spec):
    """(species, step, defaults) for the atlas and every lineage, from config and lineage_spec."""
    round2 = config["annotation_round2"]
    steps = []
    for species, lineages in config["lineages"].items():
        v8_release = f"{V8_RELEASE_DIR}/{V8_RELEASE_OBJECT[species]}"
        steps.append((species, "atlas", dict(
            object_path=f"data/interim/latent_{species}.h5ad", cluster_col=spec["atlas"][species]["res"],
            decisions=round2["atlas_decisions"][species], v8_release=v8_release,
            v8_decisions=f"resources/annotations/{species}_lineage.csv", july_object=None, fit="July whole-tissue fit")))
        for lineage in lineages:
            cluster_col = spec["species"][species]["lineages"][lineage].get("res", spec["defaults"]["res"])
            if lineage in round2["lineages"].get(species, []):
                defaults = dict(object_path=f"{round2['outdir']}/latent_{species}_{lineage}.h5ad",
                                decisions=f"resources/annotations/round2/{species}_{lineage}.csv", v8_decisions=None,
                                july_object=f"data/interim/latent_{species}_{lineage}.h5ad", fit="round-2 lineage re-fit")
            else:
                table = f"resources/annotations/{species}_{lineage}.csv"
                defaults = dict(object_path=f"data/interim/latent_{species}_{lineage}.h5ad", decisions=table,
                                v8_decisions=table, july_object=None, fit="July lineage fit")
            steps.append((species, lineage, dict(cluster_col=cluster_col, v8_release=v8_release, **defaults)))
    return steps


def build(species, step, defaults):
    level = "whole-tissue Leiden 2.5" if step == "atlas" else "Leiden 3.0"
    step_title = "whole-tissue atlas" if step == "atlas" else step + (", round 2" if defaults["july_object"] else "")
    values = {"species": species, "step": step, **{k: v for k, v in defaults.items() if k != "fit"}}
    parameters = new_code_cell(PARAMETERS.format(**{k: python_literal(v) for k, v in values.items()}))
    parameters.metadata["tags"] = ["parameters"]
    cells = [
        new_markdown_cell(INTRO.format(Species=species.capitalize(), step_title=step_title, level=level, **defaults)),
        parameters,
        new_code_cell(SETUP),
        new_markdown_cell("## Decisions and inputs\nChecks that every cluster has one allowed label and every removed cluster a recorded reason."),
        new_code_cell(LOAD),
        new_markdown_cell("## Review dot plot\nRows: clusters grouped by label, removed clusters in the red EXCLUDE block. "
                          "Amber: the fixed checks of each flagged or removed cluster. Red: the recorded exclusion reason and its source."),
        new_code_cell(DOTPLOT),
        new_markdown_cell("## Removed clusters\nThis step's rows of the removed-cluster ledger: QC, kit verdict, final reason and "
                          "where the nuclei were in v8."),
        new_code_cell(REMOVED),
        new_markdown_cell("Recorded evidence behind each removal, and the v8 wording where a reason was reworded."),
        new_code_cell(EVIDENCE),
        new_markdown_cell("## Flags on kept clusters"),
        new_code_cell(FLAGS),
        new_markdown_cell("## Where each cluster's nuclei were in v8\nv8 status: the nucleus's label in the v8 extended release "
                          "(`(extended only)`: supplementary population), or `not released`."
                          + (" Origin: the July cluster of the nucleus, or the atlas cluster it was reinstated from." if defaults["july_object"] else "")),
        new_code_cell(V8_ORIGIN),
        new_markdown_cell("## UMAP\nThe step's own embedding (`X_umap_scVI`)."),
        new_code_cell(UMAP),
        new_markdown_cell("## Label counts"),
        new_code_cell(COUNTS),
        new_markdown_cell("## Files written"),
        new_code_cell(FILES),
    ]
    metadata = {"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                "language_info": {"name": "python"}}
    return new_notebook(cells=cells, metadata=metadata)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    parser.add_argument("--spec", type=Path, default=Path("resources/lineage_spec.yaml"))
    parser.add_argument("--out-dir", type=Path, default=Path("notebooks/review"))
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    spec = yaml.safe_load(args.spec.read_text())

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for species, step, defaults in review_steps(config, spec):
        path = args.out_dir / f"{species}_{step}_review.ipynb"
        nbformat.write(build(species, step, defaults), path)
        print(path)


if __name__ == "__main__":
    main()
