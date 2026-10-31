import os

import yaml as _review_yaml

# 99b_review_notebooks -- review notebooks of the annotation decisions; OFF the critical path on purpose.
#
# One papermill notebook per annotation step: the whole-tissue atlas and every lineage of each species. Each
# shows the review dot plot (final labels, removed clusters in a red EXCLUDE block with their recorded reason),
# the kit's fixed per-nucleus checks and where every nucleus was in v8, and writes the step's rows of the
# removed-cluster ledger. Nothing in rule all depends on these rules, like 99_reports:
#
#     snakemake --profile profiles/slurm review_notebooks
#
# Templates: notebooks/review/{species}_{step}_review.ipynb (workflow/scripts/review/make_review_notebook.py).
# Code: workflow/scripts/review/. Recorded reasons and decision provenance: resources/review/.
# Needs papermill and pdfunite (poppler-utils) on PATH.

REVIEW_VERSION = "v9"                         # the decision set reviewed (vocabulary and tables v9)
REVIEW_DIR = "results/annotation_review"
REVIEW_NOTEBOOK_DIR = "results/notebooks/review"
REVIEW_KIT = ["annotation_review", "review_dotplot", "review_step", "review_ledger"]

# The v8 release: the ledger's origin column says where each removed nucleus was in v8. Pinned on purpose
# (config deliverables.release moves on with every release cut) and passed as a parameter, not an input: the
# release is immutable, and as an input it would tie this target to the 04g release rules, which re-run
# whenever the vocabulary is newer than the release.
V8_RELEASE_DIR = "data/deliverables/2026_10_29b"
V8_RELEASE_OBJECT = {"human": "Cevrim_Human_All_Cells_extended.h5ad", "mouse": "Cevrim_XMens_All_Cells_extended.h5ad"}

# Page order of the all-steps dot plots: atlas first, then the lineages as in the v8 plots.
REVIEW_LINEAGE_ORDER = ["stromal", "epithelial", "vascular", "perivascular", "perimetrium", "smc", "immune"]

_REVIEW_SPEC = _review_yaml.safe_load(open("resources/lineage_spec.yaml"))
_REVIEW_ROUND2 = config["annotation_round2"]


def review_steps(species):
    """Annotation steps of one species: the atlas, then its lineages in page order."""
    return ["atlas"] + sorted(config["lineages"][species], key=REVIEW_LINEAGE_ORDER.index)


def _is_round2(wildcards):
    return wildcards.step in _REVIEW_ROUND2["lineages"].get(wildcards.species, [])


def review_object(wildcards):
    """The clustering the step's decisions label: whole-tissue fit, round-2 re-fit or July lineage fit."""
    if wildcards.step == "atlas":
        return f"data/interim/latent_{wildcards.species}.h5ad"
    if _is_round2(wildcards):
        return f"{_REVIEW_ROUND2['outdir']}/latent_{wildcards.species}_{wildcards.step}.h5ad"
    return f"data/interim/latent_{wildcards.species}_{wildcards.step}.h5ad"


def review_decisions(wildcards):
    """The step's decision table: round-2 atlas table, round-2 lineage table or the July lineage table."""
    if wildcards.step == "atlas":
        return _REVIEW_ROUND2["atlas_decisions"][wildcards.species]
    if _is_round2(wildcards):
        return f"resources/annotations/round2/{wildcards.species}_{wildcards.step}.csv"
    return f"resources/annotations/{wildcards.species}_{wildcards.step}.csv"


def review_v8_decisions(wildcards):
    """v8 table of the same clustering (the July atlas table; a July-fit lineage's own table); none for round 2."""
    if wildcards.step == "atlas":
        return f"resources/annotations/{wildcards.species}_lineage.csv"
    if _is_round2(wildcards):
        return []
    return review_decisions(wildcards)


def review_july_object(wildcards):
    """July fit of a round-2 lineage: where its nuclei were before the re-fit."""
    return f"data/interim/latent_{wildcards.species}_{wildcards.step}.h5ad" if _is_round2(wildcards) else []


def review_cluster_column(wildcards):
    """The step's annotation resolution (lineage_spec res)."""
    if wildcards.step == "atlas":
        return _REVIEW_SPEC["atlas"][wildcards.species]["res"]
    lineage = _REVIEW_SPEC["species"][wildcards.species]["lineages"][wildcards.step]
    return lineage.get("res", _REVIEW_SPEC["defaults"]["res"])


def review_outputs(pattern):
    """pattern expanded over every species and step."""
    return [pattern.format(species=sp, step=step) for sp in config["lineages"] for step in review_steps(sp)]


# Renders one step's review notebook with papermill; the notebook writes the step's dot plot, cluster table,
# marker statistics and removed-cluster rows. Paths reach the notebook relative to the working directory.
rule render_review_notebook:
    input:
        notebook="notebooks/review/{species}_{step}_review.ipynb",
        adata=review_object,
        decisions=review_decisions,
        v8_decisions=review_v8_decisions,
        july_object=review_july_object,
        vocabulary="resources/annotation_vocabulary.yaml",
        panels="resources/review/panels_v9.json",
        reasons="resources/review/exclusion_reasons_v9.csv",
        steps="resources/review/decision_steps_v9.csv",
        kit=[f"workflow/scripts/review/{module}.py" for module in REVIEW_KIT],
    output:
        notebook=f"{REVIEW_NOTEBOOK_DIR}/{{species}}_{{step}}_review.ipynb",
        dotplot_pdf=f"{REVIEW_DIR}/{{species}}_{{step}}_dotplot.pdf",
        dotplot_png=f"{REVIEW_DIR}/{{species}}_{{step}}_dotplot.png",
        dotplot_svg=f"{REVIEW_DIR}/{{species}}_{{step}}_dotplot.svg",
        clusters=f"{REVIEW_DIR}/{{species}}_{{step}}_clusters.csv",
        marker_stats=f"{REVIEW_DIR}/{{species}}_{{step}}_marker_stats.csv",
        removed=f"{REVIEW_DIR}/{{species}}_{{step}}_removed_clusters.csv",
    params:
        repo=os.getcwd(),   # the working directory, so a copy of the repository renders against itself
        cluster_column=review_cluster_column,
        v8_decisions=lambda wildcards: review_v8_decisions(wildcards) or "None",
        july_object=lambda wildcards: review_july_object(wildcards) or "None",
        v8_release=lambda wildcards: f"{V8_RELEASE_DIR}/{V8_RELEASE_OBJECT[wildcards.species]}",
        out_dir=REVIEW_DIR,
    log:
        "logs/99b_review_{species}_{step}.log",
    wildcard_constraints:
        step="atlas|" + "|".join(REVIEW_LINEAGE_ORDER),
    resources:
        mem_mb=lambda wildcards: 48000 if wildcards.step in ("atlas", "stromal") else 16000,
        runtime=90,
    shell:
        """
        papermill {input.notebook} {output.notebook} \
            -p REPO {params.repo} \
            -p SPECIES {wildcards.species} \
            -p STEP {wildcards.step} \
            -p OBJECT_PATH {input.adata} \
            -p CLUSTER_COLUMN {params.cluster_column} \
            -p DECISIONS_PATH {input.decisions} \
            -p VOCABULARY_PATH {input.vocabulary} \
            -p PANELS_PATH {input.panels} \
            -p REASONS_PATH {input.reasons} \
            -p STEPS_PATH {input.steps} \
            -p V8_RELEASE_PATH {params.v8_release} \
            -p V8_DECISIONS_PATH {params.v8_decisions} \
            -p JULY_OBJECT_PATH {params.july_object} \
            -p OUT_DIR {params.out_dir} \
            --cwd {params.repo} \
            --log-output \
            > {log} 2>&1
        """


# One PDF per species with one page per step: the all-steps review dot plots used for approval.
rule review_dotplots_all_steps:
    input:
        lambda wildcards: [f"{REVIEW_DIR}/{wildcards.species}_{step}_dotplot.pdf" for step in review_steps(wildcards.species)],
    output:
        f"{REVIEW_DIR}/all_steps_dotplots_{{species}}_{REVIEW_VERSION}.pdf",
    log:
        "logs/99b_review_dotplots_all_steps_{species}.log",
    shell:
        "pdfunite {input} {output} > {log} 2>&1"


# The removed-cluster ledger: every EXCLUDE cluster of every step with its QC, kit verdict, final reason,
# evidence, v8 origin and decision provenance, plus a short summary.
rule removed_cluster_ledger:
    input:
        removed=review_outputs(REVIEW_DIR + "/{species}_{step}_removed_clusters.csv"),
    output:
        csv=f"{REVIEW_DIR}/removed_cluster_ledger_{REVIEW_VERSION}.csv",
        md=f"{REVIEW_DIR}/removed_cluster_ledger_{REVIEW_VERSION}.md",
    params:
        version=REVIEW_VERSION,
    log:
        "logs/99b_removed_cluster_ledger.log",
    script:
        "../scripts/review/removed_cluster_ledger.py"


# Target: the 13 rendered notebooks, both all-steps dot plot PDFs and the removed-cluster ledger.
rule review_notebooks:
    input:
        review_outputs(REVIEW_NOTEBOOK_DIR + "/{species}_{step}_review.ipynb"),
        expand(f"{REVIEW_DIR}/all_steps_dotplots_{{species}}_{REVIEW_VERSION}.pdf", species=list(config["lineages"])),
        rules.removed_cluster_ledger.output,
