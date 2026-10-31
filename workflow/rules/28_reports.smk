import os

# 28_reports -- OFF the critical path on purpose.
#
# Executes the public de-novo notebooks with papermill to refresh figures and
# produce an executed copy. Nothing in rule all depends on this, so a broken
# dotplot can never stall the pipeline:
#
#     snakemake --profile workflow/profiles/slurm reports
#
# The pipeline and the notebooks share workflow/scripts/atlas/apply_annotation.py, so
# they cannot drift -- this only re-renders the human-readable record.

REPORT_NOTEBOOKS = [f"04_annotate_{sp}_denovo" for sp in config["lineages"]] + [
    f"04c_annotate_{sp}_{lin}_denovo"
    for sp, lins in config["lineages"].items() for lin in lins
]

# Notebooks that also write an annotated lineage object, tables and figures. Their output
# directory is declared, so Snakemake clears it before a re-run (the notebooks refuse to
# overwrite an existing annotated object); its path and the decisions table reach the
# notebook as papermill parameters.
ANNOTATION_NOTEBOOKS = {
    "04c_human_epithelial_layers": {
        "decisions": "resources/annotations/human_epithelial.csv",
        "object": "annotated_human_epithelial.h5ad",
    },
    "04c_mouse_epithelial_compartments": {
        "decisions": "resources/annotations/mouse_epithelial.csv",
        "object": "annotated_mouse_epithelial.h5ad",
    },
}

rule render_notebook:
    input:
        nb="notebooks/public/{name}.ipynb",
    output:
        nb="results/notebooks/{name}.executed.ipynb",
    params:
        repo=config.get("repo", "."),
    log:
        "logs/99_report_{name}.log",
    shell:
        "papermill {input.nb} {output.nb} -p REPO {params.repo} "
        "--log-output --cwd {params.repo} > {log} 2>&1"

rule render_annotation_notebook:
    input:
        nb="notebooks/public/{name}.ipynb",
        decisions=lambda wc: ANNOTATION_NOTEBOOKS[wc.name]["decisions"],
    output:
        nb="results/notebooks/{name}.executed.ipynb",
        outdir=directory("results/annotation/{name}"),
    wildcard_constraints:
        name="|".join(ANNOTATION_NOTEBOOKS),
    params:
        repo=config.get("repo", "."),
        object=lambda wc: ANNOTATION_NOTEBOOKS[wc.name]["object"],
    log:
        "logs/99_report_{name}.log",
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        "papermill {input.nb} {output.nb} -p REPO {params.repo} "
        "-p DECISIONS_PATH {input.decisions} -p OUT_DIR {output.outdir} "
        "-p OUT_OBJECT_PATH {output.outdir}/{params.object} "
        "--log-output --cwd {params.repo} > {log} 2>&1"

# The annotation notebooks match render_notebook too; their own rule declares the outputs.
ruleorder: render_annotation_notebook > render_notebook

rule reports:
    input:
        expand("results/notebooks/{name}.executed.ipynb",
               name=[n for n in [*REPORT_NOTEBOOKS, *ANNOTATION_NOTEBOOKS]
                     if os.path.exists(f"notebooks/public/{n}.ipynb")]),
