# 28_reports -- OFF the critical path on purpose.
#
# Executes the two epithelial annotation notebooks of notebooks/public with papermill and keeps an executed copy.
# Nothing in rule all depends on this, so a broken notebook can never stall the pipeline:
#
#     snakemake --profile workflow/profiles/slurm reports
#
# The notebooks apply their decision tables themselves (not workflow/scripts/atlas/apply_annotation.py, which the
# release uses) and check them against tissue position; they render against the working directory, which the shell
# passes ("$PWD"), so a copy of the repository writes into itself and no recorded parameter names a path.

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
    log:
        "logs/99_report_{name}.log",
    shell:
        "papermill {input.nb} {output.nb} -p REPO \"$PWD\" "
        "--log-output --cwd \"$PWD\" > {log} 2>&1"

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
        object=lambda wc: ANNOTATION_NOTEBOOKS[wc.name]["object"],
    log:
        "logs/99_report_{name}.log",
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        "papermill {input.nb} {output.nb} -p REPO \"$PWD\" "
        "-p DECISIONS_PATH {input.decisions} -p OUT_DIR {output.outdir} "
        "-p OUT_OBJECT_PATH {output.outdir}/{params.object} "
        "--log-output --cwd \"$PWD\" > {log} 2>&1"

# The annotation notebooks match render_notebook too; their own rule declares the outputs.
ruleorder: render_annotation_notebook > render_notebook

rule reports:
    input:
        expand("results/notebooks/{name}.executed.ipynb", name=list(ANNOTATION_NOTEBOOKS)),
