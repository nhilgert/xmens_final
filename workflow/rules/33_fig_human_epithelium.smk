# workflow/rules/33_fig_human_epithelium.smk
# Fig. S26B-C of the revision: the human epithelium split into inner and outer layers, drawn from the release of
# record (config deliverables.version) to sit under the S26 grid (S26A: rule panels_human_tissue_s26,
# 32_fig_human_tissue.smk):
#   Fig. S26B   the S26 sections with the epithelial nuclei coloured inner / outer / ciliated, two backgrounds:
#                 grey    every other nucleus in the grey of the S26 grid
#                 decfib  decidual cells and endometrial fibroblasts, the classes the layer rule counts, in pale
#                         tints of their Fig. 5B colours; every other nucleus grey
#   Fig. S26C   UMAP of the epithelial nuclei, the two options of Fig. 3C:
#                 lineage_fit  the epithelial lineage fit on which the layers were decided (the fit is not deposited)
#                 atlas        the release's whole-tissue embedding, zoomed to its epithelial island
#   Fig. S26BC  S26B and S26C side by side, each background with each UMAP option
#
#   panels_human_epithelium_lineage_umap   lineage-fit UMAP of the released epithelial nuclei; the fit's sha256 is
#                                          checked against the release's own record of it (uns epithelial_reinstatement)
#   panels_human_epithelium_{s26b,s26c,s26bc}  one rule per layout, all drawn by panel_s26bc_epithelium.py
#
# Inputs: the plot table of rule panels_human_tissue_extract (results/figures/panels/human_tissue/data/<release>/), and
# panel_common.py and the S26 script, whose sections, display transform, spatial scale, point size and map geometry
# S26B reuses. Outputs, next to the S26 files:
# results/figures/panels/human_tissue/panel_<S26B_*|S26C_*|S26BC_*>_<tag>.{pdf,svg,png,provenance.json}. Only the
# release of record is drawn; the files of earlier tags stay as they are. Target: panels_human_epithelium_all.
#
# The lineage fit is ancient(): a figure never refits a lineage because an mtime moved. The panel scripts import
# panel_common.py and the S26 script from their own directory, which PYTHONPATH names because Python's safe-path
# mode (PYTHONSAFEPATH) drops it from sys.path. They run under the panels interpreter (config panels.python) with
# caches off HOME (PANEL_ENV, 91_panels_versions.smk).

include: "00_panels.smk"

PANELS_HE_RELEASES = [PANEL_LIVE_TAG]
PANELS_HE_FIT = "data/interim/latent_human_epithelial.h5ad"
PANELS_HE_OUT = "results/figures/panels/human_tissue"
PANELS_HE_TABLE = PANELS_HE_OUT + "/data/{he_release}"  # written by rule panels_human_tissue_extract
PANELS_HE_SCRIPTS = "workflow/scripts/figures/panels/human_tissue"
PANELS_HE_LOGS = "logs/panels_human_epithelium"

PANELS_HE_BACKGROUNDS = ["grey", "decfib"]
PANELS_HE_EMBEDDINGS = ["lineage_fit", "atlas"]
PANELS_HE_EXTENSIONS = [".pdf", ".svg", ".png", ".provenance.json"]

# Every input a panel rule reads besides the lineage table: the plot table, the panel script and the two
# modules it builds on, so a change to the S26 grid's geometry or transform redraws S26B.
PANELS_HE_PANEL_INPUTS = {
    "cells": f"{PANELS_HE_TABLE}/cells.parquet",
    "uns": f"{PANELS_HE_TABLE}/uns.json",
    "extract": f"{PANELS_HE_TABLE}/provenance.json",
    "script": f"{PANELS_HE_SCRIPTS}/panel_s26bc_epithelium.py",
    "common": f"{PANELS_HE_SCRIPTS}/panel_common.py",
    "s26": f"{PANELS_HE_SCRIPTS}/panel_s26_spatial_highlight.py",
}

PANELS_HE_SHELL = """
    {params.env}
    PYTHONPATH={params.module_dir} {params.py} -u {input.script} \\
        --cells {input.cells} \\
        --uns {input.uns} \\
        --extract-provenance {input.extract} \\
        --release {wildcards.he_release} \\
        {params.layout_args} \\
        --out-stem {params.stem} \\
        > {log} 2>&1
    """


wildcard_constraints:
    he_release="|".join(PANELS_HE_RELEASES),
    he_background="|".join(PANELS_HE_BACKGROUNDS),
    he_embedding="|".join(PANELS_HE_EMBEDDINGS),


def panels_he_lineage_inputs(wildcards):
    """Input function: the lineage-fit table and its provenance, for the lineage_fit UMAP option only."""
    if wildcards.he_embedding != "lineage_fit":
        return {}
    table_dir = PANELS_HE_TABLE.format(he_release=wildcards.he_release)
    return {
        "lineage_table": f"{table_dir}/lineage_umap.parquet",
        "lineage_provenance": f"{table_dir}/lineage_umap.provenance.json",
    }


def panels_he_umap_args(wildcards, input):
    """Command-line arguments of one UMAP option."""
    args = f"--embedding {wildcards.he_embedding}"
    if wildcards.he_embedding == "lineage_fit":
        args += f" --lineage-table {input.lineage_table} --lineage-provenance {input.lineage_provenance}"
    return args


# Every S26B, S26C and S26BC file of the release of record.
rule panels_human_epithelium_all:
    input:
        expand(
            f"{PANELS_HE_OUT}/panel_S26B_{{background}}_{{release}}{{ext}}",
            background=PANELS_HE_BACKGROUNDS,
            release=PANELS_HE_RELEASES,
            ext=PANELS_HE_EXTENSIONS,
        ),
        expand(
            f"{PANELS_HE_OUT}/panel_S26C_{{embedding}}_{{release}}{{ext}}",
            embedding=PANELS_HE_EMBEDDINGS,
            release=PANELS_HE_RELEASES,
            ext=PANELS_HE_EXTENSIONS,
        ),
        expand(
            f"{PANELS_HE_OUT}/panel_S26BC_{{background}}_{{embedding}}_{{release}}{{ext}}",
            background=PANELS_HE_BACKGROUNDS,
            embedding=PANELS_HE_EMBEDDINGS,
            release=PANELS_HE_RELEASES,
            ext=PANELS_HE_EXTENSIONS,
        ),


# Lineage-fit UMAP coordinates of the released epithelial nuclei, read once for the S26C lineage_fit option.
rule panels_human_epithelium_lineage_umap:
    input:
        fit=ancient(PANELS_HE_FIT),
        cells=f"{PANELS_HE_TABLE}/cells.parquet",
        uns=f"{PANELS_HE_TABLE}/uns.json",
        script=f"{PANELS_HE_SCRIPTS}/extract_epithelial_lineage_umap.py",
    output:
        table=f"{PANELS_HE_TABLE}/lineage_umap.parquet",
        provenance=f"{PANELS_HE_TABLE}/lineage_umap.provenance.json",
    params:
        py=PANEL_PY,
        env=PANEL_ENV,
    log:
        f"{PANELS_HE_LOGS}/lineage_umap_{{he_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        """
        {params.env}
        {params.py} -u {input.script} \\
            --fit {input.fit} \\
            --cells {input.cells} \\
            --uns {input.uns} \\
            --release {wildcards.he_release} \\
            --out-table {output.table} \\
            --out-provenance {output.provenance} \\
            > {log} 2>&1
        """


# Fig. S26B: the S26 sections with the epithelium coloured by layer, one file set per background.
rule panels_human_epithelium_s26b:
    input:
        **PANELS_HE_PANEL_INPUTS,
    output:
        multiext(f"{PANELS_HE_OUT}/panel_S26B_{{he_background}}_{{he_release}}", *PANELS_HE_EXTENSIONS),
    params:
        layout_args=lambda wildcards: f"--layout spatial --background {wildcards.he_background}",
        stem=lambda wildcards: f"{PANELS_HE_OUT}/panel_S26B_{wildcards.he_background}_{wildcards.he_release}",
        module_dir=PANELS_HE_SCRIPTS,
        py=PANEL_PY,
        env=PANEL_ENV,
    log:
        f"{PANELS_HE_LOGS}/s26b_{{he_background}}_{{he_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        PANELS_HE_SHELL


# Fig. S26C: UMAP of the epithelial nuclei, one file set per embedding option.
rule panels_human_epithelium_s26c:
    input:
        unpack(panels_he_lineage_inputs),
        **PANELS_HE_PANEL_INPUTS,
    output:
        multiext(f"{PANELS_HE_OUT}/panel_S26C_{{he_embedding}}_{{he_release}}", *PANELS_HE_EXTENSIONS),
    params:
        layout_args=lambda wildcards, input: f"--layout umap {panels_he_umap_args(wildcards, input)}",
        stem=lambda wildcards: f"{PANELS_HE_OUT}/panel_S26C_{wildcards.he_embedding}_{wildcards.he_release}",
        module_dir=PANELS_HE_SCRIPTS,
        py=PANEL_PY,
        env=PANEL_ENV,
    log:
        f"{PANELS_HE_LOGS}/s26c_{{he_embedding}}_{{he_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        PANELS_HE_SHELL


# Fig. S26B and S26C side by side, every background with every UMAP option.
rule panels_human_epithelium_s26bc:
    input:
        unpack(panels_he_lineage_inputs),
        **PANELS_HE_PANEL_INPUTS,
    output:
        multiext(
            f"{PANELS_HE_OUT}/panel_S26BC_{{he_background}}_{{he_embedding}}_{{he_release}}",
            *PANELS_HE_EXTENSIONS,
        ),
    params:
        layout_args=lambda wildcards, input: (
            f"--layout both --background {wildcards.he_background} {panels_he_umap_args(wildcards, input)}"
        ),
        stem=lambda wildcards: (
            f"{PANELS_HE_OUT}/panel_S26BC_{wildcards.he_background}_{wildcards.he_embedding}_{wildcards.he_release}"
        ),
        module_dir=PANELS_HE_SCRIPTS,
        py=PANEL_PY,
        env=PANEL_ENV,
    log:
        f"{PANELS_HE_LOGS}/s26bc_{{he_background}}_{{he_embedding}}_{{he_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        PANELS_HE_SHELL
