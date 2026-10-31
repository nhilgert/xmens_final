# workflow/rules/22_mi_pseudotime.smk
# Mutual information between gene expression and the harmonized pseudotime (Methods: MI).
#
# One arm = one set of nuclei scored against one pseudotime column:
#   [mi_control_input] -> mi_prep -> mi_worker x n_chunks -> mi_merge -> <outdir>/<arm>/mi_all.csv
# Main arms read a release object directly; specificity-control arms (mouse) first build
# an input object with pseudotime re-ranked within each sample. Arms, source objects and
# MI settings live in config["mi_pseudotime"]; scripts are workflow/scripts/mi/.
import re

MI = config["mi_pseudotime"]
MI_ROOT = MI["outdir"]
MI_ARMS = MI["arms"]
MI_CONTROL_ARMS = [arm for arm, spec in MI_ARMS.items() if "control" in spec]


wildcard_constraints:
    mi_arm="|".join(re.escape(arm) for arm in MI_ARMS),
    control_arm="|".join(re.escape(arm) for arm in MI_CONTROL_ARMS) or "^$",
    chunk=r"\d{3}",


def mi_source_h5ad(wc):
    """Object prep reads: the release object, or the built input for a control arm."""
    spec = MI_ARMS[wc.mi_arm]
    if "control" in spec:
        return f"{MI_ROOT}/{wc.mi_arm}/input.h5ad"
    return MI["objects"][spec["object"]]


def mi_prep_pt_key(wc):
    """Control inputs carry the within-sample rank as pt_norm."""
    spec = MI_ARMS[wc.mi_arm]
    return "pt_norm" if "control" in spec else spec["pt_key"]


def mi_chunk_csvs(wc):
    """All per-chunk result files of one arm."""
    return [f"{MI_ROOT}/{wc.mi_arm}/results/mi_{c:03d}.csv" for c in range(MI["n_chunks"])]


# Input object for a specificity-control arm: the uninduced nuclei, or a seeded
# size-matched draw of induced nuclei, with pseudotime ranked within each sample.
rule mi_control_input:
    input:
        h5ad=lambda wc: MI["objects"][MI_ARMS[wc.control_arm]["object"]],
        script="workflow/scripts/mi/build_mi_control_input.py",
    output:
        h5ad=f"{MI_ROOT}/{{control_arm}}/input.h5ad",
        meta=f"{MI_ROOT}/{{control_arm}}/build_meta.json",
    params:
        python=MI["python"],
        pt_key=lambda wc: MI_ARMS[wc.control_arm]["pt_key"],
        cell_set=lambda wc: MI_ARMS[wc.control_arm]["control"],
        seed=MI["control_seed"],
    log:
        "logs/mi_dpt/{control_arm}_input.log",
    threads: 1
    resources:
        mem_mb=48000,
        runtime=30,
    shell:
        r"""
        {params.python} {input.script} \
            --h5ad {input.h5ad} \
            --pt-key {params.pt_key} \
            --cell-set {params.cell_set} \
            --seed {params.seed} \
            --out {output.h5ad} \
            > {log} 2>&1
        """


# Expression matrix, pseudotime vector, gene list and manifest shared by all workers of an arm.
rule mi_prep:
    input:
        h5ad=mi_source_h5ad,
        script="workflow/scripts/mi/mi_pseudotime.py",
    output:
        expr=f"{MI_ROOT}/{{mi_arm}}/expr.npz",
        pt=f"{MI_ROOT}/{{mi_arm}}/pseudotime.npy",
        genes=f"{MI_ROOT}/{{mi_arm}}/genes.npy",
        manifest=f"{MI_ROOT}/{{mi_arm}}/manifest.json",
    params:
        python=MI["python"],
        outdir=lambda wc: f"{MI_ROOT}/{wc.mi_arm}",
        pt_key=mi_prep_pt_key,
        layer=MI["layer"],
        min_counts=MI["min_counts"],
        min_detect_frac=MI["min_detect_frac"],
        perms=MI["perms"],
        subsample_cells=MI["subsample_cells"],
        chunks=MI["n_chunks"],
        seed=MI["seed"],
    log:
        "logs/mi_dpt/{mi_arm}_prep.log",
    threads: 1
    resources:
        mem_mb=64000,
        runtime=60,
    shell:
        r"""
        {params.python} {input.script} prep \
            --h5ad {input.h5ad} \
            --pt-key {params.pt_key} \
            --drop-missing-pt \
            --layer {params.layer} \
            --min-counts {params.min_counts} \
            --min-detect-frac {params.min_detect_frac} \
            --perms {params.perms} \
            --subsample-cells {params.subsample_cells} \
            --chunks {params.chunks} \
            --seed {params.seed} \
            --outdir {params.outdir} \
            > {log} 2>&1
        """


# MI and permutation p-values for one chunk of genes (one SLURM job per chunk).
rule mi_worker:
    input:
        manifest=f"{MI_ROOT}/{{mi_arm}}/manifest.json",
        script="workflow/scripts/mi/mi_pseudotime.py",
    output:
        csv=f"{MI_ROOT}/{{mi_arm}}/results/mi_{{chunk}}.csv",
    params:
        python=MI["python"],
        outdir=lambda wc: f"{MI_ROOT}/{wc.mi_arm}",
        chunk=lambda wc: int(wc.chunk),
        perms=MI["perms"],
        subsample_cells=MI["subsample_cells"],
    log:
        "logs/mi_dpt/{mi_arm}_worker_{chunk}.log",
    threads: 1
    resources:
        mem_mb=32000,
        runtime=lambda wc: MI_ARMS[wc.mi_arm]["worker_runtime_min"],
        cpus_per_task=2,
    shell:
        r"""
        {params.python} {input.script} worker \
            --outdir {params.outdir} \
            --chunk {params.chunk} \
            --perms {params.perms} \
            --subsample-cells {params.subsample_cells} \
            > {log} 2>&1
        """


# Per-arm table of MI, permutation p and Benjamini-Hochberg q over all genes.
rule mi_merge:
    input:
        chunks=mi_chunk_csvs,
        script="workflow/scripts/mi/mi_pseudotime.py",
    output:
        mi_all=f"{MI_ROOT}/{{mi_arm}}/mi_all.csv",
    params:
        python=MI["python"],
        outdir=lambda wc: f"{MI_ROOT}/{wc.mi_arm}",
    log:
        "logs/mi_dpt/{mi_arm}_merge.log",
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        r"""
        {params.python} {input.script} merge \
            --outdir {params.outdir} \
            > {log} 2>&1
        """


# Target: merged MI tables for every configured arm.
rule mi_all:
    input:
        expand(f"{MI_ROOT}/{{mi_arm}}/mi_all.csv", mi_arm=list(MI_ARMS)),
