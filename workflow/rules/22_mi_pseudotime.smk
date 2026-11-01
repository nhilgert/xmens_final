# workflow/rules/22_mi_pseudotime.smk
# Mutual information between gene expression and the harmonized pseudotime (Methods: MI).
#
# One arm = one set of nuclei scored against one pseudotime column:
#   [mi_control_input | mi_record_input] -> mi_prep -> mi_worker x n_chunks -> mi_merge -> <outdir>/<arm>/mi_all.csv
# Main arms read a fibroblast/decidual object directly. Specificity-control arms (mouse) first build an input object
# with pseudotime re-ranked within each sample. Pseudotime-of-record arms (human, decision G) first build an input
# object whose obs[pt_key] holds the pseudotime of record named in config["pseudotime_of_record"]. Arms, source
# objects and MI settings live in config["mi_pseudotime"]; scripts are workflow/scripts/mi/.
#
# The chunk count is a parameter of the result, not only of the parallelization: each chunk seeds its own
# permutation generator (seed * 100000 + chunk) and shares its permutations across its genes (mi_pseudotime.py,
# worker), so p and q values depend on it. An arm's n_chunks overrides the global one; the manifest of each arm
# records the count it ran with.
#
# The main arms of 22-23 Oct were computed by the sbatch arrays in their arm directories (worker.sbatch) and adopted
# by Snakemake on 26 Oct; their worker logs are empty. The control arms ran under these rules.
import re

MI = config["mi_pseudotime"]
MI_ROOT = MI["outdir"]
MI_ARMS = MI["arms"]
MI_CONTROL_ARMS = [arm for arm, spec in MI_ARMS.items() if "control" in spec]
MI_RECORD_ARMS = [arm for arm, spec in MI_ARMS.items() if "pseudotime_of_record" in spec]


wildcard_constraints:
    mi_arm="|".join(re.escape(arm) for arm in MI_ARMS),
    control_arm="|".join(re.escape(arm) for arm in MI_CONTROL_ARMS) or "^$",
    record_arm="|".join(re.escape(arm) for arm in MI_RECORD_ARMS) or "^$",
    chunk=r"\d{3,4}",  # mi_{chunk:03d}.csv has four digits from chunk 1000 on


def mi_n_chunks(arm):
    """Chunk count of one arm: its own n_chunks, else the global one."""
    return int(MI_ARMS[arm].get("n_chunks", MI["n_chunks"]))


def mi_pseudotime_of_record(arm):
    """The config pseudotime_of_record entry (path, column, sha256) that a record arm reads."""
    return config["pseudotime_of_record"][MI_ARMS[arm]["pseudotime_of_record"]]


def mi_arm_source(arm):
    """Object an arm scores: the fibroblast/decidual object, or the built input of a control or record arm."""
    spec = MI_ARMS[arm]
    if "control" in spec or "pseudotime_of_record" in spec:
        return f"{MI_ROOT}/{arm}/input.h5ad"
    return MI["objects"][spec["object"]]


def mi_source_h5ad(wc):
    """Object prep reads (see mi_arm_source)."""
    return mi_arm_source(wc.mi_arm)


def mi_prep_pt_key(wc):
    """Control inputs carry the within-sample rank as pt_norm."""
    spec = MI_ARMS[wc.mi_arm]
    return "pt_norm" if "control" in spec else spec["pt_key"]


def mi_chunk_csvs(wc):
    """All per-chunk result files of one arm."""
    return [f"{MI_ROOT}/{wc.mi_arm}/results/mi_{c:03d}.csv" for c in range(mi_n_chunks(wc.mi_arm))]


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


# Input object for a pseudotime-of-record arm: the fibroblast/decidual object with obs[pt_key] replaced by the
# record's column, matched by nucleus name; the record's sha256 is checked and both nucleus sets must be equal.
rule mi_record_input:
    input:
        h5ad=lambda wc: MI["objects"][MI_ARMS[wc.record_arm]["object"]],
        table=lambda wc: mi_pseudotime_of_record(wc.record_arm)["path"],
        script="workflow/scripts/mi/build_mi_record_input.py",
    output:
        h5ad=f"{MI_ROOT}/{{record_arm}}/input.h5ad",
        meta=f"{MI_ROOT}/{{record_arm}}/build_meta.json",
    params:
        python=MI["python"],
        pt_key=lambda wc: MI_ARMS[wc.record_arm]["pt_key"],
        column=lambda wc: mi_pseudotime_of_record(wc.record_arm)["column"],
        sha256=lambda wc: mi_pseudotime_of_record(wc.record_arm)["sha256"],
    log:
        "logs/mi_dpt/{record_arm}_input.log",
    threads: 1
    resources:
        mem_mb=48000,
        runtime=30,
    shell:
        r"""
        {params.python} {input.script} \
            --h5ad {input.h5ad} \
            --table {input.table} \
            --column {params.column} \
            --expect-sha256 {params.sha256} \
            --pt-key {params.pt_key} \
            --out {output.h5ad} \
            --meta {output.meta} \
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
        chunks=lambda wc: mi_n_chunks(wc.mi_arm),
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
        mem_mb=lambda wc: MI_ARMS[wc.mi_arm].get("worker_mem_mb", 32000),
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


# Cross-species overlap of the MI arms and the Fig. 5F gene selection, settings in config["mi"]["overlap"].
# Extracted from notebooks/mi_dpt/mi_dpt_overlap_harmonized_nall.ipynb (sections 1-7 and 9) and from
# crossspecies_gene_selection.py, which were run by hand into results/mi_dpt/overlap_harmonized_nall.
MI_OVERLAP = config["mi"]["overlap"]
OVERLAP_DIR = MI_OVERLAP["outdir"]
OVERLAP_ARMS = [MI_OVERLAP["mouse_arm"], MI_OVERLAP["human_arm"], *MI_OVERLAP["control_arms"]]
OVERLAP_BUILT_ARMS = [arm for arm in OVERLAP_ARMS if arm in MI_CONTROL_ARMS or arm in MI_RECORD_ARMS]
HEATMAP_DIR = f"{OVERLAP_DIR}/crossspecies_heatmap_selection"


# Methods overlap (N, K, n, k, hypergeometric p), the Fig. S16 C-D percentile sweep and the notebook's overlap tables.
rule mi_overlap:
    input:
        mi=expand(f"{MI_ROOT}/{{arm}}/mi_all.csv", arm=OVERLAP_ARMS),
        manifests=expand(f"{MI_ROOT}/{{arm}}/manifest.json", arm=OVERLAP_ARMS),
        builds=expand(f"{MI_ROOT}/{{arm}}/build_meta.json", arm=OVERLAP_BUILT_ARMS),
        orthologs=MI_OVERLAP["orthologs"],
        script="workflow/scripts/mi/mi_overlap.py",
    output:
        pairs=f"{OVERLAP_DIR}/ortholog_pairs_one2one.csv",
        summary=f"{OVERLAP_DIR}/overlap_summary.csv",
        null=f"{OVERLAP_DIR}/overlap_permutation_null.csv",
        cutoff_sweep=f"{OVERLAP_DIR}/overlap_cutoff_sweep.csv",
        shared=f"{OVERLAP_DIR}/shared_top_decile_genes.json",
        specificity=f"{OVERLAP_DIR}/specificity_shared_universe.csv",
        methods=f"{OVERLAP_DIR}/methods_overlap.json",
        percentile_sweep=f"{OVERLAP_DIR}/overlap_percentile_sweep.csv",
        report=f"{OVERLAP_DIR}/overlap_report.json",
    params:
        python=MI["python"],
        mi_dir=MI_ROOT,
        mouse_arm=MI_OVERLAP["mouse_arm"],
        human_arm=MI_OVERLAP["human_arm"],
        control_arms=" ".join(MI_OVERLAP["control_arms"]),
        top_fraction=MI_OVERLAP["top_fraction"],
        n_perm=MI_OVERLAP["n_perm"],
        seed=config["seed"],
        q_thresh=config["mi"]["q_thresh"],
        top_percentile=config["mi"]["top_percentile"],
        outdir=OVERLAP_DIR,
    log:
        "logs/mi_dpt/overlap.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        r"""
        {params.python} {input.script} \
            --mi-dir {params.mi_dir} \
            --mouse-arm {params.mouse_arm} \
            --human-arm {params.human_arm} \
            --control-arms {params.control_arms} \
            --orthologs {input.orthologs} \
            --top-fraction {params.top_fraction} \
            --n-perm {params.n_perm} \
            --seed {params.seed} \
            --q-thresh {params.q_thresh} \
            --top-percentile {params.top_percentile} \
            --outdir {params.outdir} \
            > {log} 2>&1
        """


# Fig. 5F gene selection on the nuclei and pseudotime the two main arms scored (the human arm: pseudotime of record).
rule mi_heatmap_selection:
    input:
        mi=expand(f"{MI_ROOT}/{{arm}}/mi_all.csv", arm=[MI_OVERLAP["mouse_arm"], MI_OVERLAP["human_arm"]]),
        pairs=f"{OVERLAP_DIR}/ortholog_pairs_one2one.csv",
        mouse_h5ad=mi_arm_source(MI_OVERLAP["mouse_arm"]),
        human_h5ad=mi_arm_source(MI_OVERLAP["human_arm"]),
        panel_genes=MI_OVERLAP["heatmap"]["panel_genes"],
        script="workflow/scripts/mi/crossspecies_gene_selection.py",
    output:
        percell_mouse=f"{HEATMAP_DIR}/percell_mouse_pool.tsv.gz",
        percell_human=f"{HEATMAP_DIR}/percell_human_pool.tsv.gz",
        concordance=f"{HEATMAP_DIR}/conserved_profiles_pool_concordance.csv",
        selected=f"{HEATMAP_DIR}/conserved_profiles_selected_genes.csv",
        summary=f"{HEATMAP_DIR}/selection_summary.json",
    params:
        python=MI["python"],
        mi_dir=MI_ROOT,
        mouse_arm=MI_OVERLAP["mouse_arm"],
        human_arm=MI_OVERLAP["human_arm"],
        mouse_pt_key=MI_ARMS[MI_OVERLAP["mouse_arm"]]["pt_key"],
        human_pt_key=MI_ARMS[MI_OVERLAP["human_arm"]]["pt_key"],
        top_fraction=MI_OVERLAP["top_fraction"],
        rmin=MI_OVERLAP["heatmap"]["rmin"],
        n_early=MI_OVERLAP["heatmap"]["n_early"],
        n_late=MI_OVERLAP["heatmap"]["n_late"],
        outdir=HEATMAP_DIR,
    log:
        "logs/mi_dpt/heatmap_selection.log",
    threads: 1
    resources:
        mem_mb=48000,
        runtime=60,
    shell:
        r"""
        {params.python} {input.script} \
            --mi-dir {params.mi_dir} \
            --mouse-arm {params.mouse_arm} \
            --human-arm {params.human_arm} \
            --pairs {input.pairs} \
            --mouse-h5ad {input.mouse_h5ad} \
            --mouse-pt-key {params.mouse_pt_key} \
            --human-h5ad {input.human_h5ad} \
            --human-pt-key {params.human_pt_key} \
            --panel-genes {input.panel_genes} \
            --top-fraction {params.top_fraction} \
            --rmin {params.rmin} \
            --n-early {params.n_early} \
            --n-late {params.n_late} \
            --outdir {params.outdir} \
            > {log} 2>&1
        """


# Target: overlap tables, Methods overlap numbers and the Fig. 5F selection.
rule mi_overlap_all:
    input:
        rules.mi_overlap.output,
        rules.mi_heatmap_selection.output,
