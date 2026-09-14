# workflow/rules/17_liana_pseudotime.smk
# =============================================================================
# LIANA along the fib/dec decidualization trajectory (added 2026-09-11).
#
# Scope (user decision): fib/dec lineage only, rank-normalised DPT as the axis. The 2026_08_28
# fib/dec deliverable objects carry the coordinate per cell, so no join is needed:
#   mouse dpt_pseudotime_rank_norm (43,197 cells) | human dpt_pseudotime_rank_normalized (47,533)
#
# The pseudotime axis is quantile-binned (equal counts) and the BINS are passed as
# --cell-type-key, so rank_aggregate scores every (source bin -> target bin) context:
#   diagonal      = within-stage signalling
#   upper/lower   = earlier->later vs later->earlier along the trajectory (a directional
#                   read-out from the standard scoring; no Granger machinery involved)
# --pool-timepoints collapses obs['timepoint'] so pseudotime is the only grouping axis.
#
# Non-spatial by construction (--skip-spatial): scores come from rank_aggregate on all fib/dec
# nuclei, with magnitude_rank < 0.05 (Dimitrov et al. 2024) as in the main deliverable. Only the
# no_commot scripts take --pt-key; the bivariate variant is not wired here (its bandwidth table
# and coordinate exclusions are keyed by timepoint).
#
# Entry point: `snakemake --profile profiles/slurm_noconda liana_pt_all`
# =============================================================================
import os

LPT_CFG   = config["liana"]
LPT       = LPT_CFG["pseudotime"]
LPT_RUNS  = LPT["runs"]                                  # species -> {h5ad, pt_key}
LPT_BINS  = LPT.get("bins", 5)
LPT_OUT   = LPT.get("outdir", "results/liana_pseudotime")
LPT_SDIR  = LPT_CFG.get("scripts_dir", "workflow/scripts/liana")
LPT_CACHE = LPT_CFG.get("cache_dir", os.path.expanduser("~/.cache"))

wildcard_constraints:
    lpt_species = "|".join(LPT_RUNS),

rule liana_pt:
    input:
        h5ad   = lambda wc: ancient(LPT_RUNS[wc.lpt_species]["h5ad"]),
        script = lambda wc: f"{LPT_SDIR}/liana_pipeline_no_commot_{wc.lpt_species}.py",
        helpers = f"{LPT_SDIR}/plot_helpers.py",
    output:
        done = f"{LPT_OUT}/{{lpt_species}}/.liana_pt_complete",
    params:
        py      = LPT_CFG["python"],
        ptkey   = lambda wc: LPT_RUNS[wc.lpt_species]["pt_key"],
        bins    = LPT_BINS,
        ep      = LPT_CFG.get("expr_prop", 0.1),
        cnts    = LPT_CFG.get("counts_layer", "counts"),
        tsum    = LPT_CFG.get("target_sum", 10000),
        rcol    = LPT_CFG.get("rank_column", "magnitude_rank"),
        rcut    = LPT_CFG.get("rank_cutoff", 0.05),
        cache   = LPT_CACHE,
        outbase = lambda wc, output: os.path.dirname(output.done),
    threads: 4
    resources:
        mem_mb = 64000, runtime = 240,
    log:
        f"logs/liana/pseudotime_{{lpt_species}}.log",
    shell:
        r"""
        export XDG_CACHE_HOME={params.cache} NUMBA_CACHE_DIR={params.cache}/numba MPLCONFIGDIR={params.cache}/mpl
        mkdir -p $XDG_CACHE_HOME $NUMBA_CACHE_DIR $MPLCONFIGDIR {params.outbase} $(dirname {log})
        {params.py} {input.script} \
            --h5ad {input.h5ad} \
            --pt-key {params.ptkey} --pt-bins {params.bins} --pool-timepoints \
            --cell-type-key pt_bin \
            --expr-prop {params.ep} --counts-layer {params.cnts} --target-sum {params.tsum} \
            --rank-column {params.rcol} --rank-cutoff {params.rcut} --skip-spatial \
            --output-base {params.outbase} --timepoints all > {log} 2>&1
        touch {output.done}
        """

rule liana_pt_all:
    input:
        expand(f"{LPT_OUT}/{{sp}}/.liana_pt_complete", sp=list(LPT_RUNS)),
