# workflow/rules/14_liana.smk
# =============================================================================
# LIANA cell-cell communication, per species x variant.
#   no_commot : rank_aggregate only when skip_spatial_no_commot (default) -> the NON-SPATIAL deliverable
#               (liana_res_full/significant.csv, top_lr_pairs.csv); with skip off it also runs Moran's I.
#   bivariate : full spatial arm (spatial subset, Moran's I, LIANA bivariate Moran's R + local cosine);
#               feeds the spatial figures; run by targeting results/liana/bivariate/{species}/.liana_complete
#
# Scripts were ported 2026-09-08 from
#   /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/liana_pipeline_final
# into workflow/scripts/liana/ with four changes (see the scripts' provenance notes):
#   * --output-base is required (the legacy hardcoded /n/eddy_lab output root is dead)
#   * --h5ad is required (legacy default pointed at /n/eddy_lab)
#   * layers['log1p'] is recomputed from layers[counts_layer] with normalize_total(target_sum)
#   * mouse timepoint strings follow the 2026_08_28 deliverable vocabulary
#     (1_dpi_early ... 4_dpi_bleeding; June inputs used 1dpi_early ... 4dpi_menstrual)
# The June 2026 manual outputs (liana_pipeline_final/outputs*) are NOT touched.
#
# Output declared to Snakemake is a per-(variant,species) sentinel; the real tables land
# next to it: results/liana/{variant}/{species}/{timepoint}/{cell_type_key}/tables/*.csv
#
# config.yaml `liana:` block: python, scripts_dir, cache_dir, expr_prop, cell_type_key,
#   cell_type_map, counts_layer, target_sum, rank_column, rank_cutoff, top_n_lr, moran_i_min, skip_spatial_no_commot, variants,
#   runs.{species}.{h5ad,timepoints}
# Entry point: `snakemake --profile profiles/slurm_noconda liana_ccc_all`
# NOTE: rule liana (20_liana.smk) is an unrelated legacy stub; rule all still requests its
#   results/tables/liana_{species}.csv, so liana_ccc_all must be targeted explicitly.
# =============================================================================
import os

LIANA_CFG   = config["liana"]
LIANA_PY    = LIANA_CFG["python"]
LIANA_SDIR  = LIANA_CFG.get("scripts_dir", "workflow/scripts/liana")
LIANA_CACHE = LIANA_CFG.get("cache_dir", os.path.expanduser("~/.cache"))
LIANA_EP    = LIANA_CFG.get("expr_prop", 0.1)
LIANA_CKEY  = LIANA_CFG.get("cell_type_key", "cell_type")
LIANA_CTMAP = LIANA_CFG.get("cell_type_map")            # YAML -> obs cell_type_merged / cell_type_xspecies (merge_immune_xspecies port)
LIANA_CNTS  = LIANA_CFG.get("counts_layer", "counts")
LIANA_TSUM  = LIANA_CFG.get("target_sum", 10000)
LIANA_RCOL  = LIANA_CFG.get("rank_column", "magnitude_rank")  # LR discovery column (Dimitrov 2024: magnitude rank < 0.05)
LIANA_RCUT  = LIANA_CFG.get("rank_cutoff", 0.05)             # LR discovery: rank_column < cutoff
LIANA_SKIPSP = LIANA_CFG.get("skip_spatial_no_commot", True) # no_commot variant stops after rank_aggregate (non-spatial deliverable)
LIANA_IMIN  = LIANA_CFG.get("moran_i_min", 0.01)          # Moran's I floor for spatially variable interactions
LIANA_TOPN  = LIANA_CFG.get("top_n_lr", None)             # optional cap after the specificity filter
LIANA_RUNS  = LIANA_CFG["runs"]                       # species -> {h5ad, timepoints}
LIANA_VARS  = LIANA_CFG.get("variants", ["no_commot"])   # variants in the liana_ccc_all target; any variant can be targeted by sentinel
_LIANA_SUFFIX = {"no_commot": "", "bivariate": "_bivariate"}

wildcard_constraints:
    li_species = "|".join(LIANA_RUNS),
    li_variant = "no_commot|bivariate",

def _li_script(wc):
    return f"{LIANA_SDIR}/liana_pipeline_no_commot_{wc.li_species}{_LIANA_SUFFIX[wc.li_variant]}.py"

rule liana_ccc:
    input:
        h5ad    = lambda wc: ancient(LIANA_RUNS[wc.li_species]["h5ad"]),
        script  = _li_script,
        helpers = f"{LIANA_SDIR}/plot_helpers.py",
        ctmap   = [LIANA_CTMAP] if LIANA_CTMAP else [],
    output:
        done = "results/liana/{li_variant}/{li_species}/.liana_complete",
    params:
        py = LIANA_PY, ckey = LIANA_CKEY, ep = LIANA_EP, cnts = LIANA_CNTS, tsum = LIANA_TSUM,
        rcol = LIANA_RCOL, rcut = LIANA_RCUT, imin = LIANA_IMIN,
        skipsp = lambda wc: "--skip-spatial" if (wc.li_variant == "no_commot" and LIANA_SKIPSP) else "",
        ctmap  = f"--cell-type-map {LIANA_CTMAP}" if LIANA_CTMAP else "",
        topn = f"--top-n-lr {LIANA_TOPN}" if LIANA_TOPN else "",
        cache   = LIANA_CACHE,
        tps     = lambda wc: " ".join(LIANA_RUNS[wc.li_species]["timepoints"]),
        outbase = lambda wc, output: os.path.dirname(output.done),
    threads: 4
    resources:
        mem_mb = 64000, runtime = 720,   # June 2026 manual runs: mouse 74 min, human 52 min (sacct)
    log:
        "logs/liana/{li_variant}_{li_species}.log",
    shell:
        r"""
        export XDG_CACHE_HOME={params.cache} NUMBA_CACHE_DIR={params.cache}/numba MPLCONFIGDIR={params.cache}/mpl
        mkdir -p $XDG_CACHE_HOME $NUMBA_CACHE_DIR $MPLCONFIGDIR {params.outbase} $(dirname {log})
        {params.py} {input.script} \
            --h5ad {input.h5ad} \
            --cell-type-key {params.ckey} {params.ctmap} \
            --expr-prop {params.ep} \
            --counts-layer {params.cnts} \
            --target-sum {params.tsum} \
            --rank-column {params.rcol} --rank-cutoff {params.rcut} \
            --moran-i-min {params.imin} {params.topn} {params.skipsp} \
            --output-base {params.outbase} \
            --timepoints {params.tps} > {log} 2>&1
        touch {output.done}
        """

rule liana_ccc_all:
    input:
        expand("results/liana/{v}/{sp}/.liana_complete", v=LIANA_VARS, sp=list(LIANA_RUNS)),
