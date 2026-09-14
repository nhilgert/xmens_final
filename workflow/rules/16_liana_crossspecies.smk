# workflow/rules/16_liana_crossspecies.smk
# =============================================================================
# Cross-species LIANA conservation deliverable (liana_table.xlsx, liana_immune.xlsx, CSVs, figures).
# Port (2026-09-10) of liana_pipeline_final/cross_species_lr_explorer_v8-2.ipynb cells 5-31 +
# lr_analysis.py into workflow/scripts/liana_crossspecies.py. Consumes the NON-SPATIAL no_commot
# tables (results/liana/no_commot/{species}/{tp}/{cell_type_key}/tables/liana_res_full.csv) and
# the Ensembl one2one orthology in resources/orthology/. No spatial statistic enters the table.
# Entry point: `snakemake --profile profiles/slurm_noconda liana_crossspecies_all`
# NOTE: 22_crossspecies.smk (results/tables/conservation.csv) is the unrelated legacy MI-based rule.
# =============================================================================
import os

LX_CFG   = config["liana"]
LX_XS    = LX_CFG["crossspecies"]
LX_ROOT  = LX_XS.get("liana_root", "results/liana/no_commot")
LX_OUT   = LX_XS.get("outdir", "results/liana_crossspecies")
LX_CKEY  = LX_CFG.get("cell_type_key", "cell_type_merged")

rule liana_crossspecies:
    input:
        tables  = expand("results/liana/no_commot/{sp}/.liana_complete", sp=list(LX_CFG["runs"])),
        ortho   = LX_CFG["orthology"],
        ctmap   = LX_CFG["cell_type_map"],
        script  = "workflow/scripts/liana_crossspecies.py",
    output:
        xlsx    = f"{LX_OUT}/deliverables/liana_table.xlsx",
        immune  = f"{LX_OUT}/deliverables/liana_immune.xlsx",
        hits    = f"{LX_OUT}/deliverables/lr_hits.csv",
        pooled  = f"{LX_OUT}/deliverables/lr_pooled.csv",
        raw     = f"{LX_OUT}/deliverables/raw_calls.csv",
        summary = f"{LX_OUT}/tables/conservation_summary.json",
        venn    = f"{LX_OUT}/figures/fig1_conservation_venn.pdf",
    params:
        py      = LX_CFG["python"],
        root    = LX_ROOT, ckey = LX_CKEY, out = LX_OUT,
        cache   = LX_CFG.get("cache_dir", os.path.expanduser("~/.cache")),
        homtype = LX_CFG.get("orthology_homology_type", "ortholog_one2one"),
        rcol    = LX_CFG.get("rank_column", "magnitude_rank"),
        rcut    = LX_CFG.get("rank_cutoff", 0.05),
        mtps    = " ".join(LX_CFG["runs"]["mouse"]["timepoints"]),
        htps    = " ".join(LX_CFG["runs"]["human"]["timepoints"]),
    threads: 2
    resources:
        mem_mb = 16000, runtime = 60,
    log:
        "logs/liana/crossspecies.log",
    shell:
        r"""
        export XDG_CACHE_HOME={params.cache} MPLCONFIGDIR={params.cache}/mpl
        mkdir -p $MPLCONFIGDIR $(dirname {log})
        {params.py} {input.script} \
            --liana-root {params.root} --cell-type-key {params.ckey} \
            --mouse-timepoints {params.mtps} --human-timepoints {params.htps} \
            --orthology {input.ortho} --homology-type {params.homtype} \
            --cell-type-map {input.ctmap} \
            --rank-column {params.rcol} --rank-cutoff {params.rcut} \
            --outdir {params.out} > {log} 2>&1
        """

rule liana_crossspecies_all:
    input:
        rules.liana_crossspecies.output,
