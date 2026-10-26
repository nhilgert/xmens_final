# workflow/rules/16_zone_modules.smk
# =============================================================================
# Decidual-zone modules (mouse), added 2026-10-24. Per-nucleus AUCell scores for
# the 11 uterine microenvironments of Li et al. 2022 (Biol Reprod 107:529,
# Table 1, top 20 markers each; resources/signatures/aucell_modules_mouse.yaml),
# and their spatial panels and subtype heatmaps.
#
#   zone_aucell      score one run: pySCENIC AUCell on raw counts of the release
#                    object, fixed auc_threshold and seed, all genes ranked,
#                    only the nuclei of the run's timepoints (exact; see the
#                    script).
#   zone_plot_table  one per-nucleus CSV with everything the figures need:
#                    obs (sample, timepoint, cell_type, cell_subtype, x/y),
#                    embedding_outlier (config[zone_modules][embedding_outliers])
#                    and the 11 scores. Added 2026-10-24 so the figures can be
#                    redrawn from this CSV alone.
#   zone_figures     render config[zone_modules][figures] from the plot table
#                    only. Spatial panels hide embedding_outlier nuclei; heatmaps
#                    keep them (rule of 2026-10-22).
#
# One directory per run under results/zone_aucell/<run>/ and
# results/figures/zone_modules/<run>/. Add a run instead of editing one.
# Target: zone_modules_all (not in rule all).
#
# The release object and its MANIFEST are ancient(): a new release is cut
# deliberately (04g), never re-cut because an mtime moved.
# =============================================================================
import json

ZM      = config["zone_modules"]
ZM_RUNS = ZM["runs"]
# numba (pyscenic) and matplotlib caches must not land in the full HOME
ZM_ENV  = ("export XDG_CACHE_HOME={c} NUMBA_CACHE_DIR={c}/numba MPLCONFIGDIR={c}/mpl; "
           "mkdir -p {c}/numba {c}/mpl").format(c=ZM["cache_dir"])

wildcard_constraints:
    zone_run = "|".join(ZM_RUNS),


rule zone_aucell:
    input:
        h5ad     = ancient(ZM["h5ad"]),
        manifest = ancient(ZM["release_manifest"]),
        genesets = lambda wc: ZM_RUNS[wc.zone_run]["genesets"],
        script   = "workflow/scripts/zone_modules/score_zone_aucell.py",
    output:
        scores = "results/zone_aucell/{zone_run}/zone_aucell.csv",
        means  = "results/zone_aucell/{zone_run}/subtype_means.csv",
        prov   = "results/zone_aucell/{zone_run}/provenance.json",
    params:
        py         = ZM["python"],
        obj_name   = ZM["object_name"],
        thr        = lambda wc: ZM_RUNS[wc.zone_run]["auc_threshold"],
        seed       = lambda wc: ZM_RUNS[wc.zone_run]["seed"],
        norm       = lambda wc: ZM_RUNS[wc.zone_run].get("normalization", "none"),
        timepoints = lambda wc: ",".join(ZM_RUNS[wc.zone_run].get("timepoints") or []),
        outdir     = lambda wc: f"results/zone_aucell/{wc.zone_run}",
        env        = ZM_ENV,
    threads: 8
    resources:
        mem_mb = 32000, runtime = 60, cpus_per_task = 8,
    log:
        "logs/zone_modules/zone_aucell_{zone_run}.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --obj {input.h5ad} \
            --manifest {input.manifest} --object-name {params.obj_name} \
            --genesets {input.genesets} --auc-threshold {params.thr} --seed {params.seed} \
            --normalization {params.norm} --timepoints '{params.timepoints}' \
            --workers {threads} --out {params.outdir} > {log} 2>&1
        """


rule zone_plot_table:
    input:
        h5ad   = ancient(ZM["h5ad"]),
        scores = "results/zone_aucell/{zone_run}/zone_aucell.csv",
        prov   = "results/zone_aucell/{zone_run}/provenance.json",
        flags  = ZM["embedding_outliers"],
        script = "workflow/scripts/zone_modules/make_zone_plot_table.py",
    output:
        table = "results/zone_aucell/{zone_run}/zone_plot_table.csv",
        side  = "results/zone_aucell/{zone_run}/zone_plot_table.provenance.json",
    params:
        py  = ZM["python"],
        env = ZM_ENV,
    threads: 1
    resources:
        mem_mb = 8000, runtime = 30,
    log:
        "logs/zone_modules/zone_plot_table_{zone_run}.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --obj {input.h5ad} --scores {input.scores} \
            --scores-provenance {input.prov} --flag-barcodes {input.flags} --out {output.table} > {log} 2>&1
        """


rule zone_figures:
    input:
        table  = "results/zone_aucell/{zone_run}/zone_plot_table.csv",
        side   = "results/zone_aucell/{zone_run}/zone_plot_table.provenance.json",
        script = "workflow/scripts/zone_modules/fig_zone_modules.py",
        driver = "workflow/scripts/zone_modules/render_zone_figures.py",
    output:
        directory("results/figures/zone_modules/{zone_run}"),
    params:
        py   = ZM["python"],
        jobs = json.dumps(ZM["figures"]),
        env  = ZM_ENV,
    threads: 1
    resources:
        mem_mb = 16000, runtime = 60,
    log:
        "logs/zone_modules/zone_figures_{zone_run}.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.driver} --python {params.py} --script {input.script} \
            --plot-table {input.table} --jobs '{params.jobs}' --out {output} > {log} 2>&1
        """


rule zone_modules_all:
    input:
        expand("results/zone_aucell/{zone_run}/zone_aucell.csv", zone_run=ZM_RUNS),
        expand("results/zone_aucell/{zone_run}/zone_plot_table.csv", zone_run=ZM_RUNS),
        expand("results/figures/zone_modules/{zone_run}", zone_run=ZM_RUNS),
