# workflow/rules/25_zone_modules.smk
# =============================================================================
# Decidual-zone modules (mouse), added 2026-10-24. Per-nucleus AUCell scores for
# the 11 uterine microenvironments of Li et al. 2022 (Biol Reprod 107:529,
# Table 1, top 20 markers each; resources/signatures/aucell_modules_mouse.yaml),
# and their spatial panels and subtype heatmaps.
#
#   zone_aucell      score one run: pySCENIC AUCell on raw counts of the release
#                    object, fixed auc_threshold and seed, only the nuclei of the
#                    run's timepoints (exact; see the script). Genes ranked: all
#                    (default) or, with gene_filter_frac, SCENIC's rule (count > 0
#                    in >= 1% of the object's nuclei) plus, with
#                    keep_geneset_genes, the gene-set genes below it that are
#                    still detected (added 2026-10-26).
#   zone_plot_table  one per-nucleus CSV with everything the figures need:
#                    obs (sample, timepoint, cell_type, cell_subtype, x/y,
#                    embedding_outlier as the release object stores it)
#                    and the 11 scores. Added 2026-10-24 so the figures can be
#                    redrawn from this CSV alone.
#   zone_figures     render config[zone_modules][figures] from the plot table
#                    only. Spatial panels hide embedding_outlier nuclei; heatmaps
#                    keep them (rule of 2026-10-22). Spatial panels are turned by
#                    resources/figures/mouse_section_rotation.yaml, the table every
#                    mouse spatial panel reads (added 2026-10-30; was --rotate 135
#                    in the figure arguments).
#
# The release scored is the release of record (config deliverables.version and
# releases, 91_panels_versions.smk): one directory per release tag and run under
# results/zone_aucell/<tag>/<run>/ and results/figures/zone_modules/<tag>/<run>/,
# so a release cut adds a scoring beside the earlier ones. The v9 runs written
# before version tags were used here are results/zone_aucell/<run>/ (read by the
# v9 deposit) and stay as they are. Add a run instead of editing one.
# Target: zone_modules_all.
#
# The release object and its MANIFEST are ancient(): a new release is cut
# deliberately (04g), never re-cut because an mtime moved.
# =============================================================================
import json

include: "91_panels_versions.smk"

ZM      = config["zone_modules"]
ZM_RUNS = ZM["runs"]
# Release of record: its whole-tissue object and MANIFEST, and the folders of its scores and figures.
ZM_TAG      = PANEL_LIVE_TAG
ZM_H5AD     = f"{panel_release_dir(ZM_TAG)}/{ZM['object_name']}.h5ad"
ZM_MANIFEST = f"{panel_release_dir(ZM_TAG)}/MANIFEST.json"
ZM_SCORES   = f"results/zone_aucell/{ZM_TAG}/{{zone_run}}"
ZM_FIGURES  = f"results/figures/zone_modules/{ZM_TAG}/{{zone_run}}"
# Display rotation of the mouse sections: one table read by every mouse spatial panel, and its reader.
ZM_ROTATION_TABLE  = "resources/figures/mouse_section_rotation.yaml"
ZM_ROTATION_READER = "workflow/scripts/figures/section_rotation.py"
# numba (pyscenic) and matplotlib caches must not land in the full HOME
ZM_ENV  = ("export XDG_CACHE_HOME={c} NUMBA_CACHE_DIR={c}/numba MPLCONFIGDIR={c}/mpl; "
           "mkdir -p {c}/numba {c}/mpl").format(c=PANEL_CACHE)


def zm_gene_args(run):
    """score_zone_aucell.py options for the genes AUCell ranks (all genes when gene_filter_frac is absent or 0)."""
    frac = float(run.get("gene_filter_frac") or 0)
    if frac <= 0:
        assert not run.get("keep_geneset_genes"), "zone_modules: keep_geneset_genes needs gene_filter_frac"
        return ""
    args = f"--gene-filter-frac {frac}"
    if run.get("keep_geneset_genes"):
        args += f" --keep-geneset-genes --geneset-min-cells {int(run.get('geneset_min_cells', 1))}"
    return args

wildcard_constraints:
    zone_run = "|".join(ZM_RUNS),


rule zone_aucell:
    input:
        h5ad     = ancient(ZM_H5AD),
        manifest = ancient(ZM_MANIFEST),
        genesets = lambda wc: ZM_RUNS[wc.zone_run]["genesets"],
        script   = "workflow/scripts/zone_modules/score_zone_aucell.py",
    output:
        scores = f"{ZM_SCORES}/zone_aucell.csv",
        means  = f"{ZM_SCORES}/subtype_means.csv",
        genes  = f"{ZM_SCORES}/genes_ranked.txt",
        prov   = f"{ZM_SCORES}/provenance.json",
    params:
        py         = PANEL_PY,
        obj_name   = ZM["object_name"],
        thr        = lambda wc: ZM_RUNS[wc.zone_run]["auc_threshold"],
        seed       = lambda wc: ZM_RUNS[wc.zone_run]["seed"],
        norm       = lambda wc: ZM_RUNS[wc.zone_run].get("normalization", "none"),
        timepoints = lambda wc: ",".join(ZM_RUNS[wc.zone_run].get("timepoints") or []),
        genes      = lambda wc: zm_gene_args(ZM_RUNS[wc.zone_run]),
        outdir     = lambda wc: ZM_SCORES.format(zone_run=wc.zone_run),
        env        = ZM_ENV,
    threads: 8
    resources:
        mem_mb = 32000, runtime = 60, cpus_per_task = 8,
    log:
        f"logs/zone_modules/zone_aucell_{ZM_TAG}_{{zone_run}}.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --obj {input.h5ad} \
            --manifest {input.manifest} --object-name {params.obj_name} \
            --genesets {input.genesets} --auc-threshold {params.thr} --seed {params.seed} \
            --normalization {params.norm} --timepoints '{params.timepoints}' {params.genes} \
            --workers {threads} --out {params.outdir} > {log} 2>&1
        """


rule zone_plot_table:
    input:
        h5ad   = ancient(ZM_H5AD),
        scores = f"{ZM_SCORES}/zone_aucell.csv",
        prov   = f"{ZM_SCORES}/provenance.json",
        script = "workflow/scripts/zone_modules/make_zone_plot_table.py",
    output:
        table = f"{ZM_SCORES}/zone_plot_table.csv",
        side  = f"{ZM_SCORES}/zone_plot_table.provenance.json",
    params:
        py  = PANEL_PY,
        env = ZM_ENV,
    threads: 1
    resources:
        mem_mb = 8000, runtime = 30,
    log:
        f"logs/zone_modules/zone_plot_table_{ZM_TAG}_{{zone_run}}.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --obj {input.h5ad} --scores {input.scores} \
            --scores-provenance {input.prov} --out {output.table} > {log} 2>&1
        """


rule zone_figures:
    input:
        table  = f"{ZM_SCORES}/zone_plot_table.csv",
        side   = f"{ZM_SCORES}/zone_plot_table.provenance.json",
        script = "workflow/scripts/zone_modules/fig_zone_modules.py",
        driver = "workflow/scripts/zone_modules/render_zone_figures.py",
        rotation = ZM_ROTATION_TABLE,
        rotation_reader = ZM_ROTATION_READER,
    output:
        directory(ZM_FIGURES),
    params:
        py   = PANEL_PY,
        jobs = json.dumps(ZM["figures"]),
        env  = ZM_ENV,
    threads: 1
    resources:
        mem_mb = 16000, runtime = 60,
    log:
        f"logs/zone_modules/zone_figures_{ZM_TAG}_{{zone_run}}.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.driver} --python {params.py} --script {input.script} \
            --plot-table {input.table} --jobs '{params.jobs}' --section-rotation {input.rotation} \
            --out {output} > {log} 2>&1
        """


# Every run on the release of record: scores, plot table and figures.
rule zone_modules_all:
    input:
        expand(f"{ZM_SCORES}/zone_aucell.csv", zone_run=ZM_RUNS),
        expand(f"{ZM_SCORES}/zone_plot_table.csv", zone_run=ZM_RUNS),
        expand(ZM_FIGURES, zone_run=ZM_RUNS),
