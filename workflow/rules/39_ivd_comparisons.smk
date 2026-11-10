# workflow/rules/39_ivd_comparisons.smk
# =============================================================================
# Comparison with in vitro and in vivo decidualization models (Methods, "Comparison with in vitro and in vivo
# decidualization models"), ported 2026-11-02 from paper_code/ivd_comparisons (docs/PORTING.md item 3; the Stadtmauer
# consensus bundle of item 5). Target: ivd_all. Not in the default target.
#
#   ivd_consensus          the 123-gene human in vitro consensus (vote over Aghajanova 2011, Rytkonen 2019, Stadtmauer &
#                          Wagner 2022 with the authors' criteria of 22 Oct 2026) and its Ensembl 116 one-to-one mouse
#                          orthology; regression checks against the released list and orthology table
#   ivd_mouse_pseudobulk   X-Mens decidual / fibroblast pseudobulks per timepoint, with the uninduced fibroblasts (27 Oct)
#   ivd_mouse_fit          PyDESeq2 ~ timepoint + group -> xmens_broad_all_genes.csv
#   ivd_xmens_compare      X-Mens vs the consensus (through the orthology) and vs the three Wang et al. 2020 lists
#   ivd_consensus_tables   Supplementary Tables S9 and S10 (S10 with the X-Mens columns of the fit above)
#   ivd_human_pseudobulk   human fib/dec pseudobulks per donor and section (2026_09_15 object)
#   ivd_human_fit          PyDESeq2 ~ donor + group (and the sensitivity fits) -> human_decidual_vs_fibroblast_all_genes.csv
#   ivd_ulrich_download    the Ulrich et al. 2024 stromal atlas from CELLxGENE (sha256-pinned; copies a local mirror if present)
#   ivd_ulrich_config      the Ulrich DE config for this checkout, from the tracked template
#   ivd_ulrich_de          PyDESeq2 ~ donor + state for the configured contrasts (dS1-3 vs C7-Fibroblast and the others),
#                          with the acceptance counts of the 28 Oct 2026 release pinned in the config
#   ivd_workbooks          pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx and IVD_comparisons_with_genes.xlsx
#   ivd_verify             every output against resources/ivd/expected and the workbooks of record -> VERIFICATION.md
#
# The X-Mens object and the human object are inputs of record (ancient()): the fits were made on data/interim/final_mouse.h5ad
# and on the 2026_09_15 fib/dec release, not on the v9 release objects. Outputs go to results/ivd/<run>/.
# Environments: the fits run in the frozen PyDESeq2 environment (config ivd.python_fit; workflow/envs/ivd), the table and
# workbook steps in rsc_uv (ivd.python). PyDESeq2 p-values depend on the NumPy build (workflow/envs/ivd/README.md).
# =============================================================================

IVD = config["ivd"]
IVD_RUN = str(IVD["run"])
IVD_DIR = f"results/ivd/{IVD_RUN}"
IVD_LOGS = f"logs/ivd/{IVD_RUN}"
IVD_PY = IVD["python"]
IVD_PY_FIT = IVD["python_fit"]
IVD_S = "workflow/scripts/ivd"
IVD_R = "resources/ivd"
IVD_ENV = (f"export HDF5_USE_FILE_LOCKING=FALSE PYTHONDONTWRITEBYTECODE=1 XDG_CACHE_HOME={IVD['cache_dir']} "
           f"MPLCONFIGDIR={IVD['cache_dir']}/mpl NUMBA_CACHE_DIR={IVD['cache_dir']}/numba; mkdir -p {IVD['cache_dir']}/mpl {IVD['cache_dir']}/numba")
IVD_ULRICH = IVD["ulrich"]
IVD_ULRICH_H5AD = IVD_ULRICH["h5ad"]
IVD_ULRICH_RELEASE = f"{IVD_DIR}/ulrich/release"
IVD_ULRICH_CID = IVD_ULRICH["contrast"]

IVD_SIGNATURES = [f"{IVD_R}/consensus/{fn}" for fn in (
    "aghajanova2011_invitro_decidualization_nonendo.csv", "rytkonen2019_complete_signature.csv",
    "stadtmauer2022_hits.csv", "lucas2020_S2vsS3_decidualization_onset.csv")]
IVD_WANG = [f"{IVD_R}/wang/wang_{k}.csv" for k in ("ivd", "oil", "pregnancy")]
IVD_HUMAN_FITS = ["A_new_celltype_old_settings", "B_new_celltype_release", "C_new_coarse_release",
                  "D_new_celltype_release_by_section", "E_new_coarse_release_by_section"]


rule ivd_consensus:
    input:
        signatures = IVD_SIGNATURES,
        table      = f"{IVD_R}/consensus/stadtmauer_tables_S2_S3.csv",
        overlay    = f"{IVD_R}/consensus/stadtmauer_manual_overlay.csv",
        released   = f"{IVD_R}/consensus/released_consensus123.csv",
        released_orth = f"{IVD_R}/consensus/released_consensus123_mouse_orthology.csv",
        heca       = f"{IVD_R}/consensus/consensus123_heca_mapping.csv",
        components = "resources/orthology/ensembl_human_mouse_components.tsv",
        oma        = f"{IVD_R}/orthology/oma_human_mouse_ancestral_hogs.tsv",
        annotation = f"{IVD_R}/orthology/ensembl_116_homo_sapiens_gene_annotation.tsv.gz",
        alias      = f"{IVD_R}/orthology/aghajanova_alias_to_current_audit.csv",
        script     = f"{IVD_S}/consensus/build_consensus123.py",
        modules    = [f"{IVD_S}/consensus/{m}" for m in ("build_consensus.py", "rebuild_consensus123_stadtmauer_thresholds.py", "derive_consensus123_orthology.py")],
    output:
        consensus = f"{IVD_DIR}/consensus/consensus123_stadtmauer_authors_criteria.csv",
        orthology = f"{IVD_DIR}/consensus/consensus123_stadtmauer_authors_criteria_mouse_orthology.csv",
        hits      = f"{IVD_DIR}/consensus/stadtmauer2022_hits_authors_criteria.csv",
        table     = f"{IVD_DIR}/consensus/stadtmauer_tables_S2_S3_parsed_and_called.csv",
        manifest  = f"{IVD_DIR}/consensus/manifest.json",
    params:
        py = IVD_PY, env = IVD_ENV, outdir = f"{IVD_DIR}/consensus", sigdir = f"{IVD_R}/consensus",
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        f"{IVD_LOGS}/consensus.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --signatures {params.sigdir} --stadtmauer-table {input.table} --overlay {input.overlay} \
            --released-consensus {input.released} --released-orthology {input.released_orth} --heca-mapping {input.heca} \
            --components {input.components} --oma {input.oma} --annotation {input.annotation} --alias {input.alias} \
            --out {params.outdir} > {log} 2>&1
        """


rule ivd_mouse_pseudobulk:
    input:
        h5ad     = ancient(IVD["mouse_object"]),
        manifest = f"{IVD_R}/manifests/xmens_feature_manifest.tsv",
        script   = f"{IVD_S}/xmens/mouse_pseudobulk_job.py",
    output:
        counts   = f"{IVD_DIR}/xmens/pb/xmens_broad_carrier_counts.tsv",
        metadata = f"{IVD_DIR}/xmens/pb/xmens_broad_metadata.tsv",
        genes    = f"{IVD_DIR}/xmens/pb/xmens_broad_genes.tsv",
    params:
        py = IVD_PY_FIT, env = IVD_ENV, outdir = f"{IVD_DIR}/xmens/pb",
    threads: 8
    resources:
        mem_mb = 64000, runtime = 90, cpus_per_task = 8,
    log:
        f"{IVD_LOGS}/mouse_pseudobulk.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --object {input.h5ad} --manifest {input.manifest} --design with_uninduced \
            --out {params.outdir} > {log} 2>&1
        """


rule ivd_mouse_fit:
    input:
        counts   = rules.ivd_mouse_pseudobulk.output.counts,
        metadata = rules.ivd_mouse_pseudobulk.output.metadata,
        genes    = rules.ivd_mouse_pseudobulk.output.genes,
        script   = f"{IVD_S}/xmens/fit_mouse_de.py",
    output:
        table = f"{IVD_DIR}/xmens/fit/xmens_broad_all_genes.csv",
        diag  = f"{IVD_DIR}/xmens/fit/fit_diagnostics.json",
    params:
        py = IVD_PY_FIT, env = IVD_ENV, pbdir = f"{IVD_DIR}/xmens/pb", outdir = f"{IVD_DIR}/xmens/fit",
    threads: 8
    resources:
        mem_mb = 16000, runtime = 60, cpus_per_task = 8,
    log:
        f"{IVD_LOGS}/mouse_fit.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --pb-dir {params.pbdir} --out {params.outdir} --n-cpus {threads} > {log} 2>&1
        """


rule ivd_xmens_compare:
    input:
        xmens     = rules.ivd_mouse_fit.output.table,
        consensus = rules.ivd_consensus.output.consensus,
        orthology = rules.ivd_consensus.output.orthology,
        wang      = IVD_WANG,
        script    = f"{IVD_S}/xmens/compare_xmens_with_references.py",
        module    = f"{IVD_S}/xmens/release_comparison_functions.py",
    output:
        support = f"{IVD_DIR}/xmens/compare/consensus_ortholog_xmens_support_rows.csv",
        wang    = f"{IVD_DIR}/xmens/compare/wang_ivd_support_rows.csv",
        genes   = expand(f"{IVD_DIR}/xmens/compare/genes_xmens_vs_{{ref}}.csv", ref=["consensus", "oil", "pregnancy", "in_vitro"]),
        summary = f"{IVD_DIR}/xmens/compare/SUMMARY_xmens_comparisons.csv",
    params:
        py = IVD_PY, env = IVD_ENV, outdir = f"{IVD_DIR}/xmens/compare", wang_dir = f"{IVD_R}/wang",
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        f"{IVD_LOGS}/xmens_compare.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --xmens {input.xmens} --consensus {input.consensus} --orthology {input.orthology} \
            --wang-dir {params.wang_dir} --out {params.outdir} > {log} 2>&1
        """


rule ivd_consensus_tables:
    input:
        manifest = rules.ivd_consensus.output.manifest,
        support  = rules.ivd_xmens_compare.output.support,
        script   = f"{IVD_S}/consensus/make_consensus_tables.py",
        module   = f"{IVD_S}/consensus/build_outputs_stadtmauer_fix.py",
    output:
        s9  = f"{IVD_DIR}/consensus/tables/TABLE_consensus123_sources_by_study.csv",
        s10 = f"{IVD_DIR}/consensus/tables/TABLE_consensus123_per_gene_calls.csv",
    params:
        py = IVD_PY, env = IVD_ENV, consensus_dir = f"{IVD_DIR}/consensus", outdir = f"{IVD_DIR}/consensus/tables",
    threads: 1
    resources:
        mem_mb = 4000, runtime = 10,
    log:
        f"{IVD_LOGS}/consensus_tables.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --consensus-dir {params.consensus_dir} --xmens-support {input.support} --out {params.outdir} > {log} 2>&1
        """


rule ivd_human_pseudobulk:
    input:
        h5ad   = ancient(IVD["human_object"]),
        script = f"{IVD_S}/human/human_pseudobulk_job.py",
    output:
        meta   = f"{IVD_DIR}/human/pb/human_pb_meta.csv",
        counts = f"{IVD_DIR}/human/pb/human_pb_counts.csv.gz",
        detect = f"{IVD_DIR}/human/pb/human_pb_detect.csv.gz",
        var    = f"{IVD_DIR}/human/pb/human_var.csv",
    params:
        py = IVD_PY_FIT, env = IVD_ENV, outdir = f"{IVD_DIR}/human/pb",
    threads: 4
    resources:
        mem_mb = 48000, runtime = 60, cpus_per_task = 4,
    log:
        f"{IVD_LOGS}/human_pseudobulk.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --object {input.h5ad} --out {params.outdir} > {log} 2>&1
        """


rule ivd_human_fit:
    input:
        meta   = rules.ivd_human_pseudobulk.output.meta,
        counts = rules.ivd_human_pseudobulk.output.counts,
        detect = rules.ivd_human_pseudobulk.output.detect,
        var    = rules.ivd_human_pseudobulk.output.var,
        script = f"{IVD_S}/human/fit_human_de.py",
        tables = f"{IVD_S}/human/build_human_de_tables.py",
    output:
        fits  = expand(f"{IVD_DIR}/human/fits/DE_{{fit}}.csv.gz", fit=IVD_HUMAN_FITS),
        diag  = f"{IVD_DIR}/human/fits/fit_diagnostics.json",
        table = f"{IVD_DIR}/human/tables/human_decidual_vs_fibroblast_all_genes.csv",
        coarse = f"{IVD_DIR}/human/tables/human_decidual_vs_fibroblast_all_genes_coarse_refined_labels.csv",
    params:
        py = IVD_PY_FIT, env = IVD_ENV, pbdir = f"{IVD_DIR}/human/pb", fitdir = f"{IVD_DIR}/human/fits", outdir = f"{IVD_DIR}/human/tables",
    threads: 8
    resources:
        mem_mb = 16000, runtime = 60, cpus_per_task = 8,
    log:
        f"{IVD_LOGS}/human_fit.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --pb-dir {params.pbdir} --out {params.fitdir} --n-cpus {threads} > {log} 2>&1
        {params.py} -u {input.tables} --pb-dir {params.pbdir} --fits {params.fitdir} --out {params.outdir} >> {log} 2>&1
        """


rule ivd_ulrich_download:
    input:
        script = f"{IVD_S}/ulrich/download_ulrich_atlas.py",
    output:
        h5ad = IVD_ULRICH_H5AD,
    params:
        py = IVD_PY, env = IVD_ENV, url = IVD_ULRICH["url"], sha = IVD_ULRICH["sha256"],
        mirror = IVD_ULRICH.get("mirror", ""),
    threads: 1
    resources:
        mem_mb = 4000, runtime = 180,
    log:
        f"{IVD_LOGS}/ulrich_download.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --url {params.url} --sha256 {params.sha} --out {output.h5ad} \
            $( [ -n "{params.mirror}" ] && echo --mirror {params.mirror} ) > {log} 2>&1
        """


rule ivd_ulrich_config:
    input:
        template = f"{IVD_S}/ulrich/config_ulrich_dS_vs_C7fib.template.yaml",
        script   = f"{IVD_S}/ulrich/render_ulrich_config.py",
        h5ad     = IVD_ULRICH_H5AD,
    output:
        config = f"{IVD_DIR}/ulrich/config_ulrich_dS_vs_C7fib.yaml",
    params:
        py = IVD_PY, upstream = f"{IVD_S}/ulrich/upstream", py_fit = IVD_PY_FIT, release = IVD_ULRICH_RELEASE,
    threads: 1
    resources:
        mem_mb = 2000, runtime = 5,
    log:
        f"{IVD_LOGS}/ulrich_config.log",
    shell:
        r"""
        {params.py} -u {input.script} --template {input.template} --root . --ulrich-h5ad {input.h5ad} \
            --upstream {params.upstream} --python {params.py_fit} --output {params.release} --out {output.config} > {log} 2>&1
        """


rule ivd_ulrich_de:
    input:
        config   = rules.ivd_ulrich_config.output.config,
        h5ad     = IVD_ULRICH_H5AD,
        script   = f"{IVD_S}/ulrich/ulrich_human_de.py",
        pipeline = f"{IVD_S}/ulrich/ulrich_ivd_pipeline.py",
        upstream = f"{IVD_S}/ulrich/upstream/run_pipeline.py",
        resources_ = [f"{IVD_R}/consensus/released_consensus123.csv", f"{IVD_R}/manifests/consensus123_ulrich_mapping.csv",
                      f"{IVD_R}/manifests/ulrich_feature_manifest.tsv"],
    output:
        contrasts = f"{IVD_ULRICH_RELEASE}/ulrich_contrasts.csv",
        table     = f"{IVD_ULRICH_RELEASE}/ulrich_{IVD_ULRICH_CID}_all_genes.csv",
        donors    = f"{IVD_ULRICH_RELEASE}/intermediate/ulrich_{IVD_ULRICH_CID}_metadata.tsv",
        run_meta  = f"{IVD_ULRICH_RELEASE}/run_metadata.json",
        manifest  = f"{IVD_ULRICH_RELEASE}/release_manifest.json",
    params:
        py = IVD_PY_FIT, env = IVD_ENV, release = IVD_ULRICH_RELEASE,
    threads: 8
    resources:
        mem_mb = 48000, runtime = 240, cpus_per_task = 8,
    log:
        f"{IVD_LOGS}/ulrich_de.log",
    shell:
        r"""
        {params.env}
        mkdir -p $(dirname {params.release}); rm -rf {params.release} {params.release}.inprogress
        {params.py} -u {input.script} --config {input.config} --out {params.release} --n-cpus {threads} > {log} 2>&1
        """


rule ivd_workbooks:
    input:
        xmens     = rules.ivd_mouse_fit.output.table,
        slidetags = rules.ivd_human_fit.output.table,
        contrasts = rules.ivd_ulrich_de.output.contrasts,
        donors    = rules.ivd_ulrich_de.output.donors,
        config    = rules.ivd_ulrich_config.output.config,
        s10       = rules.ivd_consensus_tables.output.s10,
        compare   = rules.ivd_xmens_compare.output.summary,
        consensus = rules.ivd_consensus.output.consensus,
        wang      = IVD_WANG,
        de_script = f"{IVD_S}/workbooks/build_de_workbook.py",
        ivd_script = f"{IVD_S}/workbooks/build_ivd_workbook.py",
        modules   = [f"{IVD_S}/workbooks/build_house_style_workbooks.py", f"{IVD_S}/workbooks/add_human_comparisons_to_ivd_workbook.py"],
    output:
        de  = f"{IVD_DIR}/workbooks/pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx",
        ivd = f"{IVD_DIR}/workbooks/IVD_comparisons_with_genes.xlsx",
    params:
        py = IVD_PY, env = IVD_ENV, release = IVD_ULRICH_RELEASE, release_name = IVD_ULRICH.get("release_name", f"results/ivd/{IVD_RUN}/ulrich/release"),
        cid = IVD_ULRICH_CID, label = IVD_ULRICH["label"], house = f"{IVD_S}/workbooks/build_house_style_workbooks.py",
        consensus_dir = f"{IVD_DIR}/consensus", tables = f"{IVD_DIR}/consensus/tables", compare_dir = f"{IVD_DIR}/xmens/compare", wang_dir = f"{IVD_R}/wang",
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        f"{IVD_LOGS}/workbooks.log",
    shell:
        r"""
        {params.env}
        mkdir -p $(dirname {output.de}); rm -f {output.de} {output.ivd}
        {params.py} -u {input.de_script} --xmens {input.xmens} --slidetags {input.slidetags} \
            --hammoud-release {params.release} --hammoud-release-name {params.release_name} --hammoud-config {input.config} \
            --hammoud-donors {input.donors} --hammoud-contrast {params.cid} --hammoud-label "{params.label}" \
            --house-style {params.house} --out {output.de} > {log} 2>&1
        {params.py} -u {input.ivd_script} --consensus-dir {params.consensus_dir} --consensus-tables {params.tables} \
            --xmens-compare {params.compare_dir} --wang-dir {params.wang_dir} --slidetags {input.slidetags} \
            --hammoud-release {params.release} --hammoud-release-name {params.release_name} --hammoud-config {input.config} \
            --hammoud-donors {input.donors} --hammoud-contrast {params.cid} --hammoud-label "{params.label}" \
            --out {output.ivd} >> {log} 2>&1
        """


rule ivd_verify:
    input:
        rules.ivd_workbooks.output,
        rules.ivd_consensus_tables.output,
        rules.ivd_xmens_compare.output,
        rules.ivd_human_fit.output.table,
        rules.ivd_ulrich_de.output.table,
        script = f"{IVD_S}/verify_ivd_outputs.py",
        expected = f"{IVD_R}/expected/consensus/TABLE_consensus123_per_gene_calls.csv",
    output:
        report = f"{IVD_DIR}/verification/VERIFICATION.md",
        json   = f"{IVD_DIR}/verification/verification.json",
    params:
        py = IVD_PY, env = IVD_ENV, results = IVD_DIR, expected = f"{IVD_R}/expected", cid = IVD_ULRICH_CID,
        record = IVD.get("record_workbooks", ""), outdir = f"{IVD_DIR}/verification",
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        f"{IVD_LOGS}/verify.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} --results {params.results} --expected {params.expected} --ulrich-contrast {params.cid} \
            $( [ -n "{params.record}" ] && echo --record-workbooks {params.record} ) --out {params.outdir} > {log} 2>&1
        """


# The whole comparison: consensus, the three fits, the comparison tables, S9/S10, both workbooks and the verification.
rule ivd_all:
    input:
        rules.ivd_verify.output.report,
        rules.ivd_workbooks.output,
