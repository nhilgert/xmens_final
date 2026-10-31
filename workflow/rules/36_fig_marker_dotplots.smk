# workflow/rules/36_fig_marker_dotplots.smk
# Marker dot plots of release v9 on the genes of Table M1 v9 (Methods), and the Leiden review dot plots of
# every annotation step on the same genes (scripts in workflow/scripts/figures/panels/marker_dotplots/):
#   final-label dot plots, per species: the major classes; every fine lineage (fib/dec subtypes, epithelial,
#     endothelial, smooth muscle with perivascular cells); the immune subtypes twice, all v9 subtypes
#     (immune_v9) and the subtypes of the earlier immune plots (immune_as_before). Columns: each row's Table M1
#     genes in Table M1 order, as the Methods list them (marker_panels.plot_genes)
#   review dot plots, per species and annotation step of the review bundle (27_review_notebooks.smk):
#     rows = the step's Leiden clusters with their v9 decision; columns = for a lineage, the Table M1 rows of
#     its own classes, then the class-level genes of every other class; for the atlas, the class-level genes
#     of every class (marker_panels.review_columns); every step of a species also in one PDF for reading
#
#   panels_marker_dotplots_extract         per species: marker statistics of every v9 label and the Table M1
#                                          rule check, from the release object (sha256 checked against MANIFEST)
#   panels_marker_dotplots_panel           per species and plot: SVG, PNG (300 dpi), provenance JSON
#   panels_marker_dotplots_review_extract  per species and step: cluster table and marker statistics of every
#                                          Table M1 gene of the species; its provenance names the release folder
#   panels_marker_dotplots_review_panel    per species and step: SVG, PNG, provenance JSON
#   panels_marker_dotplots_review_pdf      per species: one PDF with a page per step, drawn from the review
#                                          tables by the same code as the SVGs (byte-reproducible)
#   panels_marker_dotplots_checks          Table M1 genes missing from an object's var or failing their rule
#
# Genes come from resources/figures/table_m1_markers_v9.csv (Table M1 v9). The release objects and the
# clusterings are ancient(): a figure never rebuilds a release or refits a lineage because an mtime moved.
# The review steps, clusterings, Leiden columns and decision tables come from the functions of
# 27_review_notebooks.smk, which the Snakefile includes earlier. Target: panels_marker_dotplots_all (not in
# rule all).

PMD = config.get("panels_marker_dotplots", {})
PMD_PY = PMD.get("python", "/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python")
PMD_CACHE = PMD.get("cache_dir", "/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch")
PMD_SCRIPTS = "workflow/scripts/figures/panels/marker_dotplots"
PMD_REVIEW_KIT = "workflow/scripts/review"
PMD_TABLES = "results/panels/marker_dotplots/v9"
PMD_FIGURES = "results/figures/panels/marker_dotplots/v9"
PMD_TABLE_M1 = "resources/figures/table_m1_markers_v9.csv"
PMD_VOCABULARY = "resources/annotation_vocabulary.yaml"

# Release v9 (vocabulary v9), pinned: config deliverables.release moves on with every release cut.
PMD_RELEASE = "data/deliverables/2026_10_30"
PMD_RELEASE_OBJECT = {"mouse": "Cevrim_XMens_All_Cells.h5ad", "human": "Cevrim_Human_All_Cells.h5ad"}

# Final-label plots of each species; marker_panels.PLOTS defines their rows and genes, and the extract rule
# fails if the two lists differ.
PMD_PLOTS = {
    "mouse": ["major_classes", "fib_dec_subtypes", "epithelial", "endothelial", "smooth_muscle_perivascular",
              "immune_v9", "immune_as_before"],
    "human": ["major_classes", "epithelial", "endothelial", "smooth_muscle_perivascular", "immune_v9",
              "immune_as_before"],
}
PMD_STEP_PATTERN = "|".join(["atlas"] + REVIEW_LINEAGE_ORDER)

# Matplotlib and numba caches default to HOME, which is full; HDF5 locking fails on lab_storage reads.
PMD_ENV = ("export XDG_CACHE_HOME={c} MPLCONFIGDIR={c}/mpl NUMBA_CACHE_DIR={c}/numba HDF5_USE_FILE_LOCKING=FALSE; "
           "mkdir -p {c}/mpl {c}/numba").format(c=PMD_CACHE)


def pmd_review_figure(species, step, suffix):
    """Path of one review dot plot output (suffix: .svg, .png or .provenance.json)."""
    return f"{PMD_FIGURES}/review/{species}_{step}_review_dotplot_v9{suffix}"


def pmd_review_tables(wildcards):
    """Review tables of every annotation step of one species, in page order."""
    return [
        f"{PMD_TABLES}/review/{wildcards.species}_{step}_review_{part}"
        for step in review_steps(wildcards.species)
        for part in ("clusters.csv", "marker_stats.csv", "extract_provenance.json")
    ]


def pmd_review_provenance():
    """Extract provenance of every annotation step of both species."""
    return [
        f"{PMD_TABLES}/review/{species}_{step}_review_extract_provenance.json"
        for species in config["lineages"]
        for step in review_steps(species)
    ]


# Marker statistics of every v9 label and the Table M1 rule check of one species: everything the final-label
# plots draw comes from these files.
rule panels_marker_dotplots_extract:
    input:
        h5ad = lambda wc: ancient(f"{PMD_RELEASE}/{PMD_RELEASE_OBJECT[wc.species]}"),
        manifest = ancient(f"{PMD_RELEASE}/MANIFEST.json"),
        table_m1 = PMD_TABLE_M1,
        script = f"{PMD_SCRIPTS}/extract_release_markers.py",
        panels = f"{PMD_SCRIPTS}/marker_panels.py",
        style = f"{PMD_SCRIPTS}/dotplot_style.py",
    output:
        stats = f"{PMD_TABLES}/{{species}}_marker_stats.csv",
        check = f"{PMD_TABLES}/{{species}}_table_m1_rule_check.csv",
        prov = f"{PMD_TABLES}/{{species}}_extract_provenance.json",
    params:
        py = PMD_PY,
        env = PMD_ENV,
        scripts = PMD_SCRIPTS,
        release = PMD_RELEASE,
        name = lambda wc: PMD_RELEASE_OBJECT[wc.species],
        plots = lambda wc: " ".join(PMD_PLOTS[wc.species]),
        outdir = PMD_TABLES,
    log:
        "logs/panels_marker_dotplots/extract_{species}.log",
    threads: 2
    resources:
        mem_mb = 32000, runtime = 30,
    shell:
        r"""
        {params.env}
        PYTHONPATH={params.scripts} {params.py} -u {input.script} \
            --species {wildcards.species} \
            --release-dir {params.release} \
            --object-name {params.name} \
            --table-m1 {input.table_m1} \
            --plots {params.plots} \
            --out-dir {params.outdir} > {log} 2>&1
        """


# One final-label dot plot of one species, drawn from that species' extract alone.
rule panels_marker_dotplots_panel:
    input:
        stats = f"{PMD_TABLES}/{{species}}_marker_stats.csv",
        check = f"{PMD_TABLES}/{{species}}_table_m1_rule_check.csv",
        prov = f"{PMD_TABLES}/{{species}}_extract_provenance.json",
        table_m1 = PMD_TABLE_M1,
        vocabulary = PMD_VOCABULARY,
        script = f"{PMD_SCRIPTS}/plot_marker_dotplot.py",
        panels = f"{PMD_SCRIPTS}/marker_panels.py",
        style = f"{PMD_SCRIPTS}/dotplot_style.py",
    output:
        multiext(f"{PMD_FIGURES}/{{species}}_{{pmd_plot}}_dotplot_v9", ".svg", ".png", ".provenance.json"),
    wildcard_constraints:
        pmd_plot = "|".join(sorted({p for plots in PMD_PLOTS.values() for p in plots})),
    params:
        py = PMD_PY,
        env = PMD_ENV,
        scripts = PMD_SCRIPTS,
        tables = PMD_TABLES,
        prefix = f"{PMD_FIGURES}/{{species}}_{{pmd_plot}}_dotplot_v9",
    log:
        "logs/panels_marker_dotplots/panel_{species}_{pmd_plot}.log",
    threads: 1
    resources:
        mem_mb = 4000, runtime = 10,
    shell:
        r"""
        {params.env}
        PYTHONPATH={params.scripts} {params.py} -u {input.script} \
            --species {wildcards.species} \
            --plot {wildcards.pmd_plot} \
            --tables-dir {params.tables} \
            --table-m1 {input.table_m1} \
            --vocabulary {input.vocabulary} \
            --out-prefix {params.prefix} > {log} 2>&1
        """


# Cluster table (decision, QC side columns) and marker statistics of one annotation step, on the clustering
# and decision table the review bundle names for it.
rule panels_marker_dotplots_review_extract:
    input:
        adata = lambda wc: ancient(review_object(wc)),
        decisions = review_decisions,
        v8_table = review_v8_decisions,
        reasons = "resources/review/exclusion_reasons_v9.csv",
        steps = "resources/review/decision_steps_v9.csv",
        table_m1 = PMD_TABLE_M1,
        manifest = ancient(f"{PMD_RELEASE}/MANIFEST.json"),
        script = f"{PMD_SCRIPTS}/extract_review_step.py",
        panels = f"{PMD_SCRIPTS}/marker_panels.py",
        style = f"{PMD_SCRIPTS}/dotplot_style.py",
        kit = [f"{PMD_REVIEW_KIT}/annotation_review.py", f"{PMD_REVIEW_KIT}/review_ledger.py"],
    output:
        clusters = f"{PMD_TABLES}/review/{{species}}_{{step}}_review_clusters.csv",
        markers = f"{PMD_TABLES}/review/{{species}}_{{step}}_review_marker_stats.csv",
        prov = f"{PMD_TABLES}/review/{{species}}_{{step}}_review_extract_provenance.json",
    wildcard_constraints:
        step = PMD_STEP_PATTERN,
    params:
        py = PMD_PY,
        env = PMD_ENV,
        path = f"{PMD_SCRIPTS}:{PMD_REVIEW_KIT}",
        cluster_column = review_cluster_column,
        v8_decisions = lambda wc: review_v8_decisions(wc) or "None",
        # Where each nucleus was in v8, for the round-2 re-fits; a parameter, as in 99b, so this target never
        # depends on the release rules.
        v8_release = lambda wc: f"{V8_RELEASE_DIR}/{V8_RELEASE_OBJECT[wc.species]}",
        # The release whose labels these decisions set (recorded in the provenance).
        release = PMD_RELEASE,
        outdir = f"{PMD_TABLES}/review",
    log:
        "logs/panels_marker_dotplots/review_extract_{species}_{step}.log",
    threads: 2
    resources:
        mem_mb = lambda wc: 48000 if wc.step in ("atlas", "stromal") else 16000,
        runtime = 45,
    shell:
        r"""
        {params.env}
        PYTHONPATH={params.path} {params.py} -u {input.script} \
            --species {wildcards.species} \
            --step {wildcards.step} \
            --object {input.adata} \
            --cluster-column {params.cluster_column} \
            --decisions {input.decisions} \
            --v8-decisions {params.v8_decisions} \
            --v8-release {params.v8_release} \
            --reasons {input.reasons} \
            --steps {input.steps} \
            --table-m1 {input.table_m1} \
            --release-dir {params.release} \
            --out-dir {params.outdir} > {log} 2>&1
        """


# The review dot plot of one step.
rule panels_marker_dotplots_review_panel:
    input:
        clusters = f"{PMD_TABLES}/review/{{species}}_{{step}}_review_clusters.csv",
        markers = f"{PMD_TABLES}/review/{{species}}_{{step}}_review_marker_stats.csv",
        prov = f"{PMD_TABLES}/review/{{species}}_{{step}}_review_extract_provenance.json",
        table_m1 = PMD_TABLE_M1,
        vocabulary = PMD_VOCABULARY,
        script = f"{PMD_SCRIPTS}/plot_review_step.py",
        figure = f"{PMD_SCRIPTS}/review_figure.py",
        panels = f"{PMD_SCRIPTS}/marker_panels.py",
        style = f"{PMD_SCRIPTS}/dotplot_style.py",
        kit = [f"{PMD_REVIEW_KIT}/review_dotplot.py", f"{PMD_REVIEW_KIT}/review_step.py"],
    output:
        multiext(f"{PMD_FIGURES}/review/{{species}}_{{step}}_review_dotplot_v9", ".svg", ".png", ".provenance.json"),
    wildcard_constraints:
        step = PMD_STEP_PATTERN,
    params:
        py = PMD_PY,
        env = PMD_ENV,
        path = f"{PMD_SCRIPTS}:{PMD_REVIEW_KIT}",
        tables = f"{PMD_TABLES}/review",
        prefix = f"{PMD_FIGURES}/review/{{species}}_{{step}}_review_dotplot_v9",
    log:
        "logs/panels_marker_dotplots/review_panel_{species}_{step}.log",
    threads: 1
    resources:
        mem_mb = 4000, runtime = 10,
    shell:
        r"""
        {params.env}
        PYTHONPATH={params.path} {params.py} -u {input.script} \
            --species {wildcards.species} \
            --step {wildcards.step} \
            --tables-dir {params.tables} \
            --table-m1 {input.table_m1} \
            --vocabulary {input.vocabulary} \
            --out-prefix {params.prefix} > {log} 2>&1
        """


# One PDF per species with one page per step (atlas, then the lineages in review order), for reading.
rule panels_marker_dotplots_review_pdf:
    input:
        tables = pmd_review_tables,
        table_m1 = PMD_TABLE_M1,
        vocabulary = PMD_VOCABULARY,
        script = f"{PMD_SCRIPTS}/plot_review_pdf.py",
        pages = f"{PMD_SCRIPTS}/plot_review_step.py",
        figure = f"{PMD_SCRIPTS}/review_figure.py",
        panels = f"{PMD_SCRIPTS}/marker_panels.py",
        style = f"{PMD_SCRIPTS}/dotplot_style.py",
        kit = [f"{PMD_REVIEW_KIT}/review_dotplot.py", f"{PMD_REVIEW_KIT}/review_step.py"],
    output:
        f"{PMD_FIGURES}/review/review_dotplots_{{species}}_v9.pdf",
    params:
        py = PMD_PY,
        env = PMD_ENV,
        path = f"{PMD_SCRIPTS}:{PMD_REVIEW_KIT}",
        steps = lambda wc: " ".join(review_steps(wc.species)),
        tables = f"{PMD_TABLES}/review",
    log:
        "logs/panels_marker_dotplots/review_pdf_{species}.log",
    threads: 1
    resources:
        mem_mb = 4000, runtime = 15,
    shell:
        r"""
        {params.env}
        PYTHONPATH={params.path} {params.py} -u {input.script} \
            --species {wildcards.species} \
            --steps {params.steps} \
            --tables-dir {params.tables} \
            --table-m1 {input.table_m1} \
            --vocabulary {input.vocabulary} \
            --out {output} > {log} 2>&1
        """


# Table M1 genes missing from the var of any object these plots read, and genes failing their rule on v9.
rule panels_marker_dotplots_checks:
    input:
        checks = expand(f"{PMD_TABLES}/{{species}}_table_m1_rule_check.csv", species=list(PMD_PLOTS)),
        release = expand(f"{PMD_TABLES}/{{species}}_extract_provenance.json", species=list(PMD_PLOTS)),
        review = pmd_review_provenance(),
        script = f"{PMD_SCRIPTS}/summarize_marker_checks.py",
    output:
        csv = f"{PMD_TABLES}/table_m1_marker_checks_v9.csv",
        md = f"{PMD_TABLES}/table_m1_marker_checks_v9.md",
    params:
        py = PMD_PY,
        env = PMD_ENV,
    log:
        "logs/panels_marker_dotplots/checks.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --rule-checks {input.checks} \
            --release-provenance {input.release} \
            --review-provenance {input.review} \
            --out-csv {output.csv} \
            --out-md {output.md} > {log} 2>&1
        """


# Every final-label plot, every review plot, the two reading PDFs and the gene checks.
rule panels_marker_dotplots_all:
    input:
        [f"{PMD_FIGURES}/{sp}_{plot}_dotplot_v9.provenance.json" for sp in PMD_PLOTS for plot in PMD_PLOTS[sp]],
        [pmd_review_figure(sp, step, ".provenance.json") for sp in config["lineages"] for step in review_steps(sp)],
        expand(f"{PMD_FIGURES}/review/review_dotplots_{{species}}_v9.pdf", species=list(config["lineages"])),
        rules.panels_marker_dotplots_checks.output,
