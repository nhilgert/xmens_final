# workflow/rules/41_fig_c7_presence.smk
# The C7-fibroblast figure of the revision (reviewer comment 1): Ulrich et al. 2024's Stro-7 "C7-fibroblasts" are
# human Leiden cluster 3 (resolution 0.5) and, through the SAMap graph, mouse outer fibroblasts (eFib1).
#
#   final objects (config/c7_presence.yaml objects; verified by sha256)
#     --(c7_leiden_labels)-->  results/c7_presence/<tag>/human_leiden05.{csv,record.json}
#   SAMap object of record (samap_figures.joint_sha256) + Ulrich Datasets S5, S6 (data/external/ulrich2024/)
#     --(c7_stro7_cosine)-->   results/c7_presence/<tag>/human_stromal_cosine.{csv.gz,record.json}
#     --(c7_tables)-->         results/c7_presence/<tag>/tables/{human.csv.gz,dotplot.csv,mouse.csv.gz,
#                              mouse_background.csv.gz,checks.json}
#     --(panels_c7_presence)--> results/figures/panels/c7_presence/FIG_c7_presence_{panels,annotated}_<tag>.{pdf,svg,png}
#                               and FIG_c7_presence_<tag>.provenance.json
#   review_c7_leiden_seeds (not in the default target) --> results/review/c7_presence/<tag>/leiden_seed_sweep_*
#
# Ports the code that made the figure outside the repository (config/c7_presence.yaml names each source). Each rule
# recomputes what the final objects store and fails on a difference: the Leiden labels must equal the labels of record
# nucleus for nucleus, the Stro-7 cosine and its mouse transfer must equal the stored columns within the tolerance.
# The final objects have no producing rule and live outside the repository, so they are ancient() inputs verified by
# content, as the SAMap figures do the joint object. <tag> (config c7_presence.tag) names the objects in every output.
# Target: panels_c7_presence_all (in rule all).

import os

include: "00_panels.smk"

configfile: "config/c7_presence.yaml"

C7 = config["c7_presence"]
C7_TAG = C7["tag"]
C7_OBJ = C7["objects"]
C7_PY = PANEL_PY
C7_CACHE = PANEL_CACHE
C7_SCRIPTS = "workflow/scripts/figures/panels/c7_presence"
C7_COMMON = C7_SCRIPTS + "/c7_common.py"
C7_OUT = f"results/c7_presence/{C7_TAG}"
C7_FIG = C7.get("out_dir", "results/figures/panels/c7_presence")
C7_REVIEW = f"results/review/c7_presence/{C7_TAG}"
C7_LOGS = f"logs/c7_presence/{C7_TAG}"
# The SAMap object of record (29_fig_samap.smk, included before this file), unless config c7_presence names another.
C7_JOINT = C7.get("joint", SF_JOINT)
C7_JOINT_SHA256 = C7.get("joint_sha256", SF_JOINT_SHA256)
C7_ROTATION_TABLE = "resources/figures/mouse_section_rotation.yaml"
C7_ROTATION_READER = "workflow/scripts/figures/section_rotation.py"
C7_VARIANTS = ["panels", "annotated"]
C7_FILES = ["pdf", "svg", "png"]

# Matplotlib and numba caches default to HOME, which is full; HDF5 locking fails on lab_storage reads. The scripts
# import their shared module from the repository, where no bytecode cache should be written.
C7_ENV = (
    "export XDG_CACHE_HOME={c} MPLCONFIGDIR={c}/mpl NUMBA_CACHE_DIR={c}/numba HDF5_USE_FILE_LOCKING=FALSE "
    "PYTHONDONTWRITEBYTECODE=1; "
    "mkdir -p {c}/mpl {c}/numba"
).format(c=C7_CACHE)


def c7_object(name):
    return ancient(C7_OBJ[name]["path"])


# Leiden 0.5 of the human fibroblast/decidual nuclei, recomputed with the recipe of the labels of record.
rule c7_leiden_labels:
    input:
        obj=c7_object("human_fibdec"),
        labels_of_record=C7["leiden"]["labels_of_record"],
        script=C7_SCRIPTS + "/leiden_labels.py",
        common=C7_COMMON,
    output:
        labels=C7_OUT + "/human_leiden05.csv",
        record=C7_OUT + "/human_leiden05.record.json",
    params:
        py=C7_PY,
        env=C7_ENV,
        sha=C7_OBJ["human_fibdec"]["sha256"],
        L=C7["leiden"],
    log:
        C7_LOGS + "/c7_leiden_labels.log",
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --object {input.obj} --object-sha256 {params.sha} \
            --labels-of-record {input.labels_of_record} \
            --use-rep {params.L[use_rep]} --n-neighbors {params.L[n_neighbors]} --resolution {params.L[resolution]} \
            --flavor {params.L[flavor]} --n-iterations {params.L[n_iterations]} --seed {params.L[seed]} \
            --labels {output.labels} --record {output.record} \
            > {log} 2>&1
        """


# Contrastive cosine of each human nucleus to the 11 stromal subtype centroids of Ulrich et al. 2024.
rule c7_stro7_cosine:
    input:
        joint=ancient(C7_JOINT),
        s5=C7["ulrich"]["s5"]["path"],
        s6=C7["ulrich"]["s6"]["path"],
        reference=c7_object("human_fibdec"),
        script=C7_SCRIPTS + "/stro7_cosine.py",
        common=C7_COMMON,
    output:
        cosine=C7_OUT + "/human_stromal_cosine.csv.gz",
        record=C7_OUT + "/human_stromal_cosine.record.json",
    params:
        py=C7_PY,
        env=C7_ENV,
        joint_sha=C7_JOINT_SHA256,
        s5_sha=C7["ulrich"]["s5"]["sha256"],
        s6_sha=C7["ulrich"]["s6"]["sha256"],
        ref_sha=C7_OBJ["human_fibdec"]["sha256"],
        fib=",".join(str(n) for n in C7["cosine"]["fibroblast_subtypes"]),
        ctrl=",".join(str(n) for n in C7["cosine"]["control_subtypes"]),
        K=C7["cosine"],
    log:
        C7_LOGS + "/c7_stro7_cosine.log",
    threads: 1
    resources:
        mem_mb=48000,
        runtime=60,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --joint {input.joint} --joint-sha256 {params.joint_sha} \
            --s5 {input.s5} --s5-sha256 {params.s5_sha} --s6 {input.s6} --s6-sha256 {params.s6_sha} \
            --reference {input.reference} --reference-sha256 {params.ref_sha} \
            --fibroblast {params.fib} --controls {params.ctrl} --target {params.K[target]} \
            --n-panel {params.K[n_panel]} --n-markers {params.K[n_markers]} --tolerance {params.K[tolerance]} \
            --cosine {output.cosine} --record {output.record} \
            > {log} 2>&1
        """


# Plot tables: human nuclei, dot plot, mouse cells with their SAMap-transferred values, mouse tissue background.
rule c7_tables:
    input:
        human=c7_object("human_fibdec"),
        mouse=c7_object("mouse_fibdec"),
        mouse_all=c7_object("mouse_allcells"),
        joint=ancient(C7_JOINT),
        labels=rules.c7_leiden_labels.output.labels,
        cosine=rules.c7_stro7_cosine.output.cosine,
        script=C7_SCRIPTS + "/c7_tables.py",
        common=C7_COMMON,
    output:
        human=C7_OUT + "/tables/human.csv.gz",
        dotplot=C7_OUT + "/tables/dotplot.csv",
        mouse=C7_OUT + "/tables/mouse.csv.gz",
        background=C7_OUT + "/tables/mouse_background.csv.gz",
        checks=C7_OUT + "/tables/checks.json",
    params:
        py=C7_PY,
        env=C7_ENV,
        h_sha=C7_OBJ["human_fibdec"]["sha256"],
        m_sha=C7_OBJ["mouse_fibdec"]["sha256"],
        a_sha=C7_OBJ["mouse_allcells"]["sha256"],
        joint_sha=C7_JOINT_SHA256,
        genes=",".join(C7["dotplot_genes"]),
        K=C7["cosine"],
    log:
        C7_LOGS + "/c7_tables.log",
    threads: 1
    resources:
        mem_mb=32000,
        runtime=60,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --human-obj {input.human} --human-sha256 {params.h_sha} \
            --mouse-obj {input.mouse} --mouse-sha256 {params.m_sha} \
            --mouse-all-obj {input.mouse_all} --mouse-all-sha256 {params.a_sha} \
            --joint {input.joint} --joint-sha256 {params.joint_sha} \
            --labels {input.labels} --cosine {input.cosine} --genes {params.genes} \
            --target {params.K[target]} --tolerance {params.K[tolerance]} \
            --human-table {output.human} --dotplot-table {output.dotplot} --mouse-table {output.mouse} \
            --background-table {output.background} --checks {output.checks} \
            > {log} 2>&1
        """


# The figure, panels only and annotated, each as PDF, SVG (live text, rasterized dots) and 300-dpi PNG.
rule panels_c7_presence:
    input:
        human=rules.c7_tables.output.human,
        dotplot=rules.c7_tables.output.dotplot,
        mouse=rules.c7_tables.output.mouse,
        background=rules.c7_tables.output.background,
        colors=C7["mouse_population_colors"],
        crops=C7["mouse_section_crops"],
        rotation=C7_ROTATION_TABLE,
        rotation_reader=C7_ROTATION_READER,
        script=C7_SCRIPTS + "/c7_presence_panel.py",
        common=C7_COMMON,
    output:
        figs=expand(C7_FIG + "/FIG_c7_presence_{variant}_" + C7_TAG + ".{ext}", variant=C7_VARIANTS, ext=C7_FILES),
        provenance=C7_FIG + "/FIG_c7_presence_" + C7_TAG + ".provenance.json",
    params:
        py=C7_PY,
        env=C7_ENV,
        prefix=C7_FIG + "/FIG_c7_presence",
        tag=C7_TAG,
        joint_name=os.path.basename(C7_JOINT),
    log:
        C7_LOGS + "/panels_c7_presence.log",
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --human-table {input.human} --dotplot-table {input.dotplot} --mouse-table {input.mouse} \
            --background-table {input.background} --mouse-colors {input.colors} --mouse-crops {input.crops} \
            --section-rotation {input.rotation} --joint-name {params.joint_name} \
            --panels-pdf {params.prefix}_panels_{params.tag}.pdf --panels-svg {params.prefix}_panels_{params.tag}.svg \
            --panels-png {params.prefix}_panels_{params.tag}.png \
            --annotated-pdf {params.prefix}_annotated_{params.tag}.pdf \
            --annotated-svg {params.prefix}_annotated_{params.tag}.svg \
            --annotated-png {params.prefix}_annotated_{params.tag}.png \
            --provenance {output.provenance} \
            > {log} 2>&1
        """


rule panels_c7_presence_all:
    input:
        rules.panels_c7_presence.output,


# Review, not in the default target: Leiden seeds on the stored 15-neighbour graph and on the 30-neighbour recipe graph.
rule review_c7_leiden_seeds:
    input:
        obj=c7_object("human_fibdec"),
        labels=rules.c7_leiden_labels.output.labels,
        cosine=rules.c7_stro7_cosine.output.cosine,
        script=C7_SCRIPTS + "/leiden_seed_sweep.py",
        common=C7_COMMON,
    output:
        runs=C7_REVIEW + "/leiden_seed_sweep_runs.csv",
        labels=C7_REVIEW + "/leiden_seed_sweep_labels.csv.gz",
        summary=C7_REVIEW + "/leiden_seed_sweep_summary.json",
    params:
        py=C7_PY,
        env=C7_ENV,
        sha=C7_OBJ["human_fibdec"]["sha256"],
        seeds=C7.get("review_seeds", 20),
    log:
        C7_LOGS + "/review_c7_leiden_seeds.log",
    threads: 1
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --object {input.obj} --object-sha256 {params.sha} --labels {input.labels} --cosine {input.cosine} \
            --seeds {params.seeds} --runs {output.runs} --labels-out {output.labels} --summary {output.summary} \
            > {log} 2>&1
        """
