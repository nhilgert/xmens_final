# workflow/rules/40_fig_scenic_params.smk
# Supplementary Fig. S28: the two SCENIC cutoffs behind the floor-6 consensus regulons used in Fig. 5E.
# Ported 3 Nov 2026 from the Claude Science figure script regulon_recurrence_cutoffs_fig.py (artifact 7bf91dd3).
#
#   s28_regulon_edges   floor-2 consensus pickle -> long TF-target edge table, per species   (pyscenic env)
#   s28_perrun_sizes    per-run cisTarget pickles -> per-TF median regulon size in single runs (pyscenic env)
#   fig_s28_cutoffs     panels A-D -> results/figures/panels/scenic_params/figS28_regulon_recurrence_cutoffs
#                       (.svg, .png, .provenance.json; panel python and cache of config/panels.yaml)
#
# The inputs are outputs of the SCENIC lane under config scenic.rev (scenic_seeds/). That lane predates this DAG
# (scenic.fit: false), so its files are read with ancient(). Target: fig_s28_all (not in rule all).

include: "00_panels.smk"

S28_OUT = "results/figures/panels/scenic_params"
S28_SCRIPTS = "workflow/scripts/figures/panels/scenic_params"
S28_EDGE_PY = "workflow/scripts/scenic/s28_consensus_edges.py"
S28_PERRUN_PY = "workflow/scripts/scenic/s28_perrun_sizes.py"
S28_ORTHO = "resources/orthology/ensembl_human_mouse_components.tsv"
S28_SPECIES = ["mouse", "human"]
S28_STEM = f"{S28_OUT}/figS28_regulon_recurrence_cutoffs"

# The figure is drawn at the floor-6 consensus, the floor of record (scenic.min_occurrence).
if SCENIC_MIN_OCC != 6:
    raise ValueError(f"scenic.min_occurrence is {SCENIC_MIN_OCC}; Fig. S28 is drawn at the floor-6 consensus")


def s28_floor2_dir(sp):
    """Consensus folder without a floor suffix: the floor-2 consensus, whose edges feed panels C and D."""
    return f"{_sc_root(sp)}/consensus_{SCENIC_TAG}_{SCENIC_MODE}"


rule s28_regulon_edges:
    input:
        pkl=lambda wc: ancient(f"{s28_floor2_dir(wc.s28sp)}/consensus_regulons.pkl"),
        script=S28_EDGE_PY,
    output:
        tsv=f"{S28_OUT}/regulon_edges/{{s28sp}}_consensus_regulons_long.tsv",
    wildcard_constraints:
        s28sp="mouse|human",
    params:
        py=SCENIC_PY,
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    log:
        "logs/scenic_params/regulon_edges_{s28sp}.log",
    shell:
        r"""
        {params.py} {input.script} --pkl {input.pkl} --out {output.tsv} > {log} 2>&1
        """


rule s28_perrun_sizes:
    input:
        regs=[ancient(f"{_sc_reg_dir(sp)}/regulons_seed{s}.pkl") for sp in S28_SPECIES for s in SCENIC_SEEDS],
        script=S28_PERRUN_PY,
    output:
        js=f"{S28_OUT}/perrun_sizes.json",
    params:
        py=SCENIC_PY,
        lane=f"{SCENIC_REV}/scenic_seeds",
        regdir=f"regulons_{SCENIC_TAG}",
        n=len(SCENIC_SEEDS),
    threads: 1
    resources:
        mem_mb=32000,
        runtime=60,
    log:
        "logs/scenic_params/perrun_sizes.log",
    shell:
        r"""
        {params.py} {input.script} --lane {params.lane} --regulon-dir {params.regdir} \
            --n-expected {params.n} --out-json {output.js} > {log} 2>&1
        """


rule fig_s28_cutoffs:
    input:
        rep_mouse=ancient(f"{_sc_cons_dir('mouse')}/recurrence_report.tsv"),
        rep_human=ancient(f"{_sc_cons_dir('human')}/recurrence_report.tsv"),
        long_mouse=f"{S28_OUT}/regulon_edges/mouse_consensus_regulons_long.tsv",
        long_human=f"{S28_OUT}/regulon_edges/human_consensus_regulons_long.tsv",
        perrun=f"{S28_OUT}/perrun_sizes.json",
        orthologs=S28_ORTHO,
        script=f"{S28_SCRIPTS}/plot_regulon_recurrence_cutoffs.py",
    output:
        multiext(S28_STEM, ".svg", ".png", ".provenance.json"),
    params:
        py=PANEL_PY,
        env=PANEL_ENV,
        out=S28_STEM,
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    log:
        "logs/scenic_params/fig_s28_cutoffs.log",
    shell:
        r"""
        {params.env}; {params.py} {input.script} \
            --report-mouse {input.rep_mouse} --report-human {input.rep_human} \
            --long-mouse {input.long_mouse} --long-human {input.long_human} \
            --perrun {input.perrun} --orthologs {input.orthologs} --out {params.out} > {log} 2>&1
        """


rule fig_s28_all:
    input:
        multiext(S28_STEM, ".svg", ".png", ".provenance.json"),
