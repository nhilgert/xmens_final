# workflow/rules/97_methods_numbers.smk
# Pipeline-derived numbers of the computational Methods for the release of record (config deliverables.release):
#   methods_numbers_facts, _facts2, _small_reads (read-only collection) -> methods_numbers_build
#   -> results/methods_numbers/<release>/registry/pipeline_numbers.csv and pipeline_numbers_flags.md
# The registry (workflow/scripts/methods_numbers, README.md there) recomputes every number from the release objects,
# the MI arms and overlap, SCENIC, SAMap and in vitro outputs, and flags the ones the draft text disagrees with. The
# scripts take the repository root from their own location and the release from config.yaml, so a run reads the
# checkout it belongs to. Outputs go to a registry/ folder so the hand-built runs of results/methods_numbers/<release>/
# (submit.sh) are never rewritten. Target: methods_numbers_all.
MN = config["methods_numbers"]
MN_RELEASE = config["deliverables"]["release"]
MN_DIR = f"results/methods_numbers/{MN_RELEASE}/registry"
MN_LOGS = f"logs/methods_numbers/{MN_RELEASE}"
MN_SCRIPTS = "workflow/scripts/methods_numbers"
MN_OVERLAP = config["mi"]["overlap"]["outdir"]
MN_MI_ROOT = config["mi_pseudotime"]["outdir"]
# MI arms the registry reads: the two main arms (11_mi_pseudotime.smk) and two superseded arms quoted in its notes.
MN_MI_ARMS = [
    config["mi"]["overlap"]["mouse_arm"],
    config["mi"]["overlap"]["human_arm"],
    *MN["legacy_mi_arms"],
]
# HDF5 locking off for reads from lab_storage and netscratch; caches outside the full home directory.
MN_ENV = (
    f"export PYTHONDONTWRITEBYTECODE=1 HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR={MN['cache_dir']}/mpl "
    f"NUMBA_CACHE_DIR={MN['cache_dir']}/numba XDG_CACHE_HOME={MN['cache_dir']}"
)


# Facts from the release objects, QC logs, annotation tables, MI runs and overlap, SCENIC and in vitro tables.
rule methods_numbers_facts:
    input:
        manifest=f"data/deliverables/{MN_RELEASE}/MANIFEST.json",
        mi=[f"{MN_MI_ROOT}/{arm}/mi_all.csv" for arm in MN_MI_ARMS],
        overlap_report=f"{MN_OVERLAP}/overlap_report.json",
        cutoff_sweep=f"{MN_OVERLAP}/overlap_cutoff_sweep.csv",
        pairs=f"{MN_OVERLAP}/ortholog_pairs_one2one.csv",
        selection=f"{MN_OVERLAP}/crossspecies_heatmap_selection/selection_summary.json",
        selected=f"{MN_OVERLAP}/crossspecies_heatmap_selection/conserved_profiles_selected_genes.csv",
        script=f"{MN_SCRIPTS}/collect_facts.py",
        release_objects=f"{MN_SCRIPTS}/release_objects.py",
    output:
        facts=directory(f"{MN_DIR}/facts"),
    params:
        python=MN["python"],
        env=MN_ENV,
    log:
        f"{MN_LOGS}/facts.log",
    threads: 2
    resources:
        mem_mb=64000,
        runtime=45,
    shell:
        r"""
        {params.env}
        {params.python} -u {input.script} {output.facts} > {log} 2>&1
        """


# SAMap joint-graph statistics, the in vitro DE release, supplementary tables, reference and CellBender settings.
rule methods_numbers_facts2:
    input:
        manifest=f"data/deliverables/{MN_RELEASE}/MANIFEST.json",
        script=f"{MN_SCRIPTS}/collect_facts2.py",
        release_objects=f"{MN_SCRIPTS}/release_objects.py",
    output:
        facts=directory(f"{MN_DIR}/facts2"),
    params:
        python=MN["python"],
        env=MN_ENV,
    log:
        f"{MN_LOGS}/facts2.log",
    threads: 2
    resources:
        mem_mb=80000,
        runtime=60,
    shell:
        r"""
        {params.env}
        {params.python} -u {input.script} {output.facts} > {log} 2>&1
        """


# Small text sources: SCENIC records, run logs and the scripts whose settings the Methods quote.
rule methods_numbers_small_reads:
    input:
        script=f"{MN_SCRIPTS}/small_reads.py",
    output:
        json=f"{MN_DIR}/small_reads.json",
    params:
        python=MN["python"],
        env=MN_ENV,
    log:
        f"{MN_LOGS}/small_reads.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    shell:
        r"""
        {params.env}
        {params.python} -u {input.script} {output.json} > {log} 2>&1
        """


# Every registry row evaluated against the draft text; the draft is copied next to the table it was compared with.
rule methods_numbers_build:
    input:
        facts=f"{MN_DIR}/facts",
        facts2=f"{MN_DIR}/facts2",
        small_reads=f"{MN_DIR}/small_reads.json",
        draft=MN["draft"],
        artifact_manifest=f"{MN['artifact_inputs']}/MANIFEST.json",
        proteome_check=MN["proteome_check"],
        extra_manifest=f"{MN['extra_inputs']}/MANIFEST.json",
        ptnull="results/samap_figures/ptnull/ptnull_statistics_table.csv",
        script=f"{MN_SCRIPTS}/build_pipeline_numbers.py",
        modules=expand(f"{MN_SCRIPTS}/{{name}}.py", name=["registry_core", "registry_rows", "registry_methods"]),
    output:
        csv=f"{MN_DIR}/pipeline_numbers.csv",
        flags=f"{MN_DIR}/pipeline_numbers_flags.md",
        draft=f"{MN_DIR}/draft.txt",
    params:
        python=MN["python"],
        env=MN_ENV,
        artifact_inputs=MN["artifact_inputs"],
        extra_inputs=MN["extra_inputs"],
    log:
        f"{MN_LOGS}/build.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=20,
    shell:
        r"""
        {params.env}
        cp {input.draft} {output.draft}
        {params.python} -u {input.script} \
            --facts {input.facts} \
            --facts2 {input.facts2} \
            --small-reads {input.small_reads} \
            --artifact-inputs {params.artifact_inputs} \
            --proteome-check {input.proteome_check} \
            --draft {output.draft} \
            --out {output.csv} \
            --report {output.flags} \
            --extra-inputs {params.extra_inputs} \
            > {log} 2>&1
        """


# Target: the Methods numbers of the release of record (requested for rule all at integration).
rule methods_numbers_all:
    input:
        rules.methods_numbers_build.output,
