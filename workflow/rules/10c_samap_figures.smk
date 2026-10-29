# workflow/rules/10c_samap_figures.smk
# SAMap pseudotime figures and the permutation control for the pseudotime transfer
# (Methods: "Cross-species alignment and transfer of categorical and continuous values (SAMap)").
#
#   runs_raw/<arm>/samap_transferred_pt_harmonized.h5ad     (10b: samap_raw_attach_pt, variant _harmonized)
#     --(samap_graph_block)-->   results/samap_figures/graph_block/C_hs_mm.npz (+ names, meta.json)
#     --(samap_figure_tables)--> results/samap_figures/tables/{joint.parquet, mouse_cells.csv, human_cells.csv}
#     --(samap_edge_share)-->    results/samap_figures/edge_share/samap_compartment_edge_share.csv
#     --(samap_ptnull)-->        results/samap_figures/ptnull/{null_replicates, example_transfers,
#                                                              ptnull_statistics_table}.csv
#     --(samap_supp_figure)-->   results/figures/samap_supp/     SAMap supplementary figure, panels a-g
#     --(samap_main_panels)-->   results/figures/samap_main/     Fig. 6 transferred-pseudotime UMAPs and sections
#
# Labels, coordinates and embedding_outlier flags come from the release whose pseudotime was transferred
# (config samap_figures.release; later releases are relabel-only). Display rules as in Fig. 6: fixed 0-1
# viridis, values as stored, embedding outliers not drawn; every statistic uses all nuclei.
# methods_numbers reads results/samap_figures/ptnull/ptnull_statistics_table.csv.
# Target: samap_figures_all (not in rule all). The release objects and the SAMap object are ancient():
# neither is rebuilt because an mtime moved.

SF        = config.get("samap_figures", {})
SF_PY     = SF.get("python", "/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python")
SF_CACHE  = SF.get("cache_dir", "/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch")
SF_ARM    = SF.get("arm", "vanilla_raw")
SF_PT     = SF.get("pt_variant", "_harmonized")
SF_JOINT  = f"{SR_ROOT}/runs_raw/{SF_ARM}/samap_transferred_pt{SF_PT}.h5ad"
SF_TABLES = "data/deliverables/{release}/plot_tables/{{species}}_fibdec_plot_table.parquet".format(
    release=SF.get("release", "2026_10_22_flagged_harmonized"))

# Matplotlib and numba caches default to HOME, which is full; HDF5 locking fails on lab_storage reads.
SF_ENV = ("export XDG_CACHE_HOME={c} MPLCONFIGDIR={c}/mpl NUMBA_CACHE_DIR={c}/numba HDF5_USE_FILE_LOCKING=FALSE; "
          "mkdir -p {c}/mpl {c}/numba").format(c=SF_CACHE)

assert SF_PT in SR_PT_VARIANTS, f"samap_figures.pt_variant {SF_PT!r} is not a samap_raw pt variant"


# Human x mouse block of the joint SAMap graph: the kernel of the pseudotime transfer, held fixed under permutation.
rule samap_graph_block:
    input:
        joint  = ancient(SF_JOINT),
        script = "workflow/scripts/samap/export_graph_block.py",
    output:
        npz  = "results/samap_figures/graph_block/C_hs_mm.npz",
        hs   = "results/samap_figures/graph_block/hs_names.txt",
        mm   = "results/samap_figures/graph_block/mm_names.txt",
        meta = "results/samap_figures/graph_block/meta.json",
    params:
        py     = SF_PY,
        env    = SF_ENV,
        outdir = "results/samap_figures/graph_block",
    threads: 1
    resources:
        mem_mb = 16000, runtime = 30,
    log:
        "logs/samap_figures/graph_block.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --joint {input.joint} \
            --out-dir {params.outdir} \
            > {log} 2>&1
        """


# One row per nucleus with everything the figures and the permutation control read: joint UMAP and transfer,
# release labels, coordinates, pseudotime and embedding_outlier flags.
rule samap_figure_tables:
    input:
        joint  = ancient(SF_JOINT),
        mouse  = ancient(SF_TABLES.format(species="mouse")),
        human  = ancient(SF_TABLES.format(species="human")),
        script = "workflow/scripts/samap/export_samap_tables.py",
    output:
        joint = "results/samap_figures/tables/joint.parquet",
        mouse = "results/samap_figures/tables/mouse_cells.csv",
        human = "results/samap_figures/tables/human_cells.csv",
    params:
        py     = SF_PY,
        env    = SF_ENV,
        outdir = "results/samap_figures/tables",
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        "logs/samap_figures/figure_tables.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --joint {input.joint} \
            --mouse-table {input.mouse} \
            --human-table {input.human} \
            --out-dir {params.outdir} \
            > {log} 2>&1
        """


# Share of cross-species edge weight landing in each compartment (supplement panel d; methods: 88/82/85/92%).
rule samap_edge_share:
    input:
        joint  = ancient(SF_JOINT),
        mouse  = "results/samap_figures/tables/mouse_cells.csv",
        human  = "results/samap_figures/tables/human_cells.csv",
        script = "workflow/scripts/samap/samap_edge_share.py",
    output:
        figure   = "results/samap_figures/edge_share/samap_compartment_edge_share.csv",
        variants = "results/samap_figures/edge_share/samap_compartment_edge_share_all.csv",
    params:
        py     = SF_PY,
        env    = SF_ENV,
        outdir = "results/samap_figures/edge_share",
    threads: 1
    resources:
        mem_mb = 16000, runtime = 30,
    log:
        "logs/samap_figures/edge_share.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --joint {input.joint} \
            --human-cells {input.human} \
            --mouse-cells {input.mouse} \
            --out-dir {params.outdir} \
            > {log} 2>&1
        """


# Permutation control: mouse pseudotime shuffled across all mouse nuclei or within compartment before transfer
# (supplement panels f-g; methods paragraph on the permutation control).
rule samap_ptnull:
    input:
        npz    = "results/samap_figures/graph_block/C_hs_mm.npz",
        hs     = "results/samap_figures/graph_block/hs_names.txt",
        mm     = "results/samap_figures/graph_block/mm_names.txt",
        mouse  = "results/samap_figures/tables/mouse_cells.csv",
        human  = "results/samap_figures/tables/human_cells.csv",
        script = "workflow/scripts/samap/ptnull_permutation.py",
    output:
        replicates = "results/samap_figures/ptnull/null_replicates.csv",
        examples   = "results/samap_figures/ptnull/example_transfers.csv",
        summary    = "results/samap_figures/ptnull/null_summary.json",
        table      = "results/samap_figures/ptnull/ptnull_statistics_table.csv",
    params:
        py       = SF_PY,
        env      = SF_ENV,
        graphdir = "results/samap_figures/graph_block",
        outdir   = "results/samap_figures/ptnull",
        seed     = SF.get("ptnull_seed", 20261022),
        n_perm   = SF.get("ptnull_n_perm", 1000),
    threads: 1
    resources:
        mem_mb = 16000, runtime = 60,
    log:
        "logs/samap_figures/ptnull.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --graph-dir {params.graphdir} \
            --mouse-cells {input.mouse} \
            --human-cells {input.human} \
            --seed {params.seed} \
            --n-perm {params.n_perm} \
            --out-dir {params.outdir} \
            > {log} 2>&1
        """


# SAMap supplementary figure: transfer on the joint UMAP (a-c), edge shares (d), transferred vs human
# pseudotime (e), and the permutation control in space (f) and as statistics (g).
rule samap_supp_figure:
    input:
        joint      = "results/samap_figures/tables/joint.parquet",
        mouse      = "results/samap_figures/tables/mouse_cells.csv",
        human      = "results/samap_figures/tables/human_cells.csv",
        edge_share = "results/samap_figures/edge_share/samap_compartment_edge_share.csv",
        replicates = "results/samap_figures/ptnull/null_replicates.csv",
        examples   = "results/samap_figures/ptnull/example_transfers.csv",
        script     = "workflow/scripts/samap/plot_sfig_samap.py",
    output:
        directory("results/figures/samap_supp"),
    params:
        py  = SF_PY,
        env = SF_ENV,
    threads: 1
    resources:
        mem_mb = 8000, runtime = 30,
    log:
        "logs/samap_figures/supp_figure.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --joint {input.joint} \
            --mouse-cells {input.mouse} \
            --human-cells {input.human} \
            --edge-share {input.edge_share} \
            --null-replicates {input.replicates} \
            --example-transfers {input.examples} \
            --out {output} \
            > {log} 2>&1
        """


# Fig. 6 panels on black: transferred pseudotime on the joint UMAP and on the two human sections.
rule samap_main_panels:
    input:
        joint  = "results/samap_figures/tables/joint.parquet",
        human  = ancient(SF_TABLES.format(species="human")),
        script = "workflow/scripts/samap/plot_samap_main_dark.py",
    output:
        directory("results/figures/samap_main"),
    params:
        py   = SF_PY,
        env  = SF_ENV,
        date = SF.get("main_panel_date", "20261027"),
    threads: 1
    resources:
        mem_mb = 8000, runtime = 30,
    log:
        "logs/samap_figures/main_panels.log",
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --joint {input.joint} \
            --human {input.human} \
            --date {params.date} \
            --out {output} \
            > {log} 2>&1
        """


rule samap_figures_all:
    input:
        "results/samap_figures/edge_share/samap_compartment_edge_share.csv",
        "results/samap_figures/ptnull/ptnull_statistics_table.csv",
        "results/figures/samap_supp",
        "results/figures/samap_main",
