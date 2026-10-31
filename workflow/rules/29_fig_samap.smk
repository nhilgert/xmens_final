# workflow/rules/29_fig_samap.smk
# SAMap pseudotime figures and the permutation control for the pseudotime transfer
# (Methods: "Cross-species alignment and transfer of categorical and continuous values (SAMap)").
#
#   runs_raw/<arm>/samap_transferred_pt_harmonized.h5ad     (10b: samap_raw_attach_pt, variant _harmonized)
#     --(samap_graph_block)-->   results/samap_figures/graph_block/C_hs_mm.npz (+ names, meta.json)
#     --(samap_figure_tables)--> results/samap_figures/tables/{joint.parquet, mouse_cells.csv, human_cells.csv}
#     --(samap_edge_share)-->    results/samap_figures/edge_share/samap_compartment_edge_share.csv
#     --(samap_ptnull)-->        results/samap_figures/ptnull/{null_replicates, example_transfers,
#                                                              ptnull_statistics_table}.csv
#     --(samap_supp_panel)-->    results/figures/samap_supp/panels/  one panel of the SAMap supplement (a-g)
#     --(samap_supp_figure)-->   results/figures/samap_supp/     SAMap supplementary figure, panels a-g assembled
#     --(samap_fig6c, samap_fig6e)--> results/figures/samap_main/  Fig. 6C and 6E on black
#
# The figures are drawn by the notebooks in notebooks/samap, one per panel, run here with papermill; the executed
# copies, with the figure and the numbers behind it, go to results/notebooks/samap/.
#
# Labels, coordinates and embedding_outlier flags come from the release whose pseudotime was transferred
# (config samap_figures.release; later releases are relabel-only). Human pseudotime, which the figures and the
# permutation control compare the transfer with, comes from samap_figures.human_pt_table (see SF_HUMAN_PT). Display rules as in Fig. 6: fixed 0-1
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

# Human pseudotime compared with the transfer: recomputed on 30 Oct 2026 on the run4 human fibroblast/decidual
# latent (PCA graph, released nuclei), copied from netscratch with its provenance (README.md in that directory).
SF_HUMAN_PT     = SF.get("human_pt_table", "data/interim/human_pt_2026_10_30_run4_a40/cells.parquet")
SF_HUMAN_PT_COL = SF.get("human_pt_column", "pt_rank")

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
# release labels, coordinates, pseudotime, embedding_outlier flags and the human pseudotime they compare against.
rule samap_figure_tables:
    input:
        joint  = ancient(SF_JOINT),
        mouse  = ancient(SF_TABLES.format(species="mouse")),
        human    = ancient(SF_TABLES.format(species="human")),
        human_pt = SF_HUMAN_PT,
        script   = "workflow/scripts/samap/export_samap_tables.py",
    output:
        joint = "results/samap_figures/tables/joint.parquet",
        mouse = "results/samap_figures/tables/mouse_cells.csv",
        human = "results/samap_figures/tables/human_cells.csv",
    params:
        py           = SF_PY,
        env          = SF_ENV,
        outdir       = "results/samap_figures/tables",
        human_pt_col = SF_HUMAN_PT_COL,
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
            --human-pt {input.human_pt} \
            --human-pt-col {params.human_pt_col} \
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


# Figures ----
# Notebooks import their shared modules from their own directory, so papermill runs them with that as the
# working directory. Jupyter and IPython write to HOME by default, and the kernel spec launches "python" from PATH.
SF_NOTEBOOKS = "notebooks/samap"
SF_EXECUTED  = "results/notebooks/samap"
SF_BIN       = SF_PY.rsplit("/", 1)[0]
SF_NB_ENV = (SF_ENV + "; export PATH={b}:$PATH IPYTHONDIR={c}/ipython JUPYTER_RUNTIME_DIR={c}/jupyter_runtime "
             "JUPYTER_DATA_DIR={c}/jupyter_data JUPYTER_CONFIG_DIR={c}/jupyter_config; "
             "mkdir -p {c}/ipython {c}/jupyter_runtime {c}/jupyter_data {c}/jupyter_config").format(b=SF_BIN, c=SF_CACHE)

# Tables the figure notebooks read, by the notebooks' parameter names.
SF_FIGURE_TABLES = {
    "JOINT":             "results/samap_figures/tables/joint.parquet",
    "MOUSE_CELLS":       "results/samap_figures/tables/mouse_cells.csv",
    "HUMAN_CELLS":       "results/samap_figures/tables/human_cells.csv",
    "EDGE_SHARE":        "results/samap_figures/edge_share/samap_compartment_edge_share.csv",
    "NULL_REPLICATES":   "results/samap_figures/ptnull/null_replicates.csv",
    "EXAMPLE_TRANSFERS": "results/samap_figures/ptnull/example_transfers.csv",
}

# Supplement panels: notebook and the tables it reads.
SF_PANELS = {
    "a": ("sfig_samap_a_umap_mouse_pseudotime",       ["JOINT", "MOUSE_CELLS", "HUMAN_CELLS"]),
    "b": ("sfig_samap_b_umap_transferred_pseudotime", ["JOINT", "MOUSE_CELLS", "HUMAN_CELLS"]),
    "c": ("sfig_samap_c_umap_human_pseudotime",       ["JOINT", "MOUSE_CELLS", "HUMAN_CELLS"]),
    "d": ("sfig_samap_d_edge_share",                  ["EDGE_SHARE"]),
    "e": ("sfig_samap_e_transferred_vs_human",        ["HUMAN_CELLS"]),
    "f": ("sfig_samap_f_sections_shuffled",           ["HUMAN_CELLS", "EXAMPLE_TRANSFERS"]),
    "g": ("sfig_samap_g_permutation_statistics",      ["NULL_REPLICATES", "EXAMPLE_TRANSFERS", "HUMAN_CELLS"]),
}
SF_PANEL_DIR = "results/figures/samap_supp/panels"
SF_PANEL_FILES = [".png", "_600dpi.png", ".svg", ".provenance.json"]

SF_MAIN_DIR  = "results/figures/samap_main"
SF_MAIN_DATE = str(SF.get("main_panel_date", "20261027"))
SF_HUMAN_TABLE = SF_TABLES.format(species="human")


def sf_panel_inputs(wildcards):
    """Notebook, shared modules and tables of one supplement panel."""
    notebook, tables = SF_PANELS[wildcards.panel]
    inputs = {"notebook": f"{SF_NOTEBOOKS}/{notebook}.ipynb", "style": f"{SF_NOTEBOOKS}/sfig_samap_style.py"}
    if wildcards.panel in "abc":
        inputs["umap"] = f"{SF_NOTEBOOKS}/sfig_samap_umap.py"
    return inputs | {name: SF_FIGURE_TABLES[name] for name in tables}


def sf_panel_tables(wildcards):
    """papermill flags passing the panel's tables."""
    return " ".join(f"-p {name} {SF_FIGURE_TABLES[name]}" for name in SF_PANELS[wildcards.panel][1])


# One panel of the SAMap supplement, drawn by its notebook at its exact size in the figure's grid.
rule samap_supp_panel:
    input:
        unpack(sf_panel_inputs),
    output:
        png        = SF_PANEL_DIR + "/sfig_samap_panel_{panel}.png",
        png_600dpi = SF_PANEL_DIR + "/sfig_samap_panel_{panel}_600dpi.png",
        svg        = SF_PANEL_DIR + "/sfig_samap_panel_{panel}.svg",
        provenance = SF_PANEL_DIR + "/sfig_samap_panel_{panel}.provenance.json",
        executed   = SF_EXECUTED + "/sfig_samap_panel_{panel}.executed.ipynb",
    wildcard_constraints:
        panel = "[a-g]",
    params:
        env       = SF_NB_ENV,
        papermill = SF_BIN + "/papermill",
        cwd       = SF_NOTEBOOKS,
        out       = SF_PANEL_DIR,
        tables    = sf_panel_tables,
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        "logs/samap_figures/supp_panel_{panel}.log",
    shell:
        r"""
        {params.env}
        {params.papermill} {input.notebook} {output.executed} \
            --cwd {params.cwd} \
            -p REPO "$PWD" \
            -p NOTEBOOK {input.notebook} \
            -p OUT {params.out} \
            -p FORCE True \
            {params.tables} \
            --log-output \
            > {log} 2>&1
        """


# SAMap supplementary figure: the seven panels placed on the 180 mm grid with their letters (a-c joint UMAP,
# d edge shares, e transferred vs human pseudotime, f-g the permutation control).
rule samap_supp_figure:
    input:
        notebook = SF_NOTEBOOKS + "/sfig_samap_compose.ipynb",
        style    = SF_NOTEBOOKS + "/sfig_samap_style.py",
        panels   = expand(SF_PANEL_DIR + "/sfig_samap_panel_{panel}{suffix}", panel=SF_PANELS, suffix=SF_PANEL_FILES),
    output:
        png        = "results/figures/samap_supp/sfig_samap.png",
        png_600dpi = "results/figures/samap_supp/sfig_samap_600dpi.png",
        svg        = "results/figures/samap_supp/sfig_samap.svg",
        provenance = "results/figures/samap_supp/sfig_samap.provenance.json",
        executed   = SF_EXECUTED + "/sfig_samap_compose.executed.ipynb",
    params:
        env       = SF_NB_ENV,
        papermill = SF_BIN + "/papermill",
        cwd       = SF_NOTEBOOKS,
        panels    = SF_PANEL_DIR,
        out       = "results/figures/samap_supp",
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        "logs/samap_figures/supp_figure.log",
    shell:
        r"""
        {params.env}
        {params.papermill} {input.notebook} {output.executed} \
            --cwd {params.cwd} \
            -p REPO "$PWD" \
            -p NOTEBOOK {input.notebook} \
            -p PANELS {params.panels} \
            -p OUT {params.out} \
            -p FORCE True \
            --log-output \
            > {log} 2>&1
        """


# Fig. 6C on black: transferred pseudotime on the joint UMAP, with and without the mouse nuclei underneath.
rule samap_fig6c:
    input:
        notebook = SF_NOTEBOOKS + "/fig6c_umap_transferred_pseudotime.ipynb",
        joint    = SF_FIGURE_TABLES["JOINT"],
        human    = ancient(SF_HUMAN_TABLE),
    output:
        expand(SF_MAIN_DIR + "/{stem}.{ext}",
               stem=["umap_dark_transferred_pseudotime", "umap_dark_transferred_pseudotime_humanonly"],
               ext=["png", "svg", "pdf"]),
        provenance = SF_MAIN_DIR + "/fig6c_umap_transferred_pseudotime.provenance.json",
        executed   = SF_EXECUTED + "/fig6c_umap_transferred_pseudotime.executed.ipynb",
    params:
        env       = SF_NB_ENV,
        papermill = SF_BIN + "/papermill",
        cwd       = SF_NOTEBOOKS,
        out       = SF_MAIN_DIR,
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        "logs/samap_figures/fig6c.log",
    shell:
        r"""
        {params.env}
        {params.papermill} {input.notebook} {output.executed} \
            --cwd {params.cwd} \
            -p REPO "$PWD" \
            -p NOTEBOOK {input.notebook} \
            -p JOINT {input.joint} \
            -p HUMAN_TABLE {input.human} \
            -p OUT {params.out} \
            -p FORCE True \
            --log-output \
            > {log} 2>&1
        """


# Fig. 6E on black: transferred pseudotime in the two human sections (bare panels) and their color key.
rule samap_fig6e:
    input:
        notebook = SF_NOTEBOOKS + "/fig6e_sections_transferred_pseudotime.ipynb",
        joint    = SF_FIGURE_TABLES["JOINT"],
        human    = ancient(SF_HUMAN_TABLE),
    output:
        expand(SF_MAIN_DIR + "/" + SF_MAIN_DATE + "-Human-{phase}-transferred_pseudotime-spatial.{ext}",
               phase=["secretory", "menstrual"], ext=["png", "pdf"]),
        expand(SF_MAIN_DIR + "/" + SF_MAIN_DATE + "-spatial_key_transferred_pseudotime.{ext}",
               ext=["png", "pdf", "svg"]),
        provenance = SF_MAIN_DIR + "/fig6e_sections_transferred_pseudotime.provenance.json",
        executed   = SF_EXECUTED + "/fig6e_sections_transferred_pseudotime.executed.ipynb",
    params:
        env       = SF_NB_ENV,
        papermill = SF_BIN + "/papermill",
        cwd       = SF_NOTEBOOKS,
        out       = SF_MAIN_DIR,
        date      = SF_MAIN_DATE,
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    log:
        "logs/samap_figures/fig6e.log",
    shell:
        r"""
        {params.env}
        {params.papermill} {input.notebook} {output.executed} \
            --cwd {params.cwd} \
            -p REPO "$PWD" \
            -p NOTEBOOK {input.notebook} \
            -p JOINT {input.joint} \
            -p HUMAN_TABLE {input.human} \
            -p OUT {params.out} \
            -r DATE {params.date} \
            -p FORCE True \
            --log-output \
            > {log} 2>&1
        """


rule samap_figures_all:
    input:
        "results/samap_figures/edge_share/samap_compartment_edge_share.csv",
        "results/samap_figures/ptnull/ptnull_statistics_table.csv",
        rules.samap_supp_figure.output,
        rules.samap_fig6c.output,
        rules.samap_fig6e.output,
