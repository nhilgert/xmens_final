# workflow/rules/34_fig_cross_species.smk
# Cross-species panels of Fig. 6 that the SAMap figure rules (29_fig_samap.smk) do not draw:
#   6A  joint SAMap UMAP of human and X-Mens fibroblast/decidual nuclei, colored by species
#   6B  the same nuclei colored by cell type (eFib, dec)
#   6D  eFib and dec nuclei on the tissue: X-Mens and human, decidualized and menstrual sections, on black
#   6E  X-Mens pseudotime transferred onto the human secretory and menstrual sections, on black
#
#   results/samap_figures/<tag>/tables/joint.parquet       (10c: samap_figure_tables)
#   data/deliverables/<release>/Cevrim_{XMens,Human}_Only_Fib_Dec_Cells.h5ad
#   results/plot_tables/<release>/human_fibdec_plot_table.parquet   (10c: samap_fibdec_plot_table; 6E only)
#     --(panels_cross_species_joint_umap, panels_cross_species_spatial,
#        panels_cross_species_transferred_pseudotime)-->
#   results/figures/panels/cross_species/panel_{6A,6B,6D,6E}_<tag>.{pdf,svg,png,provenance.json}
#   results/figures/panels/cross_species/panel_6E_{secretory,menstrual,key}_<tag>.{pdf,svg,png}   (6E's parts alone)
#
# 6A and 6B use the joint UMAP and square frame of Fig. 6C; labels, embedding_outlier flags and Slide-tags positions
# come from the release objects of the release of record (config deliverables.version; <tag> in the file names,
# 91_panels_versions.smk). A release cut draws a new set beside the earlier ones; the files of release v9 drawn
# before version tags were used here are panel_*_one.* and stay as they are. The mouse row of 6D is turned by
# resources/figures/mouse_section_rotation.yaml, the table every mouse spatial panel reads; the human rows keep the
# crop, mirror and rotation of spatial_efib_dec_panel.SECTIONS.
# 6E reads the tables of rule samap_fig6e (the joint table and the human plot table of the release of record) and
# draws its nuclei, values, color map, crop and key scale, in the orientation of the main text: crop, mirror and
# rotation are those of the human row of 6D, which matches the main-text Figs. 5B, 6D and 6E. samap_fig6e draws the
# sections mirrored left-right. The 6E key is labeled without samap_fig6e's "(rank)" and "early"/"late".
# 5F and 6C are drawn elsewhere (cross-species heatmap folder; rule samap_fig6c).
# Target: panels_cross_species_all. The release objects are ancient(): they are not rebuilt because an mtime moved.
# The SAMap tables are outputs of 10c for the same release and are rebuilt with it.

include: "00_panels.smk"

PCS = config.get("panels_cross_species", {})
PCS_PY = PANEL_PY
PCS_CACHE = PANEL_CACHE
PCS_TAG = PANEL_LIVE_TAG
PCS_RELEASE = PANEL_RELEASES[PCS_TAG]
PCS_JOINT = PCS.get("joint", SF_FIGURE_TABLES["JOINT"])  # 10c, included before this file
PCS_OUT = PCS.get("out_dir", "results/figures/panels/cross_species")
PCS_SCRIPTS = "workflow/scripts/figures/panels/cross_species"
PCS_OBJECT = "data/deliverables/" + PCS_RELEASE + "/Cevrim_{species}_Only_Fib_Dec_Cells.h5ad"
PCS_PANELS = ["6A", "6B", "6D", "6E"]
PCS_FILES = ["pdf", "svg", "png", "provenance.json"]
# Display rotation of the mouse sections (6D mouse row): one table read by every mouse spatial panel, and its reader.
PCS_ROTATION_TABLE = "resources/figures/mouse_section_rotation.yaml"
PCS_ROTATION_READER = "workflow/scripts/figures/section_rotation.py"
# 6E is also written one file per section and one for the key, for assembling the figure by hand.
PCS_6E_PARTS = ["secretory", "menstrual", "key"]
PCS_6E_PART_FILES = ["pdf", "svg", "png"]

# 6E: the human table rule samap_fig6e reads (SF_HUMAN_TABLE of 10c_samap_figures.smk, included before this file),
# so the two panels draw the same nuclei.
PCS_FIG6E_HUMAN = PCS.get("fig6e_human_table", SF_HUMAN_TABLE)

# Matplotlib and numba caches default to HOME, which is full; HDF5 locking fails on lab_storage reads. The scripts
# import their shared module from the repository, where no bytecode cache should be written.
PCS_ENV = (
    "export XDG_CACHE_HOME={c} MPLCONFIGDIR={c}/mpl NUMBA_CACHE_DIR={c}/numba HDF5_USE_FILE_LOCKING=FALSE "
    "PYTHONDONTWRITEBYTECODE=1; "
    "mkdir -p {c}/mpl {c}/numba"
).format(c=PCS_CACHE)


# Fig. 6A (color = species) and 6B (color = cell type): every nucleus of the joint SAMap UMAP, in the frame of 6C.
rule panels_cross_species_joint_umap:
    input:
        joint=PCS_JOINT,
        mouse=ancient(PCS_OBJECT.format(species="XMens")),
        human=ancient(PCS_OBJECT.format(species="Human")),
        script=PCS_SCRIPTS + "/joint_umap_panel.py",
        common=PCS_SCRIPTS + "/cross_species_panels.py",
    output:
        pdf=PCS_OUT + "/panel_{panel}_" + PCS_TAG + ".pdf",
        svg=PCS_OUT + "/panel_{panel}_" + PCS_TAG + ".svg",
        png=PCS_OUT + "/panel_{panel}_" + PCS_TAG + ".png",
        provenance=PCS_OUT + "/panel_{panel}_" + PCS_TAG + ".provenance.json",
    wildcard_constraints:
        panel="6A|6B",
    params:
        py=PCS_PY,
        env=PCS_ENV,
    log:
        "logs/panels/cross_species/panel_{panel}_" + PCS_TAG + ".log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=20,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --panel {wildcards.panel} \
            --joint {input.joint} \
            --mouse-obj {input.mouse} \
            --human-obj {input.human} \
            --pdf {output.pdf} \
            --svg {output.svg} \
            --png {output.png} \
            --provenance {output.provenance} \
            > {log} 2>&1
        """


# Fig. 6D: eFib and dec nuclei of one decidualized and one menstrual section per species, on the QC'd positions.
rule panels_cross_species_spatial:
    input:
        mouse=ancient(PCS_OBJECT.format(species="XMens")),
        human=ancient(PCS_OBJECT.format(species="Human")),
        script=PCS_SCRIPTS + "/spatial_efib_dec_panel.py",
        common=PCS_SCRIPTS + "/cross_species_panels.py",
        rotation=PCS_ROTATION_TABLE,
        rotation_reader=PCS_ROTATION_READER,
    output:
        pdf=PCS_OUT + "/panel_6D_" + PCS_TAG + ".pdf",
        svg=PCS_OUT + "/panel_6D_" + PCS_TAG + ".svg",
        png=PCS_OUT + "/panel_6D_" + PCS_TAG + ".png",
        provenance=PCS_OUT + "/panel_6D_" + PCS_TAG + ".provenance.json",
    params:
        py=PCS_PY,
        env=PCS_ENV,
    log:
        "logs/panels/cross_species/panel_6D_" + PCS_TAG + ".log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=20,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --mouse-obj {input.mouse} \
            --human-obj {input.human} \
            --section-rotation {input.rotation} \
            --pdf {output.pdf} \
            --svg {output.svg} \
            --png {output.png} \
            --provenance {output.provenance} \
            > {log} 2>&1
        """


# Fig. 6E: transferred pseudotime on the human sections, in the orientation of the main text and of 6D.
rule panels_cross_species_transferred_pseudotime:
    input:
        joint=PCS_JOINT,
        human=PCS_FIG6E_HUMAN,
        script=PCS_SCRIPTS + "/sections_transferred_pseudotime_panel.py",
        sections=PCS_SCRIPTS + "/spatial_efib_dec_panel.py",
        common=PCS_SCRIPTS + "/cross_species_panels.py",
    output:
        pdf=PCS_OUT + "/panel_6E_" + PCS_TAG + ".pdf",
        svg=PCS_OUT + "/panel_6E_" + PCS_TAG + ".svg",
        png=PCS_OUT + "/panel_6E_" + PCS_TAG + ".png",
        provenance=PCS_OUT + "/panel_6E_" + PCS_TAG + ".provenance.json",
        parts=expand(PCS_OUT + "/panel_6E_{part}_" + PCS_TAG + ".{ext}", part=PCS_6E_PARTS, ext=PCS_6E_PART_FILES),
    params:
        py=PCS_PY,
        env=PCS_ENV,
        parts_prefix=PCS_OUT + "/panel_6E",
        parts_tag=PCS_TAG,
    log:
        "logs/panels/cross_species/panel_6E_" + PCS_TAG + ".log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=20,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --joint {input.joint} \
            --human-table {input.human} \
            --pdf {output.pdf} \
            --svg {output.svg} \
            --png {output.png} \
            --provenance {output.provenance} \
            --parts-prefix {params.parts_prefix} \
            --parts-tag {params.parts_tag} \
            > {log} 2>&1
        """


# Every panel of this file for the release of record, each as PDF, SVG, 300-dpi PNG and provenance record.
rule panels_cross_species_all:
    input:
        expand(PCS_OUT + "/panel_{panel}_" + PCS_TAG + ".{ext}", panel=PCS_PANELS, ext=PCS_FILES),
        expand(PCS_OUT + "/panel_6E_{part}_" + PCS_TAG + ".{ext}", part=PCS_6E_PARTS, ext=PCS_6E_PART_FILES),
