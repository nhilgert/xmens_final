# workflow/rules/30_fig_mouse_tissue.smk
# Mouse whole-tissue panels of the manuscript, one set per release version tag (config panels.versions_drawn,
# 91_panels_versions.smk) (scripts in workflow/scripts/figures/panels/mouse_tissue/):
#   Fig. 3A   UMAP of the major cell classes      Fig. 3B   class abundance per stage
#   Fig. 3C   UMAP of the immune subtypes         Fig. 3D   immune subtype abundance per stage
#   Fig. S12A Slide-tags sample table             Fig. S12B class marker dot plot
#   Fig. S12C log counts per nucleus              Fig. S12D log detected genes per nucleus
#   Fig. S13  spatial map of each major class at each stage (sections turned by
#             resources/figures/mouse_section_rotation.yaml, the table of every mouse spatial panel)
#   Two options for Fig. 3A/3C, drawn alongside the panels above:
#   3Ainset   Fig. 3A with the 3C immune UMAP as an inset, joined to the atlas's immune islands
#   3Catlas   Fig. 3C drawn on the whole-tissue embedding, zoomed to its immune islands
#
#   panels_mouse_tissue_extract  per version: per-nucleus plot table, marker statistics and 10X
#                                lanes, from the release object (sha256 checked against its
#                                MANIFEST), the immune lineage fit and config/samples.tsv
#   panels_mouse_tissue_panel    per panel and version: PDF, SVG, PNG (300 dpi) and a provenance
#                                JSON, drawn from the extract; the immune subtypes' names, colors and
#                                order come from the vocabulary and resources/figures/immune_display.yaml
#
# Each tag is drawn from its release folder (config releases) and the round-2 immune lineage fit of
# config annotation_round2, the fit of every release since v9. Outputs carry the tag in folder and file
# name, so a new release adds files and an earlier tag's panels stay as drawn. Target:
# panels_mouse_tissue_all. The release objects and fits are ancient(): a figure never rebuilds a
# release or refits a lineage because an mtime moved.

include: "91_panels_versions.smk"

PMT_PY      = PANEL_PY
PMT_SCRIPTS = "workflow/scripts/figures/panels/mouse_tissue"
PMT_TABLES  = "results/panels/mouse_tissue/{pmt_version}"
PMT_FIGURES = "results/figures/panels/mouse_tissue/{pmt_version}"

PMT_ENV = PANEL_ENV

# Immune lineage fit whose UMAP Fig. 3C draws, with the sha256 config pins for it.
PMT_IMMUNE_FIT = f"{config['annotation_round2']['outdir']}/latent_mouse_immune.h5ad"
PMT_IMMUNE_FIT_SHA256 = config["annotation_round2"]["latent_sha256"]["mouse_immune"]

# Display rotation of the mouse sections: one table read by every mouse spatial panel, and its reader.
# Spatial panels list both as inputs, so a changed angle redraws them.
PMT_ROTATION_TABLE  = "resources/figures/mouse_section_rotation.yaml"
PMT_ROTATION_READER = "workflow/scripts/figures/section_rotation.py"
PMT_ROTATED_PANELS  = ["S13"]

# Panel -> (script, its own arguments). Fig. S13 draws uninduced section B, the section the
# manuscript shows (the uninduced array carries sections A and B), each section turned by the table.
PMT_PANELS = {
    "3A":   ("plot_umap.py", "--level class"),
    "3B":   ("plot_abundance.py", "--level class"),
    "3C":   ("plot_umap.py", "--level immune"),
    "3D":   ("plot_abundance.py", "--level immune"),
    "S12A": ("plot_s12_table.py", ""),
    "S12B": ("plot_s12_dotplot.py", ""),
    "S12C": ("plot_s12_violins.py", "--metric total_counts"),
    "S12D": ("plot_s12_violins.py", "--metric n_genes"),
    "S13":  ("plot_s13_spatial.py", f"--uninduced-section B --section-rotation {PMT_ROTATION_TABLE}"),
    "3Ainset": ("plot_umap_immune_inset.py", ""),
    "3Catlas": ("plot_umap.py", "--level immune_atlas"),
}
# Scripts a panel script imports besides the style module, so that changing them redraws the panel.
PMT_HELPERS = {"3Ainset": ["plot_umap.py"]}


def pmt_release_file(name):
    """Input function: a file of the version's release directory, ancient()."""
    return lambda wc: ancient(f"{panel_release_dir(wc.pmt_version)}/{name}")


# Plot tables of one version: every number the panels draw comes from these four files.
rule panels_mouse_tissue_extract:
    input:
        h5ad     = pmt_release_file("Cevrim_XMens_All_Cells.h5ad"),
        manifest = pmt_release_file("MANIFEST.json"),
        fit      = ancient(PMT_IMMUNE_FIT),
        samples  = "config/samples.tsv",
        script   = f"{PMT_SCRIPTS}/extract_mouse_tissue_tables.py",
        style    = f"{PMT_SCRIPTS}/mouse_tissue_style.py",
    output:
        cells   = f"{PMT_TABLES}/cells.csv.gz",
        markers = f"{PMT_TABLES}/marker_stats.csv",
        lanes   = f"{PMT_TABLES}/lanes.csv",
        prov    = f"{PMT_TABLES}/extract_provenance.json",
    wildcard_constraints:
        pmt_version = "|".join(PANEL_TAGS),
    params:
        py      = PMT_PY,
        env     = PMT_ENV,
        scripts = PMT_SCRIPTS,
        release = lambda wc: panel_release_dir(wc.pmt_version),
        fit_sha = PMT_IMMUNE_FIT_SHA256,
        outdir  = PMT_TABLES,
    log:
        "logs/panels_mouse_tissue/extract_{pmt_version}.log",
    threads: 2
    resources:
        mem_mb = 16000, runtime = 30,
    shell:
        r"""
        {params.env}
        PYTHONPATH={params.scripts} {params.py} -u {input.script} \
            --release-dir {params.release} \
            --object-name Cevrim_XMens_All_Cells \
            --immune-fit {input.fit} \
            --immune-fit-sha256 '{params.fit_sha}' \
            --samples {input.samples} \
            --version {wildcards.pmt_version} \
            --out-dir {params.outdir} > {log} 2>&1
        """


# One panel of one version, drawn from that version's plot tables only.
rule panels_mouse_tissue_panel:
    input:
        cells   = f"{PMT_TABLES}/cells.csv.gz",
        markers = f"{PMT_TABLES}/marker_stats.csv",
        lanes   = f"{PMT_TABLES}/lanes.csv",
        prov    = f"{PMT_TABLES}/extract_provenance.json",
        script  = lambda wc: f"{PMT_SCRIPTS}/{PMT_PANELS[wc.pmt_panel][0]}",
        helpers = lambda wc: [f"{PMT_SCRIPTS}/{h}" for h in PMT_HELPERS.get(wc.pmt_panel, [])],
        rotation = lambda wc: [PMT_ROTATION_TABLE, PMT_ROTATION_READER] if wc.pmt_panel in PMT_ROTATED_PANELS else [],
        style   = f"{PMT_SCRIPTS}/mouse_tissue_style.py",
        vocabulary = "resources/annotation_vocabulary.yaml",
        display = config["panels"]["immune_display"],
    output:
        multiext(
            f"{PMT_FIGURES}/panel_{{pmt_panel}}_{{pmt_version}}",
            ".pdf",
            ".svg",
            ".png",
            ".provenance.json",
        ),
    wildcard_constraints:
        pmt_version = "|".join(PANEL_TAGS),
        pmt_panel   = "|".join(PMT_PANELS),
    params:
        py      = PMT_PY,
        env     = PMT_ENV,
        scripts = PMT_SCRIPTS,
        args    = lambda wc: PMT_PANELS[wc.pmt_panel][1],
        seed    = config["seed"],
        tables  = PMT_TABLES,
        prefix  = f"{PMT_FIGURES}/panel_{{pmt_panel}}_{{pmt_version}}",
    log:
        "logs/panels_mouse_tissue/panel_{pmt_panel}_{pmt_version}.log",
    threads: 1
    resources:
        mem_mb = 8000, runtime = 20,
    shell:
        r"""
        {params.env}
        PYTHONPATH={params.scripts} {params.py} -u {input.script} \
            --tables-dir {params.tables} \
            --version {wildcards.pmt_version} \
            --seed {params.seed} \
            --vocabulary {input.vocabulary} \
            --immune-display {input.display} \
            {params.args} \
            --out-prefix {params.prefix} > {log} 2>&1
        """


# Every mouse whole-tissue panel of every drawn version (frozen versions: missing files only).
rule panels_mouse_tissue_all:
    input:
        panel_targets(
            "results/figures/panels/mouse_tissue/{tag}/panel_{p}_{tag}.provenance.json",
            p=list(PMT_PANELS),
        ),
