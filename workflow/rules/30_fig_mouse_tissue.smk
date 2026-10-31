# workflow/rules/30_fig_mouse_tissue.smk
# Mouse whole-tissue panels of the manuscript, drawn from pinned releases, one set per annotation
# vocabulary (scripts in workflow/scripts/figures/panels/mouse_tissue/):
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
#                                JSON, drawn from the extract alone
#
# v8 = release 2026_10_29b (vocabulary v8) with the July immune fit; v9 = release 2026_10_30
# (vocabulary v9: perivascular cells, lymphatic endothelium, neutrophils, plasma cells) with the
# round-2 immune fit. Target: panels_mouse_tissue_all (not in rule all). The release objects and
# fits are ancient(): a figure never rebuilds a release or refits a lineage because an mtime moved.

PMT         = config.get("panels_mouse_tissue", {})
PMT_PY      = PMT.get("python", "/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python")
PMT_CACHE   = PMT.get("cache_dir", "/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch")
PMT_SCRIPTS = "workflow/scripts/figures/panels/mouse_tissue"
PMT_TABLES  = "results/panels/mouse_tissue/{pmt_version}"
PMT_FIGURES = "results/figures/panels/mouse_tissue/{pmt_version}"

# Matplotlib and numba caches default to HOME, which is full; HDF5 locking fails on lab_storage reads.
PMT_ENV = ("export XDG_CACHE_HOME={c} MPLCONFIGDIR={c}/mpl NUMBA_CACHE_DIR={c}/numba HDF5_USE_FILE_LOCKING=FALSE; "
           "mkdir -p {c}/mpl {c}/numba").format(c=PMT_CACHE)

# The release and immune fit behind each version. The July fit has no checksum pinned in config;
# its sha256 is recorded in the extract provenance instead.
PMT_VERSIONS = {
    "v8": {
        "release": "data/deliverables/2026_10_29b",
        "immune_fit": "data/interim/latent_mouse_immune.h5ad",
        "immune_fit_sha256": "",
    },
    "v9": {
        "release": "data/deliverables/2026_10_30",
        "immune_fit": "data/interim/round2_2026_10_29/latent_mouse_immune.h5ad",
        "immune_fit_sha256": config["annotation_round2"]["latent_sha256"]["mouse_immune"],
    },
}

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
    return lambda wc: ancient(f"{PMT_VERSIONS[wc.pmt_version]['release']}/{name}")


# Plot tables of one version: every number the panels draw comes from these four files.
rule panels_mouse_tissue_extract:
    input:
        h5ad     = pmt_release_file("Cevrim_XMens_All_Cells.h5ad"),
        manifest = pmt_release_file("MANIFEST.json"),
        fit      = lambda wc: ancient(PMT_VERSIONS[wc.pmt_version]["immune_fit"]),
        samples  = "config/samples.tsv",
        script   = f"{PMT_SCRIPTS}/extract_mouse_tissue_tables.py",
        style    = f"{PMT_SCRIPTS}/mouse_tissue_style.py",
    output:
        cells   = f"{PMT_TABLES}/cells.csv.gz",
        markers = f"{PMT_TABLES}/marker_stats.csv",
        lanes   = f"{PMT_TABLES}/lanes.csv",
        prov    = f"{PMT_TABLES}/extract_provenance.json",
    wildcard_constraints:
        pmt_version = "|".join(PMT_VERSIONS),
    params:
        py      = PMT_PY,
        env     = PMT_ENV,
        scripts = PMT_SCRIPTS,
        release = lambda wc: PMT_VERSIONS[wc.pmt_version]["release"],
        fit_sha = lambda wc: PMT_VERSIONS[wc.pmt_version]["immune_fit_sha256"],
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
    output:
        multiext(
            f"{PMT_FIGURES}/panel_{{pmt_panel}}_{{pmt_version}}",
            ".pdf",
            ".svg",
            ".png",
            ".provenance.json",
        ),
    wildcard_constraints:
        pmt_version = "|".join(PMT_VERSIONS),
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
            {params.args} \
            --out-prefix {params.prefix} > {log} 2>&1
        """


# Every mouse whole-tissue panel in both versions.
rule panels_mouse_tissue_all:
    input:
        expand(
            "results/figures/panels/mouse_tissue/{v}/panel_{p}_{v}.provenance.json",
            v=PMT_VERSIONS,
            p=PMT_PANELS,
        ),
