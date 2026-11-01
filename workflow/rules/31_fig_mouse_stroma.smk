# workflow/rules/31_fig_mouse_stroma.smk
# Mouse stroma supplementary panels, added 2026-10-30: Figs. S15 (Cxcl12, one version per drawn release tag),
# S16 (fibroblast/decidual subtypes on the sections), S17 (subtype UMAP, abundance, markers) and S22
# (Notum on the sections, these three for the release of record only), drawn from the objects named in
# config/panels_mouse_stroma.yaml of the releases of config releases (91_panels_versions.smk).
#
#   panels_mouse_stroma_table  one per-nucleus plot table per object and release (obs, UMAP, genes of X);
#                              the object's sha256 is checked against the release MANIFEST
#   panels_mouse_stroma_s15a   Cxcl12 dot plot per major class and stage
#   panels_mouse_stroma_s15b   Cxcl12 on the 1.5 dpi early and late sections
#   panels_mouse_stroma_s16    each subtype on the six sections
#   panels_mouse_stroma_s17    A UMAP, B subtype abundance per stage, C marker dot plot
#   panels_mouse_stroma_s17aatlas  S17A variant: the same nuclei at their whole-tissue UMAP coordinates
#   panels_mouse_stroma_s22    Notum on the six sections
#   S15B, S16 and S22 turn each section by resources/figures/mouse_section_rotation.yaml (the table of every
#   mouse spatial panel).
#
# Each panel is results/figures/panels/mouse_stroma/panel_<panel>_<tag>.{pdf,svg,png} plus a
# .provenance.json (inputs with sha256, filters, counts); S17Aatlas has no PDF and its SVG and PNG are
# byte-reproducible. Target: panels_mouse_stroma_all.
# The release objects and MANIFESTs are ancient(): a new release is cut on purpose (04g), never re-run
# because an mtime moved.

include: "00_panels.smk"

configfile: "config/panels_mouse_stroma.yaml"

PMS = config["panels_mouse_stroma"]
PMS_OUT = "results/figures/panels/mouse_stroma"
PMS_TABLES = f"{PMS_OUT}/tables"
PMS_SCRIPTS = "workflow/scripts/figures/panels/mouse_stroma"
# Display rotation of the mouse sections: one table read by every mouse spatial panel, and its reader.
# The spatial rules (S15B, S16, S22) list both as inputs, so a changed angle redraws them.
PMS_ROTATION = {
    "rotation": "resources/figures/mouse_section_rotation.yaml",
    "rotation_reader": "workflow/scripts/figures/section_rotation.py",
}
PMS_EXTS = [".pdf", ".svg", ".png", ".provenance.json"]
PMS_ATLAS_EXTS = [".svg", ".png", ".provenance.json"]  # S17Aatlas: SVG for assembly, PNG preview

# The plot scripts import panel_style from their own directory, so PYTHONPATH names it (this also works
# when PYTHONSAFEPATH keeps the script directory off sys.path); matplotlib's cache must not land in HOME.
PMS_ENV = (
    "export XDG_CACHE_HOME={cache} MPLCONFIGDIR={cache}/mpl HDF5_USE_FILE_LOCKING=FALSE PYTHONPATH={scripts}; "
    "mkdir -p {cache}/mpl"
).format(cache=PANEL_CACHE, scripts=PMS_SCRIPTS)

wildcard_constraints:
    pms_object = "|".join(PMS["objects"]),
    pms_release = "|".join(PANEL_TAGS),


def pms_h5ad(wc):
    """Release object behind one plot table (ancient: only a new release triggers a re-run)."""
    name = PMS["objects"][wc.pms_object]["name"]
    return ancient(f"{panel_release_dir(wc.pms_release)}/{name}.h5ad")


def pms_manifest(wc):
    """MANIFEST of the release behind one plot table (ancient, as the object)."""
    return ancient(f"{panel_release_dir(wc.pms_release)}/MANIFEST.json")


def pms_table(obj, release):
    """Plot table and its provenance sidecar, as named inputs."""
    return {
        "table": f"{PMS_TABLES}/{obj}_{release}.csv.gz",
        "table_prov": f"{PMS_TABLES}/{obj}_{release}.provenance.json",
    }


# Per-nucleus plot table of one object and release: everything the panel scripts read.
rule panels_mouse_stroma_table:
    input:
        h5ad = pms_h5ad,
        manifest = pms_manifest,
        script = f"{PMS_SCRIPTS}/extract_panel_table.py",
    output:
        table = f"{PMS_TABLES}/{{pms_object}}_{{pms_release}}.csv.gz",
        prov = f"{PMS_TABLES}/{{pms_object}}_{{pms_release}}.provenance.json",
    params:
        py = PANEL_PY,
        env = PMS_ENV,
        name = lambda wc: PMS["objects"][wc.pms_object]["name"],
        obs = lambda wc: ",".join(PMS["objects"][wc.pms_object]["obs"]),
        obsm = lambda wc: ",".join(PMS["objects"][wc.pms_object]["obsm"]),
        genes = lambda wc: ",".join(PMS["objects"][wc.pms_object]["genes"]),
    log:
        "logs/panels_mouse_stroma/table_{pms_object}_{pms_release}.log",
    threads: 1
    resources:
        mem_mb = 24000,
        runtime = 30,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --h5ad {input.h5ad} \
            --manifest {input.manifest} \
            --object-name {params.name} \
            --obs-cols '{params.obs}' \
            --obsm '{params.obsm}' \
            --genes '{params.genes}' \
            --out {output.table} \
            --prov {output.prov} \
            > {log} 2>&1
        """


# Fig. S15A: Cxcl12 detection and mean per major class and stage, one version per release.
rule panels_mouse_stroma_s15a:
    input:
        table = f"{PMS_TABLES}/alltissue_{{pms_release}}.csv.gz",
        table_prov = f"{PMS_TABLES}/alltissue_{{pms_release}}.provenance.json",
        script = f"{PMS_SCRIPTS}/plot_gene_dotplot_by_class.py",
        style = f"{PMS_SCRIPTS}/panel_style.py",
    output:
        multiext(f"{PMS_OUT}/panel_S15A_{{pms_release}}", *PMS_EXTS),
    params:
        py = PANEL_PY,
        env = PMS_ENV,
        gene = PMS["s15"]["gene"],
        classes = ",".join(PMS["s15"]["classes"]),
        stem = lambda wc: f"{PMS_OUT}/panel_S15A_{wc.pms_release}",
    log:
        "logs/panels_mouse_stroma/s15a_{pms_release}.log",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --table {input.table} \
            --table-prov {input.table_prov} \
            --gene {params.gene} \
            --classes {params.classes} \
            --panel S15A \
            --out {params.stem} \
            > {log} 2>&1
        """


# Fig. S15B: Cxcl12 on the 1.5 dpi early and late sections, one version per release.
rule panels_mouse_stroma_s15b:
    input:
        table = f"{PMS_TABLES}/alltissue_{{pms_release}}.csv.gz",
        table_prov = f"{PMS_TABLES}/alltissue_{{pms_release}}.provenance.json",
        script = f"{PMS_SCRIPTS}/plot_gene_spatial.py",
        style = f"{PMS_SCRIPTS}/panel_style.py",
        **PMS_ROTATION,
    output:
        multiext(f"{PMS_OUT}/panel_S15B_{{pms_release}}", *PMS_EXTS),
    params:
        py = PANEL_PY,
        env = PMS_ENV,
        gene = PMS["s15"]["gene"],
        stages = ",".join(PMS["s15"]["spatial_stages"]),
        stem = lambda wc: f"{PMS_OUT}/panel_S15B_{wc.pms_release}",
    log:
        "logs/panels_mouse_stroma/s15b_{pms_release}.log",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --table {input.table} \
            --table-prov {input.table_prov} \
            --gene {params.gene} \
            --stages {params.stages} \
            --ncols 2 \
            --panel S15B \
            --section-rotation {input.rotation} \
            --out {params.stem} \
            > {log} 2>&1
        """


# Fig. S16: each fibroblast/decidual subtype on the section of each stage (fib/dec object of the release of record).
rule panels_mouse_stroma_s16:
    input:
        **pms_table("fibdec", PANEL_LIVE_TAG),
        script = f"{PMS_SCRIPTS}/plot_subtype_spatial_grid.py",
        style = f"{PMS_SCRIPTS}/panel_style.py",
        **PMS_ROTATION,
    output:
        multiext(f"{PMS_OUT}/panel_S16_{PANEL_LIVE_TAG}", *PMS_EXTS),
    params:
        py = PANEL_PY,
        env = PMS_ENV,
        stem = f"{PMS_OUT}/panel_S16_{PANEL_LIVE_TAG}",
    log:
        "logs/panels_mouse_stroma/s16.log",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --table {input.table} \
            --table-prov {input.table_prov} \
            --panel S16 \
            --section-rotation {input.rotation} \
            --out {params.stem} \
            > {log} 2>&1
        """


# Fig. S17A-C: subtype UMAP, subtype abundance per stage and marker dot plot (fib/dec object of the release of
# record).
rule panels_mouse_stroma_s17:
    input:
        **pms_table("fibdec", PANEL_LIVE_TAG),
        script = f"{PMS_SCRIPTS}/plot_fibdec_subtypes.py",
        style = f"{PMS_SCRIPTS}/panel_style.py",
    output:
        multiext(f"{PMS_OUT}/panel_S17A_{PANEL_LIVE_TAG}", *PMS_EXTS),
        multiext(f"{PMS_OUT}/panel_S17B_{PANEL_LIVE_TAG}", *PMS_EXTS),
        multiext(f"{PMS_OUT}/panel_S17C_{PANEL_LIVE_TAG}", *PMS_EXTS),
    params:
        py = PANEL_PY,
        env = PMS_ENV,
        genes = ",".join(PMS["objects"]["fibdec"]["genes"]),
        stem_a = f"{PMS_OUT}/panel_S17A_{PANEL_LIVE_TAG}",
        stem_b = f"{PMS_OUT}/panel_S17B_{PANEL_LIVE_TAG}",
        stem_c = f"{PMS_OUT}/panel_S17C_{PANEL_LIVE_TAG}",
    log:
        "logs/panels_mouse_stroma/s17.log",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --table {input.table} \
            --table-prov {input.table_prov} \
            --genes {params.genes} \
            --out-a {params.stem_a} \
            --out-b {params.stem_b} \
            --out-c {params.stem_c} \
            > {log} 2>&1
        """


# Fig. S17Aatlas: the S17A nuclei at their whole-tissue UMAP coordinates (obsm X_umap_scVI of the release's
# whole-tissue object, matched by obs_names), coloured and labelled as S17A.
rule panels_mouse_stroma_s17aatlas:
    input:
        **pms_table("fibdec", PANEL_LIVE_TAG),
        atlas_table = f"{PMS_TABLES}/alltissue_umap_{PANEL_LIVE_TAG}.csv.gz",
        atlas_table_prov = f"{PMS_TABLES}/alltissue_umap_{PANEL_LIVE_TAG}.provenance.json",
        script = f"{PMS_SCRIPTS}/plot_fibdec_subtypes_atlas.py",
        style = f"{PMS_SCRIPTS}/panel_style.py",
    output:
        multiext(f"{PMS_OUT}/panel_S17Aatlas_{PANEL_LIVE_TAG}", *PMS_ATLAS_EXTS),
    params:
        py = PANEL_PY,
        env = PMS_ENV,
        stem = f"{PMS_OUT}/panel_S17Aatlas_{PANEL_LIVE_TAG}",
    log:
        "logs/panels_mouse_stroma/s17aatlas.log",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --table {input.table} \
            --table-prov {input.table_prov} \
            --atlas-table {input.atlas_table} \
            --atlas-table-prov {input.atlas_table_prov} \
            --out {params.stem} \
            > {log} 2>&1
        """


# Fig. S22: Notum on the section of each stage (whole-tissue object of the release of record).
rule panels_mouse_stroma_s22:
    input:
        **pms_table("alltissue", PANEL_LIVE_TAG),
        script = f"{PMS_SCRIPTS}/plot_gene_spatial.py",
        style = f"{PMS_SCRIPTS}/panel_style.py",
        **PMS_ROTATION,
    output:
        multiext(f"{PMS_OUT}/panel_S22_{PANEL_LIVE_TAG}", *PMS_EXTS),
    params:
        py = PANEL_PY,
        env = PMS_ENV,
        gene = PMS["s22"]["gene"],
        stages = ",".join(PMS["s22"]["stages"]),
        ncols = PMS["s22"]["ncols"],
        stem = f"{PMS_OUT}/panel_S22_{PANEL_LIVE_TAG}",
    log:
        "logs/panels_mouse_stroma/s22.log",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 15,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --table {input.table} \
            --table-prov {input.table_prov} \
            --gene {params.gene} \
            --stages {params.stages} \
            --ncols {params.ncols} \
            --panel S22 \
            --section-rotation {input.rotation} \
            --out {params.stem} \
            > {log} 2>&1
        """


# Every mouse stroma panel: S15A/B per drawn tag (frozen tags: missing files only); S16, S17A-C, S17Aatlas and
# S22 of the release of record.
rule panels_mouse_stroma_all:
    input:
        panel_targets(f"{PMS_OUT}/panel_{{panel}}_{{tag}}{{ext}}", panel=["S15A", "S15B"], ext=PMS_EXTS),
        expand(f"{PMS_OUT}/panel_{{panel}}_{PANEL_LIVE_TAG}{{ext}}",
               panel=["S16", "S17A", "S17B", "S17C", "S22"], ext=PMS_EXTS),
        multiext(f"{PMS_OUT}/panel_S17Aatlas_{PANEL_LIVE_TAG}", *PMS_ATLAS_EXTS),
