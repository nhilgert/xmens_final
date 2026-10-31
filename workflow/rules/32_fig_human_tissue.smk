# workflow/rules/32_fig_human_tissue.smk
# Human whole-tissue panels of the revision, drawn from release v8 (2026_10_29b) and v9 (2026_10_30):
#   Fig. 5B   spatial mapping of all cell classes, secretory and menstrual sections, on black
#   Fig. 5C   UMAP of all nuclei by cell class, direct labels
#   Fig. S25  A sample table, B class marker dot plot, C counts and D genes per nucleus (violins)
#   Fig. S26  one spatial map per cell class (rows) and sample (columns), the class in magenta
#
#   panels_human_tissue_extract          per-nucleus plot table of one release object (obs, scVI UMAP,
#                                        counts and genes per nucleus, marker log1p CP10k); the counts
#                                        layer is streamed in row blocks, the object is never rewritten
#   panels_human_tissue_{5b,5c,s25,s26}  one rule per panel, reading only that plot table
#
# Outputs: results/figures/panels/human_tissue/panel_<panel>_<v8|v9>.{pdf,svg,png,provenance.json}.
# Target: panels_human_tissue_all (not in rule all).
#
# The release objects and manifests are ancient(): a release is cut deliberately (04g), never re-cut
# because an mtime moved. The panel scripts import panel_common.py from their own directory, which
# PYTHONPATH names because Python's safe-path mode (PYTHONSAFEPATH) drops it from sys.path.

PANELS_HT_RELEASES = {
    "v8": "data/deliverables/2026_10_29b",
    "v9": "data/deliverables/2026_10_30",
}
PANELS_HT_OUT = "results/figures/panels/human_tissue"
PANELS_HT_SCRIPTS = "workflow/scripts/figures/panels/human_tissue"
PANELS_HT_LOGS = "logs/panels_human_tissue"
PANELS_HT_TABLE = PANELS_HT_OUT + "/data/{ht_release}"

# S25B genes in column order, as in the published dot plot. Its column "LEFTY" is LEFTY2, the
# decidual marker; LEFTY1 is barely detected in any class.
PANELS_HT_MARKERS = [
    "EPCAM", "KLF5", "COL5A1", "COL6A3", "DCN", "LEFTY2", "CEBPB", "IGFBP1", "PTPRC",
    "STK17B", "PECAM1", "VWF", "MCAM", "CSPG4", "ACTA2", "MYH11", "PGR", "ESR1",
]

PANELS_HT_PANELS = ["5B", "5C", "S25A", "S25B", "S25C", "S25D", "S26"]
PANELS_HT_EXTENSIONS = ["pdf", "svg", "png", "provenance.json"]


wildcard_constraints:
    ht_release="|".join(PANELS_HT_RELEASES),


def panels_ht_release_file(name):
    """Input function: file `name` in the deliverables directory of release {ht_release}."""
    return lambda wildcards: f"{PANELS_HT_RELEASES[wildcards.ht_release]}/{name}"


def panels_ht_files(panel):
    """The four output files of one panel, with {ht_release} left as a wildcard."""
    return [f"{PANELS_HT_OUT}/panel_{panel}_{{ht_release}}.{ext}" for ext in PANELS_HT_EXTENSIONS]


# Every panel of this track, both releases.
rule panels_human_tissue_all:
    input:
        expand(
            f"{PANELS_HT_OUT}/panel_{{panel}}_{{release}}.{{ext}}",
            panel=PANELS_HT_PANELS,
            release=list(PANELS_HT_RELEASES),
            ext=PANELS_HT_EXTENSIONS,
        ),


# Plot table of one release: everything the four panel rules draw, so the 2 GB object is read once.
rule panels_human_tissue_extract:
    input:
        h5ad=ancient(panels_ht_release_file("Cevrim_Human_All_Cells.h5ad")),
        manifest=ancient(panels_ht_release_file("MANIFEST.json")),
        script=f"{PANELS_HT_SCRIPTS}/extract_human_tissue.py",
    output:
        cells=f"{PANELS_HT_TABLE}/cells.parquet",
        uns=f"{PANELS_HT_TABLE}/uns.json",
        provenance=f"{PANELS_HT_TABLE}/provenance.json",
    params:
        genes=",".join(PANELS_HT_MARKERS),
    log:
        f"{PANELS_HT_LOGS}/extract_{{ht_release}}.log",
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        """
        python -u {input.script} \
            --h5ad {input.h5ad} \
            --manifest {input.manifest} \
            --release {wildcards.ht_release} \
            --genes {params.genes} \
            --out-cells {output.cells} \
            --out-uns {output.uns} \
            --out-provenance {output.provenance} \
            > {log} 2>&1
        """


# Fig. 5B: spatial mapping of all cell classes in both sections.
rule panels_human_tissue_5b:
    input:
        cells=f"{PANELS_HT_TABLE}/cells.parquet",
        extract=f"{PANELS_HT_TABLE}/provenance.json",
        script=f"{PANELS_HT_SCRIPTS}/panel_5b_spatial_classes.py",
        common=f"{PANELS_HT_SCRIPTS}/panel_common.py",
    output:
        panels_ht_files("5B"),
    params:
        stem=lambda wildcards: f"{PANELS_HT_OUT}/panel_5B_{wildcards.ht_release}",
        module_dir=PANELS_HT_SCRIPTS,
    log:
        f"{PANELS_HT_LOGS}/5b_{{ht_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        """
        PYTHONPATH={params.module_dir} python -u {input.script} \
            --cells {input.cells} \
            --extract-provenance {input.extract} \
            --release {wildcards.ht_release} \
            --out-stem {params.stem} \
            > {log} 2>&1
        """


# Fig. 5C: UMAP of all nuclei by cell class.
rule panels_human_tissue_5c:
    input:
        cells=f"{PANELS_HT_TABLE}/cells.parquet",
        extract=f"{PANELS_HT_TABLE}/provenance.json",
        script=f"{PANELS_HT_SCRIPTS}/panel_5c_umap_classes.py",
        common=f"{PANELS_HT_SCRIPTS}/panel_common.py",
    output:
        panels_ht_files("5C"),
    params:
        stem=lambda wildcards: f"{PANELS_HT_OUT}/panel_5C_{wildcards.ht_release}",
        module_dir=PANELS_HT_SCRIPTS,
    log:
        f"{PANELS_HT_LOGS}/5c_{{ht_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        """
        PYTHONPATH={params.module_dir} python -u {input.script} \
            --cells {input.cells} \
            --extract-provenance {input.extract} \
            --release {wildcards.ht_release} \
            --out-stem {params.stem} \
            > {log} 2>&1
        """


# Fig. S25A-D: sample table (lanes from the sample sheet), marker dot plot, QC violins.
rule panels_human_tissue_s25:
    input:
        cells=f"{PANELS_HT_TABLE}/cells.parquet",
        extract=f"{PANELS_HT_TABLE}/provenance.json",
        samples="config/samples.tsv",
        script=f"{PANELS_HT_SCRIPTS}/panel_s25_qc.py",
        common=f"{PANELS_HT_SCRIPTS}/panel_common.py",
    output:
        panels_ht_files("S25A"),
        panels_ht_files("S25B"),
        panels_ht_files("S25C"),
        panels_ht_files("S25D"),
    params:
        genes=",".join(PANELS_HT_MARKERS),
        out_dir=PANELS_HT_OUT,
        module_dir=PANELS_HT_SCRIPTS,
    log:
        f"{PANELS_HT_LOGS}/s25_{{ht_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        """
        PYTHONPATH={params.module_dir} python -u {input.script} \
            --cells {input.cells} \
            --extract-provenance {input.extract} \
            --samples {input.samples} \
            --genes {params.genes} \
            --release {wildcards.ht_release} \
            --out-dir {params.out_dir} \
            > {log} 2>&1
        """


# Fig. S26: per-class spatial maps, the class in magenta over all other nuclei.
rule panels_human_tissue_s26:
    input:
        cells=f"{PANELS_HT_TABLE}/cells.parquet",
        extract=f"{PANELS_HT_TABLE}/provenance.json",
        script=f"{PANELS_HT_SCRIPTS}/panel_s26_spatial_highlight.py",
        common=f"{PANELS_HT_SCRIPTS}/panel_common.py",
    output:
        panels_ht_files("S26"),
    params:
        stem=lambda wildcards: f"{PANELS_HT_OUT}/panel_S26_{wildcards.ht_release}",
        module_dir=PANELS_HT_SCRIPTS,
    log:
        f"{PANELS_HT_LOGS}/s26_{{ht_release}}.log",
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        """
        PYTHONPATH={params.module_dir} python -u {input.script} \
            --cells {input.cells} \
            --extract-provenance {input.extract} \
            --release {wildcards.ht_release} \
            --out-stem {params.stem} \
            > {log} 2>&1
        """
