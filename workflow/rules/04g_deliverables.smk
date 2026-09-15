# 04g_deliverables -- finalize_deliverables: promote the canonical source
# objects into an immutable, dated, sha256-manifested release under
# data/deliverables/<release>/.
#
# This is the pipeline's OUTBOUND boundary: everything downstream of the spine
# (SAMap, SCENIC, LIANA, pt-per-timepoint, the papermill notebooks) reads
# data/deliverables/<release>/<name>.h5ad, never a hand-copied file. Cutting a
# release is deliberate: bump config[deliverables][release], update sources,
# then
#     snakemake --rerun-triggers mtime finalize_deliverables
#
# Sources are ancient(): a release is cut on purpose, never re-cut because a
# source mtime moved. Releases are immutable -- the script refuses to overwrite
# a release file whose sha256 differs. NOT in rule all.

DELIV     = config["deliverables"]
DELIV_REL = DELIV["release"]
DELIV_DIR = f"data/deliverables/{DELIV_REL}"
DELIV_SRC = dict(DELIV["sources"])            # name -> absolute source path

rule finalize_deliverables:
    input:
        sources=[ancient(p) for p in DELIV_SRC.values()],
        vocab="resources/annotation_vocabulary.yaml",
    output:
        manifest=f"{DELIV_DIR}/MANIFEST.json",
        objects=[f"{DELIV_DIR}/{name}.h5ad" for name in DELIV_SRC],
    params:
        sources=DELIV_SRC,
        release=DELIV_REL,
        note=DELIV.get("note", ""),
    resources:
        mem_mb=32000,
        runtime=120,
    log:
        f"logs/04g_finalize_deliverables_{DELIV_REL}.log",
    script:
        "../scripts/deliverables/finalize_deliverables.py"
