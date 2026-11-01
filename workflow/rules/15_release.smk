# 15_release -- finalize_deliverables: promote the canonical source
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
#
# When config release_reembed is set, one migrated all-cells object is rebuilt with
# nuclei the July whole-tissue fit never saw and re-embedded before promotion
# (reinstate_release_cells -> reembed_release -> finish_reembedded_release); its
# deliverables source is then the finished object. Cut such a release with
#     snakemake finalize_deliverables --rerun-triggers mtime --allowed-rules \
#         migrate_release_vocabulary reinstate_release_cells reembed_release \
#         finish_reembedded_release finalize_deliverables
# on an A40 node (reembed_release is a GPU rule; see workflow/profiles/slurm).
#
# When config release_rebuild is set instead (2026_10_29b, vocabulary v8), each species' all-cells object of
# the previous release loses the clusters newly EXCLUDE in its decision tables, gains the supplementary
# populations (resources/annotations/supplementary_populations.csv) as extended-only nuclei, is re-embedded on
# all of these nuclei and is written twice: the primary object (the paper's nuclei) and the extended object
# (obs in_primary); the fib/dec objects take the new whole-tissue embedding and outlier flags
# (rebuild_release_cells -> reembed_release -> finish_rebuilt_release -> refresh_fibdec_release). Cut with
#     snakemake finalize_deliverables --rerun-triggers mtime --allowed-rules \
#         rebuild_release_cells reembed_release finish_rebuilt_release refresh_fibdec_release \
#         finalize_deliverables

import sys

DELIV     = config["deliverables"]
DELIV_REL = RELEASE                           # the release of record (00_release.smk)
DELIV_DIR = RELEASE_DIR
# name -> source: the finished objects in the release's 04g2 folder
DELIV_SRC = {name: f"{RELEASE_BUILD_DIR}/{name}.h5ad" for name in DELIV["objects"]}
MIGRATE   = config["vocab_migration"]
SPLITS    = MIGRATE.get("subtype_splits", {})
REEMBED   = config.get("release_reembed")
REBUILD   = config.get("release_rebuild")
if REEMBED and REBUILD:
    raise ValueError("config: set release_reembed or release_rebuild, not both")


def release_note() -> str:
    """MANIFEST note of the release of record (config deliverables.notes) with its placeholders filled in."""
    return DELIV["notes"].get(RELEASE_VERSION, DELIV["notes"]["default"]).format(
        from_release=RR2_FROM_RELEASE,
        embedding=release_embedding_summary(),
        release=RELEASE,
    )


# Relabels the frozen source release to the current annotation vocabulary
# (resources/annotation_vocabulary.yaml legacy_aliases), then applies any
# cell-level subtype splits (config vocab_migration.subtype_splits), so a
# vocabulary change reaches a new release through this rule and
# finalize_deliverables alone.
# The vocabulary and the decision CSVs are also inputs of 04/04c/04d, so run it
# targeted, never as part of a full run (that would re-annotate against
# data/interim/latent_*.h5ad):
#     snakemake finalize_deliverables \
#         --allowed-rules migrate_release_vocabulary finalize_deliverables
rule migrate_release_vocabulary:
    input:
        objects=[ancient(f"data/deliverables/{MIGRATE['from_release']}/{name}.h5ad")
                 for name in DELIV_SRC],
        vocab="resources/annotation_vocabulary.yaml",
        # Split inputs: the clustering is ancient() (pinned by sha256 in the script),
        # the decision CSVs are not -- editing a call re-cuts the release on purpose.
        split_clusters=[ancient(s["clusters"]) for s in SPLITS.values()],
        split_decisions=[s["decisions"] for s in SPLITS.values()],
    output:
        objects=[f"{MIGRATE['outdir']}/{name}.h5ad" for name in DELIV_SRC],
        report=f"results/tables/vocab_migration_{DELIV_REL}.tsv",
    params:
        splits=SPLITS,
    log:
        f"logs/04g_migrate_release_vocabulary_{DELIV_REL}.log",
    resources:
        mem_mb=8000,   # measured 2.6 GB peak over the four 2026_09_15 objects (one at a time)
        runtime=30,    # measured 3 min
    script:
        "../scripts/deliverables/migrate_vocabulary.py"


if REEMBED:
    REEMBED_DIR = REEMBED["outdir"]
    REEMBED_BASE = f"{MIGRATE['outdir']}/{REEMBED['object']}.h5ad"

    # Labels the object's epithelium from the decision CSV and appends the nuclei of the
    # reinstated clusters (counts from the lineage object, coordinates from the spatial table);
    # drops the July embedding. Checks in the script docstring.
    rule reinstate_release_cells:
        input:
            base=REEMBED_BASE,
            lineage=ancient(REEMBED["clusters"]),      # pinned by sha256 in the script
            spatial=ancient(REEMBED["spatial"]),
            decisions=REEMBED["decisions"],
            vocab="resources/annotation_vocabulary.yaml",
        output:
            h5ad=f"{REEMBED_DIR}/reinstated_human.h5ad",
            report=f"results/tables/reinstated_cells_{DELIV_REL}.tsv",
        params:
            lineage_sha256=REEMBED["clusters_sha256"],
            cluster_column=REEMBED["cluster_column"],
            reinstated_clusters=REEMBED["reinstated_clusters"],
        log:
            f"logs/04g_reinstate_release_cells_{DELIV_REL}.log",
        resources:
            mem_mb=64000,
            runtime=60,
        script:
            "../scripts/deliverables/reinstate_cells.py"

    # The stage-02 fit (integrate.py, config scvi/leiden settings) on the rebuilt object.
    # GPU rule, A40 only (workflow/profiles/slurm set-resources).
    rule reembed_release:
        input:
            h5ad=f"{REEMBED_DIR}/reinstated_{{species}}.h5ad",
        output:
            h5ad=f"{REEMBED_DIR}/reembedded_{{species}}.h5ad",
        params:
            counts_layer="counts",          # integrate.py reads snakemake.params.counts_layer
            fit=SCVI_FIT_SETTINGS,
        log:
            f"logs/04g_reembed_release_{{species}}_{DELIV_REL}.log",
        resources:
            mem_mb=96000,
            runtime=240,
        script:
            "../scripts/atlas/integrate.py"

    # Drops the Leiden ladder, recomputes embedding_outlier on the new embedding, and checks
    # that only the epithelial split and the outlier flags changed for the release's nuclei.
    rule finish_reembedded_release:
        input:
            h5ad=f"{REEMBED_DIR}/reembedded_human.h5ad",
            base=REEMBED_BASE,
            vocab="resources/annotation_vocabulary.yaml",
        output:
            h5ad=f"{REEMBED_DIR}/{REEMBED['object']}.h5ad",
            outliers=f"results/tables/embedding_outliers_human_{DELIV_REL}.csv",
        params:
            release=DELIV_REL,
        log:
            f"logs/04g_finish_reembedded_release_{DELIV_REL}.log",
        resources:
            mem_mb=96000,
            runtime=60,
        script:
            "../scripts/deliverables/finish_reembed.py"


if REBUILD:
    RB_DIR = REBUILD["outdir"]
    RB_FROM = f"data/deliverables/{REBUILD['from_release']}"
    RB_SPECIES = REBUILD["species"]                        # species -> settings
    RB_TAG = {s: v["tag"] for s, v in RB_SPECIES.items()}  # species -> object tag (Human, XMens)
    RB_SPECIES_OF_TAG = {tag: s for s, tag in RB_TAG.items()}

    RB_SPECIES_RE = "|".join(RB_SPECIES)
    RB_TAG_RE = "|".join(RB_TAG.values())

    def _rb(species):
        return RB_SPECIES[species]

    # Drops the clusters newly EXCLUDE in v8 and appends the supplementary populations (in_primary False);
    # drops the previous embedding. Checks in the script docstring.
    rule rebuild_release_cells:
        wildcard_constraints:
            species=RB_SPECIES_RE,
        input:
            base=lambda wc: ancient(f"{RB_FROM}/Cevrim_{RB_TAG[wc.species]}_All_Cells.h5ad"),
            atlas=lambda wc: ancient(_rb(wc.species)["atlas"]),          # pinned by sha256 in the script
            lineage_decisions=lambda wc: _rb(wc.species)["lineage_decisions"],
            dropped_clusters=lambda wc: [ancient(d["clusters"]) for d in _rb(wc.species)["dropped"]],
            dropped_decisions=lambda wc: [d["decisions"] for d in _rb(wc.species)["dropped"]],
            supplementary=REBUILD["supplementary"],
            vocab="resources/annotation_vocabulary.yaml",
        output:
            h5ad=f"{RB_DIR}/rebuilt_{{species}}.h5ad",
            report=f"results/tables/rebuilt_cells_{{species}}_{DELIV_REL}.tsv",
        params:
            species=lambda wc: wc.species,
            settings=lambda wc: _rb(wc.species),
        log:
            f"logs/04g_rebuild_release_cells_{{species}}_{DELIV_REL}.log",
        resources:
            mem_mb=96000,
            runtime=60,
        script:
            "../scripts/deliverables/rebuild_cells.py"

    # The stage-02 fit (integrate.py, config scvi/leiden settings) on all nuclei of the rebuilt object.
    # GPU rule, A40 only (workflow/profiles/slurm set-resources).
    rule reembed_release:
        wildcard_constraints:
            species=RB_SPECIES_RE,
        input:
            h5ad=f"{RB_DIR}/rebuilt_{{species}}.h5ad",
        output:
            h5ad=f"{RB_DIR}/reembedded_{{species}}.h5ad",
        params:
            counts_layer="counts",          # integrate.py reads snakemake.params.counts_layer
            fit=SCVI_FIT_SETTINGS,
        log:
            f"logs/04g_reembed_release_{{species}}_{DELIV_REL}.log",
        resources:
            mem_mb=96000,
            runtime=240,
        script:
            "../scripts/atlas/integrate.py"

    # Drops the Leiden ladder, recomputes embedding_outlier over all nuclei, writes the extended and the
    # primary object and checks the primary one against the previous release.
    rule finish_rebuilt_release:
        wildcard_constraints:
            tag=RB_TAG_RE,
        input:
            h5ad=lambda wc: f"{RB_DIR}/reembedded_{RB_SPECIES_OF_TAG[wc.tag]}.h5ad",
            rebuilt=lambda wc: f"{RB_DIR}/rebuilt_{RB_SPECIES_OF_TAG[wc.tag]}.h5ad",
            base=lambda wc: ancient(f"{RB_FROM}/Cevrim_{wc.tag}_All_Cells.h5ad"),
            vocab="resources/annotation_vocabulary.yaml",
        output:
            primary=f"{RB_DIR}/Cevrim_{{tag}}_All_Cells.h5ad",
            extended=f"{RB_DIR}/Cevrim_{{tag}}_All_Cells_extended.h5ad",
            outliers=f"results/tables/embedding_outliers_{{tag}}_{DELIV_REL}.csv",
        params:
            species=lambda wc: RB_SPECIES_OF_TAG[wc.tag],
            release=DELIV_REL,
        log:
            f"logs/04g_finish_rebuilt_release_{{tag}}_{DELIV_REL}.log",
        resources:
            mem_mb=128000,
            runtime=60,
        script:
            "../scripts/deliverables/finish_rebuild.py"

    # The fib/dec object with the new whole-tissue embedding and outlier flags; nothing else changes.
    rule refresh_fibdec_release:
        wildcard_constraints:
            tag=RB_TAG_RE,
        input:
            fibdec=lambda wc: ancient(f"{RB_FROM}/Cevrim_{wc.tag}_Only_Fib_Dec_Cells.h5ad"),
            primary=f"{RB_DIR}/Cevrim_{{tag}}_All_Cells.h5ad",
            vocab="resources/annotation_vocabulary.yaml",
        output:
            h5ad=f"{RB_DIR}/Cevrim_{{tag}}_Only_Fib_Dec_Cells.h5ad",
        params:
            species=lambda wc: RB_SPECIES_OF_TAG[wc.tag],
            release=DELIV_REL,
            obsm_map=lambda wc: _rb(RB_SPECIES_OF_TAG[wc.tag])["fibdec_embedding"],
        log:
            f"logs/04g_refresh_fibdec_release_{{tag}}_{DELIV_REL}.log",
        resources:
            mem_mb=64000,
            runtime=30,
        script:
            "../scripts/deliverables/refresh_fibdec.py"


# The release objects and manifest of the release of record: built by finalize_deliverables while the release is not
# cut, plain inputs of the DAG once it is (00_release.smk).
RELEASE_FILES = [RELEASE_MANIFEST] + [f"{DELIV_DIR}/{name}.h5ad" for name in DELIV_SRC]

# Everything rule all asks for the release of record (workflow/Snakefile); later stages of this file add theirs.
RELEASE_TARGETS = list(RELEASE_FILES)

if not RELEASE_IS_CUT:

    rule finalize_deliverables:
        input:
            sources=[ancient(p) for p in DELIV_SRC.values()],
            vocab="resources/annotation_vocabulary.yaml",
            helpers=["workflow/scripts/atlas/apply_annotation.py"],     # imported by the script
        output:
            # protected(): Snakemake deletes a job's existing outputs before running it, so the sha256 guard in
            # finalize_deliverables.py never sees an older release file. Write-protected outputs make Snakemake
            # refuse the job instead (set deliverables.version and release to cut a new release).
            manifest=protected(f"{DELIV_DIR}/MANIFEST.json"),
            objects=[protected(f"{DELIV_DIR}/{name}.h5ad") for name in DELIV_SRC],
        params:
            sources=DELIV_SRC,
            release=DELIV_REL,
            note=lambda wildcards: release_note(),
        resources:
            mem_mb=32000,
            runtime=120,
        log:
            f"logs/04g_finalize_deliverables_{DELIV_REL}.log",
        script:
            "../scripts/deliverables/finalize_deliverables.py"


# Post-cut check: every object of the release of record against the annotation vocabulary (workflow/scripts/tools/check_annotations.py;
# release roles allow the legacy label columns the releases keep). Runs after the cut, and on a cut release whose
# check is missing; a failed schema check or count gate fails the job.
RELEASE_OBJECT_SPECIES = {"Human": "human", "XMens": "mouse"}      # object tag -> species
RELEASE_CHECKS = [f"results/release_checks/{RELEASE}/{name}.json" for name in DELIV["objects"]]
RELEASE_TARGETS += RELEASE_CHECKS


# The report records the check at the cut: a later vocabulary or checker edit (ancient inputs) does not re-check a
# cut release against rules made for the next one.
rule check_release_object:
    input:
        h5ad=f"{DELIV_DIR}/{{name}}.h5ad",
        manifest=RELEASE_MANIFEST,
        vocab=ancient("resources/annotation_vocabulary.yaml"),
        spec=ancient("resources/lineage_spec.yaml"),
        script=ancient("workflow/scripts/tools/check_annotations.py"),
    output:
        report=f"results/release_checks/{RELEASE}/{{name}}.json",
    params:
        python=sys.executable,
        species=lambda w: RELEASE_OBJECT_SPECIES[w.name.split("_")[1]],
        role=lambda w: "fibdec_deposit" if w.name.endswith("_Only_Fib_Dec_Cells") else "all_deposit",
    wildcard_constraints:
        name="|".join(DELIV["objects"]),
    log:
        f"logs/04g_check_release_object_{RELEASE}_{{name}}.log",
    resources:
        mem_mb=16000,
        runtime=30,
    shell:
        r"""
        {params.python} {input.script} \
            --repo . \
            --species {params.species} \
            --h5ad {input.h5ad} \
            --role {params.role} \
            --report {output.report} \
            > {log} 2>&1
        """


# Release notes of the release of record, counted from its objects, its manifest, the previous release and the
# post-cut checks (release_notes.py); they hold no typed count.
RELEASE_NOTES = f"results/release_notes/RELEASE_NOTES_{RELEASE_VERSION}.md"
RELEASE_TARGETS.append(RELEASE_NOTES)


def previous_release_objects() -> list:
    """The objects of the previous release (00_release.smk), or none for the first release."""
    if not PREVIOUS_RELEASE:
        return []
    return [ancient(f"data/deliverables/{PREVIOUS_RELEASE}/{name}.h5ad") for name in DELIV["objects"]]


rule release_notes:
    input:
        manifest=RELEASE_MANIFEST,
        objects=[f"{DELIV_DIR}/{name}.h5ad" for name in DELIV["objects"]],
        previous=previous_release_objects(),
        checks=RELEASE_CHECKS,
    output:
        notes=RELEASE_NOTES,
    params:
        release=RELEASE,
        version=RELEASE_VERSION,
        previous_release=PREVIOUS_RELEASE or "",
        round2_lineages=config["annotation_round2"]["lineages"],
        species_of_tag=RELEASE_OBJECT_SPECIES,
    log:
        f"logs/04g_release_notes_{RELEASE}.log",
    resources:
        mem_mb=32000,
        runtime=30,
    script:
        "../scripts/deliverables/release_notes.py"
