# 16_release_round2 -- a release rebuilt with round 2 of the lineage split (04b2, 04c2).
#
# Each species' all-cells object of the base release (release_round2.from_version) loses the nuclei of the round-2
# lineages and gains the kept nuclei of their round-2 annotation (04c2 reannotate_lineage_round2), every nucleus takes
# the round-2 spatial QC (04c2 spatial_qc_round2), the whole-tissue scVI fit (integrate.py, A40, seed 0) is redone on
# these nuclei and the outlier flags recomputed; the fib/dec objects take the new embedding, flags and spatial QC.
# Every population is in the primary objects, so there are no extended objects. Outputs go to the release's own
# folder (RELEASE_BUILD_DIR, 00_release.smk) and to release-named tables; finalize_deliverables (04g) promotes them.
# reembed_release_round2 is a GPU rule (A40 only; workflow/profiles/slurm). Checks are in each script's docstring.

import sys

RR2 = config.get("release_round2")

# Defined only while the release of record is not cut (00_release.smk): a cut release is never rebuilt.
if RR2 and not RELEASE_IS_CUT:
    R2 = config["annotation_round2"]
    RR2_DIR = RELEASE_BUILD_DIR
    RR2_FROM_RELEASE = RELEASES[str(RR2["from_version"])]
    RR2_FROM = f"data/deliverables/{RR2_FROM_RELEASE}"
    RR2_SPECIES = RR2["species"]                              # species -> settings
    RR2_TAG = {s: v["tag"] for s, v in RR2_SPECIES.items()}   # species -> object tag (Human, XMens)
    RR2_SPECIES_OF_TAG = {tag: s for s, tag in RR2_TAG.items()}
    RR2_TAG_RE = "|".join(RR2_TAG.values())

    # Whole-tissue embedding of each species (release_round2.embedding). auto: reused from the previous release when
    # the species' nuclei are exactly that release's (a relabel only), re-fitted on an A40 otherwise; refit: always
    # re-fitted. Decided here, before any rule runs, so the GPU job is not scheduled at all when the embedding is
    # reused; finish_release_round2 checks the rebuilt nuclei against the previous release before reusing anything.
    sys.path.insert(0, "workflow/scripts/deliverables")
    from embedding_reuse import predicted_release_nuclei, same_nuclei

    RR2_EMBEDDING_MODE = str(RR2.get("embedding", "auto"))
    if RR2_EMBEDDING_MODE not in ("auto", "refit"):
        raise ValueError(f"config release_round2.embedding is auto or refit, not {RR2_EMBEDDING_MODE!r}")

    def _all_cells(release: str, species: str) -> str:
        return f"data/deliverables/{release}/Cevrim_{RR2_TAG[species]}_All_Cells.h5ad"

    def _embedding_plan(species: str) -> str:
        if RR2_EMBEDDING_MODE == "refit" or PREVIOUS_RELEASE is None:
            return "refit"
        lineages = [
            (f"{R2['outdir']}/latent_{species}_{lineage}.h5ad",
             f"{R2['decisions_dir']}/{species}_{lineage}.csv",
             round2_lineage_res(species, lineage))
            for lineage in R2["lineages"][species]
        ]
        predicted = predicted_release_nuclei(_all_cells(RR2_FROM_RELEASE, species), lineages)
        return "reuse" if same_nuclei(predicted, _all_cells(PREVIOUS_RELEASE, species)) else "refit"

    RR2_EMBEDDING = {species: _embedding_plan(species) for species in RR2_SPECIES}
    RR2_REFIT = [species for species, plan in RR2_EMBEDDING.items() if plan == "refit"]
    print(f"04g2 {RELEASE}: whole-tissue embedding {RR2_EMBEDDING} (previous release {PREVIOUS_RELEASE})",
          file=sys.stderr)

    def release_embedding_summary() -> str:
        """The whole-tissue embedding of each species, for the MANIFEST note."""
        return "; ".join(
            f"{species} reused from {PREVIOUS_RELEASE} (same nuclei)" if plan == "reuse"
            else f"{species} re-fitted on all nuclei (integrate.py, A40)"
            for species, plan in RR2_EMBEDDING.items()
        )

    rule rebuild_release_round2:
        input:
            base=lambda w: ancient(f"{RR2_FROM}/Cevrim_{RR2_TAG[w.species]}_All_Cells.h5ad"),
            atlas=lambda w: ancient(f"data/interim/latent_{w.species}.h5ad"),
            latents=lambda w: [
                ancient(f"{R2['outdir']}/latent_{w.species}_{l}.h5ad") for l in R2["lineages"][w.species]
            ],
            lineages=lambda w: [
                f"{R2_RELEASE_DIR}/reannotated/{w.species}_{l}.h5ad" for l in R2["lineages"][w.species]
            ],
            spatial_qc=f"{R2_RELEASE_DIR}/qc_spatial_{{species}}.h5ad",
            legacy=RR2["legacy_names"],
            vocab="resources/annotation_vocabulary.yaml",
            helpers=["workflow/scripts/atlas/apply_annotation.py"],     # imported by the script
        output:
            h5ad=f"{RR2_DIR}/rebuilt_{{species}}.h5ad",
            report=f"results/tables/round2_release_cells_{{species}}_{DELIV_REL}.tsv",
        params:
            species=lambda w: w.species,
            lineages=lambda w: R2["lineages"][w.species],
            atlas_sha256=lambda w: R2["atlas_sha256"][w.species],
            latent_sha256=lambda w: {l: R2["latent_sha256"][f"{w.species}_{l}"] for l in R2["lineages"][w.species]},
            release=DELIV_REL,
        log:
            f"logs/04g2_rebuild_release_round2_{{species}}_{DELIV_REL}.log",
        resources:
            mem_mb=128000,
            runtime=90,
        script:
            "../scripts/deliverables/rebuild_round2.py"

    # The stage-02 fit (integrate.py, config scvi/leiden settings) on all nuclei of the rebuilt object, for the
    # species whose embedding is re-fitted.
    if RR2_REFIT:

        rule reembed_release_round2:
            wildcard_constraints:
                species="|".join(RR2_REFIT),
            input:
                h5ad=f"{RR2_DIR}/rebuilt_{{species}}.h5ad",
            output:
                h5ad=f"{RR2_DIR}/reembedded_{{species}}.h5ad",
            params:
                counts_layer="counts",          # integrate.py reads snakemake.params.counts_layer
                fit=SCVI_FIT_SETTINGS,
            log:
                f"logs/04g2_reembed_release_round2_{{species}}_{DELIV_REL}.log",
            resources:
                mem_mb=96000,
                runtime=240,
            script:
                "../scripts/atlas/integrate.py"

    def _finish_embedding_input(wildcards):
        """The re-fitted object, or the previous release's object whose embedding is reused."""
        species = RR2_SPECIES_OF_TAG[wildcards.tag]
        if RR2_EMBEDDING[species] == "refit":
            return {"reembedded": f"{RR2_DIR}/reembedded_{species}.h5ad"}
        return {"previous": ancient(_all_cells(PREVIOUS_RELEASE, species))}

    rule finish_release_round2:
        wildcard_constraints:
            tag=RR2_TAG_RE,
        input:
            unpack(_finish_embedding_input),
            rebuilt=lambda w: f"{RR2_DIR}/rebuilt_{RR2_SPECIES_OF_TAG[w.tag]}.h5ad",
            base=lambda w: ancient(f"{RR2_FROM}/Cevrim_{w.tag}_All_Cells.h5ad"),
            vocab="resources/annotation_vocabulary.yaml",
            helpers=[                                                   # imported by the script
                "workflow/scripts/atlas/apply_annotation.py",
                "workflow/scripts/deliverables/embedding_outliers.py",
                "workflow/scripts/deliverables/embedding_reuse.py",
            ],
        output:
            primary=f"{RR2_DIR}/Cevrim_{{tag}}_All_Cells.h5ad",
            outliers=f"results/tables/round2_embedding_outliers_{{tag}}_{DELIV_REL}.csv",
        params:
            species=lambda w: RR2_SPECIES_OF_TAG[w.tag],
            release=DELIV_REL,
            base_release=RR2_FROM_RELEASE,
            embedding=lambda w: RR2_EMBEDDING[RR2_SPECIES_OF_TAG[w.tag]],
            previous_release=str(PREVIOUS_RELEASE),
            fit=SCVI_FIT_SETTINGS,          # a reused fit must have been made with these settings
        log:
            f"logs/04g2_finish_release_round2_{{tag}}_{DELIV_REL}.log",
        resources:
            mem_mb=128000,
            runtime=60,
        script:
            "../scripts/deliverables/finish_round2.py"

    def _refresh_previous_fibdec(wildcards):
        """With a reused embedding, the previous release's fib/dec object, against which the copy is checked."""
        species = RR2_SPECIES_OF_TAG[wildcards.tag]
        if RR2_EMBEDDING[species] == "refit":
            return []
        return [ancient(f"data/deliverables/{PREVIOUS_RELEASE}/Cevrim_{wildcards.tag}_Only_Fib_Dec_Cells.h5ad")]

    def _pseudotime_of_record(wildcards) -> dict:
        """The species' pseudotime of record (config pseudotime_of_record), or {} when the release's own stands."""
        return dict(config.get("pseudotime_of_record", {}).get(RR2_SPECIES_OF_TAG[wildcards.tag], {}))

    def _pseudotime_files(wildcards) -> list:
        record = _pseudotime_of_record(wildcards)
        return [record["path"], record["run_record"], record["script"]] if record else []

    # Decision G: a species with a pseudotime of record takes it here, so the release and its deposit copies agree.
    rule refresh_fibdec_round2:
        wildcard_constraints:
            tag=RR2_TAG_RE,
        input:
            fibdec=lambda w: ancient(f"{RR2_FROM}/Cevrim_{w.tag}_Only_Fib_Dec_Cells.h5ad"),
            primary=f"{RR2_DIR}/Cevrim_{{tag}}_All_Cells.h5ad",
            previous=_refresh_previous_fibdec,
            pseudotime=_pseudotime_files,
            vocab="resources/annotation_vocabulary.yaml",
            helpers=[                                                   # imported by the script
                "workflow/scripts/atlas/apply_annotation.py",
                "workflow/scripts/deliverables/pseudotime_of_record.py",
            ],
        output:
            h5ad=f"{RR2_DIR}/Cevrim_{{tag}}_Only_Fib_Dec_Cells.h5ad",
        params:
            species=lambda w: RR2_SPECIES_OF_TAG[w.tag],
            release=DELIV_REL,
            obsm_map=lambda w: RR2_SPECIES[RR2_SPECIES_OF_TAG[w.tag]]["fibdec_embedding"],
            embedding=lambda w: RR2_EMBEDDING[RR2_SPECIES_OF_TAG[w.tag]],
            previous_release=str(PREVIOUS_RELEASE),
            pseudotime=_pseudotime_of_record,
        log:
            f"logs/04g2_refresh_fibdec_round2_{{tag}}_{DELIV_REL}.log",
        resources:
            mem_mb=64000,
            runtime=30,
        script:
            "../scripts/deliverables/refresh_fibdec_round2.py"
