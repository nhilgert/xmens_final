# 04g2_release_round2 -- release 2026_10_30 (vocabulary v9): round 2 of the lineage split reaches the deliverables.
#
# Each species' all-cells object of from_release loses the nuclei of the round-2 lineages and gains the kept nuclei
# of their round-2 annotation (04c2 reannotate_lineage_round2), every nucleus takes the round-2 spatial QC
# (04c2 spatial_qc_round2), the whole-tissue scVI fit (integrate.py, A40, seed 0) is redone on these nuclei and the
# outlier flags recomputed; the fib/dec objects take the new embedding, flags and spatial QC. Every population is in
# the primary objects, so there are no extended objects. Cut with
#     snakemake finalize_deliverables --rerun-triggers mtime --allowed-rules \
#         reannotate_lineage_round2 atlas_kept_coords_round2 spatial_qc_round2 rebuild_release_round2 \
#         reembed_release_round2 finish_release_round2 refresh_fibdec_round2 finalize_deliverables
# reembed_release_round2 is a GPU rule (A40 only; profiles/slurm). Checks are in each script's docstring.

RR2 = config.get("release_round2")

if RR2:
    R2 = config["annotation_round2"]
    RR2_DIR = RR2["outdir"]
    RR2_FROM = f"data/deliverables/{RR2['from_release']}"
    RR2_SPECIES = RR2["species"]                              # species -> settings
    RR2_TAG = {s: v["tag"] for s, v in RR2_SPECIES.items()}   # species -> object tag (Human, XMens)
    RR2_SPECIES_OF_TAG = {tag: s for s, tag in RR2_TAG.items()}
    RR2_TAG_RE = "|".join(RR2_TAG.values())

    rule rebuild_release_round2:
        input:
            base=lambda w: ancient(f"{RR2_FROM}/Cevrim_{RR2_TAG[w.species]}_All_Cells.h5ad"),
            atlas=lambda w: ancient(f"data/interim/latent_{w.species}.h5ad"),
            latents=lambda w: [ancient(f"{R2['outdir']}/latent_{w.species}_{l}.h5ad") for l in R2["lineages"][w.species]],
            lineages=lambda w: [f"{R2['outdir']}/reannotated/{w.species}_{l}.h5ad" for l in R2["lineages"][w.species]],
            spatial_qc=f"{R2['outdir']}/qc_spatial_{{species}}.h5ad",
            legacy=RR2["legacy_names"],
            vocab="resources/annotation_vocabulary.yaml",
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

    # The stage-02 fit (integrate.py, config scvi/leiden settings) on all nuclei of the rebuilt object.
    rule reembed_release_round2:
        input:
            h5ad=f"{RR2_DIR}/rebuilt_{{species}}.h5ad",
        output:
            h5ad=f"{RR2_DIR}/reembedded_{{species}}.h5ad",
        params:
            counts_layer="counts",          # integrate.py reads snakemake.params.counts_layer
        log:
            f"logs/04g2_reembed_release_round2_{{species}}_{DELIV_REL}.log",
        resources:
            mem_mb=96000,
            runtime=240,
        script:
            "../scripts/atlas/integrate.py"

    rule finish_release_round2:
        wildcard_constraints:
            tag=RR2_TAG_RE,
        input:
            h5ad=lambda w: f"{RR2_DIR}/reembedded_{RR2_SPECIES_OF_TAG[w.tag]}.h5ad",
            rebuilt=lambda w: f"{RR2_DIR}/rebuilt_{RR2_SPECIES_OF_TAG[w.tag]}.h5ad",
            base=lambda w: ancient(f"{RR2_FROM}/Cevrim_{w.tag}_All_Cells.h5ad"),
            vocab="resources/annotation_vocabulary.yaml",
        output:
            primary=f"{RR2_DIR}/Cevrim_{{tag}}_All_Cells.h5ad",
            outliers=f"results/tables/round2_embedding_outliers_{{tag}}_{DELIV_REL}.csv",
        params:
            species=lambda w: RR2_SPECIES_OF_TAG[w.tag],
            release=DELIV_REL,
        log:
            f"logs/04g2_finish_release_round2_{{tag}}_{DELIV_REL}.log",
        resources:
            mem_mb=128000,
            runtime=60,
        script:
            "../scripts/deliverables/finish_round2.py"

    rule refresh_fibdec_round2:
        wildcard_constraints:
            tag=RR2_TAG_RE,
        input:
            fibdec=lambda w: ancient(f"{RR2_FROM}/Cevrim_{w.tag}_Only_Fib_Dec_Cells.h5ad"),
            primary=f"{RR2_DIR}/Cevrim_{{tag}}_All_Cells.h5ad",
            vocab="resources/annotation_vocabulary.yaml",
        output:
            h5ad=f"{RR2_DIR}/Cevrim_{{tag}}_Only_Fib_Dec_Cells.h5ad",
        params:
            species=lambda w: RR2_SPECIES_OF_TAG[w.tag],
            release=DELIV_REL,
            obsm_map=lambda w: RR2_SPECIES[RR2_SPECIES_OF_TAG[w.tag]]["fibdec_embedding"],
        log:
            f"logs/04g2_refresh_fibdec_round2_{{tag}}_{DELIV_REL}.log",
        resources:
            mem_mb=64000,
            runtime=30,
        script:
            "../scripts/deliverables/refresh_fibdec_round2.py"
