# 14_annotate_round2 -- round 2 of the lineage annotation (2026-10-30).
#
# The four lineages that 04b2 re-integrated (mouse immune, vascular and smc; human immune) are annotated as every
# lineage is in 04c: apply_annotation.py applies the committed round-2 table
# (annotation_round2.decisions_dir/{species}_{lineage}.csv) at the lineage_spec resolution and writes the lineage
# object without its EXCLUDE clusters. The round-2 fits are inputs (ancient): they reproduce only on A40s, and
# 04g2 rebuild_release_round2 checks them against annotation_round2.latent_sha256.
#
# Spatial QC (04f) is re-run on the nuclei that the round-2 atlas tables keep, so the reinstated populations take
# part in the kNN-distance statistics as every other nucleus did. The production script spatial_qc.py runs
# unchanged on a coordinates-only object; species without a 04f configuration (human) take the overrides in
# annotation_round2.spatial_qc_overrides (human: exact duplicate coordinates only, which reproduces the QC of the
# human release of 2026-09-02).

R2 = config["annotation_round2"]

import yaml as _yaml
_R2C_SPEC = _yaml.safe_load(open("resources/lineage_spec.yaml"))


def _round2_lineage_res(wildcards):
    """The resolution this lineage is annotated at (lineage_spec res, leiden_scVI_3.0)."""
    lineage = _R2C_SPEC["species"][wildcards.species]["lineages"][wildcards.lineage]
    return lineage.get("res", _R2C_SPEC["defaults"]["res"])


def _round2_spatial_qc_params(wildcards):
    """config spatial_qc with the species' round-2 overrides applied."""
    params = dict(config["spatial_qc"])
    params.update(R2.get("spatial_qc_overrides", {}).get(wildcards.species, {}))
    return params


rule reannotate_lineage_round2:
    input:
        h5ad=ancient(f"{R2['outdir']}/latent_{{species}}_{{lineage}}.h5ad"),
        csv=f"{R2['decisions_dir']}/{{species}}_{{lineage}}.csv",
        spec="resources/lineage_spec.yaml",
        vocab="resources/annotation_vocabulary.yaml",
    output:
        h5ad=f"{R2['outdir']}/reannotated/{{species}}_{{lineage}}.h5ad",
    params:
        level="fine",
        species=lambda w: w.species,
        lineage=lambda w: w.lineage,
        res=_round2_lineage_res,
        repo=config.get("repo", "."),
    log:
        "logs/04c2_reannotate_round2_{species}_{lineage}.log",
    resources:
        mem_mb=32000,
        runtime=30,
    script:
        "../scripts/atlas/apply_annotation.py"


# The atlas-kept nuclei of the round-2 atlas tables: coordinates, sample and timepoint only.
rule atlas_kept_coords_round2:
    input:
        atlas=ancient("data/interim/latent_{species}.h5ad"),
        decisions=lambda w: R2["atlas_decisions"][w.species],
        spec="resources/lineage_spec.yaml",
    output:
        h5ad=f"{R2['outdir']}/atlas_kept_coords_{{species}}.h5ad",
    params:
        species=lambda w: w.species,
        atlas_sha256=lambda w: R2["atlas_sha256"][w.species],
        sqc=_round2_spatial_qc_params,
    log:
        "logs/04c2_atlas_kept_coords_round2_{species}.log",
    resources:
        mem_mb=64000,
        runtime=30,
    script:
        "../scripts/deliverables/atlas_kept_coords.py"


rule spatial_qc_round2:
    input:
        h5ad=f"{R2['outdir']}/atlas_kept_coords_{{species}}.h5ad",
    output:
        h5ad=f"{R2['outdir']}/qc_spatial_{{species}}.h5ad",
        tally="results/tables/spatial_qc_{species}_round2.csv",
    params:
        sqc=_round2_spatial_qc_params,
    log:
        "logs/04c2_spatial_qc_round2_{species}.log",
    resources:
        mem_mb=16000,
        runtime=30,
    script:
        "../scripts/atlas/spatial_qc.py"
