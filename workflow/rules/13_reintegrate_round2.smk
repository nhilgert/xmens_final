# 13_reintegrate_round2 -- round 2 of the lineage split and re-integration (2026-10-29).
#
# Four whole-tissue clusters were EXCLUDE at the atlas step in July although each is an identifiable cell
# population: mouse 5 (lymphatic endothelium), 7 (pericytes) and 19 (plasma cells), and human 30 (B cells).
# Round 2 assigns them to lineages (config annotation_round2.atlas_decisions: mouse 5 -> vascular, 7 -> smc,
# 19 -> immune; human 30 -> immune) and, for the four lineages that gain nuclei, does what stages 04 and 04b
# did in July: write the lineage subset, then fit scVI on it (reintegrate.py) and cluster at the annotation
# resolution only (lineage_spec res, leiden_scVI_3.0), not the full July ladder. Every other lineage keeps its
# July subset and fit, because its nuclei are the same under both tables.
#
# Like 04b, these rules exist only while config annotation_round2.fit is true: the fits are stochastic and
# reproduce only within one GPU class (A40, see workflow/profiles/slurm), so once made they are inputs. To refit:
#   snakemake --config 'annotation_round2={fit: true}' <round-2 latents>

R2 = config["annotation_round2"]

if R2.get("fit", False):

    import yaml as _yaml
    _R2_SPEC = _yaml.safe_load(open("resources/lineage_spec.yaml"))
    _ATLAS_SPEC = _R2_SPEC["atlas"]

    def _annotation_resolution(wildcards):
        """[3.0] for leiden_scVI_3.0: the resolution this lineage is annotated at (lineage_spec res)."""
        lineage = _R2_SPEC["species"][wildcards.species]["lineages"][wildcards.lineage]
        res = lineage.get("res", _R2_SPEC["defaults"]["res"])
        return [float(res.rsplit("_", 1)[1])]

    for _sp, _lins in R2["lineages"].items():

        rule:
            name:
                f"split_lineages_round2_{_sp}"
            input:
                atlas=f"data/interim/latent_{_sp}.h5ad",
                decisions=R2["atlas_decisions"][_sp],
                july_decisions=f"resources/annotations/{_sp}_lineage.csv",
                july_latents=[f"data/interim/latent_{_sp}_{_l}.h5ad" for _l in _lins],
            output:
                subsets=[f"{R2['outdir']}/lineages/{_sp}_{_l}.h5ad" for _l in _lins],
                report=f"{R2['outdir']}/split_{_sp}.tsv",
            params:
                species=_sp,
                lineages=_lins,
                res=_ATLAS_SPEC[_sp]["res"],           # the atlas clustering, as in 04_annotate
                slug=_ATLAS_SPEC[_sp]["labels"],
                atlas_sha256=R2["atlas_sha256"][_sp],
            log:
                f"logs/04b2_split_lineages_round2_{_sp}.log"
            resources:
                mem_mb=96000,
                runtime=60,
            script:
                "../scripts/atlas/split_lineages_round2.py"

    # GPU rule: workflow/profiles/slurm routes it to an A40, like reembed_release.
    rule reintegrate_lineage_round2:
        input:
            h5ad=f"{R2['outdir']}/lineages/{{species}}_{{lineage}}.h5ad",
        output:
            h5ad=f"{R2['outdir']}/latent_{{species}}_{{lineage}}.h5ad",
        params:
            counts_layer="counts",          # same raw-counts layer as 04b
            resolutions=_annotation_resolution,
        log:
            "logs/04b2_reintegrate_round2_{species}_{lineage}.log",
        resources:
            mem_mb=48000,
            runtime=240,
        script:
            "../scripts/atlas/reintegrate.py"
