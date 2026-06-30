# Stage 00a — Slide-tags spatial coordinates (per-cell x,y)
#
# The Slide-tags mapping pipeline is vendored as a pinned submodule at
# external/slidetags (nhilgert/slide-tags @ v0.1-fasrc). It is a separate
# pipeline: own conda env (`slidetags`) + Julia project, and it submits its
# own SLURM jobs via slidetag_pipe.sh. We do NOT execute it from inside this
# DAG (that would nest cluster-job submission and tangle env systems). The
# boundary between the two pipelines is the per-sample coords CSV.
#
# Run the fork separately (see external/slidetags/FASRC_SETUP.md) to produce a
# barcode,x,y CSV per sample; this rule normalizes it into the xmens layout.

# Where the fork writes its per-sample coords. Override in config/config.yaml
# (key: slidetags_coords_src) once you confirm export_coords.R's output path.
SLIDETAGS_COORDS_SRC = config.get(
    "slidetags_coords_src",
    "external/slidetags/data/{sample}/coords/{sample}.csv",
)


rule slidetags:
    input:
        fork_coords=lambda w: SLIDETAGS_COORDS_SRC.format(sample=w.sample),
    output:
        coords="data/interim/coords/{sample}.csv",
    log:
        "logs/00a_slidetags_{sample}.log",
    script:
        "../scripts/slidetags.py"
