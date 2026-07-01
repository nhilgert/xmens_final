# Stage 00 - CellBender ambient-RNA removal (UPSTREAM; not in the default DAG).
#
# Documents/reproduces the ambient-RNA correction that produced the pipeline's
# inputs. DELIBERATELY NOT wired into the main DAG: the default pipeline starts
# from the already-corrected .h5 files (ingest reads xmens_inputs/), so the
# validated downstream cell counts stay fixed. CellBender is a stochastic
# variational model -- re-running changes the corrected counts and would shift
# every downstream number. This exists for completeness / provenance /
# regenerating inputs, run ONLY on explicit request, e.g.:
#
#   snakemake --profile profiles/slurm \
#     preprocessing/cellbender/1_dpi_early_mck1/1_dpi_early_mck1_cellbender.h5
#
# Parameters match Cevrim et al. Methods:
#   CellBender remove-background (v0.2.2), --total-droplets-included <N>,
#   DEFAULT number of epochs (Methods: "default number of epochs" -> no --epochs),
#   --cuda.
#
# Input is the CellRanger RAW matrix under the velocity CellRanger outputs
# (raw_feature_bc_matrix.h5 -- the full droplet set, incl. empties, that
# CellBender needs to model ambient). The lane directory == the sample's
# coords_id (mck_1, mck_8a, ...), so we map sample -> coords_id.

rule cellbender:
    input:
        raw=lambda w: config["cellbender"]["cellranger_raw"].format(
            lane=samples.loc[w.sample, "coords_id"]
        ),
    output:
        h5="preprocessing/cellbender/{sample}/{sample}_cellbender.h5",
    params:
        total_droplets=lambda w: config["cellbender"]["total_droplets_included"],
        env=lambda w: config["cellbender"]["env"],
        ld_lib=lambda w: config["cellbender"].get("ld_library_path", ""),
    log:
        "logs/00_cellbender_{sample}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.h5})"

        # CellBender lives in its own conda env. Source conda + activate
        # explicitly so this works in a bare SLURM job shell.
        source /n/sw/Miniforge3-25.3.1-0/etc/profile.d/conda.sh
        conda activate {params.env}
        if [ -n "{params.ld_lib}" ]; then
            export LD_LIBRARY_PATH="{params.ld_lib}:${{LD_LIBRARY_PATH:-}}"
        fi

        cellbender remove-background \
            --input {input.raw} \
            --output {output.h5} \
            --total-droplets-included {params.total_droplets} \
            --cuda \
            > {log} 2>&1
        """
