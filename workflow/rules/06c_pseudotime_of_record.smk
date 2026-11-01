# workflow/rules/06c_pseudotime_of_record.smk
# Human fibroblast/decidual pseudotime of record (decision G: the refit of 30 Oct 2026) as rules:
#
#   fib/dec object --pt_record_latent (reintegrate.py, A40)--> stromal scVI latent of the fib/dec nuclei
#     --pt_record_pseudotime (harmonized_pseudotime.py --graph pca --cells released)--> cells.parquet, provenance.json
#
# The record was made by hand on netscratch (README.md in its folder) and is the file config
# pseudotime_of_record.human names (path, column, sha256); the MI arm of record and the deposit read it from there.
# These rules write exactly those files, so the DAG shows how they were made. They are defined only with
# config pseudotime.record.fit: true (default false), like scvi.retrain: scVI on a GPU reproduces only within one GPU
# class, and the run of record trained on an A40. The outputs are protected(); with the record folder write-protected,
# Snakemake refuses to replace them. A new fit goes to a new folder: point pseudotime_of_record.human.path and
# pseudotime.record.human.latent_dir there, run pt_record_all, then pin the new sha256.
# The run of record trained on human_stromal_fibdec_input.h5ad (1,532,822,298 bytes, builder not saved) holding the
# 47,533 nuclei of pseudotime.record.human.object; pt_record_latent reads that object directly.
# Target: pt_record_all (not in rule all).
import os

PTR = config["pseudotime"]["record"]

if PTR["fit"]:
    PTR_HUMAN = PTR["human"]
    PTR_RECORD = config["pseudotime_of_record"]["human"]["path"]
    PTR_RECORD_DIR = os.path.dirname(PTR_RECORD)
    PTR_LATENT = f"{PTR_HUMAN['latent_dir']}/latent_human_stromal_fibdec.h5ad"

    # Stromal scVI latent of the fib/dec nuclei: HVGs reselected on these nuclei, settings from config scvi (3,000
    # HVGs, 30 latent dimensions, 2 layers, batch key sample, seed 0). reintegrate.py runs at least one Leiden
    # resolution; X_scVI does not depend on it (the run of record skipped the ladder, reintegrate_diff.txt).
    rule pt_record_latent:
        input:
            h5ad=PTR_HUMAN["object"],
        output:
            h5ad=protected(f"{PTR_HUMAN['latent_dir']}/latent_{{species}}_{{lineage}}_fibdec.h5ad"),
        wildcard_constraints:
            species="human",
            lineage="stromal",
        params:
            counts_layer="counts",
            resolutions=[PTR_HUMAN["leiden_resolution"]],
        log:
            "logs/pt_record/latent_{species}_{lineage}.log",
        resources:
            mem_mb=48000,
            runtime=240,
        script:
            "../scripts/atlas/reintegrate.py"

    # Harmonized pseudotime on that latent: moscot TemporalProblem (secretory -> menstrual), CellRank RealTimeKernel
    # with a same-day PCA graph, 15 diffusion components, DPT; released nuclei. CPU, about 4 minutes. The script is an
    # ancient() input so that installing or reformatting it does not replace the record.
    rule pt_record_pseudotime:
        input:
            h5ad=PTR_HUMAN["object"],
            latent=PTR_LATENT,
            script=ancient("workflow/scripts/pseudotime/harmonized_pseudotime.py"),
        output:
            cells=protected(PTR_RECORD),
            provenance=protected(f"{PTR_RECORD_DIR}/provenance.json"),
        params:
            python=PTR["python"],
            outdir=PTR_RECORD_DIR,
            latent_sha256=PTR_HUMAN["latent_sha256"],
        log:
            "logs/pt_record/pseudotime_human.log",
        resources:
            mem_mb=96000,
            runtime=60,
        shell:
            r"""
            echo "{params.latent_sha256}  {input.latent}" | sha256sum --check - > {log} 2>&1
            {params.python} {input.script} \
                --species human \
                --graph pca \
                --cells released \
                --object {input.h5ad} \
                --latent {input.latent} \
                --out {params.outdir} \
                >> {log} 2>&1
            """

    rule pt_record_all:
        input:
            rules.pt_record_pseudotime.output,
