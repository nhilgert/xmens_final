# workflow/rules/23_spatial_gp.smk
# =============================================================================
# Spatial GP fit of pseudotime + gradient field (single species, native pt —
# the cross-species transfer arm is DEFERRED; see spatial_pt_gp.py docstring).
#
# config.yaml (runs of record):
#   spatial_gp:
#     python: <rsc_uv python>          # torch + gpytorch (+ geosketch)
#     runs:
#       mouse_dpt_1dpi_late:
#         h5ad: archive/pseudotime_mouse_v2_induced_no_apop/2026_08_28_mouse_fibdec_pseudotime_rooted.h5ad
#         pt_key: dpt_pseudotime       # the 28 Aug DPT (root idx 6766), not pt_harmonized
#         timepoint: 1_dpi_late
# The object is outside the pipeline: notebooks/root_selection_mouse.ipynb (cell 14) wrote its dpt_pseudotime.
# =============================================================================
import re

SGP_PY   = config["spatial_gp"]["python"]
SGP_RUNS = config["spatial_gp"]["runs"]   # run_name -> {h5ad, pt_key, timepoint}

wildcard_constraints:
    sgp_run = "|".join(re.escape(r) for r in SGP_RUNS),

rule spatial_gp:
    input:
        h5ad   = lambda wc: SGP_RUNS[wc.sgp_run]["h5ad"],
        script = "workflow/scripts/pseudotime/spatial_pt_gp.py",
    output:
        summary = "results/spatial_gp/{sgp_run}/gp_fit_summary.json",
        field   = "results/spatial_gp/{sgp_run}/gp_field.npz",
    params:
        py     = SGP_PY,
        pt_key = lambda wc: SGP_RUNS[wc.sgp_run]["pt_key"],
        tp     = lambda wc: SGP_RUNS[wc.sgp_run].get("timepoint"),
        tp_arg = lambda wc: (f"--timepoint {SGP_RUNS[wc.sgp_run]['timepoint']}"
                             if SGP_RUNS[wc.sgp_run].get("timepoint") else ""),
        outdir = lambda wc: f"results/spatial_gp/{wc.sgp_run}",
        seed   = config["seed"],           # geosketch, permutations and the GP fit
    threads: 4
    resources:
        mem_mb = 32000, runtime = 120,     # CPU exact GP at 5k sketch
    log:
        "logs/spatial_gp/{sgp_run}.log",
    shell:
        r"""
        {params.py} {input.script} \
            --h5ad {input.h5ad} --pt-key {params.pt_key} {params.tp_arg} \
            --seed {params.seed} \
            --outdir {params.outdir} > {log} 2>&1
        """

rule spatial_gp_all:
    input:
        expand("results/spatial_gp/{r}/gp_fit_summary.json", r=list(SGP_RUNS)),
