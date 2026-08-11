# workflow/rules/10_samap.smk        (supersedes 10_samap_transfer.smk)
# =============================================================================
# SAMap arm, now HOOKED UP to the annotation pipeline:
#   refined_{human,mouse}.h5ad --(prep_samap_input)--> inputs/h5ads_samap/*.h5ad
#     --(run_samap)--> runs/<run>/samap.pkl --(samap_label_stain)--> results/.../*.png
# Env is STOCK sc-samap v3. Everything hangs off SAMAP_ROOT; rules run from the xmens_final root.
#
# ── The bridge (prep_samap_input) ────────────────────────────────────────────
# The old inputs/h5ads_samap/*.h5ad were built (Apr) from a pre-July annotation and are STALE.
# prep_samap_input regenerates them from your CURRENT refined_* objects: subset fib/dec on the
# coarse cell-type column, rename to the shared schema (cell_type_annotations /
# cell_subtype_annotations), drop counts/.raw (log data already in .X), split by timepoint.
# The two refined objects use different column names -> _REFINED below (config, not hardcoding).
#
# ── refined_* is ancient() ON PURPOSE ────────────────────────────────────────
# It isolates the SAMap arm from the annotation DAG: `snakemake samap_all` will NOT pull in or
# rebuild any 04* rule (important while the annotation-side config paths are still un-migrated).
# Trade-off: a future re-annotation won't auto-cascade here — re-trigger by deleting the stale
# SAMap input h5ad (or --touch refined_*). Fine, since annotation is done.
#
# ── First run: delete the stale artifacts so the DAG rebuilds them fresh ──────
#   rm $SAMAP_ROOT/inputs/h5ads_samap/{hs_fib_dec_menstrual,hs_fib_dec_secretory,mm_fib_dec_4_dpi_bleeding,mm_fib_dec_1_dpi_late}.h5ad
#   rm $SAMAP_ROOT/runs/{menstrual_x_4dpi_hmmer,secretory_x_1dpi_late_hmmer}/samap.pkl
# Then `snakemake samap_all` runs: 4x prep -> 2x run_samap -> 2x stain. (The old v3-unreadable
# pkls MUST go regardless — anndata 0.12 can't open them.)
#
# ENV — config[samap][python] = stock-v3 venv (.../envs/samap/.venv/bin/python), ABSOLUTE.
# run_samap.py + prep_samap_input.py + transfer_label_staining.py live in workflow/scripts/.
#
# config.yaml:
#   samap:
#     root:        /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate
#     coords_root: /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/crossspecies_radial_gp
#     python:      /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/envs/samap/.venv/bin/python
#     n_iter:      3
#     pairs: ["menstrual_x_4dpi_hmmer", "secretory_x_1dpi_late_hmmer"]
# Snakefile:  include: "rules/10_samap.smk"
# =============================================================================
import csv, re

SAMAP_ROOT  = config["samap"]["root"]
COORDS_ROOT = config["samap"]["coords_root"]
SAMAP_PY    = config["samap"]["python"]
N_ITER      = config["samap"].get("n_iter", 3)
SAMAP_PAIRS = config["samap"]["pairs"]
PT_H5AD     = config["samap"].get("pt_h5ad", {})   # run_name -> mouse per-timepoint pseudotime h5ad (optional)

# manifest = source of truth for run_name -> {hs_h5ad, mm_h5ad, homology}
_manifest = {}
with open(f"{SAMAP_ROOT}/slurm/pairs.csv", newline="") as fh:
    for row in csv.DictReader(fh):
        _manifest[row["run_name"].strip()] = {k: (v.strip() if isinstance(v, str) else v)
                                              for k, v in row.items()}
_missing = [r for r in SAMAP_PAIRS if r not in _manifest]
assert not _missing, f"run_name(s) not in pairs.csv: {_missing}"
assert all(_manifest[r]["homology"] == "hmmer" for r in SAMAP_PAIRS), \
    "paper deliverables are scoped to hmmer; a requested pair is blast"
PAIRS = {r: _manifest[r] for r in SAMAP_PAIRS}

# per-species source object + canonical annotation columns + expected fib/dec count
_REFINED = {
    "hs": dict(path="data/interim/refined_human.h5ad", coarse="cell_type", fine="cell_subtype", expect=47533),
    "mm": dict(path="data/interim/refined_mouse.h5ad", coarse="cell_type", fine="cell_subtype", expect=43197),
}

# the SAMap input basenames the two pairings need, e.g. {hs_fib_dec_secretory, mm_fib_dec_4_dpi_bleeding, ...}
_NEEDED_INPUTS = set()
for _r in PAIRS:
    _NEEDED_INPUTS |= {PAIRS[_r]["hs_h5ad"], PAIRS[_r]["mm_h5ad"]}

def _parse_input_basename(bn):
    # "<sp>_fib_dec_<timepoint>" -> (sp, timepoint)   (timepoint slug taken from pairs.csv, NOT run_name)
    for sp in ("hs", "mm"):
        pre = f"{sp}_fib_dec_"
        if bn.startswith(pre):
            return sp, bn[len(pre):]
    raise ValueError(f"unexpected samap input basename: {bn}")

def _phase_stage(run):
    assert run.endswith("_hmmer"), run
    phase, _, stage = run[:-len("_hmmer")].partition("_x_")
    return phase, stage

def _hs_coords(run): return f"{COORDS_ROOT}/isodepth_{_phase_stage(run)[0]}.npz"
def _mm_coords(run): return f"{COORDS_ROOT}/transfer_all/radius_{_phase_stage(run)[1]}.npz"
def _maps_dir(run):  return f"{SAMAP_ROOT}/homology/hsmm_{PAIRS[run]['homology']}"
def _hs_input(run):  return f"{SAMAP_ROOT}/inputs/h5ads_samap/{PAIRS[run]['hs_h5ad']}.h5ad"
def _mm_input(run):  return f"{SAMAP_ROOT}/inputs/h5ads_samap/{PAIRS[run]['mm_h5ad']}.h5ad"

wildcard_constraints:
    run      = "|".join(re.escape(r) for r in PAIRS),
    basename = "|".join(re.escape(b) for b in _NEEDED_INPUTS),

# -----------------------------------------------------------------------------
# BRIDGE — refined_{species} -> one per-timepoint SAMap input. Cheap.
rule prep_samap_input:
    input:
        # ancient() -> isolates samap_all from the annotation DAG (no 04* rebuilds). refined_*
        # is used as-is; delete a stale output to force regen after a re-annotation.
        refined = lambda wc: ancient(_REFINED[_parse_input_basename(wc.basename)[0]]["path"]),
        script  = "workflow/scripts/prep_samap_input.py",   # non-ancient: cheap, re-run on edit
    output:
        h5ad = f"{SAMAP_ROOT}/inputs/h5ads_samap/{{basename}}.h5ad",
    params:
        py      = SAMAP_PY,
        species = lambda wc: _parse_input_basename(wc.basename)[0],
        tp      = lambda wc: _parse_input_basename(wc.basename)[1],
        coarse  = lambda wc: _REFINED[_parse_input_basename(wc.basename)[0]]["coarse"],
        fine    = lambda wc: _REFINED[_parse_input_basename(wc.basename)[0]]["fine"],
        expect  = lambda wc: _REFINED[_parse_input_basename(wc.basename)[0]]["expect"],
    log:
        "logs/prep_samap_input/{basename}.log",
    shell:
        r"""
        module load gcc/14.2.0-fasrc01
        {params.py} {input.script} \
            --refined       {input.refined} \
            --species       {params.species} \
            --coarse-col    {params.coarse} \
            --fine-col      {params.fine} \
            --timepoint     {params.tp} \
            --expect-fibdec {params.expect} \
            --out           {output.h5ad} \
            > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# FIT — fresh runs/<run>/samap.pkl (stock v3). hs/mm inputs are now PRODUCED by prep_samap_input
# (no longer ancient external files), so SAMap sits downstream of the annotation.
rule run_samap:
    input:
        hs         = lambda wc: _hs_input(wc.run),               # <- prep_samap_input output
        mm         = lambda wc: _mm_input(wc.run),               # <- prep_samap_input output
        map_hm     = lambda wc: ancient(f"{_maps_dir(wc.run)}/hsmm/hs_to_mm.txt"),
        map_mh     = lambda wc: ancient(f"{_maps_dir(wc.run)}/hsmm/mm_to_hs.txt"),
        fit_script = ancient("workflow/scripts/run_samap.py"),  # ancient: fit is expensive, don't re-fit on edit
    output:
        pkl = f"{SAMAP_ROOT}/runs/{{run}}/samap.pkl",           # meta.json + umap left undeclared (crash trail / best-effort)
    params:
        py       = SAMAP_PY,
        maps_dir = lambda wc: _maps_dir(wc.run),
        out_dir  = lambda wc: f"{SAMAP_ROOT}/runs/{wc.run}",
        n_iter   = N_ITER,
    threads: 8
    resources:
        mem_mb  = 128000,   # TUNE to eddy; first fit also pays a one-time numba JIT cost
        runtime = 360,      # TUNE: minutes
    log:
        f"{SAMAP_ROOT}/slurm/logs/run_samap_{{run}}.log",
    shell:
        r"""
        module load gcc/14.2.0-fasrc01
        {params.py} {input.fit_script} \
            --hs-h5ad  {input.hs} \
            --mm-h5ad  {input.mm} \
            --maps-dir {params.maps_dir} \
            --out-dir  {params.out_dir} \
            --n-iter   {params.n_iter} \
            --ncpus    {threads} \
            > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# CONSUMER — coarse dec/fib label-transfer staining. Reads the fresh pkl.
rule samap_label_stain:
    input:
        pkl       = f"{SAMAP_ROOT}/runs/{{run}}/samap.pkl",
        hs_coords = lambda wc: ancient(_hs_coords(wc.run)),     # external isodepth/radius npz, frozen
        mm_coords = lambda wc: ancient(_mm_coords(wc.run)),
        script    = "workflow/scripts/transfer_label_staining.py",
    output:
        png   = "results/figures/samap_transfer/{run}_labelstain.png",
        stats = "results/tables/samap_transfer/{run}_labelstain.json",
    params:
        py  = SAMAP_PY,     # same stock-v3 env (unpickling the pkl needs samap importable)
        run = "{run}",
    log:
        "logs/samap_label_stain/{run}.log",
    shell:
        # If pickle.load fails on the v3 pkl, switch the loader in transfer_label_staining.py to
        # samap's load_samap. Verify on the first fit.
        r"""
        module load gcc/14.2.0-fasrc01
        {params.py} {input.script} \
            --run-name  {params.run} \
            --samap-pkl {input.pkl} \
            --hs-coords {input.hs_coords} \
            --mm-coords {input.mm_coords} \
            --out-png   {output.png} \
            --out-stats {output.stats} \
            > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# DELIVERABLE OBJECT — the joint post-integration object with transferred labels attached,
# for a collaborator to open and plot (joint UMAP for both species, or per-species spatial).
# Same transfer computation as the stain (make_transferred_object.py imports it). Pseudotime
# later just adds another obs column here — no new rule.
rule samap_transferred_object:
    input:
        pkl       = f"{SAMAP_ROOT}/runs/{{run}}/samap.pkl",
        hs_coords = lambda wc: ancient(_hs_coords(wc.run)),
        mm_coords = lambda wc: ancient(_mm_coords(wc.run)),
        pt        = lambda wc: ancient(PT_H5AD[wc.run]),             # mouse per-timepoint pseudotime object
        script    = "workflow/scripts/make_transferred_object.py",
        core      = "workflow/scripts/transfer_label_staining.py",   # imported by the script
    output:
        h5ad = "results/samap_transfer/{run}_transferred.h5ad",
    params:
        py  = SAMAP_PY,     # stock-v3 env (unpickling the pkl needs samap importable)
        run = "{run}",
    log:
        "logs/samap_transferred_object/{run}.log",
    shell:
        r"""
        module load gcc/14.2.0-fasrc01
        {params.py} {input.script} \
            --samap-pkl     {input.pkl} \
            --run-name      {params.run} \
            --hs-coords     {input.hs_coords} \
            --mm-coords     {input.mm_coords} \
            --transfer-fine \
            --pt-h5ad       {input.pt} \
            --out           {output.h5ad} \
            > {log} 2>&1
        """

# -----------------------------------------------------------------------------
# Paper target. Full fresh chain: 4x prep_samap_input -> 2x run_samap -> 2x {stain, object}.
rule samap_all:
    input:
        expand("results/figures/samap_transfer/{run}_labelstain.png", run=list(PAIRS)),
        expand("results/tables/samap_transfer/{run}_labelstain.json", run=list(PAIRS)),
        expand("results/samap_transfer/{run}_transferred.h5ad", run=list(PAIRS)),
