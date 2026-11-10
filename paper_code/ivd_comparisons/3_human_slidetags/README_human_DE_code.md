# Human decidual vs fibroblast DE: code

This code reruns the human decidual-vs-fibroblast pseudobulk DE on the current human fib/dec object
(`xmens_final/data/deliverables/2026_09_15/Cevrim_Human_Only_Fib_Dec_Cells.h5ad`). Results and methods are in
`README_human_DE_rerun.md`. Primary fit: `cell_type` labels (confirmed 25 October), one pseudobulk per donor × arm,
`~ donor + group`, PyDESeq2 with the settings shared with the X-Mens and Hammoud fits: pseudobulks need at least 20 nuclei, and a gene is kept
if CPM ≥ 1 in at least half of either arm's pseudobulks, which with two donors means in either donor.

## Steps, in order

| step | script | where | output |
|---|---|---|---|
| 1 | `cluster_jobs/human_pseudobulk_job.py` | Slurm job (reads the object read-only) | `human_pb_meta.csv`, `human_pb_counts.csv.gz`, `human_pb_detect.csv.gz`, `human_var.csv`, `human_old_vs_new.json` |
| 2 | `cluster_jobs/old_object_and_mouse_labels_job.py` | Slurm job; only needed for the regression check | `old_pb_meta.csv`, `old_pb_counts.csv.gz`, `old_object_and_mouse_labels.json` |
| 3 | `fit_human_de.py` | Slurm job, release environment (`rsc_uv_py31111`) | `DE_<fit>.csv.gz` for six fits, `detection_pct_by_donor_and_arm.csv.gz`, `fit_diagnostics.json` |
| 4 | `build_human_de_tables.py` | same job | `human_decidual_vs_fibroblast_all_genes.csv` (primary), `..._coarse_refined_labels.csv` (sensitivity) |
| 5 | `compare_downstream.py` | same job | `SUMMARY_human_tables_downstream.csv`, `SUMMARY_human_table_agreement.csv`, `TABLE_marker_genes_by_version.csv`, `TABLE_consensus123_vs_human_by_version.csv` |

Diagnostics only, not needed to reproduce any table: `cluster_jobs/inspect_human_objects_job.py` (compares the four
human objects) and `cluster_jobs/relabel_check_job.py` (shows the 2 September relabel changed `cell_type` but not
`coarse_refined`).

## Commands

Steps 1 and 2 have their object paths written in the script and write into the working directory. They ran with this
header (step 1 used 4 CPUs, 48 GB and 1 hour; step 2 used 2 CPUs, 40 GB and 45 minutes):

```bash
#!/bin/bash
#SBATCH --partition=eddy
#SBATCH --account=eddy_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=60
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
export XDG_CACHE_HOME=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch
export NUMBA_CACHE_DIR=$XDG_CACHE_HOME/numba MPLCONFIGDIR=$XDG_CACHE_HOME/mpl TMPDIR=/n/lab_storage/eddy_lab/users/nhilgert/tmp
mkdir -p $XDG_CACHE_HOME $NUMBA_CACHE_DIR $MPLCONFIGDIR $TMPDIR
export HDF5_USE_FILE_LOCKING=FALSE PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python human_pseudobulk_job.py
```

Steps 3 to 5 (PB = step 1 outputs, OLDPB = step 2 outputs) run as one Slurm job in the release environment, the one
the X-Mens and Hammoud fits used. Same header as above, with 8 CPUs, 16 GB and 45 minutes:

```bash
PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111/bin/python
$PY fit_human_de.py --pb-dir PB --old-pb-dir OLDPB --old-table our_human_invivo_DE.csv --n-cpus 8 --out fits
$PY build_human_de_tables.py --pb-dir PB --fits fits --out tables
$PY compare_downstream.py --fits fits --var PB/human_var.csv --old-table our_human_invivo_DE.csv \
    --xmens xmens_broad_all_genes_with_uninduced.csv --orthology ortholog_pairs_one2one.csv \
    --consensus consensus123_stadtmauer_authors_criteria.csv --wang-dir WANG --out compare
```

Leave out `--old-pb-dir` and `--old-table` in step 3 to skip the regression fit.

## Other inputs

| file | source |
|---|---|
| `our_human_invivo_DE.csv` | the old notebook's output, `decidualization/revisions/in_vitro_genes_repo/` |
| `xmens_broad_all_genes_with_uninduced.csv` | the primary X-Mens fit of 27 October, which includes the uninduced fibroblasts (`xmens_with_uninduced_20261027.zip`, made by the mouse DE code). The 19 October release's `xmens_broad_all_genes.csv` (`release_20261019_v6_dS_vs_C7fib_luteal_tables.zip`) also works; the summary row says which was used. |
| `consensus123_stadtmauer_authors_criteria.csv` | corrected consensus list, `consensus123_stadtmauer_fix_bundle.zip` (`job_outputs/`) |
| `wang_oil.csv`, `wang_pregnancy.csv`, `wang_ivd.csv` | `cagri_ivd_end_to_end_20260815/resources/` |
| `ortholog_pairs_one2one.csv` | `xmens_final/results/mi_dpt/overlap/`: 16,275 human–mouse symbol pairs from Ensembl 116 Compara (`xmens_final/resources/orthology/ensembl_human_mouse_components.tsv`; see the README there). Keeps `ortholog_one2one`, protein-coding genes, and drops the two pairs whose mouse symbol (Tsc22d4) belongs to two different mouse genes. Re-derived from the Ensembl table on 26 October. The consensus scoring also uses Ensembl 116, but matches by Ensembl ID. |

## Environment

- **Steps 3 to 5:** the cluster's release environment `rsc_uv_py31111`: Python 3.11.11, NumPy 2.2.0, SciPy 1.17.1,
  pandas 2.3.3, PyDESeq2 0.5.4, scikit-learn 1.8.0, formulaic 1.2.2. The project environment `rsc_uv` (Python 3.11.15,
  the same versions of these packages) gives identical output. Until 29 October these steps ran outside the cluster
  with NumPy 2.4.6; see "Cluster rerun" below for what that changed.
- **Steps 1 and 2:** the cluster's `rsc_uv` environment (anndata, h5py, scipy).
- **CPUs:** `fit_human_de.py` defaults to one CPU; 1 and 8 CPUs give identical output (checked 29 October).

## Verification (26 October)

Steps 3 to 5 were rerun from the saved step 1 and 2 outputs and compared with the saved results:

- **Exactly equal:** all six fits, the detection table and `fit_diagnostics.json`.
- **Byte-identical:** both final tables, `SUMMARY_human_tables_downstream.csv`, `SUMMARY_human_table_agreement.csv` and `TABLE_consensus123_vs_human_by_version.csv`.
- **Deliberate change:** `TABLE_marker_genes_by_version.csv` now keeps full-precision p-values. The earlier copy had every value rounded to three decimals, so small p-values showed as 0 (PRL in the primary fit: 4.2 × 10⁻¹⁴). The values agree after rounding.
- **Detection counts:** they agree with an independently computed per-donor table (8 genes × 4 groups, within 0.005
  percentage points). No pseudobulk has more detecting nuclei than UMIs for any gene.
- **Results note:** every number in `README_human_DE_rerun.md` was re-derived from these outputs on 26 October.
- **Readability rewrite:** the rewritten steps 3 to 5 reproduce all 14 output files above exactly (compressed files
  compared after decompression). The four rewritten cluster jobs were rerun on the cluster and reproduce all ten of
  their output files exactly, except that `human_objects_report.json` now lists its column comparisons in sorted
  order; its content is identical.

## Cluster rerun (29 October)

Steps 3 to 5 were rerun in the release environment with the 20-nucleus minimum, from the saved step 1 and 2 outputs
(checksums match):

- **Identical across runs:** 8 vs 1 CPU, and `rsc_uv_py31111` vs `rsc_uv`, for every file the fit step writes (the table and comparison steps ran once, from those fits).
- **Regression fit (R0):** now reproduces the old notebook's table exactly: all 3,713 calls, fold changes within
  9 × 10⁻¹⁶. The 26 October run outside the cluster gave 3,711 calls and 66 genes with different p-values, so that
  mismatch came from the environment, not the code.
- **Primary fit (B):** 3,198 calls instead of 3,195. NME5, AC009320.1 and AC008403.3 move from adjusted p 0.0500–0.0501
  to 0.0496–0.0499. No call is lost and no gene changes direction; fold changes differ by at most 0.005 log2.
- **Fragile p-values:** 18 highly expressed genes in B change standard error by more than 1%, up to about 13-fold in
  either direction; none changes call. With one residual degree of freedom a gene's dispersion is barely constrained,
  so use fold changes and per-donor directions, and don't quote the two-donor p-values.
- **Unchanged:** the calls of fits C, D and E, and for fit B the consensus, X-Mens and Wang comparisons. Fit A (old
  settings) gains 10 calls and loses 5.

## Where to look when reviewing

Each script opens with a docstring saying what it does, what it reads and what it writes. The places most worth a
second pair of eyes:

- **`fit_human_de.py`:**
  - `fit_decidual_vs_fibroblast()`: the design, the contrast (fibroblast is the reference, so a positive log2FC means
    higher in decidual) and the two groups of settings;
  - `FITS`: what each fit compares;
  - `genes_passing_release_filter()`: with two donors a gene is kept if it reaches 1 CPM in either donor's pseudobulk
    of either arm. That is permissive by design: each arm is checked on its own, so genes switched on in only one arm
    aren't dropped.
- **`cluster_jobs/human_pseudobulk_job.py`:**
  - `ARM`: which label values count as decidual and as fibroblast;
  - `add_pseudobulk()`: the sums over `layers/counts` and the detection counts (nuclei with at least one UMI);
  - the old vs new comparison, which matches nuclei by sample and barcode.
- **`compare_downstream.py`:**
  - `RENAMED_CONSENSUS_GENES`: five consensus genes renamed since publication;
  - `direction_agreement()`: the chance formula;
  - `reference_lists()`: the symbol-based orthology join. Wang's gene symbols are matched as published, so a few
    renamed mouse genes drop out; every human table is affected equally.
- **`build_human_de_tables.py`, `per_donor_log2_ratios()`:** a 0.5 CPM pseudocount. The ratios are filled in only for
  genes that pass the filter.

## Changes from the session scripts

These scripts replace `fit_human_variants.py` and `compare_human_versions.py`, which only ran inside the analysis
session. The final-table step had not been saved as a script at all. The logic is unchanged, and every input and
output is now a command-line path, so each script can become one Snakemake rule. The Ulrich comparisons were removed
on 25 October. On 26 October the nuclei minimums (25 for the release's settings, 10 for the old notebook's) were made
explicit in the code. No pseudobulk is near either, and rerunning steps 3 to 5 left all 14 output files unchanged. On
29 October the release minimum was lowered to 20, the Hammoud release's cell floor, so all three fits share it (it
excludes nothing), and the `--n-cpus` help text was updated after the 8-CPU check; no computation changed.

All seven scripts were then rewritten for readability: descriptive names, one statement per line, short docstrings
and comments, constants named at the top. The computations are the same, with three small changes to the cluster
jobs:

- they check that a count layer holds integers before converting it;
- detection is counted directly as UMI > 0, rather than from the matrix's stored entries, which differs only if the
  matrix stores explicit zeros;
- the object inspection lists shared columns in sorted order, so its report no longer depends on Python's set order.

Two lines of unused code were dropped from the pseudobulk job.
