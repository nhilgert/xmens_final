# Code behind the DE and IVD workbooks (snapshot, 31 October 2026)

This folder gathers all the code behind the two workbooks on Drive, and the notebook that draws the IVD comparison bar plots.
The workbooks are now in `xmens_figures/2026-10-30_revision_panels/DESeq2_DE_and_IVD_comparisons/`; the folder was moved
there on 30 October, so comments in the scripts still give its old location, `xmens_figures/DESeq2_DE_and_IVD_comparisons/`.

| workbook | SHA-256 of the file these scripts build |
|---|---|
| `pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx` | `b6584047b5b4615a56fb2dc8c505324e41be94ae142d1ec7fe75a3003772dae5` |
| `IVD_comparisons_with_genes.xlsx` | `24b3a50898037a809f9166d9ed5383419dc3d22fcf905a9981cb906797ce9ce9` |

The DE workbook on Drive is that file. The IVD workbook on Drive was re-saved on 30 October at 13:20, apparently from Google
Sheets (SHA-256 `a17ed39d4488a967396034f6bac373734f53f644dfd0ed51feb26b76ee25de43`). Its values, formulas and cell formatting
are identical to the file above; the re-save stored the computed Summary counts and reset the width of Summary column C.

Every file is a copy of a saved Claude Science artifact or of a cluster file. `MANIFEST.csv` gives each file's source and SHA-256.
None of this code is in the Snakemake repository (`xmens_final`) yet.

## Run order

| step | folder | what it does | where it ran |
|---|---|---|---|
| 1 | `1_consensus_list/` | rebuilds the 123-gene human in vitro consensus list with the Stadtmauer criteria, and its mouse orthology (`run_rebuild_job.sh` runs `rebuild_consensus123_stadtmauer_thresholds.py`) | cluster, 22 Oct |
| 2 | `2_xmens/` | X-Mens pseudobulks, PyDESeq2 fit, comparisons with Wang et al. and the consensus list, and the X-Mens sheets of the IVD workbook | cluster, 27 Oct |
| 3 | `3_human_slidetags/` | Slide-tags pseudobulks, PyDESeq2 fits, tables and comparisons | cluster, 29 Oct |
| 4 | `4_hammoud/` | Ulrich et al. 2024 atlas pseudobulks and PyDESeq2 fits (`ulrich_ivd_pipeline.py` with the config). `jobs/record.sh` records the counts; `jobs/verify_template.sh` reruns, checks them and writes the release | cluster, 28 Oct |
| 5 | `5_workbooks/`, `run_workbooks.sh` | the two workbooks, from the tables of steps 1 to 4; then `5_workbooks/ivd_comparison_bars.ipynb` draws the IVD comparison bar plots from the IVD workbook | anywhere |

Steps 2 and 3 have their own READMEs with commands, inputs and checks. Step 1's README lists its inputs and checks; full input
paths and checksums are in `1_consensus_list/job_outputs/manifest.json`.

## Where each part of the workbooks comes from

DE workbook: all of it from step 5 (`build_de_workbook.py`), using the tables of steps 2 to 4.

IVD workbook:

| part | made by |
|---|---|
| sheets "Slide-tags vs consensus" and "Hammoud vs consensus"; Summary rows 9-10; Notes rows 12-15 | step 5, `add_human_comparisons_to_ivd_workbook.py` |
| the four X-Mens sheets; Summary rows 2-5; Notes row 11 | step 2, `update_ivd_workbook.py` |
| sheets "consensus 123 genes", "oil vs pregnancy", "in vitro vs oil" and "in vitro vs pregnancy"; Summary rows 6-8; Notes rows 2-10; the layout | code run inside the analysis session, 19-23 October (`6_not_yet_scripts/`); no script yet |

## Inputs for step 5

`bash run_workbooks.sh INPUT_DIR OUTPUT_DIR` expects these files in INPUT_DIR:

| file | saved artifact (version) | path inside |
|---|---|---|
| `xmens_broad_all_genes_with_uninduced.csv` | `xmens_with_uninduced_20261027.zip` (217436c4-0f7d-4faa-bc49-ce5599b249a6) | `xmens_with_uninduced_20261027/` |
| `human_decidual_vs_fibroblast_all_genes.csv` | `human_DE_rerun_2026_09_15_object.zip` (d83977d9-a85b-4e92-830e-a89a6fafbbd3) | `out/` |
| `hammoud_release/` | `release_20261028_v7_dS_vs_C7fib_luteal_follicular_tables.zip` (1884c383-0b99-4cfe-96c9-d1373a440c51) | the whole archive |
| `IVD_comparisons_with_genes_without_human.xlsx` | `IVD_comparisons_with_genes.xlsx`, version 8 (fb72892e-6433-4fbe-8159-e1f1ce30477e), identical to version 5 | |
| `consensus123_stadtmauer_authors_criteria.csv`, `consensus123_stadtmauer_authors_criteria_mouse_orthology.csv` | `consensus123_stadtmauer_fix_bundle.zip` (2e82e361-4e61-48f3-8abe-79572352cc8f) | `consensus123_stadtmauer_fix/job_outputs/` |

The Hammoud tables are also on the cluster, byte-identical, in
`/n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular/`.
There the pseudobulk tables are in `intermediate/` rather than `intermediate_pseudobulk_tables/`, and the config is not in the
folder (use `4_hammoud/config_ulrich_dS_vs_C7fib_luteal_follicular.yaml`).

## Software

- All fits (steps 2 to 4) ran on the cluster in `/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111`: Python 3.11.11,
  NumPy 2.2.0, PyDESeq2 0.5.4. The environment of each pseudobulk job is in its step's README.
- Step 5 needs only pandas, numpy, openpyxl and PyYAML; the bar-plot notebook needs pandas, openpyxl and matplotlib.

## Code these steps call that is not in this folder

- Step 4 imports `run_pipeline.py` from `/n/lab_storage/eddy_lab/users/nhilgert/cagri_ivd_end_to_end_20260815`, a local git
  repository with no remote (last commit 45755c0, 1 September 2026), after checking its SHA-256 against the value in the config.
  It also reads that repository's `resources/` folder and two objects:
  `/n/netscratch/eddy_lab/Everyone/nhilgert/ulrich2024/ulrich2024_20donor_stromal.h5ad` and
  `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/final_mouse.h5ad`.
- Step 1 calls `build_consensus.py` (in `human_ivd_repro_final_defense_20260812`) and `derive_consensus123_orthology.py` (in
  `cagri_ivd_end_to_end_20260815`); their paths and checksums are in `1_consensus_list/job_outputs/manifest.json`.
- Steps 2 and 3 read the X-Mens and Slide-tags objects and the Wang et al. lists; their READMEs give the paths.

## Checks run on 29 October 2026 for this snapshot

- `run_workbooks.sh`, on inputs taken from the saved artifacts above, reproduces both Drive workbooks cell for cell, with
  identical formatting.
- `2_xmens/update_ivd_workbook.py`, applied to version 4 of the IVD workbook, reproduces version 5 cell for cell.
- `4_hammoud/`: the pipeline script and config match the SHA-256 values the release recorded, and the feature manifest and
  mapping table match the values pinned in the config. All 75 of the release's validation checks pass, and the release's
  tables on the cluster are byte-identical to the saved archive.
- `4_hammoud/resources/ulrich_feature_manifest.tsv` had been kept only in a temporary job folder on the cluster; this folder
  is now its saved copy.

`verify_bundle.py` reruns the first three checks and prints PASS or FAIL for each (8 checks; all passed on 29 October):

```
python verify_bundle.py --inputs INPUT_DIR --reference-dir DRIVE_COPIES \
    --ivd-v4 IVD_v4.xlsx --ivd-v5 IVD_v5.xlsx --xmens-comparisons XMENS_COMPARISONS_DIR
```

- `--inputs`: the step-5 input folder. `--reference-dir`: a folder holding the two workbooks from Drive.
- `--ivd-v4` and `--ivd-v5`: versions 4 (cac1e20c-1a14-4f2c-8286-39f204aedfba) and 5
  (87fe6ea1-ef1e-408c-a16e-8e4326e1695b) of the `IVD_comparisons_with_genes.xlsx` artifact.
- `--xmens-comparisons`: `xmens_comparisons_with_uninduced.zip` (c63eb41a-aa06-4285-9416-bee528863c4e), unpacked.
- As a test that it can fail: given the 28 October IVD workbook as the reference, it reports the 49 cells the
  29 October Slide-tags rerun changed (48 adjusted p-values and one fold change).

The cluster comparison of the release tables was a separate command; it is not in `verify_bundle.py`.
Earlier checks are in each step's README.

## Checks run on 31 October 2026 for this snapshot

- Every file except this README, `MANIFEST.csv` and the notebook is byte-identical to the 29 October snapshot
  (`DE_IVD_code_2026-10-29.zip`).
- `verify_bundle.py`, with the two workbooks as they are on Drive today as the reference: 7 of 8 checks pass. The DE
  workbook passes. The IVD workbook check finds 0 of 39,011 cells different and the same fonts, fills, number formats,
  wrapping, filters and frozen panes, and fails only because the Drive file is no longer byte-identical to the 29 October
  file (the re-save described at the top).
- `5_workbooks/ivd_comparison_bars.ipynb`, run on the IVD workbook from Drive, writes the same `bar_values.csv` and
  byte-identical SVGs as on 30 October. Those figures are on Drive in
  `xmens_figures/2026-10-30_revision_panels/IVD_comparison_bars/`.

## Not yet in a script

`6_not_yet_scripts/` holds the code, as recorded by the analysis session, that built the Wang comparison workbook (19 October)
and the first four versions of the IVD workbook (20 to 23 October). It is there for reading: it points at files in the
Claude Science store and does not run elsewhere. The sheets it made have not changed since. Turning it into a script is the
one step left before the whole IVD workbook can be rebuilt from the fit tables.
