# X-Mens decidual vs fibroblast DE: code

This package makes the X-Mens decidual vs fibroblast differential expression and the comparisons built on it.

**Primary fit since 27 October 2026.** The uninduced sample's fibroblasts are included as a sixth timepoint, as
decided on 27 October. The fit of the release of 19 October 2026 (`release_20261019_v6_dS_vs_C7fib_luteal`), which used
the five induced timepoints only, stays available as `--design induced_only` and is reproduced byte for byte.

The release made the mouse DE with two files that also do other work: `ulrich_ivd_pipeline.py` (the fit and the
results table) and `prepare_xmens_broad()` with its helpers in the upstream `cagri_ivd_end_to_end_20260815/run_pipeline.py`
(the pseudobulks). These scripts contain only the mouse DE and its comparisons.

## Steps, in order

| step | script | where | reads | writes |
|---|---|---|---|---|
| 1 | `cluster_jobs/mouse_pseudobulk_job.py` | cluster (reads the 2 GB object) | `final_mouse.h5ad`, `xmens_feature_manifest.tsv` | `xmens_broad_carrier_counts.tsv`, `xmens_broad_metadata.tsv`, `xmens_broad_genes.tsv` |
| 2 | `fit_mouse_de.py` | anywhere with PyDESeq2 0.5.4 (the release's environment for exact numbers) | step 1's three files | `xmens_broad_all_genes.csv`, `fit_diagnostics.json` |
| 3 | `compare_xmens_with_references.py` (with `release_comparison_functions.py`) | anywhere | step 2's table, the consensus list and its orthology, Wang et al.'s three lists | support rows, the genes compared, `SUMMARY_xmens_comparisons.csv` |
| 4 | `update_ivd_workbook.py` | anywhere | the previous `IVD_comparisons_with_genes.xlsx`, step 3's folder | the new workbook |

```bash
PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111/bin/python   # the release's environment
$PY cluster_jobs/mouse_pseudobulk_job.py --out pb          # Slurm job; 8 CPUs and 64 GB were allocated
$PY fit_mouse_de.py --pb-dir pb --out fit --n-cpus 8
$PY compare_xmens_with_references.py --xmens fit/xmens_broad_all_genes.csv \
    --consensus consensus123_stadtmauer_authors_criteria.csv \
    --orthology consensus123_stadtmauer_authors_criteria_mouse_orthology.csv --wang-dir WANG --out compare
$PY update_ivd_workbook.py --previous IVD_comparisons_with_genes.xlsx --comparisons compare \
    --note "revision note" --out IVD_comparisons_with_genes_new.xlsx
```

Add `--design induced_only` to step 1 for the release's fit. Every script takes its inputs and outputs as
command-line paths, so each can become one Snakemake rule. The saved results of 27 October name step 2's table
`xmens_broad_all_genes_with_uninduced.csv`, to keep it apart from the release's file of the same format.

## The analysis

- **Nuclei:** every nucleus labelled "decidual cells" or "endometrial fibroblasts" in `coarse_refined`, so all decidual
  subtypes (decidual_apoptotic included) and all fibroblast subtypes, from the five induced timepoints, plus the
  uninduced sample's fibroblasts (9,270 nuclei). The uninduced sample's 64 decidual-labelled nuclei are left out.
- **Pseudobulks:** one per timepoint and arm, summing UMI counts over the timepoint's samples: 11 pseudobulks from
  43,133 nuclei (22,855 decidual, 20,278 fibroblast; 1,170 to 9,270 per pseudobulk). The uninduced timepoint has a
  fibroblast pseudobulk only. Each pseudobulk needs at least 20 nuclei, the cell floor of the human and Hammoud fits (25 until 29 October 2026; the smallest pseudobulk has 1,170, so the change excludes nothing).
- **Genes:** the feature manifest maps the object's 33,698 features to Ensembl 116 genes. It excludes 465 features,
  all absent from its gene annotation; each of the other 33,233 is its own gene. (The code would sum features that
  share a gene, but none do.) 64 of the genes have no symbol in the manifest, and their symbol is left empty.
  16,784 genes pass the expression filter: CPM ≥ 1 in at least half (rounded up) of either arm's pseudobulks, that is
  3 of the 5 decidual or 3 of the 6 fibroblast pseudobulks.
- **Fit:** PyDESeq2, `~ timepoint + group`, Wald test of decidual vs fibroblast (positive log2 fold change = higher in
  decidual). 11 pseudobulks and 7 coefficients leave 4 residual degrees of freedom. The uninduced pseudobulk has its own
  timepoint coefficient, which fits it exactly, so it does not enter the decidual vs fibroblast estimate directly; it
  changes the size factors and the expression filter, and the dispersions slightly. Cook's outlier refit and filtering
  on; no independent filtering, so adjusted p-values are Benjamini-Hochberg over the 16,784 fitted genes. apeglm
  shrinkage gives `log2FoldChange_shrunk`, for display only.
- **Calls:** adjusted p < 0.05 and |unshrunk log2 fold change| > 1: 4,703 genes, 2,520 up and 2,183 down.
- **Release design** (`--design induced_only`): 10 pseudobulks from 33,863 nuclei, 16,525 genes fitted (3 of 5 in either
  arm), 4,626 calls (2,504 up, 2,122 down), 4 residual degrees of freedom.
- **Comparisons (step 3):** a gene counts as called in X-Mens when it passed the filter, is estimable, and meets the
  call rule. The consensus genes are matched through strict one-to-one Ensembl 116 orthologs and Wang's in vitro list
  through current symbols or source feature IDs, both with the release's own functions (copied verbatim into
  `release_comparison_functions.py`). Wang's oil and pregnancy lists are matched by symbol among X-Mens genes with a
  unique symbol. Expected by chance is p_up(reference)·p_up(X-Mens) + p_down(reference)·p_down(X-Mens) over the genes
  compared.

## Inputs

| file | where | checked |
|---|---|---|
| `final_mouse.h5ad` | `xmens_final/data/interim/` (57,033 nuclei × 33,698 features; `layers/counts`) | SHA-256 must match the release's |
| `xmens_feature_manifest.tsv` | `cagri_ivd_end_to_end_20260815/resources/` | SHA-256 must match the release's |
| `consensus123_stadtmauer_authors_criteria.csv`, `consensus123_stadtmauer_authors_criteria_mouse_orthology.csv` | `consensus123_stadtmauer_fix_bundle.zip`, `job_outputs/` (the list corrected on 22 October) | the orthology must give 114 one-to-one genes |
| `wang_oil.csv`, `wang_pregnancy.csv`, `wang_ivd.csv` | `cagri_ivd_end_to_end_20260815/resources/` | 4,294, 5,277 and 1,501 genes |

## Environment

The release ran in `/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111`: Python 3.11.11, PyDESeq2 0.5.4,
numpy 2.2.0, pandas 2.3.3, scipy 1.17.1, anndata 0.12.11. The fit of 27 October ran there too.

## Verification

**Primary fit, 27 October 2026.** A Slurm job in the release's environment (job `f9a74573`) ran steps 1 and 2 with both
designs:
- **Release design:** the rebuilt pseudobulks (`xmens_broad_carrier_counts.tsv`, `xmens_broad_metadata.tsv`) and
  `xmens_broad_all_genes.csv` are byte-identical to the release's files.
- **With uninduced:** the table is identical on 8 CPUs and on 1 CPU.
- **Against the trial run of 24 October** (session code, another environment): its per-gene table covers the 16,525
  genes testable in the release. For all of them the calls agree, and fold changes differ by at most 0.017.

The changes from the release's fit are listed in `NOTE_xmens_with_uninduced.md` and
`TABLE_calls_changed_vs_release_v6.csv` (in `xmens_with_uninduced_20261027.zip`): no called gene changes direction.

**Steps 3 and 4, given the release's table:**
- **Step 3 support rows.** Step 3 reproduces the release's `wang_ivd_support_rows.csv` and the corrected consensus
  support rows of 22 October. Every status and text column is identical, and the numbers agree to the precision of
  the written table: the largest relative difference is 3 × 10⁻¹³, because the release built these tables from
  unrounded values in memory.
- **Step 3 gene lists.** Its four gene lists equal the workbook's four X-Mens gene sheets cell for cell, and its
  summary reproduces the workbook's four X-Mens Summary rows.
- **Step 4.** Given the release's comparisons, step 4 rebuilds the workbook of 24 October with every cell unchanged.
  The only addition is the revision note.

**Earlier check of the release design (also 27 October):** the release's own fitting code, rerun on the same
pseudobulks in the same environment, gave the byte-identical table, so that job's environment matched the
release's.

**The environment matters, not only the code.** On another machine (Python 3.11.15, numpy 2.4.6 instead of 2.2.0,
the same PyDESeq2 0.5.4), from the same pseudobulks:
- **Same genes:** both designs fit and call the same genes (4,626 and 4,703).
- **Fold changes** differ by at most 0.002 and 0.003 log2.
- **Standard errors** differ by more than 1% for a few genes (13 and 6), which moves their p-values. The largest
  changes are *Snd1* in the release design (adjusted p 1.5 × 10⁻⁷ vs 6.7 × 10⁻⁵) and *Plce1* with uninduced
  (2.3 × 10⁻³⁴ vs 8.4 × 10⁻⁵⁹).

To reproduce the reported numbers exactly, use the release's environment.

## Where to look when reviewing

- **`mouse_pseudobulk_job.py`:**
  - `DESIGNS`, `ARM`, `INDUCED_TIMEPOINTS`: which nuclei go into which pseudobulk in each design;
  - `map_features_to_genes()`: how features become genes. Excluded features are dropped and features of one gene
    would be summed (none share a gene in this object); a gene's symbol is left empty when the manifest gives none
    (64 genes) or when its features disagree (none);
  - `sum_pseudobulks()`: two sparse matrix products (nuclei into pseudobulks, features into genes), followed by a check
    that every UMI is counted exactly once.
- **`fit_mouse_de.py`:**
  - `load_pseudobulks()`: checks that every induced timepoint has one pseudobulk of each arm and the uninduced
    timepoint, if present, a fibroblast one;
  - `genes_passing_filter()`: the expression filter, applied to each arm separately;
  - `fit_decidual_vs_fibroblast()`: the design, the contrast, and the orientation of the coefficient that the shrinkage
    uses (PyDESeq2 makes "decidual" the reference level, so that coefficient's sign is flipped);
  - `build_results_table()`: genes that were not fitted stay in the table with empty statistics, a missing adjusted
    p-value counts as not significant, and `source_id` names the design.
- **`compare_xmens_with_references.py`:** `compare_wang_by_symbol()` (oil and pregnancy) and the call rule
  `CALL_RULE`. `release_comparison_functions.py` holds unedited copies of the release's functions, each marked with its
  line number in the upstream `run_pipeline.py` (SHA-256 `1c7696a2…`, the version the release recorded).
- **`update_ivd_workbook.py`:** `SHEETS` maps each gene sheet to its Summary row; the checks at the end of `main()`
  guarantee that nothing else in the workbook changes.

## Differences from the release code

The computation is the same; these are the changes made for readability.

- **Pinned numbers:** the release's config pins the expected counts (nuclei selected and per pseudobulk, genes fitted,
  calls) and the run fails if any changes. Here the two input files are pinned by checksum instead, and the outputs are
  checked by comparing them with the release's files (see Verification).
- **Pseudobulks:** step 1 builds only the pseudobulks its design uses. The release built twelve, including an
  uninduced pair, and then dropped that pair. It also read the counts matrix 4,096 nuclei at a time; step 1 reads the
  whole counts layer (about 115 million stored values) and sums with sparse matrix products. Integer sums are exact
  either way.
- **Bookkeeping:** no parquet copy of the table; no run manifest of output checksums (Snakemake can keep that record).
- **Comparison tables:** step 3 builds the release's two comparison tables with the release's own functions, but from
  the consensus list corrected on 22 October; the release used the list from before the correction.
- **Not included:** the decidual-zone outputs.
