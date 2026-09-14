# Verification notebooks — recompute and check every major claim

Six notebooks. Each recomputes a claim from source data and prints `[PASS]` / `[MISMATCH]` against
the value recorded during the session, ending with a `check_summary()` table.

**Run them in the `rsc_uv` environment on the cluster:**
`/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python -m jupyter lab`

| notebook | claim | key recorded values |
|---|---|---|
| `01_isodepth_construction` | isodepth is anchored in tissue geometry; which cells are FIT vs SMOOTHED | n_seeded 20,342 / n_unseeded 22,855; uninduced fit exists (4,155 cells) but was NOT seeded; induced blob rho 0.447–0.947, median 0.893 |
| `02_isodepth_vs_pseudotime` | the two coordinates agree; cross-species transfer by timepoint | human 0.671 / mouse 0.690 overall; transfer 0.591 overall, 4dpi 0.207, uninduced 0.711 |
| `03_mutual_information` | the MI estimator; batch stratification; MI is not redundant with correlation | stratified-vs-pooled mouse 0.927 / human 0.816; rho(MI,\|rho\|) 0.433/0.469; Jaccard 0.223 |
| `04_conservation_and_controls` | cross-species enrichment + the specificity control + the 37-gene gain set | fold 2.90 (isodepth) vs 2.91 (pseudotime); uninduced 2.581 vs induced 2.377; set B n=147, 2.33x, p=1.5e-24 |
| `05_samap_integration_readonly` | how the integrations were built (**no refitting**) | CP10K in `.X`, `.raw` deleted, two-sided guard; stock/vanilla; specificity peaks 1dpi_late–2dpi |
| `06_rewiring_and_drivers` | rewiring not heterochrony; what drives cross-species similarity | dpeak median 0, 33% displaced, two-axis agreement 0.732; dec HK 10.6% vs fib 0.7% |

## Errors these notebooks are designed to catch

Every one of these was an error I made in prose during the session; each now has an explicit check:

- **WNT5A is a LOSS gene** (gain −0.0223, z −1.91), not part of the cross-method convergence set.
  Its SAMap pair is `hs_WNT5A;mm_Wnt4` — a *different mouse gene* than the 1:1 `Wnt5a` used for MI.
- **Anchor validation is rho 0.447–0.947**, not 0.82–0.95. The 0.447 blob
  (`1_dpi_late_mck12b` blob 1, 867 cells) is unexplained and must not be dropped from the range.
- **dpeak tails are asymmetric**: −9 -> 17 genes, +9 -> 36. Not "36 and 36".
- **4dpi is degraded (rho 0.207 in the full object), not inverted.** The negative values came only
  from isolated single-timepoint fits, and only in the secretory reference.
- **Transfer directions**: in `coordinate_transfer.json`, `iso_to_mouse`=0.792 is human isodepth ->
  mouse scored vs mouse *pseudotime* (cross-scalar); `dpt_to_human`=0.655 is the reverse. Neither is
  the same quantity as the vanilla object's 0.389 (transferred isodepth vs native isodepth).
- **Regulator consistency**: only RUNX1 is 6/6 human and 7/7 mouse. WT1 is 2/6 human, pct 0.920 —
  do not cite it as transferable.
- **Stratified-vs-pooled MI agreement differs by species**: mouse 0.927, human 0.816 (0.727 on the
  matched 9,550-gene universe). Do not quote a single "~0.93" for both.
- **Per-run housekeeping does NOT track mapping specificity** — the compartment contrast
  (dec 10.6% vs fib 0.7%) is the durable finding, not any timepoint ordering.

## Known-open, do not quote yet

Three cross-species transfer numbers do not reconcile — 0.389 (vanilla stored column, mm->hs,
like-for-like), 0.591 (full-object per-timepoint, hs->mm), 0.688 (mean of 12 per-timepoint fits).
Different arms, coordinate vintages, and cell subsets. Resolve before putting a transfer number in
the manuscript.

## Traps in the file layout

- Mouse deliverable is `2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad` — **"XMens", not "Mouse"**.
  A `Mouse_*` glob returns empty and has killed a job.
- April-dated SAMap pickles under `samap_integrate/runs/` are **unreadable**
  (anndata `KeyError: '_adata_ref'`). Only the 2026-08-31 pickles load.
- SAMap alignment scores track dataset size — do not read timecourse trends off them.
