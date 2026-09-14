> **SUPERSEDED (2026-09-01) by `HAMMOUD_HANDOFF.md`.** This file is stale: the mouse arm HAS run
> (2026-08-30) and the combined notebook `hammoud_atlas_mapping_ALL.ipynb` is the deliverable.
> Some numbers below do not reproduce from the stored artifacts (see the drift report in the
> combined notebook, §2.1 and §2.7). Kept for context only.

# Context handoff: Hammoud fibroblast centroid mapping

## What this is

I'm a postdoc in the Eddy lab (Harvard). We have human single-**nucleus** RNA-seq of
endometrium (fibroblast/decidual cells, Slide-tags spatial) and are answering two reviewer
comments on a paper:

1. Do our data recapitulate biology the reviewer already knows (C7+ fibroblasts etc.)?
2. "What do you mean you only sequenced 2 humans?"

The plan my PI signed off on: map the centroids of **all fibroblast subtypes** from Ulrich
et al. 2024 (Hammoud lab, PNAS 121:e2404775121, full-thickness uterus, single-**cell**,
20 donors, 39 subtypes) onto our human data, then later onto our XMens (mouse) data.
A previous Claude session built and ran this. It is **done for the human arm**; the mouse
arm is not started.

## Access

- `ssh harvard-eddy` works (FASRC login node).
- Python env: `source /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/activate`
  (scanpy 1.11.5, anndata 0.12.11, torch, igraph/leidenalg, nbconvert — all present).
- Compute: SLURM, use `-p eddy`. Everything here runs in ~3 min on 8 cores / 96 GB.
  Don't run heavy jobs on the login node.
- Working dir: `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud`

## Files

```
xmens_final/hammoud/
  hammoud_fib_centroid_mapping.ipynb   <- executed, figures inline. START HERE.
  pnas.2404775121.sd05.xlsx            <- atlas centroids (sheet "5B Subtype centroids")
  pnas.2404775121.sd06.xlsx            <- atlas DE markers (sheet "6C Stro markers")
  fib_mapping/
    hammoud_fib_mapping.py             <- core module, dataset-agnostic, REUSE THIS
    run_mapping.py                     <- main pipeline, writes all cached results
    resolution_stability.py            <- mapping recomputed at every leiden resolution
    why_no_match.py                    <- marker detectability / why 3 subtypes don't map
    make_notebook.py                   <- generates the .ipynb (edit here, not the ipynb)
    submit*.sh                         <- sbatch wrappers
  fib_centroid_results/                <- 39 output files
    TABLE_correspondence_summary.csv   <- the deliverable
    fig1..fig7 *.png
```

Data objects (note: `netscratch` is scratch, could be purged):
`/n/netscratch/eddy_lab/Everyone/nhilgert/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad`
(47,533 nuclei × 36,601 genes) and the matching `XMens_Only_Fib_Dec_Cells.h5ad`
(43,197 cells) for the mouse arm.

## Results (already obtained, don't redo)

Six Hammoud fibroblast subtypes: Stro-5 RUNX3+, Stro-7 C7, Stro-8 dS1, Stro-9 dS2,
Stro-10 dS3, Stro-11 progenitor. (Stro-1/2/3 PV, Stro-4 vascular progenitor, Stro-6 uSMC
are used as negative controls — our object is fibroblast/decidual only.)

- **Confident:** Stro-7 C7 → our cluster 3 (score 0.343, margin 0.236); Stro-10 dS3 →
  cluster 1 (0.389 / 0.172). Both reciprocal, gene-bootstrap 1.00, held-out-gene
  validation 1.00, reproduced in all 6 sections. Stro-7 is reciprocal at *every*
  clustering resolution 0.3–2.0.
- **Probable:** Stro-8 dS1 → cluster 5. Stro-9 dS2 maps to the same cluster and is
  **not resolved** from dS1 (its defining genes are immediate-early/stress genes, damped
  in nuclei vs dissociated cells).
- **No counterpart:** Stro-5 (defined by immune genes; our object excludes immune cells,
  ~1% of nuclei positive) and Stro-11 (defined by cell-cycle genes; ~3% positive, MKI67
  0.006%). These are structural, not failures.
- **C7 detail:** C7 detected in 0.90% of nuclei but 16.2× enriched in cluster 3 (6.63% vs
  0.41%), which holds 58.2% of all C7+ nuclei while being 7.9% of the object.
- Negative controls band at score 0.081–0.124 and none are reciprocal; calls are
  calibrated against that band rather than a hand-picked threshold.

## Methodological decisions — please don't silently undo these

Each was made for a reason and reversing one quietly would break the result:

1. **Cluster-level centroids, NOT per-cell correlation, and NOT kNN-smoothed.** A median
   nucleus has only ~584 non-zero entries on the 3,000-gene panel (19.5%), so a per-cell
   correlation against a dense centroid mostly tracks library depth. An older notebook in
   this directory (`c7_cleaned.ipynb`) does kNN-smoothed per-cell cosine — that approach
   was deliberately dropped: it is circular with the scVI embedding it smooths over, makes
   scores autocorrelated, and borrows expression across samples. Don't reintroduce it.
2. **Headline metric is double-centered ("contrastive") cosine**; raw Spearman is shown
   only for continuity with the paper's Fig. 3A. Ours is single-nucleus, theirs is
   single-cell, so raw correlation carries a dataset-level offset — all scores bunch at
   0.36–0.58 and **raw Spearman actually gets Stro-7 wrong**, sending it to the largest
   cluster. Centering fixes it.
3. **Centroids are `log1p(mean(expm1(x)))`** — average in linear CP10K space, then log.
   Not the mean of logged values. `X` in the h5ad is verified log1p(CP10K) exactly
   (row sums 10,000).
4. **Leiden resolution is deliberately 0.5, not auto-selected.** Maximizing reciprocal
   matches picks 2.0, but for a bad reason: the extra "match" is Stro-5 at a noise-level
   0.08 while real margins collapse and Stro-7's held-out validation drops from 1.00 to
   0.44 (22 clusters split one population into twins). `resolution_stability.py` shows the
   conclusion holds at every resolution.
5. **Robustness checks resample *genes*, not cells.** Since `log1p(mean(expm1))` is
   per-gene independent, they row-index a cached centroid matrix instead of re-reducing the
   count matrix — this is validated by a unit test and is why 500 bootstraps × 6 subtypes
   is seconds, not hours.
6. **`donor` is NOT `sample`.** There are 6 samples but 2 donors: `mck5*` (menstrual,
   24,092 nuclei) and `mck6*` (secretory, 23,441). Donor and cycle phase are **perfectly
   confounded**. Six-section replication is real replication, but it cannot separate donor
   from phase, and the rebuttal says so explicitly.

## Known open items

- **Spatial depth numbers use raw y-coordinates** and need a per-timepoint rotation
  (menstrual 270°, secretory 90°, from the older `c7_cleaned.ipynb`) before the secretory
  panel is interpretable. The menstrual signal is strong regardless (cluster-3 median
  relative depth 0.17 vs 0.64 for the rest).
- **The mouse/XMens arm is not started.** The module is dataset-agnostic; the only extra
  piece needed is a human→mouse ortholog map on the gene axis (or SAMap-projected labels).
  The XMens object already has `fine_refined` labels (eFib1/2/2i/3, eFibP, dec1/2/3,
  decAp), so it needs no clustering step. There's a `SAMap` conda env available.
- Optional sharpener: the paper's Dataset S13 has phase-resolved centroids (proliferative /
  early-secretory / mid-late secretory). We only have S5 and S6 downloaded. S5B averages
  across phases, and our menstrual samples have no counterpart phase in their panel.

## How to work here

Edit `make_notebook.py` (not the `.ipynb`), regenerate, then re-execute with
`sbatch submit_nb.sh`. The notebook reads cached CSVs written by `run_mapping.py`, so
rendering is cheap and the heavy steps only rerun when you rerun that script.
