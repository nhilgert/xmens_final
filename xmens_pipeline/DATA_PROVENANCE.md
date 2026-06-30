# Data provenance — `xmens_final`

Records where the X-Mens / decidualization Slide-tags data comes from, how it was
produced, and the decisions made while wiring it into this repo. The goal is that a
future session (or another person) can run the pipeline without re-discovering any of
this. Paper: Çevrim et al. 2025, *Induction of menstruation in mice reveals the
regulation of menstrual shedding*, bioRxiv 2025.10.08.681007.

---

## 1. Sample inventory (15 samples)

9 mouse (GsD X-Mens induction timecourse) + 6 human (endometrial phases). Replicates
share a timepoint. `sample_id` is the canonical name used throughout the repo;
`coords_id` (= the original "mck" lane ID) is the join key to the coordinate files.

| sample_id              | species | timepoint        | coords_id | n_cells (raw CellBender) |
|------------------------|---------|------------------|-----------|--------------------------|
| uninduced_mck8a        | mouse   | uninduced        | mck_8a    | 15322 |
| uninduced_mck8b        | mouse   | uninduced        | mck_8b    | 11625 |
| 1_dpi_early_mck1       | mouse   | 1_dpi_early      | mck_1     | 3938  |
| 1_dpi_late_mck12a      | mouse   | 1_dpi_late       | mck_12a   | 9325  |
| 1_dpi_late_mck12b      | mouse   | 1_dpi_late       | mck_12b   | 8304  |
| 2_dpi_mck15            | mouse   | 2_dpi            | mck_15    | 7390  |
| 3_dpi_mck10a           | mouse   | 3_dpi            | mck_10a   | 12300 |
| 3_dpi_mck10b           | mouse   | 3_dpi            | mck_10b   | 10786 |
| 4_dpi_bleeding_mck11a  | mouse   | 4_dpi_bleeding   | mck_11a   | 11466 |
| secretory_mck6a        | human   | secretory        | mck_6a    | 17958 |
| secretory_mck6b        | human   | secretory        | mck_6b    | 17464 |
| secretory_mck6c        | human   | secretory        | mck_6c    | 10024 |
| menstrual_mck5a        | human   | menstrual        | mck_5a    | 17554 |
| menstrual_mck5b        | human   | menstrual        | mck_5b    | 18202 |
| menstrual_mck5c        | human   | menstrual        | mck_5c    | 16268 |

**Canonical `timepoint` strings (use these EXACTLY):** `uninduced`, `1_dpi_early`,
`1_dpi_late`, `2_dpi`, `3_dpi`, `4_dpi_bleeding` (mouse); `secretory`, `menstrual`
(human). All underscores. Note `4_dpi_bleeding` (not `4dpi_menstrual`).

**Notes / caveats:**
- `mck_5a` shares its GEX index (`SI-TT-A5`) with `mck_1` ("duplicate RNA index" in the
  McKinley sample sheet). Different SB indices; separated into distinct objects. Minor
  asterisk on the menstrual_mck5a sample if anything looks off downstream.
- `mck_12b`'s `cellbender_with_metadata.h5ad` in the backup is **truncated/corrupt**, but
  its `cellbender_output_filtered.h5` is intact (193 MB) — which is what we ingest from,
  so 12b is fine.
- `species` and `timepoint` are NOT uniformly present in the source objects (human
  objects have `timepoint = NaN`). Therefore **`samples.tsv` is the authoritative source**
  for species/timepoint; ingest writes them into `obs`.

---

## 2. Inputs: count matrices (CellBender output)

**Use:** `cellbender_output_filtered.h5` per sample (raw CellBender output, no metadata).
Read with `scanpy.read_10x_h5(...)` then `var_names_make_unique()` (gene symbols have
duplicates — emits a warning otherwise).

- `X` = CellBender-corrected counts (uint32). A `raw` layer in the `_with_metadata`
  variants is just a copy of `X`; the bare `.h5` has only the matrix.
- Barcodes carry a `-1` suffix (e.g. `GAGTTGTTCGCAACAT-1`).

**Durable location (copied from backup, verified to open):**
```
/n/eddy_lab/users/nhilgert/decidualization/xmens_inputs/{sample_id}.h5
```
Copied from (do NOT point the pipeline here — backup may be cleaned):
```
/n/eddy_lab/Lab/mckinley_backup_2025-12-18/data/andy_raw_data/{coords_id}/cellbender_output_filtered.h5
```
All 15 copies verified (open + cell-count match) on copy.

**Why ingest from the `.h5` and not Andy's `cellbender_with_metadata.h5ad`:**
uniform/reproducible (every sample built the same way by this pipeline), fixes the human
`timepoint=NaN` gap and the 12b corruption in one stroke, and makes metadata the
pipeline's responsibility (sourced from `samples.tsv`). The `_with_metadata.h5ad` files
are kept only as a *validation reference* (spot-check ingest reproduces them).

---

## 3. Inputs: spatial coordinates

There are **two coordinate sources** for the same underlying positions. Use Source 1.

### Source 1 — USE THIS (barcode-keyed, join-ready)
```
/n/eddy_lab/users/nhilgert/Slide-tags/coords/{coords_id}_coords.csv
```
- Produced by `export_coords.R` (our export from the slide-tags Seurat objects).
- Keyed by **barcode** (with `-1` suffix) → joins directly to the CellBender `.h5`
  (verified: e.g. mck_1 → 3590 overlapping barcodes, 0 after stripping suffix, i.e.
  suffixes already match — join as-is, NO suffix normalization).
- Columns include `barcode, sample, x_um, y_um, dbscan_1, dbscan_2, ...`.
- **Coordinate columns to use:** `dbscan_1` → x, `dbscan_2` → y. These are the
  DBSCAN-placed coords (`x_um_dbscan`/`y_um_dbscan`), a strict **superset** of the
  stricter `x_um`/`y_um` (e.g. mck_10a: 5148 vs 2511 placed; identical where both exist,
  median diff = 0). More permissive → keep more cells. No fallback needed.

### Source 2 — provenance only (richer, but cb_index-keyed)
```
.../andy_raw_data/{coords_id}/Positions/coords.csv
```
- The raw slide-tags positioning output. Keyed by **`cb_index`** (integer), NOT barcode —
  would need an index→barcode mapping to join, so not used for the merge.
- Has extra per-cell diagnostics not in Source 1: `SNR`, `pct.placed`, `clusters`,
  `minPts`, `eps`, KDE coords (`x_um_kde`/`y_um_kde`). Only reach for this if you ever
  want to filter coordinates by placement confidence (e.g. SNR threshold).

**Cell handling:** ~31–63% of CellBender cells are spatially placed (varies by sample).
**Keep ALL cells** through ingest/QC/clustering; coords stay as nullable `obs` columns
(`x_um_dbscan`/`y_um_dbscan`). Only the explicitly spatial rules (07+) subset to
`obs['x_um_dbscan'].notna()`. Do NOT drop unplaced cells globally — the non-spatial
backbone (scVI, clustering, annotation, pseudotime) uses every cell, matching the paper's
57,242 mouse / 69,525 human post-QC totals.

---

## 4. Reference genomes (cellranger count 7.2.0)

Counts were generated with **Cell Ranger count 7.2.0**; BCLs demuxed with mkfastq 8.0.1.
Do NOT re-run cellranger or CellBender — the existing outputs are what the paper's results
are built on, and CellBender is stochastic (re-running gives different counts). Treat
these as immutable upstream artifacts; document, don't regenerate. Stay on 7.2.0 (do not
upgrade) so counts and the custom reference stay consistent.

### Mouse — custom transgenic reference
```
/n/eddy_lab/users/nhilgert/reference/transgenic_mus   (mkref 7.2.0)
```
Standard GRCm assembly + 2 appended DREADD transgene contigs (the irreplaceable part):
- `CHRM3_Jax026220` — GqD / Chrm3* DREADD construct (+ EGFP, WPRE)
- `GsD_Berdeaux`     — GsD (Gs-DREADD) construct
Transgene FASTA/GTF should be extracted + committed to
`preprocessing/reference/transgenes/` so the custom build is reproducible-in-principle
without rebuilding (rebuild is a once-ever, documented preprocessing step, NOT a DAG rule).

### Human — stock 10x reference
```
/n/eddy_lab/users/nhilgert/reference/human/refdata-gex-GRCh38-2024-A
```
Stock, downloadable; no transgenes. Reproducible by `wget` for an outside user.

---

## 5. Ambient RNA correction (CellBender)

- **CellBender remove-background 0.2.2**, `--total-droplets-included 50000`.
- Run by collaborator (Andy) — outputs are the `cellbender_output_filtered.h5` files above.
- Provenance only; the DAG starts FROM these corrected counts (does not run CellBender).

---

## 6. QC thresholds (from the paper Methods — pin these in config.yaml)

Operates on CellBender-corrected counts, via Scanpy (v1.10.4 in the paper).
- Minimum **500** corrected UMIs/cell.
- Adaptive **3-MAD** filters on `log1p(total_counts)` and `log1p(n_genes_by_counts)`.
- Mitochondrial fraction **> 8%** excluded. (Mito prefix: mouse `mt-`, human `MT-`.)
- Doublets: **DoubletDetection 4.2** (BoostClassifier, 40 iters, voter_thresh 0.5,
  per-iter cluster p < 1e-7); remove flagged doublets.
- **NB:** the scaffold `config.yaml` shipped with `max_pct_mito: 20` — must be **8** to
  match the paper.

**Regression targets (post-QC, from paper):** 57,242 mouse / 69,525 human cells. A full
pipeline run should land near these.

---

## 7. Downstream structure (for orientation)

- Per-species concatenated objects (mouse, human) — NOT one combined object.
- Integration: scVI (scvi-tools 1.1.2) on 3000 seurat_v3 HVGs (excl mito+ribo),
  `sample_id` as batch key, 2 layers, 30 latent, NB likelihood.
- Clustering: Leiden res 0.1 (coarse), knn 100.
- Pseudotime: diffusion (knn 15, 10 diffmap comps); root = uninduced outer eFib.
- Pseudotime-associated genes: mutual information (n_neighbors 3, 1000 perms, q<0.05,
  top 90th percentile).
- SAMap cross-species pairings: mouse 1_dpi_late × human secretory;
  mouse 4_dpi_bleeding × human menstrual.

---

## 8. Pipeline integration decisions (summary)

- **Ingest (`00b`):** read `xmens_inputs/{sample_id}.h5` → `var_names_make_unique()` →
  write `species`/`timepoint` from sheet → join Source-1 coords on barcode (direct,
  no suffix strip) → `dbscan_1/2` into `obs['x_um_dbscan'/'y_um_dbscan']` (nullable) →
  keep all cells → write `data/interim/raw/{sample_id}.h5ad`.
- **`00a` (coords):** retired. Coords come in via `00b`, not a separate fork run. The
  `external/slidetags` submodule stays as documented provenance of how coords were made.
- **`01_qc`:** apply §6 thresholds + DoubletDetection on all cells.
- **Spatial rules (07+):** subset to placed cells internally.
- **`samples.tsv` columns:** `sample_id, species, timepoint, modality, coords_id, path`
  where `path = /n/eddy_lab/users/nhilgert/decidualization/xmens_inputs/{sample_id}.h5`
  and `coords_id` is the mck lane ID.
