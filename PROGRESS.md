# PROGRESS

Working notes for `xmens_final` — the reproducible Snakemake pipeline for the
cross-species (mouse X-Mens + human) Slide-tags analysis in Cevrim et al. 2025
(bioRxiv 2025.10.08.681007). This file is the "where things stand / how to
resume" note; it is not polished public docs.

Last updated: 2026-06-30

---

## STATUS

### Done and validated
The reproduction spine runs end-to-end via the SLURM profile, both species:

```
ingest -> qc -> concat -> integrate
```

- **ingest** (`00b_ingest`): reads CellBender-corrected .h5 from
  `xmens_inputs/{sample}.h5`, joins Slide-tags DBSCAN coords, writes
  `data/interim/raw/{sample}.h5ad` (all cells kept).
- **qc** (`01_qc`): sequential filtering matching the QC notebook exactly —
  500-count floor FIRST, then 3-MAD outliers computed on the floored object
  (log10 scale, raw MAD), then mito >8%, then DoubletDetection 4.2. Writes
  `data/interim/qc/{sample}.h5ad`.
- **concat** (`01c_concat`): per-species concatenation, all cells kept, obs
  trimmed to a whitelist, adds `cell_barcode` (‑1 stripped, for joining to
  collaborator objects). Writes `data/interim/qc_passed_{species}.h5ad`.
- **integrate** (`02_integrate`): scVI (2 layers, 30 latent, NB) on 3000
  seurat_v3 HVGs (mito/ribo excluded), batch_key=sample; neighbors k=100;
  Leiden at 14 resolutions (igraph, seed 0); UMAP. GPU rule. Writes
  `data/interim/latent_{species}.h5ad` with `X_scVI`, `X_umap_scVI`,
  `leiden_scVI_{0.1..3.0}`.

Final integrated objects: **mouse 68,069 cells / human 70,363 cells**
(both carry X_scVI + X_umap_scVI + spatial + all 14 leiden resolutions).

### Documented but NOT in the default DAG
- **cellbender** (`00_cellbender`): reproduces ambient-RNA removal from the
  CellRanger raw matrices under `velocity/{coords_id}/cellranger/outs/`. Params
  match Methods (`--total-droplets-included 50000`, default epochs, `--cuda`).
  Deliberately decoupled: the pipeline starts from the already-corrected inputs,
  so re-running (stochastic) won't shift the validated counts. Dry-run verified;
  not executed.

### Next up (fresh-energy tasks)
1. **Comparison notebook** (`notebooks/compare_to_collaborator.ipynb`, written,
   not run): set-difference my `latent_{species}` objects against the
   collaborator's FINAL objects on (sample, cell_barcode) to isolate exactly
   which cells he removed downstream, and fingerprint them (unplaced /
   manual-exclusion-box / kNN-outlier / lineage-cluster residual). Need to fill
   in the **human** final-object path (mouse final is
   `/n/eddy_lab/Lab/mckinley/data/2025_04_15_AmhrGsD_Final.h5ad`, 57,242 cells).
2. **`03_cluster` / annotation** — downstream of integrate.
3. **Spatial QC** — deferred until AFTER lineage separation (so the manual
   exclusion boxes can be verified in the integrated/lineage context). The
   `filter_spatially` logic (per-sample kNN outlier removal k=5, keep 95%, plus
   hardcoded 2_dpi/1_dpi_late coordinate cuts) will be folded in then.
4. Eventually: SAMap cross-species (mouse 1_dpi_late × human secretory;
   mouse 4_dpi_bleeding × human menstrual), SCENIC, LIANA, pseudotime.

### Known stubs / issues
- `07_spatial_pseudotime.smk` has a stale input path (`data/interim/coords/…`,
  nothing produces it) — fix when reaching the spatial stage (point at the QC'd
  h5ad, subset to placed cells).
- `01b_qc_report` is written but was only recently added to the include list —
  confirm it runs if you want the per-sample QC figures.
- `cuda11.8_env` is a kitchen-sink env (overlapping CUDA packages, harmless
  cupy/tensorflow warnings) — "clean up later".

---

## RUNNING (operational gotchas — these are not in any code file)

These are the hard-won facts that make the pipeline runnable. Set them each
session (the PATH export is NOT persistent).

1. **PATH** — put rsc_uv's python first (for interactive scanpy work) AND the
   conda binary on PATH (so the SLURM orchestrator's subshells can find it):
   ```bash
   export PATH="/n/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin:/n/sw/Miniforge3-25.3.1-0/bin:$PATH"
   ```
   (`conda` is a shell *function*, not a binary, until you do this — snakemake's
   `--use-conda` subshells need the real binary on PATH.)

2. **Run the orchestrator from a LOGIN node**, not from inside an `salloc`.
   Submitting SLURM jobs from within a SLURM allocation fails *silently* (the
   job vanishes, snakemake spins on "No active jobs"). Login node = clean submit.

3. **GPU rules use the `gres:` resource key in the profile, NOT `slurm_extra`.**
   The snakemake-executor-plugin-slurm reads `job.resources.get("gres")` and
   builds `--gres=…` itself via `set_gres_string`. Putting `--gres` in
   `slurm_extra` gets malformed and the job is silently rejected. The GPU rules
   also set `slurm_extra: ""` so they don't inherit the CPU-node exclude list.
   ```yaml
   integrate:
     gres: "gpu:1"
     slurm_extra: ""
   ```

4. **Env split (why there are two envs):**
   - driver/spine = `rsc_uv` (uv, py3.11): snakemake + scanpy + scVI (+cu129
     torch, GPU-capable). Runs the orchestrator and most rules.
   - DoubletDetection = `cuda11.8_env` (conda, py3.10): DoubletDetection 4.2.
     The `louvain` dep won't build under uv/pip, so it stays conda; the qc rule
     reaches it via a `conda:` directive + `--use-conda`.

### Typical commands
```bash
# from a login node, with PATH set as above:
# one target:
snakemake --profile profiles/slurm data/interim/latent_mouse.h5ad
# whole QC stage (all samples), expanded from samples.tsv:
snakemake --profile profiles/slurm \
  $(for s in $(tail -n +2 config/samples.tsv | cut -f1); do echo data/interim/qc/$s.h5ad; done)
# dry-run to check a target resolves:
snakemake -n <target>
# detach a long run:  tmux new -s xmens  (then run inside; Ctrl-b d)
```

---

## VALIDATION

QC reproduces the paper's initial-QC cell counts to within **0.1%**. Threshold
filtering is deterministic and matches the QC notebook exactly (e.g. mck_1 lands
on 3508 post-threshold, on the dot); the residual vs the paper is
DoubletDetection stochasticity.

Post-QC per-species totals (this run):

| species | mine   | paper initial-QC | diff        |
|---------|--------|------------------|-------------|
| mouse   | 68,069 | 68,146           | −77  (−0.1%)|
| human   | 70,363 | 70,250           | +113 (+0.2%)|

> NOTE: this is the paper's **initial-QC** checkpoint (68k/70k), NOT the final
> 57,242 mouse / 69,525 human. Those come after downstream lineage-reintegration
> removal (518 mouse / 692 human) + annotation cluster exclusions — steps not yet
> reproduced. The comparison notebook (see Next up) will decompose the
> ~11k-cell mouse gap between my 68k and the collaborator's final 57k.

Per-sample post-QC counts (mouse | human), from the concat logs:

```
uninduced_mck8a      11892     secretory_mck6a      13024
uninduced_mck8b       8446     secretory_mck6b      12685
1_dpi_early_mck1      3458     secretory_mck6c       7517
1_dpi_late_mck12a     7497     menstrual_mck5a      12794
1_dpi_late_mck12b     6860     menstrual_mck5b      12795
2_dpi_mck15           5014     menstrual_mck5c      11548
3_dpi_mck10a          8273
3_dpi_mck10b          8196
4_dpi_bleeding_mck11a 8433
```

(DoubletDetection is stochastic, so these wobble ~0.1–0.2% between runs. The
on-disk objects are the source of truth; freeze downstream work on them rather
than re-running qc.)
