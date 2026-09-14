# SAMap → Snakemake integration runbook  (v2 — stock SAMap v3, post-env-setup)

Wiring the SAMap arm (homology build + fit + label-transfer) into `xmens_final`.
Run everything from the `xmens_final` repo root. Env is **stock sc-samap 3.0.1** (no fork).

> **What changed from v1:** the env chat confirmed the pipeline uses *stock* SAMap v3, not the
> instrumented fork, and that the 36 old `samap.pkl` files are **unreadable** in the new env
> (anndata 0.12 vs 0.10.8). So the fit must actually run and the two paper pkls must be
> regenerated. `run_samap.py` was ported to v3 (no `iteration_callback`, no `snapshots.pkl`) and
> now lives in this repo. All interpreter/tool paths are absolute.

---

## 0. How it attaches (now hooked up to the pipeline)

The SAMap arm is a downstream *arm*; the annotation pipeline (04*) is untouched. It now attaches
via a **bridge rule** rather than frozen files:

```
refined_{human,mouse}.h5ad  (ancient -> annotation DAG stays isolated)
     -> prep_samap_input   (subset fib/dec, rename to shared schema, split by timepoint)
     -> inputs/h5ads_samap/{hs,mm}_fib_dec_<tp>.h5ad
     -> run_samap          -> runs/<run>/samap.pkl
     -> samap_label_stain  -> results/figures|tables/samap_transfer/*
```

`refined_*` is `ancient()` so `snakemake samap_all` never pulls in or rebuilds a 04* rule (safe
while the annotation-side config paths are still un-migrated). The SAMap inputs are now *derived
from your current annotation* instead of the stale April snapshot.

---

## 1. Place the files (five)

```
workflow/rules/09_homology.smk
workflow/rules/10_samap.smk
workflow/scripts/transfer_label_staining.py     # refactored (args, not hardcoded)
workflow/scripts/run_samap.py                    # STOCK-v3 port
workflow/scripts/prep_samap_input.py             # NEW: the refined_* -> SAMap-input bridge
rm workflow/rules/10_samap_transfer.smk          # superseded, if present
```
The fork `run_samap.py` in `samap_integrate/scripts/` is left untouched; the DAG uses this repo's
ported copy.

## 2. `config/config.yaml` — one block, real paths (no placeholders)

```yaml
samap:
  root:            /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate
  coords_root:     /n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/crossspecies_radial_gp
  proteomes_root:  /n/lab_storage/eddy_lab/users/nhilgert/data/proteomes
  python:          /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/envs/samap/.venv/bin/python
  saturn_python:   /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/saturn_uv/bin/python
  phmmer:          /n/lab_storage/eddy_lab/users/nhilgert/software/hmmer-3.4/src/phmmer
  n_iter:          3
  pairs: ["menstrual_x_4dpi_hmmer", "secretory_x_1dpi_late_hmmer"]
```

## 3. Snakefile — two includes

```python
include: "rules/09_homology.smk"
include: "rules/10_samap.smk"
```

## 4. Environments — three, all called by ABSOLUTE PATH

SLURM jobs may not source `.bashrc`, so no `uva`/`conda activate`/`$PATH` reliance.

| context | absolute path | used by |
|---|---|---|
| annotation | `rsc_uv` (existing) | 04* — unchanged |
| SAMap fit + stain | `.../xmens_final/envs/samap/.venv/bin/python` — **stock sc-samap 3.0.1** | `run_samap`, `samap_label_stain` |
| homology aggregate | `.../.uv/envs/saturn_uv/bin/python` | `phmmer_aggregate_*` (only if homology is rebuilt) |
| phmmer binary | `.../software/hmmer-3.4/src/phmmer` | `phmmer_chunk` (only if homology is rebuilt) |

## 5. Path fixes / pre-flight (do before the first fit)

1. **Verify the v3 SAMAP constructor** (only `.run()` was confirmed upstream). `run_samap.py`'s
   `SAMAP(...)` call carries fork-era kwargs (`f_maps`, `save_processed`):
   ```bash
   SAMAP_PY=/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/envs/samap/.venv/bin/python
   $SAMAP_PY -c "import inspect, samap; print('init', inspect.signature(samap.SAMAP.__init__)); print('run', inspect.signature(samap.SAMAP.run))"
   ```
   If `__init__` lacks `f_maps`/`save_processed`, adjust the `SAMAP(...)` line in
   `workflow/scripts/run_samap.py` (it's flagged in the docstring).
2. **phmmer on a fresh compute node** (only matters for a homology rebuild):
   `.../software/hmmer-3.4/src/phmmer -h | head -2` → `HMMER 3.4`.

---

## 6. Verify — in this order. The dry-run `reason:` field is the arbiter.

### 6.1 It parses
```bash
snakemake --profile profiles/slurm -n samap_all
```
A clean job list = config + includes + `pairs.csv` resolved. `KeyError` here = a missing
`config[samap]` key.

### 6.2 Delete the stale artifacts so the DAG rebuilds them from current data
```bash
S=/n/lab_storage/eddy_lab/users/nhilgert/decidualization/revisions/samap_integrate
# stale April inputs (built from pre-July annotation) — prep regenerates from refined_*:
rm -f $S/inputs/h5ads_samap/{hs_fib_dec_menstrual,hs_fib_dec_secretory,mm_fib_dec_4_dpi_bleeding,mm_fib_dec_1_dpi_late}.h5ad
# stale, v3-unreadable pkls — run_samap regenerates:
rm -f $S/runs/menstrual_x_4dpi_hmmer/samap.pkl $S/runs/secretory_x_1dpi_late_hmmer/samap.pkl
```

### 6.3 Dry-run again — the full fresh chain should appear
```bash
snakemake --profile profiles/slurm -n samap_all
```
Expect: **4 × `prep_samap_input` + 2 × `run_samap` + 2 × `samap_label_stain`** (+ `samap_all`).
What should still be absent: any `04*` rule (refined_* is `ancient`), any `phmmer_*`/`homology_build`.
Also confirm `snakemake -n homology_build` says **nothing to do**.

### 6.4 Run it (fit is the expensive part)
```bash
snakemake --profile profiles/slurm samap_all
```
Tune `run_samap`'s `mem_mb`/`runtime` to the `eddy` partition first. The first fit pays a one-time
numba JIT-compile cost — don't set the timeout too tight.

### 6.5 Verify the fresh pkl is readable + the labels are current
- **Stain succeeded** = the fresh v3 pkl unpickled. If the stain step errors on `pickle.load`,
  `save_samap` uses a non-plain format — switch the loader in `transfer_label_staining.py` to
  samap's `load_samap` and re-run the stain (cheap).
- `results/tables/samap_transfer/*.json` → `pkl_label_categories.hs` ==
  `["decidual cells","endometrial fibroblasts"]` (confirms current fib/dec annotation).
- `transfer_agreement.human` present with non-degenerate `dec_recall`/`fib_recall`; open the PNG.

**Schema note the bridge encodes** (why prep_samap_input is species-parametrized): the two
refined objects label the same content under different obs columns —
human `cell_type`/`subtype`, mouse `coarse_refined`/`fine_refined` — both renamed to the shared
`cell_type_annotations`/`cell_subtype_annotations`. Fib/dec counts the prep asserts against:
human 47,533, mouse 43,197. And **`{hs,mm}_cell_barcode` in the fitted joint obj**: the stain's
`join_coords` needs those; they come from SAMap prefixing the input's `cell_barcode`. If v3 doesn't
prefix the same way, check the fitted pkl's obs columns on the first fit and adjust the stain's
`HS_BC_COL`/`MM_BC_COL` constants.

---

## 7. Not done yet / loose ends

- **`transfer_pseudotime.py`** — the SAMap-side pseudotime transfer, the sibling of the stain leaf
  and part of the paper. Send its argparse + outputs and it wires the same way (fresh pkl input).
  It will need the same fork→v3 audit `run_samap.py` needed — grep it for `iteration_callback`,
  `gnnm_*`, `X_wPCA` before wiring.
- **Homology aggregate output paths** are `TODO:CONFIRM` in `09_homology.smk` — blocked on the
  `phmmer_to_samap_symbols.py` body (to match declared output to where it writes) and a fixpath of
  that script. Also verify `saturn_uv` actually runs it (a sibling env had a stale-path `pip`
  shebang post-migration). None of this blocks the paper — the maps exist and are frozen.
- **Git hygiene in `samap_integrate`.** The arm writes into that repo (runs/). Gitignore the
  generated dirs and commit the real scripts.
