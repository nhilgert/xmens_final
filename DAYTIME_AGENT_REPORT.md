# Daytime agent report — 2026-08-11

Repo: `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final`  
State inspected: 2026-08-11 EDT  
Delivery state: **uncommitted review diff**; no promotion was performed.

## Executive result

- Built a seeded, parameterized, scratch-only mouse pseudotime rule and tested 1dpi late and 4dpi end to end. Accepted correlations are **+0.332406** and **+0.449217**, respectively (`data/interim/_pt_test_20260811/{1dpi_late,4dpi}_metrics.json`).
- Produced a diagnostic-only 3dpi root decision aid. **3dpi root awaiting human choice — no root committed, no canonical 3dpi object written.** (`data/interim/_pt_test_20260811/3dpi_diagnostic_summary.json`; `workflow/scripts/compute_pseudotime.py:263-308`).
- `XMENS_PIPELINE_HANDOFF.md` is absent, so Task 2 was not applicable and no substitute document was invented.
- Classified all 106 notebooks in `NOTEBOOK_TRIAGE.md`: 5 canonical, 90 stale, 11 pipeline-depended (`NOTEBOOK_TRIAGE.md:5-27`).
- Authored but did not execute `tools/promote_simple.sh`. The elaborate 21-file promotion plan **does not include McKinley**; the simple script explicitly adds the four independent McKinley mirrors (`tools/promote_simple.sh:35-44`).
- All accepted scientific outputs are confined to `data/interim/_pt_test_20260811/`. No SAMap fit or promotion ran.

## Task 1 — parameterized pseudotime rule and scratch test

### Implementation

| File | Purpose | Evidence |
|---|---|---|
| `workflow/scripts/compute_pseudotime.py` | Canonical-label fib/dec subset; counts-based scoring/HVG/scVI; neighbors/diffmap; tuned DPT or diagnostic-only candidate evaluation; atomic scratch writes | canonical checks `:168-202`; seeded scVI `:239-253`; geometry `:255-261`; 3dpi diagnostics `:263-308`; tuned root/orientation gates `:310-359`; transfer contract `:361-378` |
| `workflow/rules/06b_timepoint_pseudotime.smk` | One tuned job per timepoint plus a separate 3dpi diagnostic rule; all outputs under the scratch root | `:1-48` |
| `config/config.yaml` | Preserves the legacy `pseudotime` stub parameters and adds isolated `pseudotime_tests` settings | `:96-144` |
| `profiles/slurm/config.yaml` | One GPU, 8 CPUs, 48 GB, 240 minutes for each scratch rule | `:53-65` |
| `workflow/Snakefile` | Includes the new rule file, but neither new rule is in `rule all` | include `:19-23`; `rule all` `:25-37` |

The new settings use a separate `pseudotime_tests` namespace rather than replacing the existing `pseudotime` mapping. This prevents the legacy global-pseudotime stub from becoming stale solely because this scratch test was added (`config/config.yaml:96-108`). The rule passes common settings plus only the current stage's mapping, so a future 4dpi change cannot invalidate 1dpi or 3dpi provenance (`workflow/rules/06b_timepoint_pseudotime.smk:6-24,34-44`).

Final code hashes:

- `workflow/scripts/compute_pseudotime.py`: `90cfa0115dcb269e2675bbb25c4ec6d14a431be18bfa2143c04a5a06c3878f5c`
- `workflow/rules/06b_timepoint_pseudotime.smk`: `581890d50b108362e52cf785e38587e0dd5c8c919b5711f37c7c9078897fcade`
- `config/config.yaml`: `9656aba3cff86ae622300e39f31d35ae3b6bc212fe33afeeebc3230b4a73e7c9`

### Tuned orientation results

| Stage | Cells | Root | `corr(DPT, decidual_score)` | Mean DPT: fibroblast → decidual | Result |
|---|---:|---|---:|---:|---|
| 1dpi_late | 10,656 (10,787 before excluding 131 `decidual_apoptotic`) | DC1 argmin; `efib_myometrial`; `ACAAGCTAGATCCGAG-1-1_dpi_late_mck12a` | **+0.3324063115** | 0.446607 → 0.492429 | PASS |
| 4dpi | 5,442 (6,141 before excluding 699 `decidual_apoptotic`) | DC1 argmin; `efib_myometrial`; `AAGATAGCACTAAACC-1-4_dpi_bleeding_mck11a` | **+0.4492171770** | 0.493544 → 0.726515 | PASS |

Sources: `data/interim/_pt_test_20260811/1dpi_late_metrics.json`, `4dpi_metrics.json`, and their subtype TSVs. Both DPT vectors are finite and span exactly `[0,1]`; each root has DPT 0. Independent recomputation from the H5ADs reproduced every JSON orientation number.

Per-subtype mean ordering (ascending DPT):

- **1dpi_late:** `efib_myometrial` .2332 → `efib_middle_induced` .4534 → `efib_middle` .4638 → `efib_proliferative` .4659 → `decidual_border` .4827 → `decidual_proliferative` .4852 → `efib_subluminal` .5104 → `decidual_mature` .5764 (`1dpi_late_subtype_pseudotime.tsv`). The aggregate fibroblast mean is below decidual; `efib_subluminal` overlaps the decidual range, so this is a continuum rather than a perfectly separated ladder.
- **4dpi:** `efib_myometrial` .2534 → `efib_middle_induced` .5159 → `efib_proliferative` .5159 → `efib_subluminal` .5947 → `decidual_border` .6656 → `efib_middle` .7079 → `decidual_proliferative` .7119 → `decidual_mature` .7797 (`4dpi_subtype_pseudotime.tsv`). Aggregate orientation passes; the small `efib_middle` group (n=24) lies within the decidual range.

#### 4dpi endpoint correction exposed by the test

The first seeded 4dpi attempt intentionally used the notebook's global DC1 argmax. It landed in `decidual_mature` (`TCACTATGTGTCCAAT-1-4_dpi_bleeding_mck11a`) and failed the mandatory orientation gate at **−0.595377**; job `38270495` exited 1 and published no declared output. Evidence remains in `.snakemake/slurm_logs/rule_compute_pseudotime/4dpi/38270495.log` and `data/interim/_pt_test_20260811/snakemake_controller.log`.

The reproducible scratch config therefore uses the opposite DC1 end (argmin), which landed in `efib_myometrial` and passed at +0.449217 (`config/config.yaml:124-135`). This is not a hidden sign flip: both the failed endpoint and the tested correction are recorded. Diffusion-component signs/endpoints are fit-dependent; the biological orientation gate remains required on future refits.

### 3dpi root decision aid — human choice required

**3dpi root awaiting human choice — candidates evaluated below. No root was selected, no sign was forced, and no 3dpi H5AD or temporary H5AD exists.** The diagnostic used 10,649 cells after excluding 24 `decidual_apoptotic` cells and evaluated the three explicit candidates in `workflow/scripts/compute_pseudotime.py:263-308`.

| Candidate | Root subtype | Correlation | Mean DPT: fibroblast → decidual | Direction |
|---|---|---:|---:|---|
| DC1 argmin | `decidual_mature` | **−0.425524** | 0.604367 → 0.486646 | fails |
| DC1 argmax | `efib_myometrial` | **+0.348331** | 0.621510 → 0.771710 | passes |
| Most-extreme cell in `efib_middle_induced` | `efib_middle_induced` | **+0.502649** | 0.388137 → 0.584026 | passes; strongest correlation |

The third candidate's subtype was selected mechanically as the fibroblast subtype with the lowest mean `decidual_score`; its cell was the largest absolute DC1 deviation from the global DC1 median within that subtype (`compute_pseudotime.py:265-285`). That rule is a decision aid, not a biological endorsement.

Subtype sequences (ascending mean DPT; values rounded to 4 decimals):

- **DC1 argmin:** `decidual_mature` .4264 → `efib_middle` .4955 → `decidual_proliferative` .4986 → `decidual_border` .5211 → `efib_proliferative` .5407 → `efib_middle_induced` .5711 → `efib_subluminal` .5997 → `efib_myometrial` .7717. This is reversed overall.
- **DC1 argmax:** `efib_myometrial` .2857 → `efib_middle_induced` .6869 → `efib_proliferative` .7039 → `efib_middle` .7329 → `decidual_border` .7453 → `decidual_proliferative` .7661 → `decidual_mature` .8085 → `efib_subluminal` .8247.
- **Extreme `efib_middle_induced`:** `efib_myometrial` .2272 → `efib_middle_induced` .4104 → `efib_proliferative` .4630 → `efib_middle` .5075 → `decidual_border` .5330 → `decidual_proliferative` .5729 → `decidual_mature` .6558 → `efib_subluminal` .6975.

Source tables: `data/interim/_pt_test_20260811/3dpi_root_candidates.tsv` and `3dpi_candidate_subtype_pseudotime.tsv`; status contract: `3dpi_diagnostic_summary.json` (`AWAITING_HUMAN_ROOT_CHOICE`, `final_h5ad_written=false`). Both positive candidates place the small `efib_subluminal` group last, so the human should consider the full ordering, not correlation alone.

### Artifact validation and execution provenance

- Accepted 1dpi H5AD: 10,656 × 33,698, exact staged subset/order/canonical annotations; finite `X_scVI` (30D), `X_diffmap` (15D), DPT, and counts-backed expression. Its bare 16-bp barcode index intentionally contains 10 duplicate keys/20 rows, preserving the known multi-sample transfer caveat (`logs/compute_1dpi_late.log`).
- Accepted 4dpi H5AD: 5,442 × 33,698 with the same gates; bare barcodes are unique.
- TSV group counts sum exactly to each H5AD; independently recalculated means/medians/counts agree to floating precision.
- Scratch root: `data/interim/_pt_test_20260811/`; no hidden output temporary files remain.

The final provenance refresh regenerated all nine artifacts from scratch under the finalized rule. Every output was **bit-for-bit identical** to the first accepted run, providing a direct seeded reproducibility test. Final hashes are:

| Artifact | SHA-256 |
|---|---|
| `adata_1dpi_late.h5ad` | `f1b4b8b3003f9ec176f0e71fe7487038c5e0ed9026641804d4a7f5b3fc4331ff` |
| `1dpi_late_metrics.json` | `5a6f27c10d68b425a505f2779906f909a5c6489d43293389fcdf82d7bb01e1b2` |
| `1dpi_late_subtype_pseudotime.tsv` | `ef899849b573d81371241a7cd9fd5a29c9e7727776b7c3fea346398e59ca766f` |
| `adata_4dpi.h5ad` | `3ba4c45e9b02cfab524df8478e9b0bf5f7f8d782069257a6dc175c5d8b8ef464` |
| `4dpi_metrics.json` | `345844d8e623618995e2299def9c996df26470a18f89f3803248d71c209b422f` |
| `4dpi_subtype_pseudotime.tsv` | `ea438cce4e0d1ddec9003a41e519eeddd508aec57edde131780148d867f0a4b7` |
| `3dpi_root_candidates.tsv` | `1511cb002b3c0e0c69fbe43e49a4e51134ab5308da46e587d1605566d3a1be4d` |
| `3dpi_candidate_subtype_pseudotime.tsv` | `c5b3e278d83c1ca03508b5aa0a9b03f6892183fe3bcf1eb78248884971558f53` |
| `3dpi_diagnostic_summary.json` | `b8752329e0bdd87cdb8f53b45348e1444296d3faa26605db5f20eb30d6e6d0a7` |

Scheduler provenance: accepted jobs `38270490` (1dpi; exit 0), `38270508` (3dpi diagnostic; exit 0), and `38277843` (corrected 4dpi; exit 0); the gate-failed 4dpi job was `38270495` (exit 1). A final three-job refresh under the provenance-isolated rule used jobs `38285174`, `38285175`, and `38285266`.

The final restricted dry-run (only the three scratch targets and the two new allowed rules) exited 0 with: `Nothing to be done (all requested files are present and up to date).` Its complete command/output is preserved at `data/interim/_pt_test_20260811/dryrun_final_noop.log`. An earlier dry-run had correctly exposed cross-timepoint parameter coupling; `PT_GLOBAL_CFG` now excludes the `timepoints` table and each job records only its own stage mapping (`workflow/rules/06b_timepoint_pseudotime.smk:6-24,34-44`).

### Cache caveat

Two early canceled controller attempts used Snakemake's default runtime-source cache before this was noticed. They left small source-only cache trees in home:

- `~/.cache/snakemake/snakemake/source-cache/snakemake-runtime-cache/tmpc8g830mr`: 29 files, 59,501 bytes, created 09:42:33–34 EDT.
- `~/.cache/snakemake/snakemake/source-cache/snakemake-runtime-cache/tmp6480245y`: 29 files, 60,188 bytes, created 09:45:57–58 EDT.

They were not cleaned or modified. No scientific script reached output in those attempts (`data/interim/_pt_test_20260811/snakemake_controller_cancelled_attempt1.log` and `snakemake_controller_homecache_attempt3.log`). All successful runs explicitly pinned XDG, Python bytecode, Matplotlib, Numba, Torch, and Snakemake runtime caches under the scratch root. This is the one exception to the requested no-home-write condition and is reported rather than concealed.

## Task 2 — handoff reconciliation

`XMENS_PIPELINE_HANDOFF.md` does not exist at repo root, anywhere in the worktree, in tracked files, or in path history as of HEAD `87fd380a8bf1e9f9d7ff618a741b99d752059d7a`. The prompt made this task conditional (“if present”), so no document was created and no corrections were inferred.

## Task 3 — notebook triage

Created `NOTEBOOK_TRIAGE.md` (179 lines; SHA-256 `b5ad832b6a5b43c89be7303e8a687409c06c9829e4fdc1ada6d560e70ef7782b`). Headline inventory (`NOTEBOOK_TRIAGE.md:5-27`):

| Classification | Count |
|---|---:|
| CANONICAL | 5 |
| STALE | 90 |
| PIPELINE-DEPENDED | 11 |
| **Total** | **106** |

No notebook was executed, edited, moved, or deleted. The report explicitly identifies the 13 current `resources/annotations/*.csv` files, interpreted with `resources/annotation_vocabulary.yaml` and `resources/lineage_spec.yaml`, as the annotation authority (`NOTEBOOK_TRIAGE.md:25-29`).

## Task 4 — simple promotion script and McKinley destination audit

Created `tools/promote_simple.sh` (80 lines, mode 0755, SHA-256 `a3e047fe38fc96c100961c0274ac37308135ca77b7cf027de1e9fd9c504c0522`). It was syntax/static checked but **never executed**.

The elaborate plan `tools/promotion_plan_xmens_annotation_v1_20260810.json` has exactly 21 mappings and **does not include McKinley**; deposit mappings target only `dataverse/objects/`. The exact independent mirror directory is `/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data`. Its four H5ADs are regular, single-link files with distinct inodes from the Dataverse copies; their current SHA-256 values match the plan's reviewed old-deposit hashes.

The simple script adds those four destinations to the 21 reviewed targets (`tools/promote_simple.sh:35-44`), checks the plan/script are committed and clean (`:20-32`), verifies all staged/live hashes (`:45-52`), prepares and verifies all 25 backups and adjacent install files before any move (`:53-70`), then atomically renames and rechecks in dependency order (`:71-80`). Backups go to the requested `data/interim/_promotion_backup_<UTC timestamp>/` path (`:14-15`). Because the script is deliberately left untracked/uncommitted, its own safety gate prevents execution pending human review.

## Final safety checks and human decisions

- **No promotion; no SAMap.** `tools/promote_simple.sh` and `tools/migration_promotion.py` were not run, and no SAMap rule/job ran.
- The protected-live verification against `data/interim/migrated/xmens_annotation_v1_20260810/manifests/protected_live_before.json` passed for all 26 protected live/deposit/SAMap objects with `changes=[]`.
- The staged mouse input remained exactly 1,976,268,956 bytes with SHA-256 `7f2dbd16d966645134f633e7cf1a0169135a569f5f406635ceaf4546879cba7e`.
- Human decision required: choose the 3dpi root using the diagnostic table; no default is encoded (`config/config.yaml:136-144`).
- Human decision required: review/commit `tools/promote_simple.sh` before it can run. The elaborate plan's McKinley omission is real.
- Human cleanup decision: the 90 stale notebooks are recommendations only; nothing was moved (`NOTEBOOK_TRIAGE.md`).
