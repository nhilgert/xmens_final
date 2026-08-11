# Annotation vocabulary migration — overnight staging report

**Run:** `xmens_annotation_v1_20260810`  
**Schema:** `xmens.annotation.v1`  
**Repository:** `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final`  
**Safety model:** staging-only; no promotion was performed.

## Top-line status

**PASS — the staging migration completed.** Snakemake finished 12/12 steps. `MIGRATION_REVIEW_READY` contains `PASS`, and `manifests/MIGRATION_STATUS.json` reports overall `PASS` with no validator failure or suspect. No promotion was performed.

The migration was implemented as an isolated Snakemake workflow rooted at:

- `data/interim/migrated/xmens_annotation_v1_20260810/`
- `dataverse/objects_migrated/xmens_annotation_v1_20260810/`
- `logs/migration/xmens_annotation_v1_20260810/`

Existing refined/final/reannotated/SAMap/Dataverse objects were treated as protected, read-only inputs. A full SHA-256 baseline covering 26 protected objects was written before the run. No object was promoted into a live path.

Final staged inventory is exact and contains no temporary files:

- 11 reannotated lineage H5ADs
- 2 refined all-cell H5ADs
- 2 final/re-embedded all-cell H5ADs
- 2 final fib/dec H5ADs
- 4 versioned Dataverse H5ADs
- audit manifests, count-gate TSVs, and the review-ready marker

Code checkpoints:

- `8c84132` — pre-migration snapshot
- `b3803e4` — staging-only annotation vocabulary migration
- `2ecd1dd` — Snakemake script-wrapper compatibility fix
- the commit containing this report — final documented handoff (created after workflow completion)

## Completed data stages

All 11 staged 04c lineage objects passed independent checks for retained cell IDs/order, unique observations and genes, a retained raw-count layer, canonical annotation columns, vocabulary validity, and absence of legacy annotation columns.

Both staged 04d merged objects passed independent checks:

- Human: 66,830 cells × 36,601 genes
- Mouse: 57,033 cells × 33,698 genes

For both, the cell set matches the staged lineage inputs exactly; gene order is identical across inputs; annotation counts match the locked vocabulary; and the merge count-gate manifest contains no count discrepancy.

Both 04e GPU jobs completed on the `eddy` partition using one A40 each:

- Mouse job `38206843`: 57,033 × 33,698, finite `X_scVI` (57,033×30) and `X_umap_scVI` (57,033×2)
- Human job `38206987`: 66,830 × 36,601, finite `X_scVI` (66,830×30) and `X_umap_scVI` (66,830×2)

Full streamed raw-count digests match refined→final for each species, and sampled rows independently reconstruct `.X` as `log1p(normalize_total(counts, 1e4))` within 3.31×10⁻⁷. The staged fib/dec subsets are exact, order-preserving masks of the final atlases: human 47,533 cells and mouse 43,197 cells. Their latent coordinates, UMAP coordinates, annotations, and count rows match the corresponding final-atlas rows exactly.

Completed semantic manifests:

- `manifests/export_human.json`: PASS
- `manifests/export_mouse.json`: PASS
- `manifests/check_annotations_human.json`: PASS for all 10 human staged/deposit objects
- `manifests/check_annotations_mouse.json`: PASS for all 11 mouse staged/deposit objects
- `manifests/validation_human.json`: PASS; failures=[], suspects=[]
- `manifests/validation_mouse.json`: PASS; failures=[], suspects=[]
- `manifests/protected_live_after.json`: PASS; 26 protected files; changes=[]
- `manifests/MIGRATION_STATUS.json`: PASS

All four staged deposits are byte-identical to their declared staged sources:

| Artifact pair | SHA-256 |
|---|---|
| Human final / all-cell deposit | `3600942422285ae3aee63adb5a35880b475775b4a0b0516ff6d80e4d769b0fa9` |
| Human fib/dec / fib-dec deposit | `44ee8c26466f1f4d86494c7c74674dfe50198c9fb209636b81ac006fc691cb42` |
| Mouse final / all-cell deposit | `7f2dbd16d966645134f633e7cf1a0169135a569f5f406635ceaf4546879cba7e` |
| Mouse fib/dec / fib-dec deposit | `13cd2730fd792186532545e88e0c71aefe0d1ea59d0ef14a67e43637629c33ed` |

## Count gates

### Human coarse types

| Gate | Expected | Actual | Status |
|---|---:|---:|---|
| decidual_cell | 21,648 | 21,648 | PASS |
| endometrial_fibroblast | 25,885 | 25,885 | PASS |
| epithelial_cell | 10,699 | 10,699 | PASS |
| immune_cell | 4,709 | 4,709 | PASS |
| smooth_muscle_cell | 2,084 | 2,084 | PASS |
| endothelial_cell | 1,334 | 1,334 | PASS |
| perivascular_cell | 471 | 471 | PASS |
| total | 66,830 | 66,830 | PASS |

### Mouse coarse types

| Gate | Expected | Actual | Status |
|---|---:|---:|---|
| decidual_cell | 22,919 | 22,919 | PASS |
| endometrial_fibroblast | 20,278 | 20,278 | PASS |
| epithelial_cell | 6,265 | 6,265 | PASS |
| immune_cell | 2,091 | 2,091 | PASS |
| smooth_muscle_cell | 1,443 | 1,443 | PASS |
| endothelial_cell | 3,487 | 3,487 | PASS |
| mesothelial_cell | 550 | 550 | PASS |
| total | 57,033 | 57,033 | PASS |

### Content tasks

| Gate | Expected | Actual | Status |
|---|---:|---:|---|
| Human endothelial total | 1,334 | 1,334 | PASS |
| Human endothelial_blood | 828 | 828 | PASS |
| Human endothelial_lymph | 506 | 506 | PASS |
| Human mural perivascular_cell | 471 | 471 | PASS |
| Mouse endothelial_blood | 3,487 | 3,487 | PASS |
| Mouse endothelial_lymph | 0 | 0 | PASS |
| Human non-ciliated epithelial blank subtype | 10,608 | 10,608 | PASS |
| Human smooth muscle blank subtype | 2,084 | 2,084 | PASS |
| Human perivascular blank subtype | 471 | 471 | PASS |
| Mouse epithelial blank subtype | 6,265 | 6,265 | PASS |
| Mouse smooth muscle blank subtype | 1,443 | 1,443 | PASS |
| Mouse mesothelial blank subtype | 550 | 550 | PASS |

No count was adjusted to make a gate pass.

## Issues and SUSPECT items

One plumbing failure occurred on the first 04d attempt. Snakemake prepends wrapper code to files invoked with `script:`, causing `from __future__ import annotations` to raise a syntax error. Human job `38206244` and mouse job `38206304` failed before writing a merged object. The unnecessary future imports were removed from all migration scripts invoked through `script:`, committed as `2ecd1dd`, and the workflow resumed from the 11 intact staged 04c objects. The corrected merge jobs completed and passed independent validation. This was a wrapper compatibility issue, not a count/content discrepancy.

**No staged object is marked SUSPECT.** Both validators report `failures=[]` and `suspects=[]`; all 28 deep bundle gates have delta zero. The protected-live full rehash job `38213182` completed successfully in 24m23s and found `changes=[]` across all 26 protected files.

The protected baseline was recorded at commit `b3803e4`; workflow execution finished at `2ecd1dd`. That one expected mid-run commit contains only the Snakemake wrapper compatibility fix described above. The repository advanced once more afterward solely to commit this report; no workflow or data job ran from that report-only commit.

## Files changed

### Vocabulary and decision records

- `resources/annotation_vocabulary.yaml` (new locked source of truth)
- `resources/lineage_spec.yaml`
- all 13 committed decision CSVs under `resources/annotations/` (human: epithelial, immune, lineage, perivascular, SMC, stromal; mouse: epithelial, immune, lineage, perimetrium, SMC, stromal, vascular)

The raw collaborator exports `human_immune_fine.csv` and `mouse_immune_fine.csv` were not edited.

### Producers, validators, and migration workflow

- `workflow/scripts/apply_annotation.py`
- `workflow/scripts/merge_reannotated.py`
- `workflow/scripts/integrate.py`
- `workflow/scripts/prep_samap_input.py`
- `workflow/scripts/check_migrated_annotations.py` (new)
- `workflow/scripts/export_migrated.py` (new)
- `workflow/scripts/validate_migration_bundle.py` (new)
- `workflow/scripts/finalize_migration.py` (new)
- `workflow/migration/Snakefile` (new, staging-only DAG)
- `workflow/rules/04_annotate.smk`
- `workflow/rules/04c_reannotate.smk`
- `workflow/rules/04d_merge.smk`
- `workflow/rules/10_samap.smk`
- `tools/check_annotations.py`
- `tools/compile_denovo.py`
- `tools/migrate_annotation_records.py` (new)
- `tools/protected_live_manifest.py` (new)
- `profiles/slurm/config.yaml`

### DRAFT notebook scaffolding

- `tools/gen_public_04c_drafts.py` (new)
- `notebooks/public/DRAFT_04c_MANIFEST.json`
- 11 `notebooks/public/DRAFT_04c_annotate_*_denovo.ipynb` scaffolds

These notebooks are generated from committed decisions and are explicitly DRAFT; they are not publication figures or new human annotation judgments.

### Handoff documentation

- `MIGRATION_OVERNIGHT_REPORT.md` (this report; committed after the staged workflow completed)

### Stale path cleanup

Known `/n/eddy_lab/` paths were replaced with `/n/lab_storage/eddy_lab/` in every touched file containing them:

- `config/config.yaml`
- `tools/compile_denovo.py`
- `tools/gen_public_04c_drafts.py`

`config/config.yaml` also had two top-level `references:` mappings; these were renamed to `genome_references:` and `annotation_references:` so the latter no longer silently replaces the former during YAML parsing. This was a configuration correctness fix, not an annotation-content change.

## Promotion checklist — human action only

Do not promote unless `MIGRATION_STATUS.json`, all object-validation reports, both annotation-check reports, and `protected_live_after.json` are PASS.

1. Stop writers/consumers of the live annotation objects and record a fresh backup or snapshot.
2. Review the two commits and this report; inspect all DRAFT notebooks separately.
3. Verify the following staged-to-live object mappings and hashes.
4. Promote in dependency order: reannotated lineage objects → refined atlases → final atlases → fib/dec subsets → Dataverse copies.
5. Rerun or explicitly force downstream consumers that are intentionally insulated by `ancient()`, especially SAMap inputs/fits/transfers.
6. Run the live validator and a Snakemake dry-run after promotion.

Exact staged source roots:

- `data/interim/migrated/xmens_annotation_v1_20260810/reannotated_migrated/`
- `data/interim/migrated/xmens_annotation_v1_20260810/refined_{human,mouse}_migrated.h5ad`
- `data/interim/migrated/xmens_annotation_v1_20260810/final_{human,mouse}_migrated.h5ad`
- `data/interim/migrated/xmens_annotation_v1_20260810/final_{human,mouse}_fib_dec_migrated.h5ad`
- `dataverse/objects_migrated/xmens_annotation_v1_20260810/`

Exact object mappings:

| Staged source (under the versioned run root) | Live destination |
|---|---|
| `reannotated_migrated/human_stromal.h5ad` | `data/interim/reannotated/human_stromal.h5ad` |
| `reannotated_migrated/human_epithelial.h5ad` | `data/interim/reannotated/human_epithelial.h5ad` |
| `reannotated_migrated/human_immune.h5ad` | `data/interim/reannotated/human_immune.h5ad` |
| `reannotated_migrated/human_perivascular.h5ad` | `data/interim/reannotated/human_perivascular.h5ad` |
| `reannotated_migrated/human_smc.h5ad` | `data/interim/reannotated/human_smc.h5ad` |
| `reannotated_migrated/mouse_stromal.h5ad` | `data/interim/reannotated/mouse_stromal.h5ad` |
| `reannotated_migrated/mouse_epithelial.h5ad` | `data/interim/reannotated/mouse_epithelial.h5ad` |
| `reannotated_migrated/mouse_immune.h5ad` | `data/interim/reannotated/mouse_immune.h5ad` |
| `reannotated_migrated/mouse_vascular.h5ad` | `data/interim/reannotated/mouse_vascular.h5ad` |
| `reannotated_migrated/mouse_smc.h5ad` | `data/interim/reannotated/mouse_smc.h5ad` |
| `reannotated_migrated/mouse_perimetrium.h5ad` | `data/interim/reannotated/mouse_perimetrium.h5ad` |
| `refined_human_migrated.h5ad` | `data/interim/refined_human.h5ad` |
| `refined_mouse_migrated.h5ad` | `data/interim/refined_mouse.h5ad` |
| `final_human_migrated.h5ad` | `data/interim/final_human.h5ad` |
| `final_mouse_migrated.h5ad` | `data/interim/final_mouse.h5ad` |
| `final_human_fib_dec_migrated.h5ad` | `data/interim/final_human_fib_dec.h5ad` |
| `final_mouse_fib_dec_migrated.h5ad` | `data/interim/final_mouse_fib_dec.h5ad` |

Deposit mappings, from `dataverse/objects_migrated/xmens_annotation_v1_20260810/` to `dataverse/objects/`, preserve each filename:

- `2026_08_05_Cevrim_Human_All_Cells.h5ad`
- `2026_08_05_Cevrim_Human_Only_Fib_Dec_Cells.h5ad`
- `2026_08_05_Cevrim_XMens_All_Cells.h5ad`
- `2026_08_05_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad`

Before promotion, back up all 21 existing destinations into a dated recovery directory. Use same-filesystem temporary names and atomic renames; do not copy directly over an H5AD. Live destinations must also be confirmed against the promotion manifest. This agent deliberately did not write them.

To discard this run after reviewing the paths, remove only these versioned staging roots (never an unversioned parent):

```bash
set -euo pipefail
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final
test "$(pwd -P)" = /n/lab_storage/eddy_lab/users/nhilgert/xmens_final
find data/interim/migrated/xmens_annotation_v1_20260810 \
     dataverse/objects_migrated/xmens_annotation_v1_20260810 \
     logs/migration/xmens_annotation_v1_20260810 -maxdepth 2 -print
rm -rf /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/migrated/xmens_annotation_v1_20260810
rm -rf /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/dataverse/objects_migrated/xmens_annotation_v1_20260810
rm -rf /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/logs/migration/xmens_annotation_v1_20260810
```

The code/CSV commits are separate from staged data; revert them only with an explicit reviewed Git operation.

## Follow-ups

- SAMap was intentionally deferred. The annotation migration does not require it, and existing SAMap objects remained protected. After promotion, regenerate SAMap into a new run directory using canonical `cell_type`/`cell_subtype` inputs and raw counts in `.X`.
- Human review/polish of all 11 DRAFT de novo notebooks.
- Integrate the two pseudotime notebooks under stable, unique, sample-aware cell IDs before rerunning transferred pseudotime; the previously audited 1-dpi object contains ten ambiguous bare-barcode pairs.
- Decide when to remove transition/legacy aliases from downstream packaging after all consumers use `cell_type` and `cell_subtype`.
