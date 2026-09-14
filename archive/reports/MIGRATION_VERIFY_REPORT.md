# Annotation migration: biological verification and promotion preparation

**Run:** `xmens_annotation_v1_20260810`

**Repository:** `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final`

**Staged run root:** `data/interim/migrated/xmens_annotation_v1_20260810/`

**Audit date:** 2026-08-11

**Disposition:** biologically and structurally ready for human review; **nothing was promoted**.

This audit was read-only for every staged and live H5AD. It submitted no GPU job, ran no SAMap fit, and changed no live object. The only authored artifacts are this report and the promotion-preparation files under `tools/`.

## 1. Biological verification of the vascular split

Expression was read from the staged final objects only:

- `data/interim/migrated/xmens_annotation_v1_20260810/final_human_migrated.h5ad`
- `data/interim/migrated/xmens_annotation_v1_20260810/final_mouse_migrated.h5ad`

In these files, `.X` is per-cell natural-log normalized expression,
`ln(1 + 10,000 × raw_count / library_total)`. The tables report the arithmetic mean of that value across the exact requested group. Detection fractions were independently checked against `layers['counts']`; `.X > 0` and raw-count detection agreed exactly.

### Human mural versus endothelial identity

| Exact staged group | n | PDGFRB | RGS5 | NOTCH3 | ACTA2 | PECAM1 | VWF | CDH5 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `perivascular_cell` | 471 | 1.5417 | 0.7735 | 1.0253 | 1.2364 | 0.1403 | 0.0551 | 0.0176 |
| `endothelial_cell` | 1,334 | 0.0419 | 0.0408 | 0.0212 | 0.0334 | 1.7425 | 1.7299 | 0.7215 |

The 471 mural cells show the expected reciprocal pattern: PDGFRB was detected in 66.88%, NOTCH3 in 52.65%, ACTA2 in 50.96%, and RGS5 in 36.73%, whereas only 1.27–7.86% expressed the endothelial markers. The 1,334 endothelial cells were PECAM1-positive in 80.28%, VWF-positive in 62.52%, and CDH5-positive in 47.98%, while each mural marker was detected in only 1.35–3.00%.

**Verdict:** the 471 cells assigned `perivascular_cell` are biologically mural/pericytic, and the 1,334 assigned `endothelial_cell` are biologically endothelial. There is no material inversion or contamination signal.

### Human blood versus lymphatic endothelial identity

| Exact staged subtype | n | PROX1 | LYVE1 | PDPN | FLT4 |
|---|---:|---:|---:|---:|---:|
| `endothelial_blood` | 828 | 0.0510 | 0.0080 | 0.0110 | 0.1789 |
| `endothelial_lymph` | 506 | 2.0318 | 0.1365 | 0.5439 | 1.0514 |

The lymphatic group was PROX1-positive in 85.77%, FLT4-positive in 65.22%, and PDPN-positive in 40.71%, versus 3.99%, 14.01%, and 0.85% in blood endothelium. LYVE1 was more modest (9.49% in lymphatic versus 0.72% in blood), but the coherent PROX1/FLT4/PDPN separation strongly supports the label. Blood endothelial cells were much more VWF-high than lymphatic cells (mean 2.5425 versus 0.4000; 85.27% versus 25.30% detected).

**Verdict:** the 828 blood and 506 lymphatic endothelial assignments match their biology. Modest LYVE1 alone is not a red flag because three independent lymphatic markers separate strongly.

### Mouse blood-only endothelial decision

The exact staged mouse group `cell_type == endothelial_cell` and `cell_subtype == endothelial_blood` contained 3,487 cells.

| Marker | Mean `.X` | Cells detected | Mean raw counts/cell |
|---|---:|---:|---:|
| Pecam1 | 2.3675 | 85.52% | 4.0731 |
| Cdh5 | 1.4245 | 62.83% | 1.9690 |
| Prox1 | 0.0169 | 1.00% | 0.0115 |
| Lyve1 | 0.0039 | 0.23% | 0.0023 |
| Pdpn | 0.0050 | 0.40% | 0.0046 |

**Verdict:** mouse endothelium is convincingly blood endothelial (Pecam1/Cdh5-high) and essentially negative for the lymphatic program. A fabricated mouse lymphatic subtype is not supported; the blood-only decision is justified.

## 2. The two edits beyond value renaming

### `workflow/scripts/integrate.py`

Exact diff against pre-migration commit `8c84132`:

```diff
@@ -50,6 +50,8 @@ SCVI_NEIGHBORS_KEY = f"neighbors_scVI_{N_NEIGHBORS}"
 SCVI_UMAP_KEY = "X_umap_scVI"

+from pathlib import Path as _Path
+_Path(logfile).parent.mkdir(parents=True, exist_ok=True)
 with open(logfile, "w") as log:
@@ -147,5 +149,13 @@ with open(logfile, "w") as log:
     say(f"final: {adata.n_obs} cells, obsm={list(adata.obsm.keys())}")
-    adata.write_h5ad(out)
+    import os
+    from pathlib import Path
+    out_path = Path(out)
+    out_path.parent.mkdir(parents=True, exist_ok=True)
+    if out_path.exists() and "migrated" in out_path.parts:
+        raise FileExistsError(f"refusing to overwrite existing object: {out}")
+    tmp_path = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp.h5ad")
+    adata.write_h5ad(tmp_path)
+    os.replace(tmp_path, out_path)
     say(f"wrote -> {out}")
```

**Verdict: `LOGIC-CHANGED (I/O/plumbing only; scientific algorithm unchanged)`.** This is not merely cosmetic: it creates missing log/output directories, refuses an existing migrated target, and changes H5AD writing from direct to temporary-plus-atomic replacement. It does **not** change HVG selection, scVI architecture/training parameters, batch handling, neighbor construction, Leiden resolutions, UMAP, or the resulting intended AnnData content.

One caveat: the no-overwrite guard applies only when a path component is literally `migrated`; it is not a universal live-output guard. That does not affect the successfully completed staged run, but it should not be described more broadly.

### `config/config.yaml` reference-key propagation

The duplicate YAML key was split cleanly:

| Key/read | Location | Resolution |
|---|---|---|
| `genome_references` definition | `config/config.yaml:12` | Exists; currently documentation/provenance only, with no tracked code read under `workflow/`, `tools/`, or `config/`. |
| `annotation_references` definition | `config/config.yaml:123` | Exists. |
| `config["annotation_references"][w.species]` | `workflow/rules/04d_merge.smk:16` | Resolves. |
| `config["annotation_references"][w.species]` | `workflow/migration/Snakefile:74` | Resolves. |
| old `config["references"]` | exhaustive tracked search | **No reads and no remaining definition.** |

The annotation reference is used only for an informational legacy-label agreement comparison during merge; it does not assign labels, filter cells, or determine the count gates. There is no dangling `config["references"]` read that would raise a KeyError.

## 3. Promotion readiness and recovery design

### Authored, not executed

- `tools/promote_migration.sh`
- `tools/rollback_migration.sh`
- `tools/migration_promotion.py` (shared fail-closed implementation)
- `tools/promotion_plan_xmens_annotation_v1_20260810.json`

The canonical inner promotion-plan SHA-256 is:

`5c0dda6041034be2a6d941efe2430b028d143eaa0b1e40b055164972d0503ad3`

It binds exactly 21 unique, ordered source/destination pairs: 11 reannotated lineage objects, two refined atlases, two final atlases, two fib/dec subsets, and four deposits. It also pins 18 status/schema/audit inputs. The staged sources total 22,766,694,636 bytes; the current 21 live destinations total 24,198,483,771 bytes. All were stable while their SHA-256 values were recorded.

The implementation is intentionally more conservative than direct copies:

1. Default promotion invocation is a read-only full preflight.
2. Execute mode requires a unique promotion ID, exact reviewed Git HEAD, and typed `PROMOTE:xmens_annotation_v1_20260810` confirmation.
3. It verifies the plan digest, exact 21-file inventory, all 18 pinned inputs, PASS statuses, current old-live SHA/size/metadata, staged SHA/size, five protected non-target SHA values, critical-file Git cleanliness, allowed committed paths since reviewed base `2cdc6a50ef2c84c9c3340a86fafff94e95929fa6`, same-filesystem placement, and free space.
4. Before touching live data, it makes and verifies independent (non-reflink) backups of all 21 old files and prepares **all** 21 promotion and 21 restore files beside their destinations. Its free-space gate also reserves another full old-object set for a manual rollback after a hard crash: about 95 GB of object copies plus a 5 GiB margin. The audited filesystem had ample capacity.
5. It rehashes the live set after preparation, then rechecks each still-pending destination immediately before an `os.replace`.
6. A durable `RECOVERY_REQUIRED` marker and fsynced journal exist throughout the 21-file apply window. Ordinary errors/signals trigger reverse-order automatic rollback. A kill/host loss leaves the marker for explicit recovery.
7. The rollback is manifest-bound and idempotent: each destination must be the reviewed old SHA, reviewed promoted SHA, or old content with repairable metadata. Any third content hash is preserved as a conflict and never overwritten.
8. Promotion and rollback write promotion-aware live manifests, preserve verified backups, and verify the five non-target protected files. No downstream rule is launched.

The wrappers and helper passed shell/Python syntax parsing. **Neither wrapper, the helper, nor any mutation path was executed.** Fault-injection and SIGKILL recovery tests therefore remain **UNVERIFIED** and are appropriate before first production use.

### Human-run sequence after review

Run only during a short maintenance window with writers and readers stopped; there is atomicity per H5AD, but no 21-file atomic transaction. Use a durable compute/login session because the preflight, backup, staging, and verification read/copy well over 100 GB.

```bash
# Read-only; expected to hash the complete reviewed set.
tools/promote_migration.sh --dry-run \
  --expected-git-head <exact-reviewed-current-HEAD>

# Mutation only after the dry-run and human review.
tools/promote_migration.sh --execute \
  --promotion-id <unique-UTC-id> \
  --expected-git-head <exact-reviewed-current-HEAD> \
  --confirm PROMOTE:xmens_annotation_v1_20260810
```

The successful promotion prints the recovery directory and the SHA-256 of its immutable `BACKUP_MANIFEST.json`. Preserve both. A rollback begins with a read-only classification and then requires the printed manifest digest:

```bash
tools/rollback_migration.sh --dry-run \
  --promotion-id <same-id> \
  --expected-backup-manifest-sha256 <printed-sha256>

tools/rollback_migration.sh --execute \
  --promotion-id <same-id> \
  --expected-backup-manifest-sha256 <printed-sha256> \
  --confirm ROLLBACK:<same-id>
```

After any rollback, atomic replacement necessarily changes inodes. Generate a newly reviewed promotion plan rather than attempting to reuse the old plan.

**Promotion blocker found:** no biological/count/schema blocker remains. Operationally, promotion should wait for human code review, committed/executable scripts, the read-only promotion preflight, a maintenance window, and—if desired—a small disposable fault-injection rehearsal. Do not promote by hand around these gates.

## 4. Downstream impact

### Critical dry-run finding

A normal main-workflow `samap_all` call is unsafe after promotion. A real read-only dry-run at the audited HEAD scheduled **28 upstream live-DAG jobs**—2 whole-atlas integrations, 2 atlas annotations, 11 GPU lineage reintegrations, 11 lineage reannotations, and 2 merges—while scheduling **zero SAMap jobs**, because the existing August 10 SAMap outputs were treated as current.

`ancient(refined_*)` suppresses the refined timestamp as a trigger for SAMap preparation; it does **not** prevent Snakemake from resolving and rebuilding the rule that produces refined. Therefore:

- **Do not run unrestricted `samap_all`.**
- **Do not run the 28 upstream jobs after promotion.** The staged refined/final objects are already the validated products.
- Existing SAMap inputs/fits/stains/transferred objects are pre-migration products and do need a deliberate refresh after promotion.

This restricted dry-run boundary was verified to produce exactly 11 intended jobs and no stage 00–04 producers:

```bash
snakemake -s workflow/Snakefile --cores 1 -n -p samap_all \
  --allowed-rules prep_samap_input run_samap samap_label_stain \
                  samap_transferred_object samap_all \
  --forcerun prep_samap_input
```

It resolves to four preparation jobs, two SAMap fits, two label stains, two transferred objects, and the aggregate target. This is a **dry-run recipe only**. The current real rules would overwrite protected August 10 products, so first give the refresh a new immutable input/run/result root (or a dedicated standalone Snakefile), then repeat the restricted dry-run before execution.

### Consumer inventory

| Consumer | Current annotation input/behavior | Post-promotion action |
|---|---|---|
| `merge_reannotated` → refined | Reads the 11 live reannotated objects plus canonical `cell_type`/`cell_subtype` | **Do not rerun.** The staged refined objects are already the validated merge result. |
| `reembed` → final | Reads refined counts/sample; carries annotations but does not use them to change the scVI algorithm | **Do not rerun.** The staged final objects are already validated. |
| SAMap prep/fit/stain/transfer | Prep reads refined canonical labels/counts; `ancient()` prevents automatic invalidation | Regenerate deliberately in a **new immutable output root**: 4 prep, 2 fit, 2 stain, 2 transfer. Never run unrestricted `samap_all`. |
| `05_de` | Reads legacy `annotated_{species}.h5ad`; implementation is touch-only | Implement and rewire to the intended final/fib-dec canonical contract before running. No force now. |
| global pseudotime / CellRank / PT MI | Rules read legacy annotated objects or their stub outputs; implementations are touch-only | Implement/rewire first. The separate manual mouse PT notebooks still use legacy `coarse_refined`/`fine_refined` and must be migrated before regeneration. |
| `22_crossspecies` | Reads MI CSVs, not H5AD annotations; implementation is touch-only | Implement upstream MI and crossspecies first; do not run now. This rule is unrelated to SAMap despite its name. |
| `23_invitro` | Reads legacy `annotated_human.h5ad`; implementation is touch-only | Implement/rewire before running. |
| `90_figures` | Declares marker/MI/conservation inputs only; omits final/PT/SAMap dependencies; implementation is touch-only | Expand dependencies and implement before running. |

A broader annotation/downstream dry-run scheduled **39 jobs**: the 30 stage 02–04e rebuilds plus 9 placeholder downstream jobs. It must not be turned into a real unrestricted run. The main `rule all` does not currently target refined, final, or SAMap products, so promotion alone does not make the main DAG current.

## 5. Final disposition

- Human vascular split: **PASS, biologically coherent.**
- Mouse blood-only endothelial assignment: **PASS, biologically justified.**
- `integrate.py`: **LOGIC-CHANGED in I/O safety only; scientific integration logic unchanged.**
- Reference-key split: **PASS; no dangling old-key read.**
- Staged/live plan binding: **PASS, 21/21 with canonical plan digest.**
- Promotion/rollback package: **written and statically parsed; never executed.**
- Live objects: **untouched.**
- GPU/SAMap execution: **none.**
- Next decision: human review of the promotion package, then read-only preflight; promotion remains a deliberate separate action.
