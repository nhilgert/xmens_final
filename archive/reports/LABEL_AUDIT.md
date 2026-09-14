# LABEL AUDIT — cell-type clusters and labels, both species
**2026-09-02 evening. Read-only: every number below was computed from the objects/tables on disk (`label_audit.json`); nothing was modified.**

## 0. How labels are produced (from the rules, not memory)

```
03_cluster        leiden ladder on X_scVI (res 0.1 … 10.0)                 -> annotated_{sp}.h5ad carries all of it
04_annotate       leiden_scVI_3.0 x resources/annotations/{sp}_lineage.csv -> coarse_refined (6 lineages incl. "stromal"), EXCLUDE drops cells
04b_reintegrate   per-lineage scVI                                          -> latent_{sp}_{lineage}.h5ad
04c_reannotate    leiden_scVI_3.0 within lineage x {sp}_{lineage}.csv       -> reannotated/{sp}_{lineage}.h5ad   (apply_annotation.py, THE single apply step)
04d_merge         concat lineages, vocabulary rollup                        -> refined_{sp}.h5ad  (cell_type / cell_subtype = canonical snake_case)
04e_reembed                                                                 -> final_{sp}.h5ad -> deliverables
```
`config.annotate_res = leiden_scVI_3.0` for both species. `annotation_references` point at the `20251008_*` submission objects — those are the reference the human stromal "vote" was taken against.

## 1. What the label columns are, and which are redundant

Partition relations were computed for every pair of label columns on 8 objects. Result:

| column | mouse | human | verdict |
|---|---|---|---|
| `cell_type` | 7 | 7 | **canonical coarse** (snake_case) |
| `cell_subtype` | 23 | 16 | **canonical fine** (snake_case); refines `cell_type` everywhere |
| `coarse_refined` | 7 | 6 | mouse: identical to `cell_type` (display-name renaming). human: **coarser** than `cell_type` — see §2b |
| `fine_refined` | 20 | 12 | **legacy** — old codes on mouse (`dec1/dec2/dec3/decAp/eFib1/eFib2/eFib2i/eFib3/eFibP`) and verbose immune names; `""` for every lineage that got no fine label |
| `subtype` (human) | 12 | | identical to `fine_refined` — duplicate |
| `lineage` (human) | | 5 | identical to `origin_lineage` — duplicate |
| `annotation_level` (human) | | 2 | records which cells got a fine label (`subtype` 6,134 vs `cell_type` 60,696). Mouse has no equivalent |
| `stro7_margin` (human fib/dec) | | ~per-cell | a numeric score living in obs as if it were a label; harmless, misnamed |

**The "cross-cutting" flag on `fine_refined` is NOT a real inconsistency.** It is the empty string: on mouse the 11,745 `""` cells are exactly epithelial + endothelial + SMC + perimetrial (6,265 + 3,487 + 1,443 + 550); on human the 60,696 `""` are exactly fib + dec + epithelial_cell + SMC + perivascular_cell. `fine_refined` is simply unfilled for those lineages, and `cell_subtype` fills them. On the mouse fib/dec object `cell_subtype` ≡ `fine_refined` one-to-one (the old-codes-vs-new-names mapping the SAMap heatmap tripped over).

**Recommendation to decide on:** keep `cell_type` + `cell_subtype` as the only label columns on deliverables; drop `fine_refined`, `subtype`, `lineage` (or keep `fine_refined` but rename its values to the new vocabulary so no code can pick up `dec2` again). `origin_lineage` is useful provenance and nests cleanly — keep.

## 2. Real inconsistencies

**2a. Mouse stromal clusters→labels are already on the new vocabulary.** `mouse_stromal.csv` (35 leiden_3.0 clusters → 9 labels) uses `decidual_boundary`, `efib_outer`, etc. The archived "04c applies the OLD label / decidual_border / efib_myometrial" note is stale — resolved. Nothing to do.

**2b. Human perivascular: three columns disagree on 471 cells.**
- `coarse_refined = "perivascular cells"` = 1,805 = `cell_type` perivascular_cell (471) + endothelial_cell (1,334)
- `lineage = "stromal"` = 48,004 = fib + dec (47,533) **+ those same 471 perivascular_cell**
- `lineage = "perivascular"` = 1,334 = only the endothelium
So the 471 `perivascular_cell` cells were clustered *inside the stromal lineage* (`human_stromal.csv` has a `perivascular_cell` label) but rolled up to the perivascular coarse type. They are the "+471 perivascular throwbacks" the pseudotime CONFIG mentions, and they are excluded from the fib/dec deliverable. Decision: is `perivascular_cell` a stromal subtype (then `coarse_refined` is wrong) or a vascular one (then `lineage` is wrong)? Mouse has no analogue — its `perimetrial_cell` is its own lineage.

**2c. Cross-species coarse vocabulary differs for the vessel compartment.** Mouse: `vascular cells` / `endothelial_cell` (+ `endothelial_blood` subtype). Human: `perivascular cells` covering endothelium (`endothelial_blood`, `endothelial_lymph`) AND `perivascular_cell`. If SAMap/any cross-species table uses `cell_type`, `endothelial_cell` matches; if it uses `coarse_refined`, mouse "vascular" vs human "perivascular" do not. Decision: one name.

**2d. Human stromal: the vote vs the dotplot markers (`resources/lineage_spec.yaml`).** All 34 within-lineage
leiden_3.0 clusters were scored on the panel. Finding first: **of the four decidual markers only LEFTY2 discriminates.**
PRL is detected in ≤5% of cells in every cluster and IGFBP1 in ≤25% (usually <10%); PDGFRA and DCN are pan-stromal
(cluster means 1.1–1.9 on both sides), so the "fibroblast" half of the panel does not separate either. LEFTY2
detection runs 0.03 → 0.70 and tracks the vote: dec-voted clusters 0.25–0.70, fib-voted 0.03–0.16.

**There are two distinct ambiguous populations, and they barely overlap** (corrected 2026-09-02; an earlier draft of
this section wrongly said the low-purity votes were mostly the LEFTY2-intermediate band — the overlap is 4 of 9
clusters, 1,453 cells, 18% of low-purity cells):

| | clusters | cells | signature |
|---|---|---|---|
| **A — LEFTY2-intermediate** | 3, 9, 10, 12, 19, 22, 25, 29 | 6,058 | LEFTY2 detection 0.20–0.35, between the two bodies; votes split 4/4. Weak decidual signal |
| **B — double-positive** | 2, 5, 8, 11, 24 | 6,619 | LEFTY2 near-decidual (0.38–0.62) **and** fibroblast/ECM markers at the top of the range: mean DCN percentile **0.92** vs **0.26** for the high-purity decidual clusters. Co-expression, not weak signal |

**B is 82% of the low-purity cells** (6,619 of 8,072) — so the reference's disagreement is mostly on cells that carry
*both* programs, not on cells with little signal. That is what a transitional/boundary state should look like:
retained ECM program alongside an acquired decidual marker. A being LEFTY2-intermediate is the other, quieter kind.
Any single-axis threshold on decidual markers separates A and misses B entirely.

Outright contradictions: **24** (154 cells, fib-voted, LEFTY2 0.51) and **11** (1,306, fib-voted, 0.38) look decidual
by marker and are both class B. Cluster 13 (471) is the perivascular set (ACTA2/PECAM1 high), consistent with §2b.

**Cluster 33 — RESOLVED 2026-09-03, stays `decidual_cell`; user decision, and the evidence agrees.** An earlier
draft of this section called it "low-count/low-quality cells, marker support in neither direction." That was wrong:
its FOXO1 detection is **0.781**, the highest of the 13 fibroblast-range clusters and above decidual-voted clusters
8 (0.733) and 21 (0.736) — FOXO1 being the panel's best decidual discriminator (AUROC 0.716). Only LEFTY2 (0.069)
and the fibroblast-panel detection (0.465) were low, and I generalised from those two to "every marker."
What it actually is, from the spatial run (`cluster_spatial_stats.csv`, `cluster33_markers.csv`):

- **98% of its placed cells sit in decidual-dominated neighbourhoods, 0% in fibroblast-dominated**; median isodepth
  **0.947, the highest of all 34 clusters** — the deepest decidual pole, in tight compact clumps, not scattered.
- **Menstrual-specific**: 770 of 781 placed cells in the three menstrual sections, 11 in secretory.
- Markers are hypoxia + stress: VEGFA, NDRG1, P4HA1, PLOD2, STC1 (HIF targets), HSPA1A/B, HSPB1, HSPH1, HSP90AA1,
  NEAT1, MALAT1, SQSTM1, MMP10, INHBA. No mitochondrial or ribosomal genes in the top 40; depth is lower
  (n_genes 1,252 vs 2,021 decidual) but not degenerate.
- Reading: decidual cells in the menstrual **breakdown/shedding zone** — LEFTY2 is off because they have exited the
  differentiation program, not failed to enter it. Plausibly the human counterpart of mouse `decidual_apoptotic`,
  which human was noted as lacking. NOT scored against a hypoxia module yet — that check would confirm or kill it.

**Trap:** `human_decidual_boundary.ipynb`'s marker-threshold call flips cluster 33 to `endometrial_fibroblast`
(it is the only outright flip). That call is known-wrong; the notebook now pins 33 as decidual with this reason.
Tables: `human_stromal_lefty2_table.csv`, `human_stromal_cluster_marker_{means,fracs}.csv`.

Decision: (i) keep the vote; (ii) re-vote on LEFTY2 with a threshold — separates A, leaves B mislabelled; (iii) call
class B a **boundary/transitional subtype** (`decidual_boundary`, matching the mouse vocabulary) and decide A on the
LEFTY2 cut — this is the option the marker evidence actually points at, and the first step toward a human fine
vocabulary. `human_decidual_boundary.ipynb` tests it: whether B sits at the fib/dec interface in space, replicates
across samples, and lands where mouse `decidual_boundary` lands. Nothing changed yet.

**2e. `decidual_proliferative` — CORRECTED 2026-09-02 evening: it IS proliferating; keep the name.** An earlier
audit reported a median G2M score of ~0 and called it "not cycling"; that was a median on a heterogeneous population.
Recomputed on the deliverable (Tirosh S/G2M, log-CP10k): **Mki67 detected in 27%, Top2a in 50%, G2M > 0.1 in 36%**
of decidual_proliferative cells vs 2–4% / 13–14% / 2% in decidual_mature and decidual_boundary (efib_proliferative: 61% / 80% / 76%).
Its top markers vs the other decidual subtypes are mitotic (Top2a, Knl1, Kif15, Smc4, Anln, Hjurp, Ect2, Kif11, Cit, Atad2).
Cycling declines over the timecourse (Mki67 detected 47% at 1dpi_early → 18% at 4dpi). Files: `cycling_by_subtype.csv`,
`decprolif_by_timepoint.csv`, `decprolif_markers_vs_decidual.csv`. **No rename.**

## 3. Cluster exclusions (the evidence-layer supplement)

Cells dropped between annotate and refined, by `EXCLUDE` rows in the tables:

| species | annotated → refined | excluded | by lineage |
|---|---|---|---|
| mouse | 58,354 → 57,033 | **1,321** | epithelial 416 · vascular 475 · SMC 326 · perimetrium 104 · stromal 0 · immune 0 |
| human | 69,146 → 66,830 | **2,316** | epithelial 1,480 · immune 357 · stromal 471* · perivascular 8 · SMC 0 |

\* human "stromal 471" is the perivascular re-homing in §2b, not an exclusion.

Tables with `EXCLUDE`: mouse lineage 9, epithelial 2, perimetrium 4, SMC 4, vascular 3; human lineage 3, epithelial 4, immune 3, perivascular 6. Where a `note` column exists the reason is recorded (mixed markers, ambient, immune contamination); the `human_immune`, `human_lineage`, `mouse_*` two-column tables have **no reason recorded**. For the Scully-style defense each exclusion needs: the markers that flagged it, replication across samples, and the doublet score. Decision: build that audit (all inputs are on disk; it is a read-only job) or not.

## 4. Housekeeping (no decision needed, just approval)

- `resources/annotations/_BROKEN_fib_metadata_15972of43052_rows.csv` — a partial legacy metadata dump (15,972 of 43,052 rows, old `cell_subtype_annotations` = `dec2` etc.). Not read by any rule. → `_attic/` with a MANIFEST line.
- `uns['display_names']` is byte-identical on mouse and human and includes `decidual_apoptotic`, which human does not have. Cosmetic.
- `annotation_level` exists on human only; either add to mouse or drop.
- `stro7_margin` on the human fib/dec object: a score, not a label — rename or move to `obsm`.

## 5. Decision list (nothing pre-decided)

1. Label columns: **decided** — `cell_type`+`cell_subtype` canonical; `coarse_refined`/`fine_refined` ignored (drop from deliverables when next re-cut).
2. Human 471 `perivascular_cell`: **decided — stromal** (a stromal subtype, not fib or dec). 'vascular' is not a lineage.
3. Coarse vessel vocabulary: moot once `coarse_refined` is dropped; `cell_type` already agrees (`endothelial_cell`).
4. Human stromal: keep vote / LEFTY2 re-vote (class A only) / `decidual_boundary` for the double-positive class B, 6,619 cells (§2d) — see `human_decidual_boundary.ipynb`.
5. ~~Rename `decidual_proliferative`~~ — closed, it is cycling (§2e).
6. Cluster exclusions: **confirmed as mixed-marker contamination** by the user. Evidence audit (markers per excluded cluster) still optional.
7. Housekeeping in §4.

Every change above flows through `resources/annotations/*.csv` + `annotation_vocabulary.yaml` and a re-run of 04c→04d→04e; the deliverables would be re-cut as new dated files (never overwritten), per standing rule.
