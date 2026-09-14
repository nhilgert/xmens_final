# ANSWERS — Hammoud open questions (for Claude Code)
2026-09-01. Everything below was verified against files on the cluster; sources named inline.

## 1. The three exports you requested — done, next to the driver JSON

All in `xmens_final/hammoud/fib_centroid_results_mouse/`:

| file | what | verification |
|---|---|---|
| `hammoud_h_cluster_assignments.csv` | 47,533 barcodes -> h0..h9 | per-cluster counts match `hammoud_gene_drivers_v5.json` n_cells for ALL 10 clusters (checked before writing) |
| `hammoud_driver_computation.py` | the EXACT executed script, extracted verbatim from the producing job's cmd.sh (job 414e5b83, run 2026-08-31 05:57, env rsc_uv) | not a reconstruction; header documents panels, centring, frac_from_top30 |
| `hammoud_h_provenance.txt` | the one-liner below | recovered from the writer job's script (a7b1239b), not memory |

## 2. h0..h9 provenance — RESOLVED (was "unknown")

```
object     data/deliverables/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad
           (post embedding-fix X_scVI = the 04b stromal latent)
neighbors  sc.pp.neighbors(use_rep='X_scVI', random_state=0)   # default n_neighbors=15
leiden     sc.tl.leiden(resolution=1.0, random_state=0, flavor='igraph', n_iterations=2)
column     leiden_1.0, written to results/human_fib_annotation/fine_cluster_assignments.csv
           by job a7b1239b (2026-08-30 ~16:09)
labels     driver job relabeled 'h' + leiden_1.0
```

**Consequences for the notebook:**
- Retire the "unsupportable to tie h-clusters to our clusters" caveat. h0..h9 IS a clustering
  of the human deliverable itself (all 47,533 cells, resolution 1.0), sibling to the 6-cluster
  leiden-0.5 arm — not an alien vocabulary.
- Your statement "no stored file maps h0..h9 to the 6-cluster solution" is now false in the
  useful direction: `fine_cluster_assignments.csv` carries BOTH `hammoud_res05` AND
  `leiden_1.0` per cell. The h <-> 6-cluster crosstab is one pivot. Add it to section 4.
- Keying cross-species driver comparisons on Stro subtype remains correct (the two species'
  query clusterings are still different granularities), but per-h statements are now supportable.

## 3. Decision 1 (atlas scope) — CONFIRMED: (a) with (b) retained

- Primary mouse table: re-argmax of the saved similarity matrix restricted to the 6 Stro
  subtypes (pure read, no refit).
- The all-39 table appears once, as the specificity finding (stromal queries do not
  preferentially match stromal atlas subtypes when offered the whole atlas).
- **eFib2i is reported UNRESOLVED** — restricted argmax flips on a 0.0017 gap
  (Stro-5 0.1321 vs Stro-10 0.1304), below any credible floor. Not either match.

## 4. Decision 2 (missing mouse analyses) — CONFIRMED: run none

Label replication / marker ablation / spatial readout as "not run", and resolution
stability as "N/A by design" (curated labels), exactly as you proposed. Rationale on
record: a rebuttal deliverable should not gain never-reviewed analyses. If a reviewer
pushes later, your cost estimates (one ~3-5 min SLURM job) are noted.

## 5. Panel provenance — your finding, now with the numbers pinned

Two legitimate cosines exist per correspondence and MUST be labeled by source:
- saved matrices: full shared-ortholog panel (15,606 genes, mouse; human 3,000-gene variance panel)
- driver JSON: atlas-derived panels (human n=3,099; mouse n=2,728 via strict 1:1 orthologs)

Example: eFib1->Stro-7 = 0.174 (matrix) and 0.243 (JSON). A cross-file verify() is a
manufactured MISMATCH — verify each number against its own source, panel stated.
**The decidual matches exist only on the variance/atlas panels (dec3 negative, dec2 ~ -0.003
on the full panel). This is a stated caveat in the rebuttal section, not a footnote.**

## 6. Trap addendum (beyond brief section 3)

- `frac_from_top30` is a RATIO (top-30 contribution sum / total cosine), not a fraction.
  It is 1.519 for decidual_apoptotic and -3.856 for human h3 because their cosines are ~0.
  Do not report it where |cosine| is small; the computation header now says this too.
- Driver gene symbols are UPPERCASE (human space) even for the mouse arm — mouse genes were
  mapped through orthologs. Say so once in section 4 or readers will think they are mouse genes.
- Driver-class calls from our audit, to state in section 4c: decidual_apoptotic->Stro-8 is
  housekeeping-carried (7/15 ribosomal/translation, cosine 0.019 -> report as NO MATCH);
  decidual_proliferative->Stro-11 and efib_proliferative->Stro-11 are cell-cycle-driven
  (state match, not identity match); efib_outer->Stro-7 is identity-driven
  (PTGIS, FAP, CCN5, EPHX1, FRZB) and is the credible three-way anchor.
- human h3->Stro-10 cosine is NEGATIVE (-0.0074): an argmax with nowhere to go. Report as
  NO MATCH alongside decidual_apoptotic.

## 7. Unchanged

Definition of done: generated notebook + `*_executed.ipynb`, 0 MISMATCH, 0 cell errors,
every headline number recomputed by a verify() against its own source file. A failing
verify() is a result to report, not a bug to hide.
