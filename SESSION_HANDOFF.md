# ▶ START HERE — startup instructions for a new session
**Updated 2026-09-01 (evening). This section is the entry point; read it before anything else.**

## Who you are working with
A grad student leading the computational analysis for a paper in revision (cross-species
decidualization: human menstrual endometrium × induced mouse). They know the biology and the
data far better than you; they do not always remember which of many objects/columns carries what,
and they will ask you to explain lineage. When they ask "where does X come from", **read the
code and answer from it** — never from recollection.

## The one rule
> **Artifacts are trustworthy. Prose about artifacts is not — including this document.**
> Before any number goes into a figure, slide, message, or manuscript: read it back from the file.

This session (Aug 28 – Sep 1) produced ~20 corrections. Every one moved toward a tidier story
than the data supported, and every one came from writing a sentence before re-reading the table.
Late in the session the previous Claude's own errors increased as its history compacted into
summaries — that is why you are reading this in a fresh conversation.

## Startup checklist (do these in order, ~15 minutes)
1. **Read this whole file** (`xmens_final/SESSION_HANDOFF.md`). §0–§8 below are the science and
   the traps; the 2026-09-01 ADDENDUM at the end of this section supersedes anything above it
   that disagrees.
2. **Confirm the cluster:** `ssh:harvard-eddy` via `host.compute.create`. Interactive
   `call_command` has a ~60 s ceiling — use `submit_job` (SLURM `-p eddy --account=eddy_lab`)
   for anything heavier. Env: `/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python`.
   Login-node `/tmp` is NOT compute-node `/tmp` — embed scripts in the job (base64) rather than
   staging them in `/tmp`.
3. **Check for live jobs** before submitting: `host.compute.ledger()` and `squeue -u nhilgert`.
   As of this writing: **no live jobs.**
4. **Run `host.artifacts(search=...)` before rebuilding anything.** Almost every analysis you
   might think of has been done; the artifact store and `xmens_final/notebooks/verify/` hold it.
   The user has explicitly said: *"please read the whole conversation"* — i.e. consult prior work
   before building.
5. **Before quoting a headline number**, open the verification notebook that recomputes it
   (`notebooks/verify/01..07`, staged inputs in `notebooks/verify/data/`), or the JSON it came from.
6. **Before touching any object**, follow the house write pattern (§ Standing rules → "house
   pattern" in the addendum). Never `rm`; back up by copy/move.
7. **A separate Claude Code agent is working on the Hammoud atlas-mapping notebook** from
   `xmens_final/hammoud/COMPILATION_BRIEF.md` + `hammoud/ANSWERS.md`. Do not duplicate that
   work; if the user relays a question from it, answer from files and append to `ANSWERS.md`.

## Where things are (quick map)
```
xmens_final/                                   the project root (lab_storage)
  SESSION_HANDOFF.md                           THIS FILE
  data/deliverables/2026_08_28_Cevrim_{XMens,Human}_{All_Cells,Only_Fib_Dec_Cells}.h5ad
                                               the 4 canonical objects (mouse is "XMens", not "Mouse")
  data/interim/annotated_{mouse,human}.h5ad    post-clustering/annotation pipeline objects
  data/interim/qc_spatial_{mouse,human}.h5ad   04f spatial-QC outputs (see addendum)
  results/isodepth/{mouse,human}_fibdec_isodepth.h5ad   session-built working objects
                                               (isodepth, dpt_pseudotime, QC coords) — NOT pipeline products
  results/tables/spatial_qc_*                  04f tallies + per-cell verdicts
  notebooks/verify/01..07_*.ipynb (+ data/)    the verification suite — source of truth
  workflow/rules/*.smk, workflow/scripts/      snakemake pipeline (00a → 99)
  hammoud/                                     atlas-mapping thread (Claude Code's domain)
/n/netscratch/eddy_lab/Everyone/nhilgert/      mirrors of the 4 deliverables + SAMap joint exports
                                               (purgeable scratch; lab_storage is the real copy)
/n/lab_storage/.../human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv
                                               the ortholog table (LONG format — see traps)
```

## How the user likes to work
- Approve before launching; back up before writing; explain lineage from code.
- Prefers per-gene, checkable claims over aggregate scores.
- Will push on any arbitrary choice ("how was k chosen?", "what's the background?") — expect it,
  and have the file open before answering.
- Wants new analyses as config-driven pipeline modules, not one-off scripts.
- Job width: *"we can go up to 100 or 200 at once, it's totally fine"* — profile is set to 200.

---

# 2026-09-01 ADDENDUM — what changed today (supersedes §1/§3.1/§6/§8 where they disagree)

## A. Spatial QC — now one vocabulary, pushed everywhere
The pipeline rule `04f_spatial_qc` ("blank, don't drop": far / degenerate / boxed → coordinates set
NaN, cells kept) **first executed today 12:13** (mouse). Everything built before today used the
raw `x_um_dbscan`/`y_um_dbscan`. Today's push (job 7a03e7bf, verified by read-back) added to all
six working objects:
- `x_um_qc`, `y_um_qc` — QC'd coordinates (NaN where flagged)
- `spatial_qc_reason` — kept / far / degenerate / boxed / unmapped
- `uns['spatial_qc_push']` — source + notes; `x_um_dbscan` preserved untouched; `.pre_qc_push.bak` beside each

| object | cells with QC'd coords |
|---|---|
| XMens All_Cells / Fib_Dec | 29,544 / 23,570 |
| Human All_Cells / Fib_Dec | 45,730 / 33,247 |
| mouse / human isodepth working objects | 23,570 / 33,247 |

- **Mouse:** full 04f. Tally: far 776, **degenerate 1,290** (1,214 of them in 1_dpi_late), boxed 76.
- **Human:** a *minimal* run (degenerate only, thresholds parsed from the 04f source at run time;
  far/boxed deferred until the slide-tags positioning rerun, per the rule's own docstring). Found
  **4** degenerate cells — human is essentially clean. Object self-identifies via `uns['spatial_qc_minimal']`.
- Netscratch mirrors of the 4 deliverables re-synced in the same job (size-checked).
- **STILL OPEN:** downstream consumers default to `x_um_dbscan` (`spatial_pt_gp.py --x-key`, any
  plotting). Repointing to `x_um_qc` is the "separate, deliberate edit" 04f reserved for the user.
  Not yet done. Ask before doing it.

**Three unrelated uses of "DBSCAN" collide in this project — keep them apart:**
1. `x_um_dbscan` — the column NAME; "dbscan" refers to the slide-tags fork's internal cell
   positioning method. It is the RAW coordinate, not a QC'd one.
2. 04f spatial QC — kNN / exact-duplicate / box rules. Not DBSCAN at all.
3. Piece detection (session jobs only) — `DBSCAN(eps=max(6*mean5NN,150), min_samples=12)`,
   pieces ≥120 cells, used inside the boundary-depth validation to separate physically separate
   tissue pieces (uninduced slides each carry TWO uterine-horn cross-sections). Never a pipeline rule.

## B. The 0.447 "outlier" is RETIRED — it was a coordinate artifact
The geometric depth validation (hull-peel distance-from-boundary vs stored isodepth, per piece)
had one unexplained failure since Aug 31: `1_dpi_late_mck12b` blob 1, n=867, ρ=0.447. Today:
- That "piece" is **867 cells at 222 unique positions in a 26 µm window, 143 cells on one bead** —
  a coordinate-decoding stack, not tissue. 04f flags it as `degenerate`.
- Rerun on `x_um_qc` (v3, job 574ff107, `boundary_depth_v3_qc_coords.json`): the fake piece never
  forms; **all 11 real pieces validate, ρ 0.797–0.964** (induced 0.797–0.945; uninduced 0.848–0.964).
- **Paper sentence upgrade:** "geometry recovers the fitted isodepth on every tissue piece
  (ρ 0.80–0.94 induced); the single prior exception was a coordinate artifact of ~1,200 stacked
  cells, detected and removed by spatial QC."
- **Caveat that stays:** the uninduced rows compare geometry against `isodepth_corrected`, which is
  itself geometry-derived — that is consistency, not independent validation. Only the induced rows
  test the GASTON fit. State this.
- **TODO:** update §1 Tier 1, §3.1, notebook `01_isodepth_construction`, and the FRIDAY doc to this.

## C. `tissue_section` — built, then retired from the pipeline plan
The user asked for per-piece coordinate frames so uninduced's two cross-sections don't plot as two
blobs. Built as `tissue_section` + `x/y_um_section` (centroid-origin) from the boundary-depth blob
table and written to the mouse isodepth object (job 927c5e14). Then, once the 04f situation was
understood, the user decided: **"stop using that s0 s1 stuff for the pipeline — just look at
objects post clustering and annotation."** So:
- The three columns were **removed** from the mouse isodepth object in the push job;
  `uns['tissue_section_provenance']` remains with `superseded` set. Backup
  `mouse_fibdec_isodepth.h5ad.pre_tissue_section.bak` still exists (plus a `.bak.partial` from a
  cancelled job — both deletable).
- If per-piece labels are ever needed, they belong in a small rule DOWNSTREAM of 04f on `x_um_qc`.
- Figures `qc_sections_{mouse,human}.png` show all pieces on QC'd coords: mouse 11 pieces
  (uninduced 2+2), human 6 single-piece sections. "p0/p1" in those figures are plot-only piece
  tags; nothing stores them.

## D. Hammoud thread — provenance resolved, exports delivered
- `hammoud/COMPILATION_BRIEF.md` (task), `hammoud/ANSWERS.md` (answers to Claude Code's Phase 0),
  and three exports in `hammoud/fib_centroid_results_mouse/`: `hammoud_h_cluster_assignments.csv`
  (47,533 cells → h0..h9, counts verified against the driver JSON), `hammoud_driver_computation.py`
  (exact executed script from job 414e5b83), `hammoud_h_provenance.txt`.
- **h0..h9 provenance (was "unknown"):** `sc.pp.neighbors(use_rep='X_scVI', random_state=0)` +
  `sc.tl.leiden(resolution=1.0, random_state=0, flavor='igraph', n_iterations=2)` on the human
  Fib_Dec deliverable (post embedding fix); column `leiden_1.0` in
  `results/human_fib_annotation/fine_cluster_assignments.csv`, which ALSO carries `hammoud_res05`.
- Decisions recorded: Stro-restricted re-argmax as primary mouse table, all-39 kept as specificity
  check, eFib2i reported UNRESOLVED; no missing mouse analyses run.
- Driver-class verdicts: decidual_apoptotic→Stro-8 and human h3→Stro-10 are NO MATCH
  (cosine 0.019 / −0.007; housekeeping-carried); the two proliferative→Stro-11 matches are
  cell-cycle STATE matches; efib_outer→Stro-7 (PTGIS, FAP, CCN5, EPHX1, FRZB) is the credible
  identity-driven three-way anchor.

## D2. In-vitro / HECA differential expression — AUDITED Aug 30, FIXES NOT APPLIED
Separate project dir: `/n/lab_storage/eddy_lab/users/nhilgert/cagri_ivd_end_to_end_20260815/`
(limma-voom pseudobulk DE of decidualized vs non-decidualized stroma; frozen, hash-locked, its own
re-freeze protocol in `PENDING_PIPELINE_FIXES.md`). Audit: `HECA_DE_audit.md` in that directory
(installed 2026-09-01; also an artifact). Nothing in the directory has changed since 2026-08-18.
- **Model is correct** (`~ dataset + state` + `duplicateCorrelation(block=person)`; the canonical
  answer to an unbalanced paired design). Do not "fix" the design.
- **Open, biology-changing:** F1 — eligibility filter (>=10 counts in >=7 pseudobulks) drops
  switch-on genes; **IGFBP1 is absent from the production table for this reason.** Fix = arm-aware
  CPM filter. F2 — CD68 top-5 in both fitters = macrophage contamination of tiny (28–106-cell)
  decidual pseudobulks; fix = immune-score pseudobulks / flag CD68-family hits.
- Lesser: F3 `voomLmFit(sample.weights=TRUE)`; F4 phase/state partial confound (secretory-only
  sensitivity run); 3 dormant guard bugs in PENDING_PIPELINE_FIXES still unapplied.
- **2026-09-01 evening: handed to Claude Code.** Brief = `DE_BRIEF.md` in that directory (4 phases,
  each gated: recon+git init → F1/F2/F3 fixes on `_v2` copies → approved cleanup with `_attic/` and
  a MANIFEST, never `rm` → re-freeze only on explicit user instruction). **You do not edit that
  repo yourself** — you are the reviewer.

## D3. Agent-to-agent protocol (you ⇄ Claude Code) — READ IF THE USER SAYS "Claude Code has reported"
There is no direct channel; both agents share the cluster filesystem. In the DE repo:
- `_agent_comms/CLAUDE_CODE_REPORT.md` — Claude Code appends; each entry ends `STATUS: WAITING`.
- `_agent_comms/CLAUDE_SCIENCE_REPLY.md` — YOU append; header `## <UTC> — Phase N: <title>`.
Procedure when the user says a report is in: (1) `cat` the report's latest entry; (2) verify its
load-bearing claims against files yourself (hashes, path map, gene fates) — do not take its
numbers on trust any more than your own; (3) append a reply: decisions, corrections, explicit
"go / stop" for the next phase; (4) tell the user "replied — tell Claude Code to read the reply."
Decisions that need the USER, not you: anything that moves a hash-locked file, the F2 option
(re-gate vs covariate vs flag), and Phase 3 re-freeze. Ask them; do not decide on their behalf.
The same protocol can be reused for the Hammoud thread (`hammoud/ANSWERS.md` was the ad-hoc
precursor); if you do, create `hammoud/_agent_comms/` the same way.

## E. House write pattern (every object write this session used it — keep using it)
1. Copy backup beside the object (`<name>.h5ad.<purpose>.bak`); **refuse if one already exists**.
2. Refuse if the columns you are adding already exist.
3. Assert 100% barcode join of your input table.
4. Fingerprint existing key columns (sums) before and after; assert unchanged.
5. Provenance in `uns[...]` with **string values only** (`json.dump`/h5ad choke on tuple keys —
   this failed a job AFTER its write once).
6. Write to `<path>.tmp`, then `os.replace`.
7. Re-open backed and verify (columns present, counts match, `uns` present).
8. Re-sync netscratch mirrors if a deliverable changed.

## F. Compute state at handoff
**No live jobs.** Jobs referenced today (all terminal, all harvested): 927c5e14 (tissue_section
write), 7a03e7bf (QC push), 574ff107 (validation v3), b6b191b7 (human minimal QC), d3f6bdcf
(cancelled — bad provenance string, object untouched), dd46f3e3 / 4328bfc2 (failed v3 attempts:
NumPy-2 `ndarray.ptp` removal; login-node `/tmp` not visible to compute nodes).

## G. Additional traps learned today
- `ndarray.ptp()` is gone in NumPy 2 → `np.ptp(arr)`. rsc_uv is on NumPy 2.
- `c.download()` returns a dict `{'local_path', 'bytes', 'sha256'}`, not a path string.
- Files written to `/tmp` on the login node do not exist on compute nodes.
- TWO 3-MAD rules exist; do not conflate them. **04f `far` (authoritative):** mean distance to
  k=5 nearest neighbours, log1p-transformed, threshold = expm1(median + 3.0×1.4826×MAD) pooled
  per TIMEPOINT, computed after `degenerate` (>1 cell at one coord, or kNN dist < 10 µm) is
  removed so stacks cannot tighten the MAD. **Earlier `spatial_outlier_fix` (Aug 29, superseded):**
  distance from per-SAMPLE centroid, linear scale, median + 3×1.4826×MAD. Radial vs local-density
  are opposite failure modes — the radial rule structurally could not catch the mck12b stack.
- The driver JSON's `frac_from_top30` is a ratio, not a fraction (1.519 / −3.856 where cosine ≈ 0).

---

# SESSION HANDOFF — cross-species decidualization (human menstrual × induced mouse)
**Written 2026-09-01. Every number below was re-read from its artifact, not recalled.**
Supersedes `FRIDAY_MEETING.md` and `HANDOFF_NEXT_SESSION.md` where they disagree.

---

# 0. READ THIS FIRST — the working rule of this project

This session produced ~15 corrections. **Roughly half were caught by the user, half by audit.
Every single one moved in the same direction: toward a tidier story than the data supported.**

The mechanism was identified and is mechanical, not moral: *the summary sentence got written before
the table was re-read.* The practical consequence:

> **The artifacts are trustworthy. Prose about the artifacts is not, including prose in this
> document. Before you put a number in a figure, slide, or manuscript, read it back from the file.**

That is why the deliverable of the session's last third is not an analysis but **seven verification
notebooks** at `xmens_final/notebooks/verify/`. Each recomputes a claim from source and prints
`[PASS]`/`[MISMATCH]` against the recorded value. Current state: **85 checks, 0 MISMATCH, 0 errors.**

| notebook | recomputes |
|---|---|
| `01_isodepth_construction` | which cells carry a *fit* vs a *smoothed* depth value |
| `02_isodepth_vs_pseudotime` | the two coordinates, cross-species transfer, pseudotime dynamic range |
| `03_mutual_information` | MI estimator, batch stratification, redundancy vs correlation |
| `04_conservation_and_controls` | conservation enrichment end-to-end with all controls |
| `05_samap_integration_readonly` | integration code + normalization (does not refit) |
| `06_rewiring_and_drivers` | Δpeak, regulator verdicts, driver lists |
| `07_samap_gene_conservation` | SAMap's own gene-gene graphs, conservation vs divergence |

`notebooks/verify/data/` holds ~27 staged input artifacts so the suite is self-contained.
**Run the suite before trusting anything here.** If it disagrees with this document, it is right.

---

# 1. THE SCIENCE — what the paper claims

## The sentence

> **Human and mouse endometrial stroma share a conserved spatial axis and a
> transcription-factor-enriched repertoire of genes carrying positional information along it. The
> two species deploy that repertoire differently. A specific 37-gene subset acquires conserved
> positional information only upon decidualization.**

Not "decidualization is conserved." That phrase is unsupportable and a reviewer with these numbers
would say so. The user's own framing, which the data backs: **a lot of reused machinery, without
strong evidence it is reused the same way — heterochrony and rewiring.** Testing showed
predominantly *rewiring*, with a minor timing component.

## Three tiers, strongest first

**TIER 1 — the axis is real and conserved. Nothing threatens this.**
- `isodepth` exists in both species (GASTON fits per sample → quantile-normalized consensus)
- A purely geometric coordinate (distance from myometrium inward) reproduces the fitted values at
  **ρ 0.797–0.945 on all 7 induced pieces** on QC'd coordinates (v3, 2026-09-01). The earlier
  0.447 "outlier" was a coordinate-decoding stack removed by spatial QC — see ADDENDUM §B.
- Correlates with pseudotime *within every section*: human 0.671 overall (secretory 0.551,
  menstrual 0.795); mouse 0.690 overall (per timepoint 0.606–0.840)
- Transfers cross-species: induced pairings **mean ρ 0.688** (n=32, range 0.468–0.822)

**TIER 2 — the parts list is shared, but it is STROMAL, not decidual.**
- Cross-species top-decile MI overlap on isodepth: **402/1,386 vs 138.5 expected = 2.90×,
  p = 3.7e-101** (universe 13,867 strict 1:1 orthologs)
- Nearly identical on pseudotime: 399/1,370 = **2.91×** — the result does not depend on the axis
- Transcription factors are enriched among conserved-positional genes: **OR 2.446, p = 8.7e-5**
  (5.22% of TFs conserved vs 2.37% of all genes)
- **BUT the specificity control kills the decidual framing** (see §3.4): baseline uninduced uterus
  scores **2.581×** against induced **2.377×**. The program is present at full strength *before*
  decidualization.

**TIER 3 — deployment differs. This is the interesting half.**
- Transcriptome-wide verdicts over 9,550 measurable ortholog pairs:
  `NOT_POSITIONAL 7,606 · HUMAN_ONLY 746 · MOUSE_ONLY 720 · CONSERVED_POSITIONAL 226 ·
  CONSERVED_UNCLEAR_SHAPE 157 · INVERTED 95`
- **Divergent categories outnumber the conserved one ~7:1.**
- Of 626 curated regulators (624 measurable): **27 CONSERVED_POSITIONAL, 4 INVERTED
  (MITF, FOXO1, RBPJ, AEBP2), 54 MOUSE_ONLY, 42 HUMAN_ONLY, 483 NOT_POSITIONAL**
- **FOXO1 — the canonical human decidualization factor — is spatially INVERTED between species.**
  Lead with this; a reviewer will find it regardless.
- Peak displacement (Δpeak) across 478 conserved genes: median 0, mean +1.06, SD 4.36,
  **41.6% within ±1 decile, 33.1% displaced ≥5**. Replicated on the pseudotime axis
  (38.9% / 31.0%), and the two axes agree on *which* genes move: **ρ 0.732, same sign 69.7%**.
  → rewiring, not a coherent timing shift. Extreme tails are asymmetric toward human-later
  (+9: 36 genes vs −9: 17), so state "predominantly rewiring with a minor heterochronic component."

## The gene-level rescue — the session's key scientific result

The aggregate specificity test found no decidualization-specific advantage. **The user rejected the
aggregate answer** ("there aren't any decidualization specific genes that are conserved though?")
and the per-gene version found them:

| set | definition | n |
|---|---|---|
| A | human top-decile ∩ induced-mouse top-decile | 221 |
| **B** | **A minus genes already positional in uninduced** | **147** |
| C | B with gain z > 2 | 115 |
| **D** | **C where uninduced MI is at the noise floor** | **37** |
| E | reverse: positional in uninduced, LOST on induction | 31 |

**Set B is 2.327× enriched (147 observed vs 63.2 expected, p = 1.5e-24).**
Set D leaders: **Cped1 (z=+12.0), Runx1 (+10.2), C1qtnf7 (+8.3), Cd63 (+7.5), Ebf1 (+7.4),
Cebpb (+7.3).**

> **The frame this yields: decidualization does not create the axis — it reorganizes which genes
> read it.** That is sharper than either "conserved" or "not specific," and it has named genes with
> z-scores behind it.

---

# 2. THE DATA ON DISK

## Deliverables (the canonical objects)
```
xmens_final/data/deliverables/
  2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad    mouse, 43,197 cells   <- NOT "Mouse_"
  2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad    human, 47,533 cells
```
**Naming trap:** the mouse file is `XMens`, not `Mouse`. A species-name glob returns empty and
killed one job.

## Coordinate objects
```
xmens_final/results/isodepth/mouse_fibdec_isodepth.h5ad
xmens_final/results/isodepth/human_fibdec_isodepth.h5ad
```
Key columns and **what each actually means**:
- `isodepth_quantile_raw` — the **real GASTON fit**. 20,342 mouse cells. **Zero uninduced cells.**
- `isodepth_consensus` — fit + expression-space label spreading. 43,197 mouse cells, incl. 9,334
  uninduced. **Uninduced values here are imputed, not measured.**
- `isodepth_uninduced_measured` / `isodepth_corrected` — added 2026-08-31, the geometric
  boundary-depth substitute (3,958 uninduced cells). Provenance in `uns['uninduced_isodepth_fix']`.
  Backup: `.pre_uninduced_fix.bak`.
- `dpt_pseudotime` — see §5 caveat, its dynamic range is tiny.
- `x_um_dbscan` / `y_um_dbscan` — QC'd spatial coordinates.

## Integration objects
```
/n/netscratch/eddy_lab/Everyone/nhilgert/
  2026_08_31_samap_vanilla_fibdec_joint.h5ad    90,730 × 70,299   <- THE PAPER ARM
  2026_08_30_samap_harmony_fibdec_joint.h5ad    identical shape   robustness supplement
  2026_09_01_samap_secretory_x_1dpi_early_joint.h5ad
  2026_09_01_samap_secretory_x_1dpi_late_joint.h5ad
  2026_09_01_samap_menstrual_x_4dpi_joint.h5ad
```
**The all-samples exports carry BOTH gene-gene graphs** — discovered late, after I wrongly claimed
they carried none:
- `varp['homology_graph']` — sequence-derived (BLAST/phmmer bit scores, 28.3–17,570, median 241)
- `varp['homology_graph_reweighted']` — after expression refinement (0.89 max, median 0.055)
- `uns['gnnm_corr']` (39,843² — a *different* gene space, unusable without its index),
  `uns['homology_gene_names_dict']`, `uns['mapping_K']=20`
- `.X` is SAM's internally processed matrix (row sums 3,372–6,348, max 10.2), **not** the CP10K input

**Transferred columns in the two all-samples exports exist on HUMAN cells only**
(`isodepth_transferred` 34,876 human / 0 mouse). Mouse can be a transfer *source*, never a target
there. The three per-pairing exports (2026-09-01) fixed this — they carry **both directions** plus
per-cell label confidence.

## Per-pairing SAMap runs
```
decidualization/revisions/samap_integrate/runs/{menstrual,secretory}_x_{uninduced,1dpi_early,
    1dpi_late,2dpi,3dpi,4dpi}_hmmer/samap.pkl      12 runs, 8.2–10 GB, written 2026-08-31, READABLE
decidualization/revisions/samap_integrate/runs/full_{blast,hmmer}/samap.pkl   April — UNREADABLE
decidualization/revisions/samap_integrate/homology/hsmm_{blast,hmmer}/         sequence maps
```
**The April all-samples pickles fail with `KeyError: '_adata_ref'` in the SAMap venv** — confirmed,
an anndata version incompatibility the pipeline comments predicted. In `rsc_uv` they fail with
`ModuleNotFoundError`, which is merely the wrong environment — **do not confuse these two**; I did.

## Other paths that are easy to lose
```
/n/lab_storage/.../human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv
```
The ortholog table lives under a **different project root**. `xmens_final/resources/orthology/` is
an empty placeholder. The table is **long-format** — see §6 trap 3.

---

# 3. THE ANALYSIS ARMS

## 3.1 Isodepth (spatial coordinate)
Per-sample GASTON fits → quantile-normalize → consensus. **Uninduced deliberately excluded from
seeding** (standing rule). Convention: **decidual-high in both species.**
- Geometric validation: ρ 0.797–0.945 on all 7 induced pieces (v3, QC'd coords); the old 0.447 was an artifact (ADDENDUM §B)
- Uninduced imputed-vs-measured: **ρ 0.651** — the imputation preserves rank order but compresses
  scale (imputed range only 0.094–0.569). This is why the imputed column is unusable for the control.
- A per-sample fit for uninduced **does exist** (`radius_uninduced.npz`, 4,155 cells). It was never
  seeded into the consensus. **Open decision:** reseed from it with an orientation gate (better) vs
  keep the geometric substitute (current).

## 3.2 Mutual information (the positional-information layer)
`mi_pseudotime.py`, worker/merge array over 200 chunks. Permutation null; **sklearn's estimator
clips at zero, so null-subtraction is an independence-bias correction only — significance comes
from the permutation.**
- **Batch stratification matters and was run:** MI computed within each sample, then n-weighted
  average = information about position *conditional on sample*. Agreement with pooled differs by
  species: **mouse ρ 0.927, human ρ 0.816** (0.727 on the matched universe). Do not quote one value
  for both — I did, and it was wrong.
- Stratification **strengthened** the cross-species result rather than merely preserving it.
- **Caveat to fix before the code becomes the paper's:** the stratified script re-normalizes
  already-normalized data. Harmless (MI is invariant under monotone transforms) but should go.
- **MI is not redundant with correlation:** rank correlation MI vs |Spearman| is only 0.433 (mouse)
  / 0.469 (human); top-decile overlap 46.8% / 56.0%; MI-unique genes are **~1.5–1.9× more curved**
  in their depth profiles (median |quadratic| 0.0578 vs 0.0303 mouse). Conserved sets from the two
  metrics overlap at **Jaccard 0.223**, with **132 conserved genes MI finds and correlation misses.**

## 3.3 Archetypes (profile-shape clustering)
Originally k=6, **hardcoded with no selection procedure** — admitted when asked. A proper sweep
(silhouette + 50-fold bootstrap ARI) gives:

| k | silhouette | bootstrap ARI | concordant | inverted | min cluster |
|---|---|---|---|---|---|
| 3 | 0.3998 | 0.987 | 2 | 1 | 84 |
| **4** | **0.433** | **0.983** | 2 | 2 | 42 |
| 5 | 0.3433 | 0.875 | 3 | 2 | 42 |

**k=4 is optimal on silhouette with stability essentially tied to k=3; ARI collapses at k=5.**
The four archetypes **are the natural sign quadrants** (up/up, down/down, up/down, down/up) — the
user spotted this, and it verified at high purity. Exceptions are genes whose correlation in one
species is ~0, where shape assignment is the defensible call.
**Critically: the concordant-vs-inverted split holds at every k tested**, so the biological
conclusion does not depend on the selection at all. Present the 2×2; use the archetype figure as
the visual.

## 3.4 The specificity control — a genuine negative result, and it redirected the project
**The question:** every conservation number was measured against a shuffled-gene null, which only
asks "more than chance?" It cannot distinguish decidualization-specific conservation from generic
stromal conservation.

**The user supplied the control from their own data:** the uninduced sample. Getting it usable took
most of an evening (§6 trap 5). Design: power-matched (induced subsampled to n=3,000 with 10
replicate draws), identical coordinate construction on both sides, shared gene universe, per-arm
permuted-coordinate check.

| arm | fold | threshold-free ρ |
|---|---|---|
| **A: uninduced** (3,958 cells, 4 blobs) | **2.581** | 0.420 |
| **B: induced, matched n** (10 reps) | **2.377** (sd 0.094, range 2.247–2.593) | 0.451 |

**Verdict: conservation is NOT decidualization-specific.** Uninduced matches or slightly exceeds
induced. I predicted a clear gap and was wrong.
Splitting by gene class (defined from the *human* side to avoid circularity) found **no class is
decidualization-specific**, though conservation is consistently weaker at the decidual pole in both
states (~2.06× vs ~3.3×).

**Report this negative result prominently. It is what makes the 37-gene claim credible rather than
assumed.** Limit to state: **n = 1 animal, 2 sections, 4 blobs** — pseudo-replicates, so a
qualitative contrast, not a CI across animals.

## 3.5 Robustness of the conserved-gene result
- **Cutoff dependence.** The top-decile threshold was inherited from the prior paper *before*
  seeing results, so it stands as pre-specified. Sweep: **1% → 11.65×, 2% → 8.13×, 5% → 3.98×,
  10% → 2.90×, 15% → 2.49×, 20% → 2.19×, 25% → 2.01×, 30% → 1.86×.** Enrichment declines
  **monotonically**. Report the whole curve; do not mine it.
- **RRHO is p-hacking if reported raw** (the user's challenge, and correct). The map's maximum is a
  max over thousands of correlated tests. With a permutation null over the whole map: observed max
  −log10 p = 300 (clipped), **null mean 32.7, null SD 6.59, null max 50.55, empirical p = 0.0050.**
  So unpermuted RRHO maxima are uninterpretable, and with the null the observed value is real.
  **The maximum's interior location (22.5–27.5%) reflects where statistical POWER peaks, not where
  agreement is strongest** — the cutoff sweep proves agreement is strongest at the strictest
  threshold. Use the **threshold-free rank correlation with a shuffled null (ρ 0.52, z 62)** as the
  primary statistic.
- **Sign concordance is uninformative.** Conserved 0.724 vs species-specific backgrounds 0.698 /
  0.712. **Profile shape separates them:** median profile r **0.674** (conserved) vs 0.495 / 0.496,
  p = 1.8e-4 / 9.9e-4. An earlier background was selected *alphabetically* and not matched on
  information content — that inflated the effect. A third planned control was ill-posed by
  construction (`matched_random: n=0`) and was **dropped, not patched**.

## 3.6 SAMap — the paper arm is VANILLA (user-confirmed: "yeah vanilla is all good")
`run_samap.py` calls stock `SAMAP(...)` + `sm.run(NUMITERS=3)`. No harmony. A replication battery
showed every qualitative conclusion holds on stock; harmony is a robustness supplement.

**Normalization (settled, and I was wrong once):** `.X` gets **CP10K, unlogged**; `.raw` is
**deleted**; guards assert row sums == 1e4 **and** `.X` max > 20. Rationale from the library source:
SAM applies `sum_norm='cell_median'` → `log2(x+1)` → `min_expression=1`. On CP10K input the
sum-norm step is a no-op, so exactly one log2 is applied, and CP10K keeps ~99.7% of nonzero entries
vs ~71% for raw counts at this depth. **Deleting `.raw` is load-bearing** — SAM prefers `.raw.X`
and would silently swap the scale.

**Gene-level conservation/divergence** (notebook 07, from the all-samples vanilla object):
- Refinement **only prunes; it adds zero edges**. Sequence graph = candidate universe.
- Of **754,031** cross-species candidate edges, **225,553 survive (29.9%)** and **528,478 are
  pruned (70.1%)**. Seven in ten sequence-similar cross-species pairs are not corroborated by
  expression.
- **Never divide reweighted/initial** — the graphs are on incommensurable scales (bit scores median
  241 vs correlations median 0.055). The ratio reads <0.5 for 100% of edges as a pure unit artifact.
- Per-gene verdicts (13,649 genes with a 1:1 ortholog + ≥1 surviving partner):
  **CONSERVED 8,425 (61.7%) · SUBSTITUTED 3,196 (23.4%) · ORTHOLOG_REJECTED 2,028 (14.9%)**
  → **~38% of assessable genes have their 1:1 ortholog displaced or rejected.**
- **41 substitutions named independently by both the all-samples graph and GenePairFinder** on the
  12 per-pairing runs: `ITGAV→Itga5` (ortholog rank 2, 12/12), `PGR→Ar` (rank 3), `HIF1A→Epas1`,
  `NFIX→Nfia`, `LSAMP→Negr1`, `MDK→Ptn`, `KDM5B→Kdm5c`, `MAGI1→Magi3`, `TXNRD1→Gsr`, `FAT4→Pcdh9`.
- Homology-model choice (all three data-independent, so none circular):
  **Ensembl strict 1:1 → fold 3.112** (9,550 links, 0% paralog) · SAMap BLAST → 1.855 (190,472
  links, 95.5% paralog, paralog-only stratum still **1.764**) · SAMap HMMER → 1.452 (228,212 links,
  97.8% paralog, paralog-only **1.408**). **Keep strict 1:1 as primary**; the paralog-tolerant
  graphs dilute rather than recover, but their paralog-only enrichment is the positive divergence
  result and belongs in a sensitivity analysis.

## 3.7 The compartment asymmetry (GenePairFinder, 13,895 pairs, 12 pairings)
| compartment | n | housekeeping | promiscuous | clean paralog % | median MI pct (h/m) |
|---|---|---|---|---|---|
| dec↔dec | 9,828 | **10.6%** | 29.1% | 49.8% | 0.919 / 0.904 |
| fib↔fib | 4,067 | **0.7%** | 22.4% | 46.5% | 0.932 / 0.931 |

**Decidual correspondence is ~15× more housekeeping-carried than fibroblast correspondence**, in
the same direction in every one of the 12 pairings. Fibroblast correspondence runs on identity genes.
- A **promiscuity filter is mandatory**: mouse genes attracting ≥4 distinct human partners are
  homology-graph artifacts (six unrelated human genes onto `Wdr33`; the UBE2 family onto `Ube2d2a`).
  Before filtering, the decidual driver list was topped by artifacts.
- **Do NOT claim per-run housekeeping tracks mapping specificity.** I did; my own table refuted it
  (the two highest-HK runs are the 1dpi_early pairings, not uninduced, and there is no monotone
  relationship). The pooled compartment contrast is the defensible statement.

## 3.8 Cross-method convergence — the strongest validation available
Two analyses sharing **no machinery** (SAMap = cross-species expression correlation; MI = mutual
information with the spatial coordinate, within sample) name the same genes:
- decidual pole: **RUNX1, CEBPB, GJA1, CD63, NAMPT**
- fibroblast pole: **MEIS2, EBF1, NFIA, RBMS3, RORA**

**NOT convergent — asserted once in error:** WNT5A. It is a **loss** gene in the MI analysis
(gain z = −1.91, `uninduced_top10=True`, `induced_top10=False`), and its SAMap pair is
`hs_WNT5A;mm_Wnt4` — a *different mouse gene* from the 1:1 `Wnt5a` the MI used. Notebook 04 asserts
this so it cannot recur.

## 3.9 Resolution — the premise was inverted all session
I repeated that "10 human decidual clusters collapse onto 1 mouse state." **Wrong.** In the
deliverables at `cell_subtype`: **human decidual = 1 subtype; mouse decidual = 3**
(`decidual_boundary`, `decidual_mature`, `decidual_proliferative`). The "10 human clusters" were
unannotated Leiden clusters computed inside a joint object — a different vocabulary.

Both sweeps then landed:
- **Mouse** decidual: 3 clusters at res 0.2 (sil 0.0842, bootstrap ARI 0.663), and both criteria
  decline **monotonically** thereafter → a continuum, and **the decidual-pole conservation gap is
  NOT an annotation artifact.** This retires a caveat I had called load-bearing.
- **Human** decidual (21,648 cells): 4 clusters at res 0.3 (sil 0.0651, ARI 0.417) — substructure
  exists but is weak, and the markers do **not** correspond to the mouse states (two look like
  stress signatures).

**Recommendation, reversed from earlier: do NOT do fine-resolution gene pairing.** Matching cluster
labels across species here would manufacture correspondence. The continuous spatial coordinate is
the honest refinement. (`cell_subtype` transfer agreement is a hard **0.000** for exactly this
reason — 2 human labels vs 9 mouse; not a bug.)

## 3.10 Cross-species transfer and the 4dpi failure
Per-pairing, both coordinates, both directions, 48 measurements, each against a shuffled-seed null:

| arm | n | mean ρ | range |
|---|---|---|---|
| induced (excl 4dpi) | 32 | **0.688** | 0.468 – 0.822 |
| uninduced | 8 | 0.492 | 0.214 – 0.699 |
| **4dpi** | 8 | **−0.048** | −0.427 – 0.354 |

**Isodepth ≈ pseudotime as a cross-species carrier** (induced: 0.695 vs 0.681). Since isodepth is
*measured* and pseudotime *inferred*, equal performance argues for **isodepth as the primary axis —
same transferability, fewer assumptions.**

**The 4dpi (shedding) failure — diagnosed precisely, and it is NOT a bug:**
- The join is perfect: ρ(radial position, native isodepth) = **−0.9799 in the export, −0.9799 in
  the source**, bit-identical.
- The mouse 4dpi coordinate is **excellent** — ρ **0.947** vs geometry, the *best* of any timepoint.
- Composition is real (79.1% decidual) and 3dpi is *higher* (81.3%) yet transfers fine.
- Cross-species edges are plentiful: **median 20** human neighbours per mouse cell, 0.05% with zero.
- **The fibroblast arm collapses.** Cross-species edge mass from mouse fibroblasts onto human
  fibroblasts: **95.8% at 1dpi_late → 26.1% at 4dpi** (going 73.9% to human *decidual* instead).
  Mouse decidual is unaffected (66.5% → 67.7%).
- This **inverts the compartment offset**: native gap (dec−fib) **+0.411 → transferred −0.074**.
  Within-compartment ρ is +0.263 / +0.277 but **pooled ρ is only +0.121** — a Simpson reversal.
  **So it is systematic error, not noise** (I called it "speckle" and "noise"; that was wrong).
- Independently confirmed at cluster level: `diag_over_off` **4.585** at secretory×1dpi_late →
  **0.567** at menstrual×4dpi (below 1 = off-diagonal exceeds diagonal).
- **Multi-sample context partially rescues it:** in the full joint object 4dpi transfers at
  **ρ 0.207** — degraded, not inverted, because a 4dpi cell can anchor through 3dpi/2dpi neighbours.
- **Conclusion:** exclude 4dpi from *single-pairing cross-species* analyses; it remains fully usable
  for *within-mouse spatial* work. And the label transfer there is **worse than a constant guess**
  (0.571 vs a 0.816 majority baseline).
- **Proposed QC metric (not yet a rule):** self-directed fibroblast edge mass. ~0.96 healthy,
  ~0.26 broken. One cheap number that would have flagged 4dpi immediately.

Best correspondence window: **1–2 dpi**, in both human references independently. Frame as a
recommendation ("use 1–2 dpi"), not a limitation.

---

# 4. STANDING RULES (user-established — follow without asking)

1. **Nothing launches without explicit approval.**
2. **Job width may go to a full single wave** — *"we can go up to 100 or 200 at once, it's totally
   fine."* The snakemake profile had `jobs: 20` hardcoded, which is why controllers kept dying on
   walltime with zero worker failures. Now 200. Backup: `config.yaml.pre_width.bak`.
3. **Back up (move, never delete) before destructive operations.**
4. **The mouse deliverable stays untouched by embedding fixes.**
5. **Apoptotic cells are dropped before constructing any new embedding.**
6. **Roots are chosen manually and visually**, never auto-picked. Wrong-orientation coordinates are
   refused or flipped **with recorded provenance** — never silently negated.
7. **Depth convention: decidual-high in both species.**
8. **Uninduced is excluded from consensus depth seeding.** Its anchor, when needed, is the user's
   choice: *"I wanted outer tissue boundary, depth from muscle inward."*
9. **New analyses become config-driven pipeline modules**, not one-off scripts.
10. **Readiness checks precede rule-writing.**
11. **Vanilla (stock) SAMap is the paper arm.**
12. **Spatial steps are no longer part of the cell-cell communication analysis** — trim that rule.
13. **On evidence quality:** *"I really don't think these numbers are always helpful? like you can
    get almost anything you want here. but at the end of the day it's about genes - what key
    regulators are conserved and which ones aren't?"* → prefer per-gene, checkable claims over
    aggregate scores. Alignment scores track dataset size; the user rejected them and was right.
14. *"You can use plain gene expression (no scenic) too."* SCENIC supplies TF *names* only.

---

# 5. THE PIPELINE

`xmens_final/workflow/` — snakemake. Rules `10_samap.smk` … `15_pt_per_timepoint.smk`.
Config keys that already collided once: use **`mi_pseudotime`** (not `mi`) and rule **`liana_ccc`**
(not `liana`). Config expanded to 12 SAMap pairings (`config.yaml.pre_12pairs.bak`).

Six targets validated by dry-run, **nothing launched**:
`pt_per_timepoint_all` (8, GPU) · `mi_all` (405) · `spatial_gp_all` (2) · `scenic_all` (151) ·
`liana_ccc_all` (2) · `samap_all` (resolves).
**Trap:** `scenic_all` wants **74 GRNBoost jobs at ~3 days each**. Review the seed plan before ever
running it — 50 seeds (1000–1049) already exist on disk.

**Not yet written: `16_conservation.smk`.** The ~20 one-off scripts from this session should become
config-driven rules (standing rule 9). This is the main reproducibility debt.

---

# 6. TRAPS — every one of these cost real time

1. **`glob.glob("~/**/x.csv", recursive=True)` walks the entire cluster home** (thousands of job
   dirs) and times out. It killed the first notebook smoke test at 1200 s. All notebook inputs are
   now staged in `notebooks/verify/data/`.
2. **`ast.parse` is not validation.** It confirms syntax and nothing more. `Q.median` resolves to
   the DataFrame *method*, not a column named `median`, and only *executing* the notebook caught it.
   Any column named `median`/`count`/`min`/`max`/`sum`/`shape` shadows a method under attribute
   access — always bracket-index.
3. **The Ensembl ortholog table is LONG-format** (`component_id | species | gene_symbol |
   homology_type | human_copy_number | mouse_copy_number | …`). A "column containing 'human'"
   heuristic silently picks **`human_copy_number`** and yields numeric garbage pairs **with no
   error**. Pivot on `component_id` where `homology_type == 'ortholog_one2one'` → 16,277 pairs.
   Assert the schema.
4. **Deliverable naming:** mouse is `XMens`, not `Mouse`.
5. **Uninduced depth has TWO layers** and conflating them wasted an hour: a per-sample fit *exists*
   (`radius_uninduced.npz`); it was *excluded from consensus seeding*. "Never got a fit" is wrong;
   "never seeded" is right.
6. **Tissue geometry:** each uninduced section contains **two separate tissue cross-sections**.
   Per-sample centroid logic fails (the centroid falls *between* the blobs) — an earlier radial
   outlier filter could not have worked. Use per-blob DBSCAN.
7. **Two different failures look alike:** `ModuleNotFoundError` on a SAMap pickle = wrong
   environment; `KeyError: '_adata_ref'` = the real anndata incompatibility.
8. **Never divide the two homology graphs** (§3.6) — incommensurable scales.
9. **A `.` in a kwarg name** is invalid Python and killed two jobs at submit.
10. **`json.dump` fails on tuple dict keys** — this killed a job *after* its object write succeeded.
    Always verify by reading the object back, never from the script's intent.
11. **Nested shell quoting** in `call_command` silently swallows loop variables; prefer separate
    calls or base64-encoded payloads.
12. **`c.call_command` has a 60 s ceiling** — use `submit_job` for anything heavier.

---

# 7. ROUTE TO THE FRIDAY MEETING

## What the meeting needs
A story with named genes, each backed by a number, and the caveats owned rather than buried. The
material exists; **what does not exist is the figure.**

## Priority order

**1. Reconcile the three transfer numbers — BLOCKING, do this first.**
`0.389` / `0.591` / `0.688` are all "cross-species isodepth transfer" from different constructions
(different arms, coordinate vintages, cell subsets, and directions). I built a factor-of-two
"discrepancy" on an inverted direction pair once. **Until reconciled, quote NONE of them.**
Also note `coordinate_transfer.json`'s `iso_to_mouse = 0.7921` and `dpt_to_human = 0.6547` are
**cross-scalar** (isodepth vs pseudotime), not like-for-like — they are not comparable to the above.

**2. Build the figure. Panel 3 is the keystone and does not exist.**
Recommended narrative: *conservation decays along the tissue axis.*
- **Panel 1** — the coordinate exists and is conserved: isodepth in both species on tissue
  coordinates, plus the geometric validation (ρ 0.80–0.94 on every induced piece, QC'd coords — ADDENDUM §B).
- **Panel 2** — the shared repertoire: cross-species MI enrichment (2.90×) with the threshold-free
  statistic and the cutoff sweep as an inset.
- **Panel 3 (KEYSTONE)** — conservation vs depth, **with the uninduced control overlaid.** This is
  the panel that carries both the positive result and the honest negative one.
- **Panel 4** — deployment differs: the 4-archetype heatmap (40 bins), genes × depth, both species
  side by side, blocked by archetype. Inversion reads as one panel being the other's mirror.
- **Panel 5** — the closing claim: the per-gene validity map, with FOXO1 called out as inverted.
Load `figure-composer` (which loads `figure-style` per panel); consider `paper-narrative` for
ordering across the whole set.

**3. Bootstrap CIs on the five headline numbers** — 2.90× enrichment, 2.581 vs 2.377 control,
2.327× set-B enrichment, 33.1% displaced, 61.7% ortholog-conserved. Currently all point estimates.

**4. `16_conservation.smk`** — convert this session's one-off scripts into pipeline rules.

**5. Add the edge-mass QC metric as a rule** (cheap, and it catches a real failure mode).

## Presenting it — what to lead with and what to own

**Lead:** the axis is real and conserved; the parts list is shared and TF-enriched; deployment
diverges for ~⅓ of positional genes; and 37 genes acquire conserved positional information
specifically upon decidualization.

**Own, do not bury:**
- Conservation is **not** decidualization-specific in aggregate (2.581 uninduced vs 2.377 induced).
  Say it plainly — it is what makes the 37-gene result credible.
- **FOXO1 is spatially inverted between species.**
- Uninduced control is **n = 1 animal, 2 sections, 4 blobs** (pseudo-replicates).
- Space and time are **confounded by design** in the mouse.
- 4dpi is excluded from cross-species work, and **why** (integration failure, not a bad coordinate).
- Alignment scores track dataset size and are not used for any conclusion.

**On the mouse model — the likely question.** Not "the mouse recapitulates human decidualization."
Instead: a model with a **defined validity map** — 9,550 genes labelled transferable / inverted /
species-specific, so someone designing a mouse experiment can look up their gene first. That is a
more useful and more defensible product than an unqualified conservation claim, and this is a
better paper than the one where everything came out conserved: it has a control that could have
killed the story and didn't.

---

# 8. STATE OF COMPUTE (as of 2026-09-01 morning — see ADDENDUM §F for the evening state)

**No live jobs.** All 19 ledger entries have ended. Every "failed"/"cancelled" row flagged on resume
is **superseded** — do not re-harvest:
- `fa535aee` — cancelled by me (bad ortholog pivot) → superseded by `8cf0096a` ✓
- `652a6220` — missing staged files → `53ab01b7` → `8cf0096a` ✓
- `5125b4cc` — staining wait-loop failed at shell level → rerun as `e1cd3505` ✓
- `6fd98648` — smoke test v1 (the recursive-glob timeout) → `05b77a1e` → `31d08a76` → `8962b10f` ✓
- `35d42e31` — **exit 1 came from `json.dump` AFTER the object write succeeded.** The fix landed;
  verified by reading the object back. Nothing to redo.
