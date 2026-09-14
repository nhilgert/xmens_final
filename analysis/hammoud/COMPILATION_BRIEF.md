# TASK BRIEF — compile ALL Hammoud/Ulrich atlas-mapping work into one verifiable notebook

**For:** Claude Code, working on `ssh harvard-eddy` (FASRC).
**Deliverable:** ONE notebook a reader can go through top to bottom and verify, covering the human
arm, the mouse arm, and the gene-driver decomposition. Nothing recomputed that already exists on
disk unless a verification cell disagrees with it.

---

## 0. THE ONE RULE THAT MATTERS

This project's cross-species work accumulated ~15 corrections, and **every one moved toward a
tidier story than the data supported.** The mechanism: a summary sentence written before the table
was re-read.

So this notebook is not a write-up. **It is a write-up with a self-check under every claim.**
Adopt the house pattern already used in `xmens_final/notebooks/verify/` — copy the helper verbatim:

```python
_checks = []
def verify(label, got, expected, tol=0.02, note=""):
    try:
        ok = (abs(float(got)-float(expected)) <= tol*max(1.0, abs(float(expected)))
              if expected is not None else None)
    except (TypeError, ValueError):
        ok = (got == expected)
    tag = "PASS " if ok else ("MISMATCH" if ok is False else "n/a  ")
    _checks.append((label, got, expected, tag))
    print(f"[{tag}] {label}: recomputed={got}  recorded={expected}" + (f"   ({note})" if note else ""))
    return ok

def check_summary():
    d = pd.DataFrame(_checks, columns=["claim","recomputed","recorded","status"])
    print(f"\n{len(d)} checks | {int((d.status.str.strip()=='MISMATCH').sum())} MISMATCH")
    return d
```

Every headline number in prose must have a `verify()` beneath it recomputing from the CSV/JSON.
End the notebook with `check_summary()`. **Target: 0 MISMATCH, 0 cell errors.**

---

## 1. WHERE EVERYTHING IS (verified inventory, 2026-09-01)

Working dir: `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud`

```
hammoud_fib_centroid_mapping.ipynb    29 cells, HUMAN ARM ONLY  <- the thing to replace
HANDOFF.md                            prior-session context (says mouse arm "not started" - STALE)
pnas.2404775121.sd05.xlsx             atlas centroids, sheet "5B Subtype centroids"
pnas.2404775121.sd06.xlsx             atlas DE markers,  sheet "6C Stro markers"

fib_mapping/
  hammoud_fib_mapping.py              CORE MODULE - reuse, do not reimplement
  run_mapping.py                      human pipeline
  run_mapping_mouse.py                MOUSE pipeline (2026-08-30 13:43)
  resolution_stability.py             mapping at every leiden resolution
  why_no_match.py                     marker detectability / why 3 subtypes cannot map
  make_notebook.py                    ** GENERATES the .ipynb - EDIT HERE, NOT THE IPYNB **
  submit*.sh                          sbatch wrappers

fib_centroid_results/                 HUMAN arm, 39 files
  TABLE_correspondence_summary.csv     <- the human deliverable
  similarity_{fib,stro}_contrastive_cosine.csv
  similarity_{fib,stro}_spearman.csv
  correspondence_{fib,stro}_by_{cluster,subtype}[_spearman].csv
  resolution_sweep.csv, resolution_stability.csv, cluster_sizes.csv
  replication_per_{sample,phase}.csv, replication_consistency.csv
  robustness_{marker_ablation,gene_bootstrap,heldout_split}.csv
  marker_sets.csv, headline_marker_fraction_positive.csv, score_by_cluster.csv
  per_cell_subtype_scores.csv (10 MB), cluster_centroids.parquet, cluster_assignments.csv
  program_detectability.csv, marker_panel{s}_detectab*.csv
  fig1..fig7 *.png

fib_centroid_results_mouse/           MOUSE arm, 7 files (2026-08-30 13:44)
  TABLE_correspondence_mouse.csv
  TABLE_correspondence_mouse_by_query.csv
  matrix_contrastive_cosine.csv
  matrix_raw_spearman.csv
  heldout_mouse.csv, gene_bootstrap_mouse.csv, meta_mouse.json

figures/  dotplot__fib_{published_markers,negative_controls}.png
```

**Core-module API — read it, do not guess (a prior attempt failed by guessing names):**
`FIB_SUBTYPES`, `SHORT`, `load_atlas_centroids(path, sheet="5B Subtype centroids")`,
`load_atlas_markers(path, sheet="6C Stro markers")`, `linear_panel_matrix(...)`,
`group_centroids(...)`, `contrastive_cosine(...)`.
**Open the file and check each signature and return type before calling it.**

### The gene-driver file — NOW EXPORTED (2026-09-01 12:54)

It previously existed only in the Claude Science artifact store, which Claude Code cannot reach.
It is now on the cluster at
**`hammoud/fib_centroid_results_mouse/hammoud_gene_drivers_v5.json`** (13,010 bytes).

Structure — 19 entries keyed `"<species>:<query_cluster>-><Stro subtype>"`, each with:
`cosine`, `n_cells`, `n_panel`, `frac_from_top30`, `top15`.
Covers **mouse**: `decidual_{apoptotic,boundary,mature,proliferative}`,
`efib_{middle,middle_induced,outer,proliferative,subluminal}` and **human**: `h0`–`h9`.

**Do not regenerate it** — that computation failed three times before it worked (see §3). Read the
exported file. Note `hammoud_fib_mapping.py` and `why_no_match.py` are on the cluster but NOT in the
artifact store, so the cluster copies are the only copies — do not overwrite them without a backup.

---

## 2. NOTEBOOK STRUCTURE TO BUILD

Generate via `make_notebook.py` (extend it; add a `--arm {human,mouse,combined}` switch or a new
`make_notebook_combined.py` alongside it). Write to
**`hammoud/hammoud_atlas_mapping_ALL.ipynb`** — a NEW filename. Do not overwrite the existing
human-arm notebook; it is a working deliverable.

```
0.  What this is, and what question each section answers
    - the two reviewer comments this addresses (see HANDOFF.md)
    - the atlas: Ulrich 2024, PNAS 121:e2404775121, full-thickness uterus, single-CELL,
      20 donors, 39 subtypes
    - our data: single-NUCLEUS + Slide-tags spatial. State this asymmetry up front; it is
      the reason contrastive (double-centred) cosine is the headline metric rather than raw
      correlation.
    - a table of every input file consumed, with its mtime, printed by the notebook itself

1.  METHOD — why contrastive cosine
    Reproduce the existing justification, and add a verify() that the raw-Spearman and
    contrastive-cosine matrices give DIFFERENT top matches for at least one cluster
    (otherwise the choice of metric is immaterial and should be stated as such).

2.  HUMAN ARM  (from fib_centroid_results/)
    2.1 resolution choice            resolution_sweep.csv, resolution_stability.csv
    2.2 correspondence matrix        similarity_fib_contrastive_cosine.csv + fig1
    2.3 the correspondence table     TABLE_correspondence_summary.csv
    2.4 sink test                    is each match real, or argmax with nowhere else to go?
    2.5 replication                  per sample, per phase, consistency
    2.6 robustness                   marker ablation, gene bootstrap, held-out split
    2.7 marker-level evidence        marker_sets.csv, headline_marker_fraction_positive.csv
    2.8 why three subtypes cannot map  program_detectability.csv, why_no_match.py output
    2.9 per-cell + spatial readout   per_cell_subtype_scores.csv, fig5

3.  MOUSE ARM  (from fib_centroid_results_mouse/)  ** NEW - not in any current notebook **
    Same skeleton as §2 for whatever exists: matrix_contrastive_cosine.csv,
    matrix_raw_spearman.csv, TABLE_correspondence_mouse{,_by_query}.csv,
    heldout_mouse.csv, gene_bootstrap_mouse.csv, meta_mouse.json.
    Be explicit about which §2 analyses have NO mouse counterpart (no replication file, no
    marker ablation, no spatial readout) and say "not run" rather than leaving a silent gap.

4.  GENE DRIVERS — "what genes drive the correspondence"  ** NEW **
    From hammoud_gene_drivers_v5.json. For each of the 19 correspondences report
    cosine, n_cells, n_panel, frac_from_top30, top15.
    Then the three questions worth asking of it:
    (a) Is the correspondence carried by the atlas's own top-30 markers, or by the wider
        panel? -> frac_from_top30, per correspondence.
    (b) Which drivers are shared between the human and mouse routes to the SAME Stro subtype?
        e.g. human h4/h5/h8 -> Stro-7 C7 vs mouse efib_outer -> Stro-7 C7. Genes appearing on
        BOTH species' top15 for the same subtype are three-way conserved (atlas <-> human <->
        mouse) and are the strongest single claim in this whole body of work.
    (c) Which correspondences are driven by cell-cycle / ribosomal / chaperone genes rather
        than identity genes? Flag them explicitly. Expect the "proliferative" matches
        (mouse decidual_proliferative -> Stro-11 progenitor, efib_proliferative -> Stro-11)
        to be mitotic-signature-driven; if so, SAY SO - it changes what the match means.

5.  CROSS-ARM SYNTHESIS
    Which Stro subtypes are recovered in both species, which in one, which in neither.
    A single table: Stro subtype x {human match, human cosine, mouse match, mouse cosine,
    shared drivers, verdict}. This is the figure-ready object.

6.  WHAT TO WRITE IN THE REBUTTAL
    Keep and update the existing section. Separate (i) supported by both arms,
    (ii) human only, (iii) mouse only, (iv) not supported.

7.  LIMITS - write these, do not soften them
    - single-cell atlas vs single-nucleus query: compositional and dropout differences
    - n = 2 human donors (this is literally one of the reviewer comments)
    - three Stro subtypes are undetectable in our data by construction, not by biology
    - cosine on centroids is a similarity, not a statistical test; the sink test and the
      bootstrap are the evidence that a match is not an artifact
    - any proliferative match driven by cell-cycle genes is a state match, not an identity match

8.  check_summary()
```

---

## 2b. AUDIT FINDINGS — read before writing section 3 or 4

I audited the mouse files and the driver JSON on 2026-09-01. Four things will mislead you if you
present them at face value.

### (a) THREE different mouse cluster vocabularies are in play. Build a crosswalk FIRST.

| file | mouse cluster names |
|---|---|
| `TABLE_correspondence_mouse_by_query.csv`, `gene_bootstrap_mouse.csv` | short codes: `dec1 dec2 dec3 decAp eFib1 eFib2i eFib3 …` |
| `hammoud_gene_drivers_v5.json` | full names: `decidual_{apoptotic,boundary,mature,proliferative}`, `efib_{middle,middle_induced,outer,proliferative,subluminal}` |
| the deliverable h5ad `cell_subtype` | `decidual_boundary`, `decidual_mature`, `decidual_proliferative`, … |

**This exact class of mismatch already caused a wrong claim earlier in this project** (old cluster
codes read instead of renamed subtypes). Emit an explicit crosswalk table as a notebook cell, with
a `verify()` that every name in every file maps to exactly one canonical label, before joining
anything across files.

### (b) The atlas scope is INCONSISTENT between mouse files. This is the serious one.

- `TABLE_correspondence_mouse.csv` has **39 rows — the whole atlas, all lineages.** Its `best_query`
  column therefore matches endothelial and glandular subtypes (`L_Endo`, `Endo-1_Artery`,
  `Endo-2_EndoMT`, `Endo-3_PCV`) to our fibroblast/decidual clusters, which is not interpretable.
- `TABLE_correspondence_mouse_by_query.csv` picks `best_atlas` over **all 39**, so e.g.
  `dec2 -> 8_OLFM4_Prolif` and `dec3 -> 5_Gland_sec-progenitor A` — **glandular/epithelial**
  subtypes, for a stromal query.
- `gene_bootstrap_mouse.csv` and the driver JSON are restricted to the **6 mappable Stro subtypes**.

**Do not present the all-39 mouse tables as the mouse result.** Either restrict them to the Stro
subtypes and recompute the argmax (cheap, no refit — just re-argmax the saved similarity matrix),
or report them only as a specificity check ("the stromal query does not preferentially match
stromal atlas subtypes when the full atlas is offered" — which is itself a finding worth stating).
Ask the user which. The human arm does not have this problem.

### (c) Cosine magnitudes are TINY. They are not correlations and must not be read as such.

Full range across all 19 correspondences: **−0.0074 to 0.3671.**

| correspondence | cosine | note |
|---|---|---|
| `human:h5 -> Stro-7_C7` | **0.3671** | strongest overall |
| `mouse:efib_proliferative -> Stro-11 progenitor` | 0.3466 | |
| `mouse:decidual_mature -> Stro-8_dS1` | 0.2553 | strongest credible decidual match |
| `mouse:efib_outer -> Stro-7_C7` | 0.2429 | the three-way C7 anchor |
| `mouse:decidual_apoptotic -> Stro-8_dS1` | **0.0192** | ~zero; see (d) |
| `human:h3 -> Stro-10_dS3` | **−0.0074** | **NEGATIVE — an argmax over nothing** |

A double-centred cosine of 0.02 is not evidence of correspondence. **Every match needs a null**
(the sink test and gene bootstrap are the existing machinery — use them, and state a floor below
which a match is reported as "no match" rather than as a weak match). `human:h3` and
`mouse:decidual_apoptotic` should almost certainly be reported as **non-matches**.

### (d) Two matches are driven by the wrong kind of gene — my earlier prediction, confirmed.

- **`mouse:decidual_apoptotic -> Stro-8_dS1`: 7 of its top 15 drivers are ribosomal / translation /
  respiratory** (`RPS19`, `SEC61G`, `SERP1`, `UQCR11`, `ROMO1`, `EIF4G2`, `TMSB4X`). Combined with
  cosine 0.0192, this is a housekeeping artifact, not a correspondence.
- **`mouse:decidual_proliferative -> Stro-11_Stromal progenitor`: drivers are mitotic**
  (`ZWINT`, `GINS2`, `TFDP1`, `NCAPH2`, `DIAPH3`, `G2E3`). This is a **cell-cycle state match, not
  an identity match.** Say so explicitly — it changes the biological claim entirely, and this
  project has already established that the mouse "proliferative" decidual population is barely
  cycling at the median, so its name is historical.
- By contrast `efib_outer -> Stro-7_C7` (0.2429) is driven by identity genes — `PTGIS`, `FAP`,
  `CCN5`, `EPHX1`, `FRZB` — and is the credible three-way anchor.

**Add a driver-class column to every correspondence: `identity` / `housekeeping` / `cell_cycle`.**
Classify by symbol prefix (`RPS|RPL|MRP|EEF|EIF|TMSB|UQCR|ROMO|SEC61|SERP|PSM|NDUF|COX|ATP5`) and a
mitotic list, and report the count out of 15. A match whose drivers are not identity genes must not
be presented as evidence of conserved cell identity.

### (e) Two smaller gotchas

- **`frac_from_top30` is NOT a fraction.** It exceeds 1 for `decidual_apoptotic` (1.519) and is
  **negative** for `human:h3` (−3.856). It is a ratio and inherits the sign of a near-zero cosine.
  Do not label it as a percentage; if you present it, define it from the code in
  `hammoud_fib_mapping.py` and exclude the cases where the denominator is ~0.
- **Driver gene symbols are UPPERCASE (human space) even for the mouse arm** — the mouse panel was
  mapped through orthologs into human symbols. That is correct, but say it, or a reader will think
  mouse genes are being reported.

---

## 3. TRAPS — each of these already cost real time. Read before writing code.

**The gene-driver computation failed THREE times before it worked. The causes:**

1. **Guessed function names** that were not in `hammoud_fib_mapping.py`.
   -> Open the module and read the actual `def` lines first.
2. **Mis-read the module's return signature** (assumed a different shape than it returns).
   -> `print(type(x))` and `print(x.shape / x.keys())` on the first call before using the result.
3. **Silently collapsed the mouse gene panel** by intersecting HUMAN gene symbols against MOUSE
   `var_names`. Human symbols are uppercase (`RUNX1`), mouse are title-case (`Runx1`); a naive
   intersection returns a tiny accidental overlap and **produces a plausible-looking wrong answer
   with no error.** -> Map species symbols explicitly through the ortholog table, and
   `assert len(panel) > <expected floor>` on every panel you build.
4. **Assumed ortholog-table column names.** See the next trap.

**The Ensembl ortholog table is LONG-format and lives under a different project root:**
```
/n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv
```
(`xmens_final/resources/orthology/` is an empty placeholder.)
Columns: `component_id | species | ensembl_gene_id | gene_symbol | biotype | homology_type |
gene_tree_id | human_copy_number | mouse_copy_number | source_release`.

**It has `human_copy_number` and `mouse_copy_number` columns.** Any "pick the column containing
'human'" heuristic selects a COPY NUMBER and yields numeric garbage pairs **without raising**.
Correct usage:
```python
o2o = orth[orth.homology_type == "ortholog_one2one"]
w = (o2o.pivot_table(index="component_id", columns="species",
                     values="gene_symbol", aggfunc="first")
        .dropna(subset=["human","mouse"]))
pairs = dict(zip(w["human"].str.upper(), w["mouse"]))    # -> 16,277 pairs
assert pairs.get("TSPAN6") == "Tspan6"
```

**Other traps:**
- **The mouse deliverable is named `XMens`, not `Mouse`:**
  `2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad`. A species-name glob returns empty and has
  already killed one job.
- **Never `glob.glob("~/**/x", recursive=True)`** — it walks thousands of cluster job dirs and
  times out. Stage inputs into one directory and read them by explicit path.
- **`ast.parse` is not validation.** It catches syntax and nothing else. A column named
  `median`/`count`/`min`/`max`/`sum`/`shape` shadows the DataFrame method under attribute access
  (`Q.median` returns the method) — always bracket-index. Only *executing* the notebook catches it.
- **`json.dump` fails on tuple dict keys** — and it fails AFTER your writes, so a job can report
  failure with its real output already on disk. Verify by reading the file back.
- **`netscratch` is scratch and can be purged.** Prefer the `xmens_final/data/deliverables/` copies.
- Environment: `/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python`
  (scanpy 1.11.5, anndata 0.12.11, nbconvert present).
- **SLURM `-p eddy`, `--account=eddy_lab`.** Everything here runs in ~3 min on 8 cores / 96 GB.
  Do not run heavy jobs on the login node. Interactive ssh commands have a ~60 s ceiling.

---

## 4. WHAT NOT TO DO

- **Do not refit or recompute anything in `fib_centroid_results/` or
  `fib_centroid_results_mouse/`.** Those are finished. Read them, verify them, present them.
  Recompute only if a `verify()` disagrees — and then report the disagreement rather than
  quietly overwriting the file.
- **Do not edit `hammoud_fib_centroid_mapping.ipynb` in place** and do not edit the generated
  notebook by hand. The generator is `make_notebook.py`.
- **Do not run the mouse arm's missing analyses** (replication, marker ablation, spatial) without
  asking. They are absent because they were never run, not because they failed.
- **Do not trust `HANDOFF.md`'s "mouse arm is not started"** — it is stale; the mouse arm ran
  2026-08-30 13:44.
- **Do not launch any long job without explicit approval** — standing rule in this project.
- **Do not report a "conserved" or "three-way" claim without checking the driver genes.** A match
  carried by ribosomal or cell-cycle genes is not evidence of shared identity.

---

## 5. DEFINITION OF DONE

1. `hammoud/hammoud_atlas_mapping_ALL.ipynb` exists, generated from a script under `fib_mapping/`.
2. Executed end to end in `rsc_uv`, with an `*_executed.ipynb` copy saved beside it.
3. `check_summary()` reports **0 MISMATCH and 0 cell errors**.
4. Every headline number in prose has a `verify()` recomputing it from source.
5. Sections 3 (mouse) and 4 (gene drivers) are new content, not present in any existing notebook.
6. Section 5 produces the single cross-arm synthesis table.
7. A short `HAMMOUD_HANDOFF.md` replaces the stale `HANDOFF.md`, stating what is done, what was
   never run, and where every file is.
8. Report back: the count of checks, any MISMATCH with your reading of it, and an explicit list of
   analyses you found MISSING for the mouse arm.

**If a `verify()` fails, that is a result, not a bug to hide. Report it.**
