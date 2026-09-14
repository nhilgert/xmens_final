# HANDOFF — for a future Claude session on the Cevrim/XMens project
Written 2026-08-31 at the close of a multi-day session. Read this before
touching anything. The user is Nick (grad student, computational lead);
paper is in revisions; Friday meetings matter.

## Current objects — WHAT IS TRUE ON DISK
Canonical deliverables (spatial-outlier-QC'd, dpt_pseudotime attached,
embeddings corrected to per-lineage stromal fits):
  xmens_final/data/deliverables/2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad   (mouse, 43,197)
  xmens_final/data/deliverables/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad   (human, 47,533)
Identical copies: /n/netscratch/eddy_lab/Everyone/nhilgert/ and
/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/. Every modified
file has a sibling .bak (.pre_spatial_fix.bak, .pre_embedding_fix.bak, etc.).
uns['spatial_outlier_fix'] documents the coordinate QC (radial 3-MAD, coords
NaN'd, cells kept; 100 mouse cells, 0 human).

Isodepth objects (consensus spatial coordinate, decidual=high convention):
  xmens_final/results/isodepth/{human,mouse}_fibdec_isodepth.h5ad
  human: GASTON npz seeds (99.3/98.7% barcode join) + kNN label-spreading.
  mouse: radius_{tp}.npz seeds, uninduced UNSEEDED (imputed only — flag it).

SAMap exports (joint objects with transfers in obs):
  2026_08_31_samap_vanilla_fibdec_joint.h5ad  <- USE THIS ONE (vanilla = paper arm)
  2026_08_30_samap_harmony_fibdec_joint.h5ad
  both at netscratch + mckinley backup/data. obs: label_transfer_fine/coarse,
  palantir_pt_transferred, isodepth_transferred, QC'd x/y_um_dbscan.

Pseudotime: mouse rooted object (pseudotime_mouse_v2_induced_no_apop/
2026_08_28_mouse_fibdec_pseudotime_rooted.h5ad) carries dpt + palantir_pt
(rho 0.917 to each other). BEST mouse pt = results/pt_variants/varB_no_efibouter
palantir_pt (rank rho +0.852). Human integrated pt is weak (0.31) — use isodepth.

## Vocabulary
Old codes = new names, SAME partition (verified 1:1): dec1=decidual_boundary,
dec2=decidual_proliferative, dec3=decidual_mature, decAp=decidual_apoptotic,
eFib1=efib_outer, eFib2=efib_middle, eFib2i=efib_middle_induced,
eFib3=efib_subluminal, eFibP=efib_proliferative. fine_refined carries old codes;
cell_subtype carries new. USE cell_subtype. Human has only coarse cell_type
(endometrial_fibroblast/decidual_cell) + the fine leiden clusters in
results/human_fib_annotation/fine_cluster_assignments.csv (not yet named — a
naming pass with the user is pending).

## Standing rules (user-established; do not violate)
1. NOTHING launches without explicit user approval; back up (move, not delete)
   before any destructive operation; temp-then-rename atomic writes.
2. Wrong-orientation pseudotime is REFUSED at write time (flip the root, never
   negate the column). Roots are chosen manually/visually by the user.
3. decidual_apoptotic dropped before constructing any new embedding.
4. SAMap input scale = CP10K in .X (guard in prep_samap_input.py).
5. Isodepth convention: decidual=high, both species.
6. Uninduced mouse excluded from isodepth seeding and from induced-only embeddings.
7. Rank MI results by MI value (q ties at the permutation floor).
8. Vanilla SAMap = primary arm (replication battery passed 2026-08-31);
   harmony/scvi are robustness supplements. Harmony code audit not needed.
9. New analyses become config-driven snakemake modules, not one-offs.

## Pipeline state (xmens_final)
Targets that dry-run clean: pt_per_timepoint_all, mi_all, spatial_gp_all,
scenic_all, liana_ccc_all, samap_all. Bare `snakemake` (rule all) FAILS on
pre-existing stage-07 missing input — known, documented in pipeline_audit.md.
Config namespaces: mi_pseudotime (NOT mi — that key was taken), rule liana_ccc
(NOT liana — collision). SCENIC seeds on disk: 50 (1000-1049). WARNING:
scenic_all wants 74 GRNBoost jobs (~3 days each) — review before running.
MI runs: mouse_dpt, mouse_isodepth merged; human_dpt, human_isodepth merging
at handoff time. MI worker subsamples 3k cells, shares permutations per chunk.

## Traps this session hit (don't repeat)
- config.yaml edits: NEVER heredoc-pipe python that opens the config for write
  (a parse error after open('w') zeroed it once — restored from
  config.yaml.pre_new_arms.bak). Pattern: parse, write temp, re-parse, replace.
- snakemake lock: a walltimed controller leaves the dir locked — `--unlock` first.
- SCP staging (c.download / submit inputs) intermittently stalls; workaround:
  base64 via call_command stdin, or job-output harvest.
- The scvi latent table (heca integrations) prefixes obs_names 'source::' —
  split on '::'. The SAMap joint objects use bare barcode-suffix names, NO
  species prefix (distinguish by 'secretory/menstrual' substring).
- Ortholog table (human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv)
  is LONG format: pivot on component_id, filter homology_type=ortholog_one2one.
  Never intersect human symbols directly against mouse var_names (F8/C3/T trap).
- mutual_info_classif null scales with class count — match address counts across
  species or subtract null mean; 100 perms floors q at ~0.07 for 4k tests (use 300).
- Old atlas mapping module lives at xmens_final/hammoud/fib_mapping/; its
  group_centroids returns (C genes-x-groups, ns) and needs >=2 columns per side.

## Analyses complete (artifacts in project store + cluster)
Embedding fixes; SAMap CP10K rerun (2 pairings); label transfer; mouse+human
pseudotime (integrated, per-timepoint arm installed, variants A/B/C); root
sweep; isodepth consensus both species (+flip to decidual-high); coordinate
transfer (rho 0.74-0.79); isodepth-pt correlation (within-tissue 0.55-0.84);
MI: mouse_dpt, mouse_isodepth, positional both species (+batch confound test:
within-sample retention 78/84%); cNMF k=9/10/11 both species + cross-species
matching (+HOG-level negative); coexpression conservation (AUROC, TF-restricted);
Hammoud both arms + gene drivers (PTGIS/RERG/SCN7A = C7 axis; mitotic
Stro-11 match); HECA bridge (both sides); fine clustering (10/19/30 clusters);
decidual_proliferative subclustering (NO cycling subcluster — name is
historical); depth heatmap; fine cluster heatmap (staircase); staining panels;
vanilla replication battery. Reports: MORNING_REPORT.md, HECA_DE_audit.md,
pipeline_audit.md, isodepth_design.md, SESSION_SUMMARY.md, FRIDAY_MEETING.md,
ISODEPTH_TRANSFER_README.md (cluster, for the other Claude session).

## Open / pending
1. human_dpt + human_isodepth MI merges (controller c9f3d022) -> then
   cross-species top-10% isodepth-MI overlap + RRHO (rrho_mi.R staged,
   smoke-tested, uncommitted in workflow/scripts/).
2. Human fine-cluster vocabulary naming (needs user).
3. Pipeline-ify the session's one-off analyses (isodepth consensus, transfer,
   staining, positional MI, coexpression, cNMF compare).
4. Cluster-removal evidence audit (needs user's removal history) + spatial QC
   evidence figure (fix applied; figure pending).
5. LIANA: user says NO spatial steps — trim rule to rank_aggregate-only path.
6. SCENIC re-run on new objects (TF orthogroup crosswalk ready).
7. h16 -> dStromal_early HECA vote unexplained.
8. Claude Code end-to-end repo audit (user wants it; after vocabulary freeze).
9. Other Claude session works on SAMap/harmony (code/samap-scvi/) and Hammoud
   (xmens_final/hammoud/) — coordinate via ISODEPTH_TRANSFER_README.md.
