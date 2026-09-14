# QUEUE STATUS — 2026-08-31 evening (updated after 4 jobs landed)

Collect with: `host.compute.create('ssh:harvard-eddy').attach_job('<id>').result()`

## STILL IN FLIGHT (2)

| job_id | what | results land |
|---|---|---|
| fdb86d01-a790-4ff9-95e4-93ede06119bd | 12 CP10K **vanilla** SAMap integrations (hmmer): 8 prep -> 12 run_samap -> 12 label_stain | `runs/<run>/samap.pkl`, `results/figures/samap_transfer/*_labelstain.png` |
| 2fc9c3b2-d56a-4292-8efe-e896a79422c7 | GenePairFinder on those 12 — self-polling (5min/40h) for pickles >1GB AND mtime>1788206400; writes JSON incrementally | job workdir `genepairs_all_runs.json` |

## HEADLINE RESULT (revised — this is the paper's claim)

Baseline stroma already carries the conserved positional program at full strength. **On top of it,
37 genes ACQUIRE conserved positional information specifically upon decidualization** (uninduced MI
at noise floor, induced top-decile, human top-decile, gain z>2): Cped1(z=12.0), Runx1(10.2),
Cd63(7.5), Ebf1(7.4), Cebpb(7.3), C1qtnf7(8.3), Txn1, Jarid2, Rbms3, Ldlrad4, Fosl2, Smoc2, Pim1,
Bnip3l, Itgb5, Rev3l, Specc1, Pard3b, Asic2, Afdn. Wider set (n=115, z>2) adds Sat1, Ror2, Nampt,
Pfkfb4, Gja1, Ctnnb1, Adamts5, Aldh1a2.
**And 31 conserved genes LOSE positional information on induction** (conserved erasure): Kcnd2,
Tmem132c, Lsamp, Palld, Synpo2, Vcan, Ptch1, Dcn, Dio2, Pdgfd, Mmp16, Prdm1.
Set B (conserved & not-baseline-positional, n=147) is **2.33x enriched, p=1.5e-24**.
=> Decidualization does not create the axis; it reorganizes which genes read it.

## OTHER NUMBERS (final)

- Specificity control (aggregate): uninduced 2.581x vs induced 2.377+/-0.094 (matched n=3958, 10 draws).
  NOT decidualization-specific in aggregate — hence the gene-level analysis above.
- Archetype split: fibroblast-pole 3.29 (unind) / 2.75 (ind); decidual-pole 2.06 / 2.03. Conservation
  ~1.6x weaker at the decidual pole in BOTH states.
- Homology: keep strict Ensembl 1:1 primary (3.11x). BLAST 1.86x, HMMER 1.45x — 95-98% paralog edges,
  dilution not recovery. Paralog-only stratum still enriched (1.76x/1.41x) -> sensitivity analysis only.
- TFs 2.4x enriched among conserved-positional (5.2% vs 2.4%, Fisher p=8.7e-5). 9550 genes:
  226 CONSERVED_POSITIONAL, 95 INVERTED, 746 HUMAN_ONLY, 720 MOUSE_ONLY.
- **Delta-peak validated on TWO axes**: isodepth median 0 / mean +1.06 / SD 4.36 / 41.6% within +/-1 /
  33.1% at >=5; pseudotime median 1 / mean +1.35 / 38.9% / 31.0%. Spearman(iso,pt)=0.732, same sign
  69.7%, both-matched 29.1%, both-displaced 22.6%. Tails asymmetric toward human-later (+9:36, -9:17),
  Wilcoxon p=3.6e-7 and 1.7e-13. Reading: predominantly REWIRING + minor human-later heterochrony.
- Boundary-depth anchor validated: distance-from-myometrium reproduces fitted isodepth at median
  rho 0.893 across induced blobs. Stored uninduced isodepth is label-spread (rho 0.651 vs measured) —
  DO NOT USE for the control.
- Regulators: only RUNX1 fully consistent (6/6 human, 7/7 mouse sections). CEBPB/RUNX2/LEF1 strong but
  partial. WT1 weak (pct_human 0.920, consistency 2/6) — do not cite as transferable. FOXO1, MITF INVERTED.

## CORRECTIONS TO EARLIER CLAIMS (important — do not re-inherit)

1. **The resolution asymmetry runs the OTHER way.** At `cell_subtype`: human decidual = 1 subtype
   (`decidual_cell`); mouse decidual = 3 (`decidual_boundary`, `_mature`, `_proliferative`). The
   "10 human clusters onto 1 mouse state" line came from SAMap's unannotated fine leiden clusters,
   NOT cell_subtype — never a matched comparison. Open question is now on the HUMAN side.
2. **Mouse decidual is NOT under-resolved.** Leiden sweep: silhouette and bootstrap ARI both decline
   monotonically from res 0.2 (3 clusters, sil 0.084, ARI 0.663) to res 3.0 (48 clusters, sil 0.018,
   ARI 0.328). No stability optimum at finer resolution. Low silhouettes throughout => the compartment
   is a CONTINUUM, not discrete states. The decidual-pole conservation gap is therefore NOT an
   annotation artifact.
3. Delta-peak tails are asymmetric (+9:36 vs -9:17), not symmetric as first stated.
4. `prep_samap_input.py` writes **CP10K** into `.X` (raw-counts version at `.raw_counts.bak`).

## TRAPS

- Mouse deliverable is `2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad` (NOT "Mouse_") — a Mouse_* glob
  returns empty and killed one job.
- April SAMap pickles (`runs/*/samap.pkl`, Apr 24 mtime) are UNREADABLE: anndata `KeyError: '_adata_ref'`.
- SAMap alignment scores track dataset size (Spearman(n_joint, score) ~ -0.9) — do not read timecourse
  trends off them. Menstrual reference oscillates; secretory is near-monotone. Not evidence.
- Backups: `config/config.yaml.pre_12pairs.bak` (pairs 2->12), `profiles/slurm_noconda/config.yaml.pre_width.bak`
  (job cap now 200).

## NEXT SESSION, IN ORDER

1. Collect fdb86d01 + 2fc9c3b2 (the dec<->dec gene pairs across all 12 pairings incl. uninduced).
2. Subcluster the HUMAN decidual compartment (mirror of the mouse sweep) — the now-correct form of the
   resolution question.
3. Build the cross-species figure. Panel 3 (conservation vs depth with the uninduced control) is the
   keystone and does not exist yet. The 37-gene gain set + 31-gene loss set is the new panel 4.
4. Bootstrap CIs on headline numbers (2.4x TF enrichment, 37/31 gene sets, 42%/33% Delta-peak).
5. Capture the MI/conservation scripts as `16_conservation.smk` (standing rule: config-driven modules).
