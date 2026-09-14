# WORKFLOW INDEX — generated 2026-09-13 (regenerate after structural changes; do not hand-edit)

Attic: superseded rule files in `workflow/rules/attic/`, orphan scripts in `workflow/scripts/attic/` (see its MANIFEST.md).

| rule file | scripts | key outputs | purpose |
|---|---|---|---|
| 00_cellbender.smk | — | — | Stage 00 - CellBender ambient-RNA removal (UPSTREAM; not in the default DAG). |
| 00a_slidetags.smk | slidetags.py | — | Stage 00a — Slide-tags spatial coordinates (per-cell x,y) |
| 00b_ingest.smk | ingest.py | — | Stage 00b — ingest CellBender-corrected counts + Slide-tags coords -> AnnData |
| 01_qc.smk | qc.py | — | Stage 01 — quality control (Cevrim et al. 2025 Methods). |
| 01b_qc_report.smk | qc_report.py | results/qc/figures/ | Stage 01b - QC diagnostic figures (decoupled from the filtering rule 01_qc). |
| 01c_concat.smk | concat.py | — | Stage 01c - concatenate per-sample QC'd objects into a per-species object. |
| 02_integrate.smk | integrate.py | — | Stage 02 - scVI integration + Leiden clustering, per species. |
| 03_cluster.smk | cluster.py | — | Stage 03 — Leiden clustering + clustree + UMAP (per species). |
| 04_annotate.smk | apply_annotation.py | results/tables/clusters/ | 04_annotate -- apply the committed atlas choices CSV. |
| 04b_reintegrate.smk | reintegrate.py | — | Stage 04b - per-lineage scVI reintegration. |
| 04c_reannotate.smk | apply_annotation.py | — | 04c_reannotate -- apply a committed per-lineage choices CSV to a reintegrated |
| 04d_merge.smk | merge_reannotated.py | results/tables/refined_ | 04d_merge -- concatenate the cleaned per-lineage objects into one per-species |
| 04e_reembed.smk | integrate.py | — | 04e_reembed -- fresh whole-tissue scVI on the refined object. Removing EXCLUDE |
| 04f_spatial_qc.smk | spatial_qc.py | results/tables/spatial_qc_ | 04f_spatial_qc -- "blank, don't drop" spatial-coordinate QC, applied AFTER |
| 05_de.smk | de.py, de_lineage.py | results/tables/de/, results/tables/markers_ | Stage 05 — marker DE (Wilcoxon). |
| 06_pseudotime.smk | pseudotime.py | — | Stage 06 — diffusion pseudotime (global, mouse fibroblasts + decidual cells). |
| 06b_timepoint_pseudotime.smk | compute_pseudotime.py | — | Scratch-only, seeded per-timepoint mouse pseudotime tests. |
| 07_spatial_pseudotime.smk | spatial_pseudotime.py | — | Stage 07 — per-sample spatial pseudotime + GPR gradient field (slide-tags sample |
| 08_cellrank.smk | cellrank.py | — | Stage 08 — CellRank fate / probabilistic layer (on the global pseudotime). |
| 09_homology.smk | phmmer_to_samap_symbols.py | — | workflow/rules/09_homology.smk        (upstream of 10_samap.smk) |
| 10_samap.smk | make_transferred_object.py, prep_samap_input.py, run_samap.py, transfer_label_staining.py | results/figures/samap_transfer/, results/samap_transfer/, results/tables/samap_transfer/ | COMPLEMENTARY to 10b_samap_raw.smk: this file is the TRANSFER arm; 10b is the va |
| 10b_samap_raw.smk | compare_samap_arms.py, make_samap_transferred_h5ad.py, prep_samap_full_input.py, run_samap_arm.py, samap_attach_pt.py, samap_label_stain_full.py, samap_transfer_pt.py | — | COMPLEMENTARY to 10_samap.smk: this file is the VANILLA/RAW arm (the paper arm); |
| 11_mi_pseudotime.smk | mi_pseudotime.py | results/mi | workflow/rules/11_mi_pseudotime.smk |
| 12_spatial_gp.smk | spatial_pt_gp.py | results/spatial_gp/ | workflow/rules/12_spatial_gp.smk |
| 13_scenic.smk | build_consensus_aucell.py, ctx_only.py, grn_only.py, scenic_make_loom.py | — | workflow/rules/13_scenic.smk |
| 14_liana.smk | — | results/liana/ | workflow/rules/14_liana.smk |
| 15_pt_per_timepoint.smk | pt_per_timepoint.py | results/pt_per_timepoint/ | workflow/rules/15_pt_per_timepoint.smk |
| 16_liana_crossspecies.smk | liana_crossspecies.py | results/liana/no_commot, results/liana/no_commot/, results/liana_crossspecies | workflow/rules/16_liana_crossspecies.smk |
| 17_liana_pseudotime.smk | — | results/liana_pseudotime | workflow/rules/17_liana_pseudotime.smk |
| 22_crossspecies.smk | crossspecies.py | results/tables/conservation.csv, results/tables/mi_human.csv, results/tables/mi_mouse.csv | Stage 22 — cross-species conservation: orthogroups + threshold-free rank compari |
| 23_invitro.smk | invitro.py | results/tables/invitro_comparison.csv | Stage 23 — in-vitro decidualization vs in-vivo comparison. |
| 24_pseudotime_notebooks.smk | run_notebook.py | results/pt_notebook_rules/ | 14 -- pseudotime notebooks as DAG rules (dependency-free papermill-style runner) |
| 90_figures.smk | figures.py | results/figures/.done, results/tables/conservation.csv, results/tables/markers_ | Stage 90 — assemble figure panels from the tables above. |
| 98_clustree.smk | export_clustree.py, ladder_stability.py | results/figures/clustree/, results/tables/clustree/, results/tables/clustree/stability_summary_all_lineages.tsv | Clustree -- OFF the critical path, like 99_reports. |
| 99_reports.smk | apply_annotation.py | results/notebooks/ | 99_reports -- OFF the critical path on purpose. |
