# WORKFLOW INDEX (generated 2026-09-14, post-reorganization; do not hand-edit)

Layers: workflow/config/resources/external = machine | data/results = artifacts (gitignored) | analysis = living human work | archive = fossil record.

| rule file | scripts | key outputs | purpose |
|---|---|---|---|
| 00_cellbender.smk | — | — | Stage 00 - CellBender ambient-RNA removal (UPSTREAM; not in the default DAG). |
| 00a_slidetags.smk | ingest_qc/slidetags.py | — | Stage 00a — Slide-tags spatial coordinates (per-cell x,y) |
| 00b_ingest.smk | ingest_qc/ingest.py | — | Stage 00b — ingest CellBender-corrected counts + Slide-tags coords -> AnnData |
| 01_qc.smk | ingest_qc/qc.py | — | Stage 01 — quality control (Cevrim et al. 2025 Methods). |
| 01b_qc_report.smk | ingest_qc/qc_report.py | results/qc/figures/ | Stage 01b - QC diagnostic figures (decoupled from the filtering rule 01_qc). |
| 01c_concat.smk | ingest_qc/concat.py | — | Stage 01c - concatenate per-sample QC'd objects into a per-species object. |
| 02_integrate.smk | atlas/integrate.py | — | Stage 02 - scVI integration + Leiden clustering, per species. |
| 03_cluster.smk | atlas/cluster.py | — | Stage 03 — Leiden clustering + clustree + UMAP (per species). |
| 04_annotate.smk | atlas/apply_annotation.py | results/tables/clusters/ | 04_annotate -- apply the committed atlas choices CSV. |
| 04b_reintegrate.smk | atlas/reintegrate.py | — | Stage 04b - per-lineage scVI reintegration. |
| 04c_reannotate.smk | atlas/apply_annotation.py | — | 04c_reannotate -- apply a committed per-lineage choices CSV to a reintegrated |
| 04d_merge.smk | atlas/merge_reannotated.py | results/tables/refined_ | 04d_merge -- concatenate the cleaned per-lineage objects into one per-species |
| 04e_reembed.smk | atlas/integrate.py | — | 04e_reembed -- fresh whole-tissue scVI on the refined object. Removing EXCLUDE |
| 04f_spatial_qc.smk | atlas/spatial_qc.py | results/tables/spatial_qc_ | 04f_spatial_qc -- "blank, don't drop" spatial-coordinate QC, applied AFTER |
| 05_de.smk | de/de.py, de/de_lineage.py | results/tables/de/, results/tables/markers_ | Stage 05 — marker DE (Wilcoxon). |
| 06_pseudotime.smk | pseudotime/pseudotime.py | — | Stage 06 — diffusion pseudotime (global, mouse fibroblasts + decidual cells). |
| 06b_timepoint_pseudotime.smk | pseudotime/compute_pseudotime.py | — | Scratch-only, seeded per-timepoint mouse pseudotime tests. |
| 07_spatial_pseudotime.smk | pseudotime/spatial_pseudotime.py | — | Stage 07 — per-sample spatial pseudotime + GPR gradient field (slide-tags sample |
| 08_cellrank.smk | pseudotime/cellrank.py | — | Stage 08 — CellRank fate / probabilistic layer (on the global pseudotime). |
| 09_homology.smk | phmmer_to_samap_symbols.py | — | workflow/rules/09_homology.smk        (upstream of 10_samap.smk) |
| 09_pt_genes_mi.smk | mi/pt_genes_mi.py | results/tables/mi_ | Stage 09 — MI pseudotime-associated genes (per species) + GAM heatmaps. |
| 10_samap.smk | samap/make_transferred_object.py, samap/prep_samap_input.py, samap/run_samap.py, samap/transfer_label_staining.py | results/figures/samap_transfer/, results/samap_transfer/, results/tables/samap_transfer/ | COMPLEMENTARY to 10b_samap_raw.smk: this file is the TRANSFER arm; 10b is the va |
| 10b_samap_raw.smk | samap/compare_samap_arms.py, samap/make_samap_transferred_h5ad.py, samap/prep_samap_full_input.py, samap/run_samap_arm.py, samap/samap_attach_pt.py, samap/samap_label_stain_full.py, samap/samap_transfer_pt.py | — | COMPLEMENTARY to 10_samap.smk: this file is the VANILLA/RAW arm (the paper arm); |
| 11_mi_pseudotime.smk | mi/mi_pseudotime.py | results/mi | workflow/rules/11_mi_pseudotime.smk |
| 12_spatial_gp.smk | pseudotime/spatial_pt_gp.py | results/spatial_gp/ | workflow/rules/12_spatial_gp.smk |
| 13_scenic.smk | build_consensus_aucell.py, ctx_only.py, grn_only.py, scenic/scenic_make_loom.py | — | workflow/rules/13_scenic.smk |
| 14_liana.smk | — | results/liana/ | workflow/rules/14_liana.smk |
| 15_pt_per_timepoint.smk | pseudotime/pt_per_timepoint.py | results/pt_per_timepoint/ | workflow/rules/15_pt_per_timepoint.smk |
| 16_liana_crossspecies.smk | liana/liana_crossspecies.py | results/liana/no_commot, results/liana/no_commot/, results/liana_crossspecies | workflow/rules/16_liana_crossspecies.smk |
| 17_liana_pseudotime.smk | — | results/liana_pseudotime | workflow/rules/17_liana_pseudotime.smk |
| 20_liana.smk | liana/liana.py | results/tables/liana_ | Stage 20 — ligand-receptor / cell-cell communication (per species). |
| 21_scenic.smk | scenic/scenic.py | results/tables/scenic_ | Stage 21 — SCENIC gene regulatory networks / TF activity (per species). |
| 22_crossspecies.smk | crossspecies/crossspecies.py | results/tables/conservation.csv, results/tables/mi_human.csv, results/tables/mi_mouse.csv | Stage 22 — cross-species conservation: orthogroups + threshold-free rank compari |
| 23_invitro.smk | crossspecies/invitro.py | results/tables/invitro_comparison.csv | Stage 23 — in-vitro decidualization vs in-vivo comparison. |
| 24_pseudotime_notebooks.smk | pseudotime/run_notebook.py | results/pt_notebook_rules/ | 14 -- pseudotime notebooks as DAG rules (dependency-free papermill-style runner) |
| 90_figures.smk | crossspecies/figures.py | results/figures/.done, results/tables/conservation.csv, results/tables/markers_ | Stage 90 — assemble figure panels from the tables above. |
| 98_clustree.smk | atlas/export_clustree.py, atlas/ladder_stability.py | results/figures/clustree/, results/tables/clustree/, results/tables/clustree/stability_summary_all_lineages.tsv | Clustree -- OFF the critical path, like 99_reports. |
| 99_reports.smk | atlas/apply_annotation.py | results/notebooks/ | 99_reports -- OFF the critical path on purpose. |
