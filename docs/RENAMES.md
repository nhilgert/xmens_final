# Rule files renamed on 31 Oct 2026

The rule files were renumbered into one sequence in pipeline order (commit after 7fea14d). Rule names inside did
not change, so targets, jobs and outputs are the same. Use this table to read older notes, the reading kit,
methods_vs_code.csv and logs that use the old names.

| Old | New |
|---|---|
| `00_cellbender.smk` | `01_cellbender.smk` |
| `00a_slidetags.smk` | `02_slidetags.smk` |
| `00b_ingest.smk` | `03_ingest.smk` |
| `01_qc.smk` | `04_qc.smk` |
| `01b_qc_report.smk` | `05_qc_report.smk` |
| `01c_concat.smk` | `06_concat.smk` |
| `02_integrate.smk` | `07_integrate.smk` |
| `04_annotate.smk` | `08_annotate.smk` |
| `04b_reintegrate.smk` | `09_reintegrate_lineages.smk` |
| `04c_reannotate.smk` | `10_annotate_lineages.smk` |
| `04e_reembed.smk` | `11_reembed.smk` |
| `04f_spatial_qc.smk` | `12_spatial_qc.smk` |
| `04b2_reintegrate_round2.smk` | `13_reintegrate_round2.smk` |
| `04c2_annotate_round2.smk` | `14_annotate_round2.smk` |
| `04g_deliverables.smk` | `15_release.smk` |
| `04g2_release_round2.smk` | `16_release_round2.smk` |
| `95_deposit_objects.smk` | `17_deposit.smk` |
| `05_de.smk` | `18_marker_de.smk` |
| `09_homology.smk` | `19_homology.smk` |
| `10_samap.smk` | `20_samap_transfer.smk` |
| `10b_samap_raw.smk` | `21_samap_raw.smk` |
| `11_mi_pseudotime.smk` | `22_mi_pseudotime.smk` |
| `12_spatial_gp.smk` | `23_spatial_gp.smk` |
| `13_scenic.smk` | `24_scenic.smk` |
| `16_zone_modules.smk` | `25_zone_modules.smk` |
| `98_clustree.smk` | `26_review_clustree.smk` |
| `99b_review_notebooks.smk` | `27_review_notebooks.smk` |
| `99_reports.smk` | `28_reports.smk` |
| `10c_samap_figures.smk` | `29_fig_samap.smk` |
| `91_panels_mouse_tissue.smk` | `30_fig_mouse_tissue.smk` |
| `91_panels_mouse_stroma.smk` | `31_fig_mouse_stroma.smk` |
| `91_panels_human_tissue.smk` | `32_fig_human_tissue.smk` |
| `91_panels_human_epithelium.smk` | `33_fig_human_epithelium.smk` |
| `91_panels_cross_species.smk` | `34_fig_cross_species.smk` |
| `91_panels_programs.smk` | `35_fig_programs.smk` |
| `91_panels_marker_dotplots.smk` | `36_fig_marker_dotplots.smk` |

Other moves: profiles/slurm -> workflow/profiles/slurm; envs/ -> workflow/envs/; tools/check_annotations.py and
tools/gen_04c_notebooks.py -> workflow/scripts/tools/; the two epithelial notebook generators,
analysis/hammoud/fib_mapping, analysis/notebooks, notebooks/mi_dpt and three notebooks -> paper_code/;
preprocessing/reference -> resources/reference; xmens_pipeline/DATA_PROVENANCE.md -> docs/;
archive/pseudotime_mouse_v2_induced_no_apop (csv, json) -> docs/provenance/mouse_pseudotime_v2/.
Archived (out of the repository): profiles/slurm_noconda, workflow/scripts/attic, .gitmodules.

Old names still appear in a few places that rules read (editing them would make Snakemake rerun jobs); they are
updated with the next real edit of each file: config/config.yaml line 10 (preprocessing/reference),
resources/lineage_spec.yaml line 3 (tools/gen_04c_notebooks.py), docstrings of workflow/scripts/ingest_qc/qc.py,
workflow/scripts/atlas/split_lineages_round2.py, workflow/scripts/atlas/spatial_qc.py (analysis/notebooks) and
workflow/scripts/zone_modules/fig_zone_modules.py line 62 (analysis/notebooks).
