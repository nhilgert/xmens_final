# Attic manifest — moved 2026-09-13 (workflow simplification; restore = mv back one level)

- **check_migrated_annotations.py** — one-time annotation migration (xmens_annotation_v1/v2, 2026-08-10/11); migration complete
- **export_migrated.py** — same migration family
- **finalize_migration.py** — same migration family
- **validate_migration_bundle.py** — same migration family
- **propagate_relabel.py** — same migration family (label propagation)
- **merge_lineages.py** — superseded by merge_reannotated.py (04d)
- **sync_qc_coords.py** — one-time spatial-coordinate sync into deliverables (2026-09-02); done
- **variant_embed.py** — pt_variants embedding experiment; superseded by settled RTK pseudotime
- **variant_embed_v2.py** — pt_variants embedding experiment; superseded by settled RTK pseudotime
- **varD_two_roots.py** — varD root experiment; superseded by settled RTK pseudotime
- **mi_stratified.py** — session MI driver; superseded by 11_mi_pseudotime.smk (mi_pseudotime.py)
- **build_mi_dpt_input.py** — session MI driver; superseded by 11_mi_pseudotime.smk + notebooks/mi_dpt
- **samap_scvi.py** — SAMap embedding experiment; superseded by 10b_samap_raw vanilla arm

## moved 2026-10-28 (rules 06b, 15 and 24 retired; their rule files are in rules/attic/, config in config/attic/)

- **compute_pseudotime.py** — 06b_timepoint_pseudotime.smk, per-timepoint mouse pseudotime tests; retired
- **pt_per_timepoint.py** — 15_pt_per_timepoint.smk, per-timepoint scVI + DPT; retired
- **run_notebook.py** — 24_pseudotime_notebooks.smk, papermill pseudotime notebooks (legacy pseudotime, never run); retired
- **liana/** — LIANA scripts; their rules moved to rules/attic/ on 2026-10-24 (LIANA is not in the paper)
