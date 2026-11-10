# workflow/scripts/ivd: the in vitro / in vivo decidualization comparison

Ported on 2 Nov 2026 from `paper_code/ivd_comparisons` (the 31 Oct 2026 code bundle; `README.md` there gives the
history). Rules: `workflow/rules/39_ivd_comparisons.smk`; config: `ivd` in `config/config.yaml`; inputs: `resources/ivd`.
Files marked "verbatim" are byte-identical to the bundle (or to the repository they came from) and were not edited.

| folder | file | origin | role |
|---|---|---|---|
| consensus | `build_consensus.py` | verbatim, human_ivd_repro_final_defense_20260812 | the vote (collapse per study; >= 2 agreeing clean calls) |
| consensus | `rebuild_consensus123_stadtmauer_thresholds.py` | verbatim, bundle 1_consensus_list | the 22 Oct correction; imported for its functions |
| consensus | `build_outputs_stadtmauer_fix.py` | verbatim, bundle 1_consensus_list | imported for `per_gene_table` (S10) and `per_study_table` (S9) |
| consensus | `derive_consensus123_orthology.py` | adapted from cagri_ivd_end_to_end_20260815 | Ensembl 116 one-to-one orthology; inputs are arguments (logic unchanged) |
| consensus | `build_consensus123.py` | new | rule `ivd_consensus`: tables from CSV instead of .docx, vote, orthology, regression checks |
| consensus | `make_consensus_tables.py` | new | rule `ivd_consensus_tables`: S9 and S10 with the current X-Mens columns |
| xmens | `mouse_pseudobulk_job.py`, `fit_mouse_de.py`, `compare_xmens_with_references.py`, `release_comparison_functions.py` | verbatim, bundle 2_xmens | X-Mens pseudobulks, PyDESeq2 fit, comparisons with the consensus and the Wang lists |
| human | `human_pseudobulk_job.py` | adapted from bundle 3_human_slidetags | object and output folder as arguments; old-object report optional |
| human | `fit_human_de.py`, `build_human_de_tables.py` | verbatim, bundle 3_human_slidetags | the five fits and the two tables |
| ulrich | `ulrich_ivd_pipeline.py` | verbatim, bundle 4_hammoud | the 28 Oct release code; imported for its functions |
| ulrich | `upstream/run_pipeline.py` | verbatim, cagri_ivd_end_to_end_20260815 (SHA-256 1c7696a2...) | the upstream helpers the pipeline hash-checks and imports |
| ulrich | `ulrich_human_de.py` | new | rule `ivd_ulrich_de`: the human contrasts only (no X-Mens refit, no zones), same functions |
| ulrich | `config_ulrich_dS_vs_C7fib.template.yaml`, `render_ulrich_config.py` | template of the release config; new | paths as placeholders, filled per checkout |
| ulrich | `download_ulrich_atlas.py` | new | rule `ivd_ulrich_download`: CELLxGENE download or local mirror, sha256-checked |
| workbooks | `build_house_style_workbooks.py`, `build_de_workbook.py`, `add_human_comparisons_to_ivd_workbook.py` | verbatim, bundle 5_workbooks | layout helpers; the DE workbook; the human comparison functions and notes |
| workbooks | `build_ivd_workbook.py` | new | rule `ivd_workbooks`: the IVD workbook from tables alone (replaces the session code of 6_not_yet_scripts) |
| . | `verify_ivd_outputs.py` | new | rule `ivd_verify`: every output against `resources/ivd/expected` and the workbooks of record |

Run: `snakemake --profile workflow/profiles/slurm ivd_all` (dry run first). Outputs: `results/ivd/<ivd.run>/`, with
`verification/VERIFICATION.md` as the last step. The two environments are `ivd.python_fit` (PyDESeq2, `workflow/envs/ivd`)
and `ivd.python` (rsc_uv).
