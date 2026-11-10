# Correction to the 123-gene in vitro decidualization consensus: Stadtmauer table genes

Run 2026-10-22, Slurm job `e8f91bb6-d6dc-4059-bb17-40242f99c734` (partition eddy, Python 3.11.11, pandas 2.3.3).

## The problem
`derive_signatures.py` (lines 98-99, 181, 186) gave every gene in Stadtmauer & Wagner 2022 Supplementary Tables S2/S3 a direction for each decidualizing treatment from the sign of (treatment TPM − max(base TPM, 3)). The FDR and fold change printed in the same tables were not used. As a result, KLF4 and TSHZ1, which do not change significantly under either treatment, received a Stadtmauer vote and entered the list (KLF4 after a hand-entered "up", `stadtmauer_manual_overlay.csv` line 12), while ADAMTS5 and CLU, significantly up with MPA + cAMP and unchanged with MPA + PGE2, were treated as conflicting and excluded.

## The fix
A tabulated gene now gets a direction for a treatment only if that treatment meets the paper's own criteria versus base medium (p. 159): FDR < 10^-6, fold change > 1.5, expressed > 3 TPM (in either condition compared). Directions for the 37 genes reported only in the Results text are unchanged; the KLF4 hand-entered row is dropped. The vote, the orthology derivation and the X-Mens scoring are the original code, imported or run unchanged.

## Checks (all passed)
1. The new table parser, applied with the old rule, reproduces all 32 table-derived rows of `signatures/stadtmauer2022_hits.csv`, and KLF4's pre-override up/down conflict.
2. The unchanged vote on the original inputs reproduces `resources/consensus123.csv` exactly.
3. The unchanged orthology script on the released list reproduces `resources/consensus123_mouse_orthology.csv` byte for byte (SHA-256 fce4ba85…).
4. The unchanged scoring (`run_pipeline.build_consensus_ortholog_support`) on the released list reproduces the verified release's support rows exactly.

## Result

| list | genes | 1:1 orthologs | called in X-Mens | same | opposite | % same | % chance |
|---|---|---|---|---|---|---|---|
| released (as built) | 123 | 114 | 49 | 19 | 30 | 38.8 | 45.2 |
| **corrected** | 123 | 114 | 50 | 20 | 30 | 40.0 | 45.6 |
| corrected, Aghajanova 14 d only | 122 | 113 | 50 | 20 | 30 | 40.0 | 45.6 |
| corrected, no text-derived Stadtmauer calls | 102 | 94 | 37 | 12 | 25 | 32.4 | 40.5 |
| corrected, TPM filter read as a mean | 123 | 114 | 50 | 20 | 30 | 40.0 | 45.6 |
| corrected, plus SOX4 from the text | 124 | 115 | 51 | 20 | 31 | 39.2 | 45.1 |

Genes changed: ADAMTS5 added (X-Mens: supported); CLU added (X-Mens: tested not called); KLF4 removed (X-Mens: tested not called); TSHZ1 removed (X-Mens: tested not called).
The list stays at 123 genes (93 up, 30 down) with 114 one-to-one orthologs. ADAMTS5 is called in X-Mens in the same direction (log2 fold change +1.86), so concordance moves from 19 of 49 to 20 of 50.

## Inputs (read only; path and first 12 characters of the SHA-256)
- Stadtmauer tables: `…/human_ivd_repro_final_defense_20260812/papers/stadtmauer/tables2.docx` (95ede8ade759), `…/human_ivd_repro_final_defense_20260812/papers/stadtmauer/tables3.docx` (dc91010b0f8d); hand-entered directions: `…/human_ivd_repro_final_defense_20260812/papers/stadtmauer/stadtmauer_manual_overlay.csv` (7c6d13e87f84)
- Signatures: `…/human_ivd_repro_final_defense_20260812/signatures/aghajanova2011_invitro_decidualization_nonendo.csv` (586d70aaa829), `…/human_ivd_repro_final_defense_20260812/signatures/rytkonen2019_complete_signature.csv` (3afdba6091f7), `…/human_ivd_repro_final_defense_20260812/signatures/stadtmauer2022_hits.csv` (cb2bbf667c11), `…/human_ivd_repro_final_defense_20260812/signatures/lucas2020_S2vsS3_decidualization_onset.csv` (23354acde776)
- Vote: `…/human_ivd_repro_final_defense_20260812/build_consensus.py` (a4759720c3c9); orthology: `…/cagri_ivd_end_to_end_20260815/derive_consensus123_orthology.py` (d4b22175d67f) with `…/cagri_ivd_end_to_end_20260815/resources/consensus123_heca_mapping.csv` (b374c377a3ef); scoring: `…/cagri_ivd_end_to_end_20260815/run_pipeline.py` (1c7696a267f0), `…/cagri_ivd_end_to_end_20260815/config.yaml` (460106546504)
- Released list and orthology: `…/cagri_ivd_end_to_end_20260815/resources/consensus123.csv` (dab23d9fb1ff), `…/cagri_ivd_end_to_end_20260815/resources/consensus123_mouse_orthology.csv` (fce4ba85d071)
- X-Mens fit and support rows: verified release 2026-10-19 v6 (84516a99f2f0, 108617ed9135)
Full paths and checksums are in `job_outputs/manifest.json`.

## Files
- `rebuild_consensus123_stadtmauer_thresholds.py`: the rebuild, all checks, and the sensitivity lists (run by `run_rebuild_job.sh`).
- `build_outputs_stadtmauer_fix.py`: per-gene and per-study tables, workbook, summary CSV and methods text from the job outputs.
- `job_outputs/`: everything the job wrote, including `consensus123_stadtmauer_authors_criteria.csv` (the corrected list, same columns as `resources/consensus123.csv`), its orthology table, the X-Mens support rows, the parsed Stadtmauer tables with both rules, `run_log.txt` and `manifest.json`.
- `outputs/`: the regenerated deliverables.

## Not changed
- Nothing in either repository was written. The pipeline still reads the released `resources/consensus123.csv`.
- Releases built on the released list keep their numbers: the Cagri/Ulrich handoff workbook (e.g. "7 of 123", "41 of 123") and the September HECA F release ("38 of 123").
- Messages already sent with 49 / 19 / 30 (39%) should be updated to 50 / 20 / 30 (40%).

## To adopt the corrected list in the pipeline (not done; needs your go-ahead)
Copy `consensus123_stadtmauer_authors_criteria.csv` and its orthology table into `resources/` under new names, point `config.yaml` at them with new checksums, extend `consensus123_heca_mapping.csv` with ADAMTS5 and CLU, and cut a new release; the expected X-Mens result is the corrected row above.
