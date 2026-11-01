# methods_numbers: the pipeline-derived numbers in the computational methods

`pipeline_numbers.csv` has one row for each number in the methods text that comes from the pipeline: the value
recomputed from the current outputs, the file and step that produce it, and the words around it in the draft
(`doc_phrase`). Rows whose recomputed value differs from the draft are flagged.

## Rebuild

The numbers are the Snakemake target `methods_numbers_all` (`workflow/rules/97_methods_numbers.smk`, config
`methods_numbers`); `workflow/scripts/methods_numbers/README.md` describes the run and its outputs. The draft it checks
is `methods_draft.txt` in this folder. It must be plain text, for example extracted from the docx or PDF;
`registry_core.normalize_draft` joins line breaks and writes superscript exponents the way the document's plain text
reads them: `10−84`, with the U+2212 minus sign and no caret. `doc_phrase`, `draft_value` and `value` all use that
form, so a value can be substituted into its phrase directly.

## Columns

The eight requested columns come first: `key`, `value` (formatted as the draft writes it), `unit`, `description`,
`source_file`, `source_step`, `doc_phrase`, `updated`. Then `value_raw` (unrounded), `draft_value` (as written in the
draft), `status` and `note`. P-values are computed on the log scale (`scipy.stats.hypergeom.logsf`), so they cannot
underflow; for those rows `value_raw` is the p-value in e-notation derived from log10 p, and `note` gives log10 p.

`status` is one of:
- `match`: the draft agrees with the recomputed value.
- `MISMATCH`: the draft differs; `note` gives context, for example what the superseded runs gave.
- `no_source`: not recomputed; `note` says why.
- `phrase_found_Nx` or `value_not_in_phrase`: the draft wording changed, so the row's phrase in `registry_rows.py`
  needs updating.

## Files

| file | what it does |
|---|---|
| `collect_facts.py` | reads the Oct 22 objects (obs, uns, gene detection), QC logs, annotation tables, MI runs, SCENIC and in vitro tables |
| `collect_facts2.py` | SAMap joint-graph statistics (edge totals, alignment scores, label transfer), the in vitro DE release, reference and CellBender settings |
| `small_reads.py` | small text sources: SCENIC records, run logs, scripts |
| `registry_rows.py` | one entry per number: key, draft phrase, display format, and how to recompute it |
| `registry_methods.py` | the recomputations that need more than one line |
| `registry_core.py` | draft normalization, phrase matching, number formatting |
| `build_pipeline_numbers.py` | evaluates every row and writes the CSV and the flag report |
| `check_proteome_release.sh` | one-off check of the Ensembl release of the SAMap proteomes (login node; needs internet) |
| `submit.sh` | the rebuild above |

To add a number, add a `Row` in `registry_rows.py`. Its phrase must occur exactly once in the draft.

## Inputs that are not pipeline outputs

`resources/methods_numbers/` holds inputs that exist only outside the pipeline. Each is a gap to close before
submission:

| file | why it is needed |
|---|---|
| `artifact_inputs/consensus123_stadtmauer_fix_bundle.zip` | the corrected in vitro consensus (2026-10-22) that the methods describe. The cluster copy, `cagri_ivd_end_to_end_20260815/resources/consensus123.csv`, is the list before the correction (it lacks ADAMTS5 and CLU and includes KLF4 and TSHZ1) and gives 49 called, 19 concordant, 54 not called, 38.8% (45.2% expected) instead of the draft's 50, 20, 53, 40% (46%). |
| `artifact_inputs/figS16_revision_adopted_cutoff_sweep.csv` | the Fig. S16 C-D cutoff sweep, made outside Snakemake. Its 90th-percentile row is checked against the MI runs at every build. Rule `mi_overlap` now writes the same table (`overlap_percentile_sweep.csv`; identical on the 22 Oct human arm). |
| `fig5f_panel_genes.csv` | the 28 genes drawn in Fig. 5F (selection on the 30 Oct human pseudotime refit, netscratch `2026_10_30_s29_5f_bigheat/xsp_selection_fibdec_latent`). Rule `mi_heatmap_selection` reports how its selection compares with them. |
| `artifact_inputs/monocle_heatmap_repro.R`, `render_crossspecies_heatmaps.py`, `disp_table_mouse.tsv.gz`, `gene_peak_order.tsv` | the pseudotime heatmaps, made outside Snakemake |
| `proteome_release_check.txt` | output of `check_proteome_release.sh` on 2026-10-26: the local proteome files carry the Ensembl 113 FTP times to the minute (release 112's files are from February 2024). There is no download record. |
| `draft_2026_10_22.txt` | the draft text the first run compared against |

`artifact_inputs/MANIFEST.json` gives each file's project artifact version and sha256.
