# methods_numbers: the pipeline-derived numbers in the computational methods

`pipeline_numbers.csv` has one row for each number in the methods text that comes from the pipeline: the value
recomputed from the current outputs, the file and step that produce it, and the words around it in the draft
(`doc_phrase`). Rows whose recomputed value differs from the draft are flagged.

## Rebuild

```bash
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final
bash workflow/scripts/methods_numbers/submit.sh results/methods_numbers/YYYY_MM_DD path/to/methods_draft.txt [EXTRA_INPUTS]
```

This submits three read-only collection jobs (about 5 minutes; 64, 80 and 8 GB) and a build job that starts when all
three succeed. The outputs are `pipeline_numbers.csv` and `pipeline_numbers_flags.md` in the run directory. The draft
must be plain text, for example extracted from the docx or PDF; `registry_core.normalize_draft` joins line breaks and
writes superscript exponents the way the document's plain text reads them: `10−84`, with the U+2212 minus sign and no
caret. `doc_phrase`, `draft_value` and `value` all use that form, so a value can be substituted into its phrase directly.

The objects are those of the release that `config/config.yaml` names under `deliverables.release`; their paths come
from `data/deliverables/<release>/MANIFEST.json` (`release_objects.py`). Cutting a new release and updating the config
is enough to point the registry at it.

`EXTRA_INPUTS` is a directory of project artifacts that no pipeline output holds yet, with a `MANIFEST.json` giving each
file's artifact version and sha256 (same format as `resources/methods_numbers/artifact_inputs/MANIFEST.json`):
`xmens_with_uninduced_20261027.zip` (the X-Mens fit with the uninduced fibroblasts: fit diagnostics, pseudobulk table
and DE table). Without it, the rows that read it are `no_source`. The permutation statistics of the SAMap pseudotime
transfer are read from `results/samap_figures/ptnull/ptnull_statistics_table.csv`, written by rule `samap_ptnull`
(`workflow/rules/10c_samap_figures.smk`). The copy used for the
2026-10-29 rebuild is `/n/netscratch/eddy_lab/Everyone/nhilgert/2026_10_29_methods_numbers/artifact_inputs`.

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
| `collect_facts.py` | reads the released objects (obs, uns, gene detection, embedding-outlier overlaps), QC logs, annotation tables, MI runs, SCENIC and in vitro tables, and the outputs of the executed v7 epithelium notebooks |
| `collect_facts2.py` | SAMap joint-graph statistics (edge totals, alignment scores, label transfer), the in vitro DE release, the 28 Oct supplementary tables (Wang summary, methods checks), reference and CellBender settings |
| `release_objects.py` | resolves the release named in `config/config.yaml` to its MANIFEST and object paths |
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
| `artifact_inputs/figS16_revision_adopted_cutoff_sweep.csv` | the Fig. S16 C-D cutoff sweep, made outside Snakemake. Its 90th-percentile row is checked against the MI runs at every build. |
| `artifact_inputs/monocle_heatmap_repro.R`, `render_crossspecies_heatmaps.py`, `disp_table_mouse.tsv.gz`, `gene_peak_order.tsv` | the pseudotime heatmaps, made outside Snakemake |
| `proteome_release_check.txt` | output of `check_proteome_release.sh` on 2026-10-26: the local proteome files carry the Ensembl 113 FTP times to the minute (release 112's files are from February 2024). There is no download record. |
| `draft_2026_10_22.txt` | the draft text the first run compared against |

Further inputs read where they live, outside the pipeline's results:

| path | why it is needed |
|---|---|
| `/n/netscratch/eddy_lab/Everyone/nhilgert/2026_10_29b_human_epithelial_layers`, `.../2026_10_29_mouse_epithelial_compartments` | executed runs of the `04c_*` epithelium notebooks (per-cluster tables, DE tables, annotated lineage objects) behind the two epithelium paragraphs; netscratch is purged, so these belong under `results/annotation/` |
| `xmens_figures/2026-10-28_supplementary_tables/csv` | `SuppTable_Wang_comparison__Summary.csv` (the X-Mens agreement with the Wang et al. lists) and `SuppTable_IVD_methods_checks.csv`; read only |
| `EXTRA_INPUTS` (above) | the 27 Oct X-Mens fit, which exists only as a project artifact and in job scratch |

`artifact_inputs/MANIFEST.json` gives each file's project artifact version and sha256.
