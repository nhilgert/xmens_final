# SAMap figures: Fig. 6C, 6E, the cell-identity panel and the SAMap supplement

Rules: `workflow/rules/29_fig_samap.smk` (target `samap_figures_all`). Each figure is drawn by its notebook in
`notebooks/samap/`, run with papermill, and writes a `.provenance.json` naming its inputs, their sha256 and the
commit it was drawn at.

## Two object sets

| set | settings | output tag | objects read |
|---|---|---|---|
| release of record, drawn by `snakemake all` | `config/config.yaml`, `samap_figures` | release version (`v9`) | fib/dec objects of `data/deliverables/<release>`; the SAMap transfer object of record, `runs_raw/vanilla_raw/samap_transferred_pt_harmonized.h5ad` (sha256 `62fb75d7…`) |
| deposited final objects | `config/samap_figures_final_objects.yaml` | `final` | the two fib/dec objects and their `embedding_outlier` tables, pinned in `resources/final_objects/MANIFEST.json`; `Samap_Final.h5ad` (sha256 `ce15b9d4…`) |

The figures of record are the `final` set (decided 3 Nov 2026). `snakemake all` stays on the release of record until
the pipeline builds the final objects itself. Then the settings of `config/samap_figures_final_objects.yaml` move into
`config.yaml`, pointing at the pipeline's outputs, and that file and `resources/final_objects/` are removed.

## Drawing the final set

The rules read `.h5ad` files, so run them in a Slurm job, from the repository root:

```bash
snakemake samap_fig6c --configfile config/samap_figures_final_objects.yaml --cores 4
snakemake samap_figures_all --configfile config/samap_figures_final_objects.yaml --cores 4
```

Outputs: `results/plot_tables/final_objects/`, `results/samap_figures/final/` (tables, graph block, edge shares,
permutation control), `results/figures/samap_main/final/` (6C, 6E, cell identity) and
`results/figures/samap_supp/final/`. Every object is checked by sha256 before anything is drawn; an edited or purged
object stops the run.

## Where each input of the final set comes from

1. **SAMap fit and joint export** (`runs_raw/vanilla_raw/samap.pkl`, 1 Sep; `samap_transferred.h5ad`, 2 Sep). Rules
   `samap_raw_prep`, `samap_raw_fit`, `samap_raw_stain` and `samap_raw_transferred` (`21_samap_raw.smk`); the outputs
   were made by sbatch scripts before the rules existed, with the same scripts, and the prep step read the netscratch
   `2026_08_28` fib/dec objects rather than the `2026_09_15` release its rule names (header of `21_samap_raw.smk`).
2. **Pseudotime transfer** (24 Oct). Rules `samap_raw_pt_transfer` and `samap_raw_attach_pt`, variant `_harmonized`:
   mouse `pt_harmonized` of `2026_10_22_flagged_harmonized`, rank-normalized, averaged over each human nucleus's
   cross-species SAMap edges. Made by hand with no log; rebuilt by the rules on 3 Nov 2026 with identical data
   (`docs/provenance/samap_transfer_rebuild_2026_11_03`).
3. **`Samap_Final.h5ad`** (2 Nov). The object of step 2 with obs cut to six columns,
   `transferred_mm_dpt_pseudotime_ranknorm` renamed `transferred_mouse_pseudotime`, and every dataset outside obs left
   byte-identical; made by `tmp/2026_11_02_samap_final_figures/slim_samap_obs.py` (report beside it). No rule yet.
4. **Fib/dec final objects** (2 Nov). Deposited by hand; on 2 Nov at 18:08 their `embedding_outlier` flags moved out
   of obs into the tables in `final_objects/embedding_outliers/` and obsm `X_diffmap_pseudotime` was renamed
   `X_diffmap` (edit record in that folder). Not built by the pipeline yet.
5. **Human pseudotime of record** (`config.yaml`, `pseudotime_of_record.human`, sha256 checked). The human object's
   obs `pseudotime` must equal it on every nucleus (`samap_figure_tables` checks this).

## Fig. 6C

Each human nucleus colored by the mouse pseudotime transferred to it (step 2), on the joint SAMap UMAP, with the
mouse nuclei in dark grey beneath (a second version leaves them out). Fixed 0-1 viridis, values as stored, never
re-ranked; every nucleus drawn, `embedding_outlier` included (1 Nov 2026). Human nuclei without a transferred value
are grey, keyed "no transferred value": 181 nuclei, 24 with no mouse neighbor in the SAMap graph and 157 whose mouse
neighbors are all decAp, the one mouse subtype without pseudotime. Square frame fitted to the whole joint embedding,
so the UMAP panels of Fig. 6 overlay; 90 × 90 mm, 600 dpi, draw-order seed 3.

The `final` set drawn by these rules on 3 Nov 2026 was compared file by file with the 2 Nov scratch render behind
the Drive figures (`docs/provenance/samap_figures_final_check_2026_11_03`). Every figure is identical but for the grey
key of 6C; two tables differ where the deposited objects do (the record lists how), and the permutation statistics
are identical. The Drive 6C SVG equals the 2 Nov render, embedded rasters included.

## Open items

1. A rule for step 3, the reduction of the transfer object to `Samap_Final.h5ad`.
2. Pipeline-built final objects (step 4); then retire `config/samap_figures_final_objects.yaml` as above.
3. The deposited objects sit on purgeable netscratch (`/n/netscratch/eddy_lab/Everyone/nhilgert/final_objects`).
4. Fig. 6A, 6B, 6D and 6E of `34_fig_cross_species.smk` still draw from the release of record.
