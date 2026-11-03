# SAMap figure set from the deposited final objects, checked against the 2 Nov render (3 Nov 2026)

`snakemake samap_figures_all --configfile config/samap_figures_final_objects.yaml` was run in a scratch clone
(`tmp/2026_11_03_fig6c_pipeline/repo`: commit `fd1859c` plus the changes committed with this record) and every output
was compared with the scratch render of 2 Nov 2026 behind the Drive figures
(`tmp/2026_11_02_samap_final_figures/repo/results`, drawn from the final objects before their 2 Nov 18:08 edit and
from the transfer object of record before `Samap_Final.h5ad` was reduced) by `compare_final_sets.py`
(`final_set_comparison.json`).

Of 71 files, 44 are identical, 21 are PDFs or provenance records (run dates and paths, not compared), and 6 differ:

| file | difference |
|---|---|
| `umap_dark_transferred_pseudotime{,_humanonly}.{png,svg}` (Fig. 6C) | the grey key only: "human, no mouse neighbor (181)" became "no transferred value (181)"; 16,800 pixels in the key's box; SVG rasters pixel-identical |
| `tables/joint.parquet` | the three transfer variants no figure reads are no longer written; mouse `sample` and `timepoint` carry the final objects' names (33,863 nuclei, e.g. `1_dpi_early` -> `1.5_dpi_early`) |
| `edge_share/samap_compartment_edge_share_all.csv` | no rows for SAMap's own labels, which `Samap_Final.h5ad` does not keep; the release-label rows are identical, as is the figure's table |

Identical: Fig. 6E sections and key, the cell-identity panels and `transferred_cell_type.csv`, S30 panels a-g and the
assembled supplement, the permutation control (`ptnull_statistics_table.csv`, replicates, examples), the edge-share
figure table, the cell tables and the graph block.

The 6C SVG on Drive (`xmens_figures/final_figures/Fig6_crossspecies_radial_patterning/`) equals the 2 Nov render
after normalizing dates and matplotlib's ids, embedded rasters included (`fig6c_drive_check.json`).

Of the 181 human nuclei without a transferred value, 24 have no mouse neighbor in the SAMap graph (graph block
`rows_without_mouse_edge`) and 157 have only mouse neighbors without pseudotime; all 907 such mouse nuclei are decAp.
