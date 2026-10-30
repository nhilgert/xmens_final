# SAMap figure notebooks

One notebook per panel of the SAMap supplementary figure and of the SAMap pseudotime panels of Fig. 6. Each reads the
tables written by the 10c rules (`workflow/rules/10c_samap_figures.smk`), draws its panel, and writes the panel files
with a provenance record: input files with sha256, the notebook's md5 and commit, the font, and the counts and
statistics the panel shows.

| figure | panel | notebook | writes, under `results/figures/` |
|---|---|---|---|
| SAMap supplement | a | `sfig_samap_a_umap_mouse_pseudotime.ipynb` | `samap_supp/panels/sfig_samap_panel_a.*` |
| SAMap supplement | b | `sfig_samap_b_umap_transferred_pseudotime.ipynb` | `samap_supp/panels/sfig_samap_panel_b.*` |
| SAMap supplement | c | `sfig_samap_c_umap_human_pseudotime.ipynb` | `samap_supp/panels/sfig_samap_panel_c.*` |
| SAMap supplement | d | `sfig_samap_d_edge_share.ipynb` | `samap_supp/panels/sfig_samap_panel_d.*` |
| SAMap supplement | e | `sfig_samap_e_transferred_vs_human.ipynb` | `samap_supp/panels/sfig_samap_panel_e.*` |
| SAMap supplement | f | `sfig_samap_f_sections_shuffled.ipynb` | `samap_supp/panels/sfig_samap_panel_f.*` |
| SAMap supplement | g | `sfig_samap_g_permutation_statistics.ipynb` | `samap_supp/panels/sfig_samap_panel_g.*` |
| SAMap supplement | a-g assembled | `sfig_samap_compose.ipynb` | `samap_supp/sfig_samap.{png,svg}`, `sfig_samap_600dpi.png` |
| Fig. 6 | C | `fig6c_umap_transferred_pseudotime.ipynb` | `samap_main/umap_dark_transferred_pseudotime*.{png,svg,pdf}` |
| Fig. 6 | E | `fig6e_sections_transferred_pseudotime.ipynb` | `samap_main/<date>-Human-<phase>-transferred_pseudotime-spatial.*`, `<date>-spatial_key_transferred_pseudotime.*` |

Code shared by the supplement notebooks sits next to them: `sfig_samap_style.py` (font, text sizes, the 12-column
grid, panel files, provenance) and `sfig_samap_umap.py` (the nuclei, draw order and frame of the UMAP panels a-c).
The two Fig. 6 notebooks are self-contained.

Fig. 6A-B (species and cell-type UMAPs) and 6D (cell types in the sections) are not drawn by code in this repository.

## Running

The pipeline runs every notebook with papermill and keeps the executed copy, with the figure and the numbers behind
it, in `results/notebooks/samap/`:

```bash
snakemake samap_figures_all --cores 4 --resources mem_mb=100000   # inside a Slurm job
```

By hand, open a notebook with this directory as the working directory and run all cells, or from the repository root:

```bash
papermill notebooks/samap/<notebook>.ipynb results/notebooks/samap/<notebook>.executed.ipynb \
    --cwd notebooks/samap -p FORCE True
```

Paths are relative to the repository root, and a notebook refuses to overwrite its outputs unless `FORCE = True`.
Away from the cluster, copy `results/samap_figures/{tables,edge_share,ptnull}` (and, for Fig. 6, the human plot
table of the release) to the same paths in a checkout, or pass their paths as parameters.

Fonts: Helvetica Neue on macOS, which the versions of record use. On the cluster the supplement falls back to Nimbus
Sans and Fig. 6 to DejaVu Sans; each provenance record names the font used.

Commit the notebooks without outputs: the provenance record notes whether a notebook matches its commit.
