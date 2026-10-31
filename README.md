# xmens_final

Snakemake pipeline for Çevrim et al., *Induction of menstruation in mice reveals the regulation of menstrual
shedding* (bioRxiv 2025.10.08.681007): single-nucleus Slide-tags spatial transcriptomics of the mouse X-Mens
induced-menstruation time course and of human endometrium (15 samples), from CellBender-corrected counts and
Slide-tags coordinates to the released objects, the figure panels and the numbers quoted in the Methods.

## Run it

```bash
# Snakemake 9.23.1, scanpy and scvi-tools in the rsc_uv environment on FASRC
export PATH=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin:$PATH
snakemake --profile profiles/slurm -n      # always look at the dry run first
snakemake --profile profiles/slurm         # build the default target
```

The default target builds the release named by `deliverables.release` in `config/config.yaml`
(`data/deliverables/<release>/`, four objects and a `MANIFEST.json` with sha256 checksums), its deposit copies
(`data/deposit/`), the figure panels, the SAMap figures and the zone modules. Release and deposit outputs are
protected: Snakemake will not overwrite them.

The stages that made the objects of record are not rerun by default. They are defined only when a config flag is
set: `atlas_annotation.rerun` (atlas annotation and spatial QC), `scvi.retrain` (scVI integration, lineage
reintegration, re-embedding) and `annotation_round2.fit` (the round-2 lineage split and fits).

## Where things are

| Path | Contents |
|---|---|
| `workflow/Snakefile` | Included rule files and the default target |
| `workflow/rules/` | One file per stage, numbered in pipeline order |
| `workflow/scripts/` | The Python and R scripts the rules call, by topic |
| `config/` | `config.yaml` (all paths and parameters), `samples.tsv`, deposit and panel settings |
| `resources/` | Curated inputs: annotation decision tables (`annotations/`, `annotations/round2/`), vocabulary, |
| | lineage spec, figure settings, review records |
| `profiles/slurm/` | SLURM settings for the FASRC eddy partition |
| `notebooks/` | `public/` and `review/` are rendered by rules 99 and 99b, `samap/` by 10c; `mi_dpt/` holds the MI overlap record |
| `tools/` | `check_annotations.py` and the generators of the public annotation notebooks |
| `analysis/hammoud/` | Code of the comparison with the Ulrich et al. atlas (cited in the Methods; not yet a rule) |
| `external/slidetags` | Pinned submodule: the Slide-tags mapping code (provenance) |
| `xmens_pipeline/DATA_PROVENANCE.md` | Sample inventory, input locations and reference genomes |

Not in git: `data/`, `results/` and `logs/` (outputs) and installed environments.

## Stages

| Rule files | Stage |
|---|---|
| `00a`-`01c` | Ingest, CellBender (upstream), QC, per-species concatenation (provenance; not rerun) |
| `02`, `04`, `04e`, `04f` | scVI integration, atlas annotation, re-embedding, spatial QC |
| `04b`, `04c` / `04b2`, `04c2` | Lineage reintegration and annotation: July / round 2 (October) |
| `04g2`, `04g`, `95` | Release: rebuild, re-embed and finish (round 2), finalize with MANIFEST, deposit copies |
| `05`, `09`, `10`, `10b`, `11`, `12`, `13`, `16` | DE, homology, SAMap, MI against pseudotime, spatial GP, SCENIC, zone modules |
| `10c`, `91_panels_*` | Figures |
| `98`, `99`, `99b` | Clustering stability, annotation notebooks, annotation review kit |

## Release of record

v9 = `data/deliverables/2026_10_30` (`Cevrim_XMens_All_Cells`, `Cevrim_XMens_Only_Fib_Dec_Cells`,
`Cevrim_Human_All_Cells`, `Cevrim_Human_Only_Fib_Dec_Cells`); deposit copies in
`data/deposit/2026_10_30c/{with,without}_new_cell_types/`.

## Reproducibility notes

- QC (`qc.py`) keeps 68,069 mouse and 70,363 human nuclei; the QC outputs of 26 July 2026 are the pinned
  starting point (a rerun on 29 Oct reproduced the threshold filters exactly; DoubletDetection calls differed for
  about 1% of nuclei).
- scVI fits repeat exactly only within one GPU class; release embeddings are fitted on A40 GPUs and record device
  and versions in `uns['scvi_fit']`.
- Still outside the rules (being brought in before submission): the Ulrich comparison (`analysis/hammoud/`), the
  30 Oct human pseudotime refit, the in vitro decidualization DE, and the Methods-numbers build.

## Citation

Çevrim et al. 2025, bioRxiv 2025.10.08.681007. License: to be added.
