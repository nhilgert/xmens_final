# xmens_final

Snakemake pipeline for Çevrim et al., *Induction of menstruation in mice reveals the regulation of menstrual
shedding* (bioRxiv 2025.10.08.681007): single-nucleus Slide-tags spatial transcriptomics of the mouse X-Mens
induced-menstruation time course and of human endometrium (15 samples), from CellBender-corrected counts and
Slide-tags coordinates to the released objects, the figure panels and the numbers quoted in the Methods.

## Run it

```bash
# Snakemake 9.23.1, scanpy and scvi-tools in the rsc_uv environment on FASRC
export PATH=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin:$PATH
snakemake --profile workflow/profiles/slurm -n      # always look at the dry run first
snakemake --profile workflow/profiles/slurm         # build the default target
```

The default target builds the release named by `deliverables.release` in `config/config.yaml`
(`data/deliverables/<release>/`, four objects and a `MANIFEST.json` with sha256 checksums), its deposit copies
(`data/deposit/`), the figure panels, the SAMap figures and the zone modules. Release and deposit outputs are
protected: Snakemake will not overwrite them.

The steps that made the objects of record are not rerun by default. They are defined only when a config flag is
set: `atlas_annotation.rerun` (atlas annotation and spatial QC), `scvi.retrain` (scVI integration, lineage
reintegration, re-embedding) and `annotation_round2.fit` (the round-2 lineage split and fits).

## Steps

One rule file per step in `workflow/rules/`, numbered in pipeline order; `workflow/Snakefile` includes them in
this order.

| Steps | What they do |
|---|---|
| 01-06 | Raw data to per-species objects: CellBender (upstream), Slide-tags coordinates, ingest, QC, QC report, concatenation (provenance; not rerun) |
| 07-12 | Atlas: scVI integration, atlas annotation, lineage reintegration and annotation (July), re-embedding, spatial QC |
| 13-17 | Round 2 (October) lineage reintegration and annotation; the release (15 defines and finalizes it with MANIFEST, 16 rebuilds its round-2 objects); deposit copies |
| 18-25 | Analyses: marker DE, homology, SAMap (transfer and raw arms), MI against pseudotime, spatial GP, SCENIC, zone modules |
| 26-28 | Review (not in the default target): clustering stability, annotation review notebooks, public annotation notebooks |
| 29-36 | Figure panels: SAMap figures, mouse tissue, mouse stroma, human tissue, human epithelium, cross-species, programs, marker dot plots |

## Where things are

| Path | Contents |
|---|---|
| `workflow/Snakefile` | The steps, in order, and the default target |
| `workflow/rules/` | One file per step (01-36) |
| `workflow/scripts/` | The Python and R scripts the rules call, by topic; `tools/` holds `check_annotations.py` |
| `workflow/envs/` | Environment specifications |
| `workflow/profiles/slurm/` | SLURM settings for the FASRC eddy partition |
| `config/` | `config.yaml` (all paths and parameters), `samples.tsv`, deposit and panel settings |
| `resources/` | Curated inputs: annotation decision tables, vocabulary, lineage spec, markers, figure settings, review records, transgene reference |
| `notebooks/` | Notebooks the rules render: `public/` (step 28), `review/` (step 27), `samap/` (step 29) |
| `paper_code/` | Methods code that no rule runs yet; it empties as each piece becomes a step (`docs/PORTING.md`) |
| `docs/` | Data provenance, the Slide-tags fork, the old-to-new rule-file names, the porting plan |

Not in git: `data/`, `results/` and `logs/` (outputs) and installed environments.

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
- Still outside the rules (being brought in before submission): `paper_code/` and the items in
  `docs/PORTING.md` (the 30 Oct human pseudotime refit, the in vitro decidualization DE, the Ulrich
  comparison, the Methods-numbers build).

## Citation

Çevrim et al. 2025, bioRxiv 2025.10.08.681007. License: to be added.
