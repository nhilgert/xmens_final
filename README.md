# xmens_final

Reproducible analysis pipeline for the X-Mens / decidualization spatial
transcriptomics study (Çevrim et al. 2025, *Induction of menstruation in mice
reveals the regulation of menstrual shedding*, bioRxiv 2025.10.08.681007).

Cross-species single-nucleus **Slide-tags** spatial transcriptomics of the
endometrium: a transgenic mouse model (GsD X-Mens, an induction timecourse) and
human endometrial samples (secretory + menstrual phases). The pipeline takes
CellBender-corrected counts + spatial coordinates through QC, integration,
clustering, annotation, pseudotime, and cross-species comparison.

> **Status:** under active development. The data-provisioning + QC stages
> (00–01) are wired and validated on all 15 samples; downstream stages are in
> progress. <!-- TODO: update as stages are completed -->

---

## Quickstart

<!-- TODO: finalize once envs/base is built and the slurm profile is confirmed -->

```bash
# 1. clone with submodules (the slide-tags fork is a pinned submodule)
git clone --recurse-submodules <repo-url> xmens_final
cd xmens_final
# if already cloned without submodules:
#   git submodule update --init --recursive

# 2. build environments    # TODO: document env build (uv for base, conda for isolates)

# 3. configure: edit config/config.yaml (paths, thresholds) and config/samples.tsv

# 4. dry-run, then run
snakemake -n                              # full DAG dry-run
snakemake --cores 4 data/interim/raw/<sample>.h5ad   # single stage/sample
# TODO: document the SLURM-profile production run once profiles/slurm is validated
```

---

## Inputs & data

The pipeline starts from **CellBender-corrected count matrices** + **Slide-tags
spatial coordinates** (it does NOT re-run cellranger/CellBender — those are
upstream, documented for provenance). See **`DATA_PROVENANCE.md`** for the full
record: the 15-sample inventory, where the inputs live, the reference genomes
(custom transgenic mouse + stock human), and the key data decisions.

- `config/samples.tsv` — the sample manifest (one row per sample). Adding a row
  adds a sample to the DAG. Authoritative source for `species` and `timepoint`.
- `config/config.yaml` — all parameters/thresholds, pinned to the paper Methods.

Raw sequencing data: available on request, to be released on publication.
<!-- TODO: add accession/DOI when data is deposited -->

---

## Pipeline structure

Built with [Snakemake](https://snakemake.readthedocs.io/). Three layers:

1. **Spine** (`workflow/rules/00–09`) — the linear must-reproduce path:
   ingest → QC → integrate → cluster → annotate → DE → pseudotime → spatial → …
2. **Analysis modules** (`workflow/rules/20+`) — LIANA, SCENIC, cross-species,
   in-vitro.
3. **`src/xmens/`** — installable shared utilities, editable-installed into envs.

```
xmens_final/
├── config/            config.yaml + samples.tsv  (drive everything)
├── workflow/
│   ├── Snakefile      includes + rule all
│   ├── rules/         one .smk per stage
│   └── scripts/       one script per rule (uses the snakemake object)
├── src/xmens/         installable shared utils
├── envs/              per-tool environments
├── external/slidetags pinned submodule: the Slide-tags mapping fork (provenance)
├── preprocessing/     documented upstream steps (reference build, etc.)
├── analysis/          exploratory notebooks (not part of the DAG)
├── data/              raw/interim  (gitignored)
└── results/           figures/tables  (gitignored)
```

<!-- TODO: a figure→rule map (which paper figure each stage produces) -->

---

## Environments

Mixed by design — uv for the modern Python spine, conda for legacy-dependency
tools:

- **Driver + spine** (snakemake + scanpy): <!-- TODO: rsc_uv now / envs/base later -->
- **`01_qc`** (DoubletDetection 4.2): conda env via the rule's `conda:` directive
  (run snakemake with `--use-conda`). <!-- TODO: replace abs path with committed env.yml -->
- **scenic / samap**: existing conda envs (see `config.yaml` `envs:`).
- **slide-tags mapping**: lives entirely inside the submodule; not run by this DAG.

<!-- TODO: commit lockfiles (uv.lock / environment.yml) for full reproducibility -->

---

## Reproducibility notes

This repo aims for **analytical + provenance reproducibility**, not bit-identity
(CellBender and scVI are stochastic; bit-identical output is not a goal). The
recipe is fully pinned; conclusions should be robust to numerical noise.

- Versions pinned in `config.yaml` (cellranger 7.2.0, CellBender 0.2.2,
  scvi-tools 1.1.2, DoubletDetection 4.2, …).
- Expected post-QC cell counts as a regression target: **57,242 mouse /
  69,525 human** (from the paper).
- Stochastic upstream steps (CellBender) are checkpointed — the pipeline consumes
  their output rather than re-running them.

---

## Citation

<!-- TODO: add CITATION.cff and full citation on publication -->
Çevrim et al. 2025, bioRxiv 2025.10.08.681007.

## License

<!-- TODO: add LICENSE -->

