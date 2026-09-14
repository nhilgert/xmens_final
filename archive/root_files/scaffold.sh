#!/usr/bin/env bash
# =============================================================================
# scaffold.sh — build the ENTIRE xmens_final repo from scratch (run once).
#
#   bash scaffold.sh
#
# Self-contained: creates the foundation (config, package, base env, real
# ingest/qc, spine rules) AND every downstream stage (rule + placeholder script),
# per-tool env stubs, SLURM profile, and a complete Snakefile.
#
# Idempotent: never overwrites an existing file (so it's safe to re-run to fill
# gaps); the only exception is the Snakefile, which is regenerated after backing
# up to Snakefile.bak. Placeholder scripts just `touch` their outputs, so the
# full pipeline is runnable immediately — replace one stub at a time with real logic.
# =============================================================================
set -euo pipefail

mk() {  # mk <path> — write heredoc from stdin only if the file is absent
  local p="$1"; mkdir -p "$(dirname "$p")"
  if [ -f "$p" ]; then echo "  skip (exists): $p"; cat >/dev/null
  else cat > "$p"; echo "  wrote: $p"; fi
}

echo ">> directories"
mkdir -p config workflow/rules workflow/scripts \
         envs/{base,scenic,liana,cellrank,samap,r} \
         resources/{orthology,markers,signatures,lr_db} \
         data/{raw,interim} results/{qc,tables,figures} logs \
         profiles/slurm analysis/{invitro_decidualization,exploratory} src/xmens

# ---------------------------------------------------------------------------
echo ">> foundation: config"
mk config/config.yaml <<'EOF'
# All parameters/thresholds live here — the run is described by this file + samples.tsv.
seed: 0
qc:
  min_genes: 200
  min_counts: 500
  max_pct_mito: 20
  # NOTE: in shedding tissue, high-mito/low-count cells can be REAL dying/shed cells.
mito_prefix:
  mouse: "mt-"
  human: "MT-"
scvi:       {n_hvg: 3000, n_latent: 30, n_layers: 2, batch_key: sample_id}
leiden:     {resolution: 0.1, knn: 100}
pseudotime: {knn: 15, n_diffmap: 10}
mi:         {n_neighbors: 3, n_perm: 1000, top_percentile: 90}
EOF

if [ -f config/samples.tsv ]; then echo "  skip (exists): config/samples.tsv"; else
  { printf 'sample_id\tspecies\ttimepoint\tmodality\tpath\n'
    printf 'mm_unind_1\tmouse\tuninduced\tslidetags\tdata/raw/mm_unind_1.h5ad\n'
    printf 'mm_3dpi_1\tmouse\t3dpi\tslidetags\tdata/raw/mm_3dpi_1.h5ad\n'
    printf 'mm_4dpi_1\tmouse\t4dpi\tslidetags\tdata/raw/mm_4dpi_1.h5ad\n'
    printf 'hs_secr_1\thuman\tsecretory\tslidetags\tdata/raw/hs_secr_1.h5ad\n'
    printf 'hs_mens_1\thuman\tmenstrual\tslidetags\tdata/raw/hs_mens_1.h5ad\n'
  } > config/samples.tsv; echo "  wrote: config/samples.tsv (edit to your real samples)"
fi

echo ">> foundation: package + envs + gitignore"
mk pyproject.toml <<'EOF'
[project]
name = "xmens"
version = "0.0.1"
description = "Shared utils for the X-Mens analysis."
requires-python = ">=3.10"
dependencies = []
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"
[tool.hatch.build.targets.wheel]
packages = ["src/xmens"]
EOF

mk envs/base/pyproject.toml <<'EOF'
[project]
name = "xmens-base-env"
version = "0.0.1"
requires-python = ">=3.10"
dependencies = ["snakemake>=8","snakemake-executor-plugin-slurm","pandas",
                "scanpy","anndata","scvi-tools==1.1.2","networkx","pybiomart","xmens"]
[tool.uv.sources]
xmens = { path = "../..", editable = true }
EOF

mk .gitignore <<'EOF'
data/raw/
data/interim/
results/figures/
*.h5ad
.venv*/
**/.venv*/
__pycache__/
.snakemake/
*.egg-info/
EOF

echo ">> foundation: src/xmens"
mk src/xmens/__init__.py <<'EOF'
"""xmens — shared utilities (io, gpu, ortho, spatial, lr, ...)."""
__version__ = "0.0.1"
EOF
mk src/xmens/io.py <<'EOF'
import pandas as pd
def read_samples(path="config/samples.tsv"):
    return pd.read_csv(path, sep="\t").set_index("sample_id", drop=False)
EOF
mk src/xmens/ortho.py <<'EOF'
"""Ensembl Compara orthogroups with many-to-many relationships (connected components).
Do NOT filter on orthology confidence by default — it deletes expanded families
(e.g. 25/26 prolactin edges are low-confidence)."""
from __future__ import annotations
import networkx as nx, pandas as pd

def load_compara_homology(host="http://www.ensembl.org",
                          dataset="mmusculus_gene_ensembl", target="hsapiens"):
    import pybiomart
    d = pybiomart.Dataset(name=dataset, host=host)
    df = d.query(attributes=["ensembl_gene_id",
        f"{target}_homolog_ensembl_gene", f"{target}_homolog_orthology_type",
        f"{target}_homolog_orthology_confidence"])
    df.columns = ["mm","hs","type","conf"]
    return df.dropna(subset=["hs"])

def orthogroups_from_homology(df, drop_low_confidence=False):
    if drop_low_confidence: df = df[df["conf"] == 1]
    g = nx.Graph(); g.add_edges_from((("mm",a),("hs",b)) for a,b in zip(df["mm"],df["hs"]))
    return list(nx.connected_components(g))

def gene_to_orthogroup(orthogroups):
    return {gene:i for i,og in enumerate(orthogroups) for gene in og}
EOF
for m in gpu spatial lr pseudotime plotting; do
  mk "src/xmens/$m.py" <<EOF
"""xmens.$m — TODO: migrate your $m utilities here."""
EOF
done

# ---------------------------------------------------------------------------
echo ">> spine: real ingest + qc"
mk workflow/scripts/ingest.py <<'EOF'
"""Stage 00b — load one matrix into AnnData, tag with sample metadata."""
import scanpy as sc
snk = snakemake  # noqa: F821
path = str(snk.input[0]); meta = snk.params.meta
adata = sc.read_10x_h5(path) if path.endswith(".h5") else sc.read_h5ad(path)
adata.var_names_make_unique()
for k, v in meta.items(): adata.obs[k] = v
adata.obs_names = [f"{meta['sample_id']}_{bc}" for bc in adata.obs_names]
adata.write_h5ad(snk.output[0])
EOF

mk workflow/scripts/qc.py <<'EOF'
"""Stage 01 — explicit, logged cell QC (one species)."""
import logging
import anndata as ad, pandas as pd, scanpy as sc
snk = snakemake  # noqa: F821
q = snk.params.qc; mito_prefix = snk.params.mito_prefix
logging.basicConfig(filename=str(snk.log[0]), level=logging.INFO, format="%(asctime)s  %(message)s")
adata = ad.concat([sc.read_h5ad(f) for f in snk.input], join="outer", merge="same")
n0 = adata.n_obs
adata.var["mito"] = adata.var_names.str.startswith(mito_prefix)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, inplace=True)
keep = ((adata.obs["n_genes_by_counts"] >= q["min_genes"]) &
        (adata.obs["total_counts"]      >= q["min_counts"]) &
        (adata.obs["pct_counts_mito"]   <= q["max_pct_mito"]))
s = (pd.DataFrame({"sample": adata.obs["sample_id"].values, "kept": keep.values})
       .groupby("sample")["kept"].agg(n_total="size", n_kept="sum"))
s["n_dropped"] = s.n_total - s.n_kept
s.to_csv(snk.output.summary)
logging.info("%s: %d -> %d (%.1f%% dropped)", snk.wildcards.species, n0, int(keep.sum()), 100*(1-keep.mean()))
adata.obs.loc[~keep].to_csv(snk.output.dropped)
adata[keep].copy().write_h5ad(snk.output.kept)
EOF

mk workflow/rules/00b_ingest.smk <<'EOF'
rule ingest:
    input:  lambda w: samples.loc[w.sample, "path"]
    output: "data/interim/raw/{sample}.h5ad"
    params: meta=lambda w: samples.loc[w.sample].to_dict()
    log:    "logs/ingest_{sample}.log"
    script: "../scripts/ingest.py"
EOF

mk workflow/rules/01_qc.smk <<'EOF'
rule qc:
    input:  species_samples
    output: kept="data/interim/qc_passed_{species}.h5ad",
            dropped="results/qc/dropped_{species}.csv",
            summary="results/qc/qc_summary_{species}.csv"
    params: qc=config["qc"], mito_prefix=lambda w: config["mito_prefix"][w.species]
    log:    "logs/qc_{species}.log"
    script: "../scripts/qc.py"
EOF

echo ">> downstream: placeholder scripts"
for s in slidetags integrate cluster annotate de pseudotime spatial_pseudotime \
         cellrank pt_genes_mi liana scenic crossspecies invitro figures; do
  mk "workflow/scripts/$s.py" <<'EOF'
"""STUB — replace with real logic. For now just creates declared outputs so the
full DAG runs. Implement: read snakemake.input, work, write snakemake.output."""
from pathlib import Path
snk = snakemake  # noqa: F821
for f in snk.output:
    p = Path(str(f)); p.parent.mkdir(parents=True, exist_ok=True); p.touch()
print("STUB:", snk.rule, "->", list(snk.output))
EOF
done

echo ">> downstream: stage rules"
mk workflow/rules/00a_slidetags.smk <<'EOF'
rule slidetags:        # raw sequencing -> per-cell (x,y)
    input:  lambda w: samples.loc[w.sample, "path"]
    output: "data/interim/coords/{sample}.csv"
    log:    "logs/slidetags_{sample}.log"
    script: "../scripts/slidetags.py"
EOF
mk workflow/rules/02_integrate.smk <<'EOF'
rule integrate:        # scVI (per species)
    input:  "data/interim/qc_passed_{species}.h5ad"
    output: "data/interim/latent_{species}.h5ad"
    params: scvi=config["scvi"]
    log:    "logs/integrate_{species}.log"
    script: "../scripts/integrate.py"
EOF
mk workflow/rules/03_cluster.smk <<'EOF'
rule cluster:          # Leiden + clustree + UMAP
    input:  "data/interim/latent_{species}.h5ad"
    output: "data/interim/clustered_{species}.h5ad"
    params: leiden=config["leiden"]
    log:    "logs/cluster_{species}.log"
    script: "../scripts/cluster.py"
EOF
mk workflow/rules/04_annotate.smk <<'EOF'
rule annotate:         # markers, lineage re-integration, MAD outlier drop
    input:  "data/interim/clustered_{species}.h5ad"
    output: "data/interim/annotated_{species}.h5ad"
    log:    "logs/annotate_{species}.log"
    script: "../scripts/annotate.py"
EOF
mk workflow/rules/05_de.smk <<'EOF'
rule de:               # marker DE (Wilcoxon)
    input:  "data/interim/annotated_{species}.h5ad"
    output: "results/tables/markers_{species}.csv"
    log:    "logs/de_{species}.log"
    script: "../scripts/de.py"
EOF
mk workflow/rules/06_pseudotime.smk <<'EOF'
rule pseudotime_global:   # diffusion pseudotime (mouse)
    input:  "data/interim/annotated_mouse.h5ad"
    output: "data/interim/pt/global_mouse.h5ad"
    params: pt=config["pseudotime"]
    log:    "logs/pseudotime_global.log"
    script: "../scripts/pseudotime.py"
EOF
mk workflow/rules/07_spatial_pseudotime.smk <<'EOF'
rule spatial_pseudotime:  # per-sample spatial pt + GPR field
    input:  coords="data/interim/coords/{sample}.csv"
    output: "data/interim/pt/spatial_{sample}.h5ad"
    log:    "logs/spatial_pt_{sample}.log"
    script: "../scripts/spatial_pseudotime.py"
EOF
mk workflow/rules/08_cellrank.smk <<'EOF'
rule cellrank:         # fate/probabilistic layer  (env: cellrank)
    input:  "data/interim/pt/global_mouse.h5ad"
    output: "data/interim/pt/cellrank_mouse.h5ad"
    log:    "logs/cellrank.log"
    script: "../scripts/cellrank.py"
EOF
mk workflow/rules/09_pt_genes_mi.smk <<'EOF'
rule pt_genes_mi:      # MI pseudotime-associated genes (per species)
    input:  "data/interim/annotated_{species}.h5ad"
    output: "results/tables/mi_{species}.csv"
    params: mi=config["mi"]
    log:    "logs/pt_genes_mi_{species}.log"
    script: "../scripts/pt_genes_mi.py"
EOF
mk workflow/rules/20_liana.smk <<'EOF'
rule liana:            # ligand-receptor  (env: liana)
    input:  "data/interim/annotated_{species}.h5ad"
    output: "results/tables/liana_{species}.csv"
    log:    "logs/liana_{species}.log"
    script: "../scripts/liana.py"
EOF
mk workflow/rules/21_scenic.smk <<'EOF'
rule scenic:           # GRN / TF activity  (env: scenic)
    input:  "data/interim/annotated_{species}.h5ad"
    output: "results/tables/scenic_{species}.csv"
    log:    "logs/scenic_{species}.log"
    script: "../scripts/scenic.py"
EOF
mk workflow/rules/22_crossspecies.smk <<'EOF'
rule crossspecies:     # orthogroups + threshold-free rank comparison
    input:  human="results/tables/mi_human.csv", mouse="results/tables/mi_mouse.csv"
    output: "results/tables/conservation.csv"
    log:    "logs/crossspecies.log"
    script: "../scripts/crossspecies.py"
EOF
mk workflow/rules/23_invitro.smk <<'EOF'
rule invitro:          # in-vitro vs in-vivo decidualization
    input:  "data/interim/annotated_human.h5ad"
    output: "results/tables/invitro_comparison.csv"
    log:    "logs/invitro.log"
    script: "../scripts/invitro.py"
EOF
mk workflow/rules/90_figures.smk <<'EOF'
rule figures:
    input:  expand("results/tables/markers_{sp}.csv", sp=SPECIES),
            expand("results/tables/mi_{sp}.csv", sp=SPECIES),
            "results/tables/conservation.csv"
    output: "results/figures/.done"
    log:    "logs/figures.log"
    script: "../scripts/figures.py"
EOF

echo ">> per-tool env stubs + profile + resources"
for tool in scenic liana cellrank samap; do
  case $tool in
    scenic) dep='"pyscenic"';; liana) dep='"liana"';;
    cellrank) dep='"cellrank"';; samap) dep='"samap"';;
  esac
  mk "envs/$tool/pyproject.toml" <<EOF
[project]
name = "xmens-$tool-env"
version = "0.0.1"
requires-python = ">=3.10"
dependencies = [$dep, "scanpy", "anndata", "xmens"]
[tool.uv.sources]
xmens = { path = "../..", editable = true }
EOF
done
mk envs/r/README.md <<'EOF'
# R env (uv is Python-only). Use renv:
#   R -e 'renv::init()'; renv::install(c("clustree","fgsea","RRHO2")); renv::snapshot()
EOF
mk profiles/slurm/config.yaml <<'EOF'
executor: slurm
jobs: 50
default-resources:
  slurm_partition: "shared"   # <- set your partition
  runtime: 120
  mem_mb: 16000
  cpus_per_task: 4
set-resources:
  integrate: {slurm_partition: "gpu", slurm_extra: "'--gres=gpu:1'"}
  cellrank:  {mem_mb: 64000}
EOF
mk resources/orthology/README.md <<'EOF'
Compara homology dump here (record Ensembl release + download date).
EOF

echo ">> Snakefile (backing up any existing to Snakefile.bak)"
[ -f workflow/Snakefile ] && cp workflow/Snakefile workflow/Snakefile.bak
cat > workflow/Snakefile <<'EOF'
import pandas as pd
configfile: "config/config.yaml"
samples = pd.read_csv("config/samples.tsv", sep="\t").set_index("sample_id", drop=False)
SPECIES   = sorted(samples["species"].unique())
SLIDETAGS = list(samples.index[samples["modality"] == "slidetags"])

def species_samples(w):
    ids = samples.index[samples["species"] == w.species]
    return expand("data/interim/raw/{sample}.h5ad", sample=ids)

for r in ["00a_slidetags","00b_ingest","01_qc","02_integrate","03_cluster",
          "04_annotate","05_de","06_pseudotime","07_spatial_pseudotime",
          "08_cellrank","09_pt_genes_mi","20_liana","21_scenic",
          "22_crossspecies","23_invitro","90_figures"]:
    include: f"rules/{r}.smk"

rule all:
    input:
        expand("data/interim/coords/{s}.csv", s=SLIDETAGS),
        expand("results/tables/markers_{sp}.csv", sp=SPECIES),
        expand("results/tables/mi_{sp}.csv", sp=SPECIES),
        expand("results/tables/liana_{sp}.csv", sp=SPECIES),
        expand("results/tables/scenic_{sp}.csv", sp=SPECIES),
        "results/tables/conservation.csv",
        "data/interim/pt/global_mouse.h5ad",
        "data/interim/pt/cellrank_mouse.h5ad",
        expand("data/interim/pt/spatial_{s}.h5ad", s=SLIDETAGS),
        "results/tables/invitro_comparison.csv",
        "results/figures/.done",
EOF
echo "  wrote: workflow/Snakefile"

echo ""
echo "Done. Next:"
echo "  cd envs/base && uv sync && cd ../..       # build the driver env"
echo "  uv run --project envs/base snakemake -n   # dry-run the full pipeline"
echo "  (add real matrices to data/raw/ per config/samples.tsv; or 'touch' them to test wiring)"
