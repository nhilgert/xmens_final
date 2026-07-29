# Clustree -- OFF the critical path, like 99_reports.
#
#     snakemake clustree                      # every object in CLUSTREE_TARGETS
#     snakemake results/figures/clustree/mouse_stromal.png
#
# Nothing in rule all depends on these: a figure is terminal, so a plotting
# failure (or a missing R install) can never stall real computation. But they
# ARE in the DAG, so they regenerate when the underlying object changes -- the
# figure can't silently go stale relative to the data.
#
# Requires R with clustree + readr + ggplot2 on PATH (module load / renv).

CLUSTREE = config.get("clustree", {})
CLUSTREE_TARGETS = CLUSTREE.get("targets", ["mouse_stromal"])
CLUSTREE_RANGE = CLUSTREE.get("range", {"start": 0.1, "stop": 0.8, "step": 0.1})


rule clustree_export:
    input:
        h5ad="data/interim/latent_{species}_{lineage}.h5ad",
    output:
        csv="results/tables/clustree/{species}_{lineage}.csv",
    params:
        start=CLUSTREE_RANGE["start"],
        stop=CLUSTREE_RANGE["stop"],
        step=CLUSTREE_RANGE["step"],
    log:
        "logs/98_clustree_export_{species}_{lineage}.log",
    script:
        "../scripts/export_clustree.py"


rule clustree_plot:
    input:
        csv="results/tables/clustree/{species}_{lineage}.csv",
    output:
        png="results/figures/clustree/{species}_{lineage}.png",
    params:
        title=lambda w: f"{w.species} {w.lineage} -- Leiden ladder",
    log:
        "logs/98_clustree_plot_{species}_{lineage}.log",
    script:
        "../scripts/clustree_plot.R"


rule clustree:
    input:
        expand("results/figures/clustree/{target}.png", target=CLUSTREE_TARGETS),
