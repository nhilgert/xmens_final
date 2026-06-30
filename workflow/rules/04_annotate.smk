# Stage 04 — marker-based annotation, lineage re-integration, MAD outlier drop.
rule annotate:
    input:  "data/interim/clustered_{species}.h5ad"
    output: "data/interim/annotated_{species}.h5ad"
    log:    "logs/annotate_{species}.log"
    script: "../scripts/annotate.py"
