# Stage 05 — marker DE (Wilcoxon).
rule de:
    input:  "data/interim/annotated_{species}.h5ad"
    output: "results/tables/markers_{species}.csv"
    log:    "logs/de_{species}.log"
    script: "../scripts/de.py"
