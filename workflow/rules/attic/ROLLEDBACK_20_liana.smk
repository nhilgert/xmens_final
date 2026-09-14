# SUPERSEDED 2026-09-13: use 14_liana.smk (+16_liana_crossspecies, 17_liana_pseudotime). Moved to rules/attic/ (workflow simplification).
# Stage 20 — ligand-receptor / cell-cell communication (per species).
rule liana:
    input:  "data/interim/annotated_{species}.h5ad"
    output: "results/tables/liana_{species}.csv"
    log:    "logs/liana_{species}.log"
    script: "../scripts/liana.py"      # uses envs/liana
