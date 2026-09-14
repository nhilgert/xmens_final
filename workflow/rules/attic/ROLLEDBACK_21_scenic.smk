# SUPERSEDED 2026-09-13: use 13_scenic.smk. Moved to rules/attic/ (workflow simplification).
# Stage 21 — SCENIC gene regulatory networks / TF activity (per species).
rule scenic:
    input:  "data/interim/annotated_{species}.h5ad"
    output: "results/tables/scenic_{species}.csv"
    log:    "logs/scenic_{species}.log"
    script: "../scripts/scenic.py"     # uses envs/scenic
