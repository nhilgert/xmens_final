# Stage 02 — scVI integration (per species).
rule integrate:
    input:  "data/interim/qc_passed_{species}.h5ad"
    output: "data/interim/latent_{species}.h5ad"
    params: scvi=config["scvi"]
    log:    "logs/integrate_{species}.log"
    script: "../scripts/integrate.py"
