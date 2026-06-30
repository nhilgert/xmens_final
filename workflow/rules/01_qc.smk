rule qc:
    input:  species_samples
    output: kept="data/interim/qc_passed_{species}.h5ad",
            dropped="results/qc/dropped_{species}.csv",
            summary="results/qc/qc_summary_{species}.csv"
    params: qc=config["qc"], mito_prefix=lambda w: config["mito_prefix"][w.species]
    log:    "logs/qc_{species}.log"
    script: "../scripts/qc.py"
