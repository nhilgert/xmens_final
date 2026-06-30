rule ingest:
    input:  lambda w: samples.loc[w.sample, "path"]
    output: "data/interim/raw/{sample}.h5ad"
    params: meta=lambda w: samples.loc[w.sample].to_dict()
    log:    "logs/ingest_{sample}.log"
    script: "../scripts/ingest.py"
