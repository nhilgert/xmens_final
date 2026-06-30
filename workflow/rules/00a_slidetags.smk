# Stage 00a — Slide-tags spatial reconstruction: raw sequencing -> per-cell (x,y).
rule slidetags:
    input:  lambda w: samples.loc[w.sample, "path"]
    output: "data/interim/coords/{sample}.csv"
    log:    "logs/slidetags_{sample}.log"
    script: "../scripts/slidetags.py"
