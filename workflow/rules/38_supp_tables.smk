# Step 38 - supplementary tables built from the pipeline's own records.
#
# One rule per table. Each writes into results/supp_tables/{version}/ so a new release gets new files and no
# earlier supplement is overwritten, and each writes a provenance JSON next to its CSV naming the source files and
# their sha256. Every rule checks its own output against the numbers written in the Methods (config supp_tables.*
# expected / methods_totals) and stops rather than shipping a table that disagrees with the paper.
#
# S5  per-sample QC cutoffs and counts  <- the per-sample QC logs of record (logs/01_qc_{sample}.log)
# S13 clusters removed during annotation <- the removed-cluster ledger (rule removed_cluster_ledger)
#
# S6 (over-clustering) is owned elsewhere and is expected to land in this file as rule supp_table_s6.
#
# Note on S5's source: qc.py prints every cutoff and per-stage count but writes only the filtered object, so the
# logs are the only record of the applied numbers. Re-running QC to regenerate them would rewrite the inputs and
# re-run DoubletDetection, which is not reproducible. When QC is next re-run, qc.py should write the same numbers
# as a JSON beside each object (docs/PORTING.md) and this rule should read those instead.

SUPP = config.get("supp_tables", {})
SUPP_DIR = f"results/supp_tables/{RELEASE_VERSION}"


rule supp_table_s5:
    input:
        logs=expand("logs/01_qc_{sample}.log", sample=list(samples.index)),
    output:
        csv=f"{SUPP_DIR}/Table_S5_qc_per_sample.csv",
        provenance=f"{SUPP_DIR}/Table_S5_qc_per_sample.provenance.json",
    params:
        methods_totals=SUPP.get("s5", {}).get("methods_totals", {}),
    log:
        "logs/supp_tables/s5.log",
    script:
        "../scripts/supp_tables/qc_per_sample.py"


rule supp_table_s13:
    input:
        ledger=f"results/annotation_review/removed_cluster_ledger_{RELEASE_VERSION}.csv",
    output:
        csv=f"{SUPP_DIR}/Table_S13_removed_clusters.csv",
        provenance=f"{SUPP_DIR}/Table_S13_removed_clusters.provenance.json",
    params:
        spec=SUPP.get("s13", {}),
    log:
        "logs/supp_tables/s13.log",
    script:
        "../scripts/supp_tables/removed_clusters.py"


# Target: every supplementary table this pipeline builds. Add new tables here as their rules are written.
rule supp_tables:
    input:
        rules.supp_table_s5.output,
        rules.supp_table_s13.output,
