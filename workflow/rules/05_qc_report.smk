# Step 05 - QC diagnostic figures (decoupled from the filtering rule 04_qc).
#
# Reads the PRE-QC raw object (to show full distributions with threshold lines)
# and the POST-QC object (to report survival), and writes per-sample figures to
# results/qc/figures/{sample}/. Mirrors the project QC notebook's plots:
#   - total_counts, n_genes_by_counts histograms (log10) with MAD bounds
#   - pct_counts_mt histogram with the fixed 8% bound
#   - doublet score distribution (if doublet calls present in the QC'd object)
#
# Filtering lives in 04_qc; this rule only visualizes. Runs in the same env as
# the spine (scanpy + matplotlib); does NOT need DoubletDetection.

# Defined only when config ingest.rerun is true (default false), like scvi.retrain: the ingest and QC stages are the
# provenance of the pipeline's inputs and are kept to document them; by default nothing rebuilds those inputs.
if config.get("ingest", {}).get("rerun", False):

    rule qc_report:
        input:
            raw="data/interim/raw/{sample}.h5ad",
            qcd="data/interim/qc/{sample}.h5ad",
        output:
            # a sentinel so Snakemake tracks completion; figures go alongside it
            done=touch("results/qc/figures/{sample}/.done"),
        params:
            species=lambda w: samples.loc[w.sample, "species"],
            outdir=lambda w: f"results/qc/figures/{w.sample}",
        log:
            "logs/01b_qc_report_{sample}.log",
        script:
            "../scripts/ingest_qc/qc_report.py"
