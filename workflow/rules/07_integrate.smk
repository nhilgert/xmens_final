# Step 07 - scVI integration + Leiden clustering, per species.
#
# Ported from integrate_and_cluster.ipynb (faithful to Cevrim et al. Methods):
#   - HVG selection (seurat_v3, 3000, batch_key=sample) on mito/ribo-EXCLUDED
#     genes, computed on the raw-counts layer
#   - scVI (n_layers=2, n_latent=30, gene_likelihood=nb) trained on the raw
#     counts of the HVG subset, with sample as batch key; early stopping
#   - neighbors (k=100) on the scVI latent
#   - Leiden at many resolutions (igraph, seed 0)
#   - UMAP on the scVI graph
#
# GPU rule: workflow/profiles/slurm routes `integrate` to a GPU node (--gres=gpu:1).
# Input is the per-species concatenated object from 06_concat.

# The settings integrate.py and reintegrate.py fit with: config seed, scvi (without the retrain switch) and leiden.
# Every rule that runs either script declares them as params, so a changed setting re-runs the fit (rerun trigger
# params) instead of being read silently from the config.
SCVI_FIT_SETTINGS = {
    "seed": config["seed"],
    "scvi": {key: value for key, value in config["scvi"].items() if key != "retrain"},
    "leiden": config["leiden"],
}

# Defined only when config scvi.retrain is true. scVI fits are deterministic only
# within one GPU class and the annotation tables refer to the July 2026 fits, so by
# default data/interim/latent_{species}.h5ad is an input the pipeline reads,
# never an output it rewrites.
if config["scvi"].get("retrain", False):

    rule integrate:
        input:
            h5ad="data/interim/qc_passed_{species}.h5ad",
        output:
            h5ad="data/interim/latent_{species}.h5ad",
        params:
            counts_layer="counts",          # this pipeline's raw-counts layer (= notebook 'raw')
            fit=SCVI_FIT_SETTINGS,
        log:
            "logs/02_integrate_{species}.log",
        script:
            "../scripts/atlas/integrate.py"
