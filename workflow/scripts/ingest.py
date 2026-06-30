"""Stage 00b — load one matrix into AnnData, tag with sample metadata."""
import scanpy as sc
snk = snakemake  # noqa: F821
path = str(snk.input[0]); meta = snk.params.meta
adata = sc.read_10x_h5(path) if path.endswith(".h5") else sc.read_h5ad(path)
adata.var_names_make_unique()
for k, v in meta.items(): adata.obs[k] = v
adata.obs_names = [f"{meta['sample_id']}_{bc}" for bc in adata.obs_names]
adata.write_h5ad(snk.output[0])
