"""Stage 01 — explicit, logged cell QC (one species)."""
import logging
import anndata as ad, pandas as pd, scanpy as sc
snk = snakemake  # noqa: F821
q = snk.params.qc; mito_prefix = snk.params.mito_prefix
logging.basicConfig(filename=str(snk.log[0]), level=logging.INFO, format="%(asctime)s  %(message)s")
adata = ad.concat([sc.read_h5ad(f) for f in snk.input], join="outer", merge="same")
n0 = adata.n_obs
adata.var["mito"] = adata.var_names.str.startswith(mito_prefix)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], percent_top=None, inplace=True)
keep = ((adata.obs["n_genes_by_counts"] >= q["min_genes"]) &
        (adata.obs["total_counts"]      >= q["min_counts"]) &
        (adata.obs["pct_counts_mito"]   <= q["max_pct_mito"]))
s = (pd.DataFrame({"sample": adata.obs["sample_id"].values, "kept": keep.values})
       .groupby("sample")["kept"].agg(n_total="size", n_kept="sum"))
s["n_dropped"] = s.n_total - s.n_kept
s.to_csv(snk.output.summary)
logging.info("%s: %d -> %d (%.1f%% dropped)", snk.wildcards.species, n0, int(keep.sum()), 100*(1-keep.mean()))
adata.obs.loc[~keep].to_csv(snk.output.dropped)
adata[keep].copy().write_h5ad(snk.output.kept)
