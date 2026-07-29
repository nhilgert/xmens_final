"""Stage 04b helper - merge lineage sub-objects into a group for joint reintegration.

Concatenates the member lineage h5ads on the cell axis. obs['lineage'] rides
along, so members stay distinguishable and split back after the joint scVI fit.
Genes are identical across members (all are cell-subsets of one annotated
object) and layers['counts'] is preserved for scVI.
"""
import anndata as ad
import scanpy as sc

inp = list(snakemake.input)            # noqa: F821  member h5ads
out = snakemake.output[0]              # noqa: F821
species = snakemake.wildcards.species  # noqa: F821
group = snakemake.wildcards.group      # noqa: F821
logfile = snakemake.log[0]             # noqa: F821

with open(logfile, "w") as log:
    def say(m):
        print(f"[merge:{species}:{group}] {m}", file=log, flush=True)

    parts = []
    for p in inp:
        a = sc.read_h5ad(p)
        say(f"{p}: {a.n_obs} cells x {a.n_vars} genes | layers={list(a.layers)}")
        parts.append(a)

    merged = ad.concat(parts, axis=0, join="outer", merge="first",
                       uns_merge="first", index_unique=None)

    if "counts" not in merged.layers:
        raise KeyError(
            "merged object lost layers['counts']; scVI needs raw counts "
            f"(members had: {[list(a.layers) for a in parts]})"
        )
    if "lineage" in merged.obs:
        say(f"lineage split: {merged.obs['lineage'].value_counts().to_dict()}")
    else:
        say("WARNING no obs['lineage'] column -- members can't be split back")

    say(f"merged -> {merged.n_obs} cells x {merged.n_vars} genes | layers={list(merged.layers)}")
    merged.write_h5ad(out)
    say(f"wrote -> {out}")
