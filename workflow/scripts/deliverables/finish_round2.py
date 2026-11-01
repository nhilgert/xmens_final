"""Finish one species' rebuilt all-cells object as a source of the release being cut.

04g2 finish_release_round2. The whole-tissue embedding comes one of two ways (params.embedding, decided in 04g2):
  refit  reembed_release_round2 ran integrate.py (the stage-02 scVI fit and its settings) on rebuild_release_round2's
         output. This step removes what integrate.py adds that the releases do not carry (the Leiden ladder and the
         neighbour-graph matrices; the parameters in uns['neighbors_scVI_100'] are kept), keeps the released .X of the
         nuclei that were in the base (the recomputation must equal it within float32 rounding, X_TOLERANCE), and
         recomputes embedding_outlier on the new embedding over all nuclei with the species' release rule
         (embedding_outliers.py).
  reuse  the nuclei are exactly those of the previous release (params.previous_release), so its whole-tissue
         embedding (obsm X_scVI and X_umap_scVI, the HVG, neighbour, UMAP and scVI-fit records, the var columns
         integrate.py adds) and its embedding_outlier flags are copied nucleus by nucleus onto the rebuilt object. The
         step fails unless the rebuilt nuclei are the previous release's; the rule recomputed on the reused embedding
         is written next to the copied flags.
Either way the flagged nuclei, with their neighbour statistics, go to a table, and the object is checked: the rebuilt
nuclei in order with their obs values (embedding_outlier set here), counts, .X, X_spatial and genes; the base nuclei
outside the round-2 lineages unchanged except embedding_outlier and the spatial QC that rebuild_release_round2
reported; counts and labels as in the vocabulary.
"""
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "atlas"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from apply_annotation import validate_canonical_annotations  # noqa: E402
from embedding_outliers import outlier_table  # noqa: E402
from embedding_reuse import fit_settings_differ, reuse_embedding  # noqa: E402

snk = snakemake  # noqa: F821 (injected by Snakemake)

SPATIAL_OBS = ("spatial_qc_reason", "x_um_qc", "y_um_qc", "tissue_section")
HVG_VAR = ("highly_variable", "highly_variable_rank", "highly_variable_nbatches")
RULE_RETENTION = "Flagged cells are retained in all analyses and omitted from embedding and spatial plots."
X_TOLERANCE = 1e-6        # see rebuild_cells.py: the July mouse .X differs from a recomputation by <= 3.4e-7

LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[finish_round2:{snk.params.species}] {msg}", file=LOG, flush=True)


def same_values(a: pd.Series, b: pd.Series) -> np.ndarray:
    """Element-wise equality, NaN matching NaN; strings compared as strings."""
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        x, y = a.to_numpy(float), b.to_numpy(float)
        return (x == y) | (np.isnan(x) & np.isnan(y))
    return (a.astype(object).where(a.notna(), "<NA>").astype(str).values
            == b.astype(object).where(b.notna(), "<NA>").astype(str).values)


def max_difference(a, b) -> float:
    a, b = sp.csr_matrix(a), sp.csr_matrix(b)
    if ((a != 0) != (b != 0)).nnz:
        return float("inf")
    diff = (a.astype(np.float64) - b.astype(np.float64)).tocsr()
    return float(np.abs(diff.data).max()) if diff.nnz else 0.0


def drop_clustering(adata):
    """The Leiden ladder and graph matrices integrate.py adds; releases carry neither."""
    ladder = [c for c in adata.obs.columns if str(c).startswith("leiden_scVI_")]
    adata.obs = adata.obs.drop(columns=ladder)
    for key in [k for k in adata.uns if str(k).startswith("leiden_scVI_")]:
        del adata.uns[key]
    graphs = list(adata.obsp.keys())
    for key in graphs:
        del adata.obsp[key]
    say(f"dropped {len(ladder)} Leiden columns and obsp {graphs}")


def keep_released_x(adata, base):
    """.X of the nuclei that were in the base, as released; the recomputation must agree within X_TOLERANCE."""
    in_base = adata.obs_names.isin(base.obs_names)
    n_base = int(in_base.sum())
    if not in_base[:n_base].all():
        raise ValueError("the base nuclei are not first in the re-embedded object")
    released = sp.csr_matrix(base.X)[base.obs_names.get_indexer(adata.obs_names[:n_base])]
    recomputed = sp.csr_matrix(adata.X)
    difference = max_difference(recomputed[:n_base], released)
    if difference > X_TOLERANCE:
        raise ValueError(f"recomputed .X differs from the released .X by {difference:.3g}")
    adata.X = sp.vstack([released, recomputed[n_base:]], format="csr").astype(recomputed.dtype)
    say(f".X: released values kept for {n_base:,} base nuclei (recomputation within {difference:.2g}); "
        f"{adata.n_obs - n_base:,} new nuclei keep the recomputed values")


def main():
    species = str(snk.params.species)
    release = str(snk.params.release)
    embedding = str(snk.params.embedding)
    vocab = yaml.safe_load(open(snk.input.vocab))
    rebuilt = ad.read_h5ad(snk.input.rebuilt)
    base = ad.read_h5ad(snk.input.base)
    say(f"rebuilt {rebuilt.n_obs:,} nuclei; base {base.n_obs:,}; whole-tissue embedding: {embedding}")

    if embedding == "refit":
        adata = ad.read_h5ad(snk.input.reembedded)
        if not adata.obs_names.equals(rebuilt.obs_names):
            raise ValueError("the re-embedded nuclei are not the rebuilt nuclei in order")
        for column in rebuilt.obs.columns:
            if not same_values(adata.obs[column], rebuilt.obs[column]).all():
                raise ValueError(f"integrate.py changed obs {column}")

        # 1-2. Releases carry no clustering; base nuclei keep their released .X.
        drop_clustering(adata)
        keep_released_x(adata, base)

        # 3. Outlier flags over all nuclei of the new embedding.
        table = outlier_table(adata.obs, adata.obsm, species=species)
        source = (f"recomputed on the {release} embedding over all nuclei of this object by finish_release_round2 "
                  "(workflow/scripts/deliverables/embedding_outliers.py)")
    elif embedding == "reuse":
        previous_release = str(snk.params.previous_release)
        previous = ad.read_h5ad(snk.input.previous, backed="r")

        # The reused fit must be the one the current settings would make (params.fit, 02_integrate.smk).
        differ = fit_settings_differ(dict(previous.uns.get("scvi_fit", {})), dict(snk.params.fit))
        if differ:
            raise ValueError(f"embedding reuse: {previous_release} was fitted with other scVI settings (recorded, "
                             f"current): {differ}; set config release_round2.embedding: refit")
        adata, flags, copied = reuse_embedding(rebuilt, previous)
        previous.file.close()
        say(f"reused from {previous_release}: obsm {copied['obsm']}, uns {copied['uns']}, var {copied['var']}, "
            f"embedding_outlier ({copied['n_flagged']} flagged)")

        # The copied flags are the release's; the rule recomputed with the current cell types is reported beside them.
        table = outlier_table(adata.obs, adata.obsm, species=species)
        table["recomputed_flag"] = table["embedding_outlier"].to_numpy(bool)
        table["embedding_outlier"] = flags
        n_differ = int((table["recomputed_flag"] != table["embedding_outlier"]).sum())
        say(f"embedding_outlier: {int(flags.sum())} flags of {previous_release} copied; the rule recomputed with the "
            f"current cell types differs for {n_differ} nuclei (column recomputed_flag of {snk.output.outliers})")
        source = (f"the flags of release {previous_release}, whose whole-tissue embedding this release reuses because "
                  "its nuclei are the same (finish_release_round2)")
    else:
        raise ValueError(f"params.embedding is refit or reuse, not {embedding!r}")
    adata.obs["embedding_outlier"] = table["embedding_outlier"].to_numpy(bool)
    flagged = table[table.embedding_outlier].drop(columns="embedding_outlier")
    extra = adata.obs.loc[flagged.index, ["cell_subtype", "sample", "timepoint", "spatial_qc_reason"]].astype(str)
    umap = pd.DataFrame(np.asarray(adata.obsm["X_umap_scVI"])[adata.obs_names.get_indexer(flagged.index)],
                        index=flagged.index, columns=["umap_1", "umap_2"]).round(3)
    flagged = pd.concat([flagged[["cell_type"]], extra, flagged[["lineage"]], umap,
                         flagged.drop(columns=["cell_type", "lineage"])], axis=1)
    flagged.index.name = "cell"
    Path(snk.output.outliers).parent.mkdir(parents=True, exist_ok=True)
    flagged.round(3).to_csv(snk.output.outliers)
    base_flags = base.obs["embedding_outlier"].astype(bool)
    still = int(base_flags.reindex(flagged.index).fillna(False).astype(bool).sum())
    say(f"embedding outliers: {len(flagged)} ({still} of the base's {int(base_flags.sum())} still flagged); by class "
        f"{flagged['class'].value_counts().to_dict()}; by cell type {flagged['cell_type'].value_counts().to_dict()}")
    rule = str(rebuilt.uns["embedding_outlier_rule"])
    if rule.count(RULE_RETENTION) != 1:
        raise ValueError("uns['embedding_outlier_rule'] does not have the expected wording")
    adata.uns["embedding_outlier_rule"] = (
        f"{rule.split(RULE_RETENTION)[0]}{RULE_RETENTION} {len(flagged)} cells flagged in this object, {source}; "
        f"flagged cells and their neighbour statistics: {snk.output.outliers}.")
    adata.obs = adata.obs[list(base.obs.columns)]

    # 4. Checks.
    for column in rebuilt.obs.columns:
        if not same_values(adata.obs[column], rebuilt.obs[column]).all():
            raise ValueError(f"obs {column} differs from the rebuilt object")
    for name, after, before in (("counts", adata.layers["counts"], rebuilt.layers["counts"]),
                                (".X", adata.X, rebuilt.X)):
        if (sp.csr_matrix(after) != sp.csr_matrix(before)).nnz:
            raise ValueError(f"{name} differs from the rebuilt object")
    if not np.array_equal(np.asarray(adata.obsm["X_spatial"]), np.asarray(rebuilt.obsm["X_spatial"]), equal_nan=True):
        raise ValueError("X_spatial differs from the rebuilt object")
    if not adata.var_names.equals(base.var_names):
        raise ValueError("genes differ from the base")
    if sorted(adata.obsm.keys()) != sorted(base.obsm.keys()):
        raise ValueError(f"obsm {sorted(adata.obsm.keys())} != base {sorted(base.obsm.keys())}")
    lineages = list(rebuilt.uns["release_round2"]["lineages"])
    outside = base.obs_names[~base.obs["origin_lineage"].astype(str).isin(lineages).to_numpy()]
    if not outside.isin(adata.obs_names).all():
        raise ValueError("base nuclei outside the round-2 lineages are missing")
    rows = adata.obs_names.get_indexer(outside)
    n_changed_qc = 0
    for column in base.obs.columns:
        if column == "embedding_outlier":
            continue
        equal = same_values(adata.obs[column].iloc[rows], base.obs.loc[outside, column])
        if column in SPATIAL_OBS:
            n_changed_qc = max(n_changed_qc, int((~equal).sum()))
        elif not equal.all():
            raise ValueError(f"obs {column} changed for {int((~equal).sum())} base nuclei outside the round-2 lineages")
    reported = int(rebuilt.uns["release_round2"]["spatial_qc"]["n_base_changed"])
    if n_changed_qc > reported:
        raise ValueError(f"{n_changed_qc} base nuclei outside the round-2 lineages changed spatial QC; rebuild "
                         f"reported {reported} base nuclei in all")
    for name, after, before in (("counts", adata.layers["counts"], base.layers["counts"]), (".X", adata.X, base.X)):
        if (sp.csr_matrix(after)[rows] != sp.csr_matrix(before)[base.obs_names.get_indexer(outside)]).nnz:
            raise ValueError(f"{name} changed for base nuclei outside the round-2 lineages")
    found = {k: int(v) for k, v in adata.obs["cell_type"].astype(str).value_counts().items()}
    expected = {k: int(v) for k, v in vocab["expected_coarse_counts"][species].items()}
    if found != expected or adata.n_obs != int(vocab["expected_coarse_counts"]["totals"][species]):
        raise ValueError(f"expected_coarse_counts: found {found} ({adata.n_obs}), vocabulary expects {expected}")
    validate_canonical_annotations(adata, vocab, context=f"{release} {species} all-cells")
    hvg_before = set(base.var_names[base.var["highly_variable"].astype(bool)])
    hvg_after = set(adata.var_names[adata.var["highly_variable"].astype(bool)])
    say(f"checks passed: {adata.n_obs:,} nuclei; {len(outside):,} base nuclei outside the round-2 lineages unchanged "
        f"except outlier flags and spatial QC (<= {n_changed_qc} changed); HVGs {len(hvg_after):,} "
        f"({snk.params.base_release} fit {len(hvg_before):,}; shared {len(hvg_after & hvg_before):,}); "
        f"cell types {found}")

    path = Path(snk.output.primary)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp.h5ad")
    adata.write_h5ad(tmp)
    os.replace(tmp, path)
    say(f"wrote {adata.n_obs:,} nuclei -> {path}")


main()
LOG.close()
