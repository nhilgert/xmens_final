"""Finish one species' re-embedded extended object as the 2026_10_29b release sources (vocabulary v8).

04g finish_rebuilt_release. reembed_release ran integrate.py (the stage-02 scVI fit and its settings) on
rebuild_release_cells' output: the primary nuclei and the supplementary populations together, so both release
objects share one embedding. This step
  1. removes what integrate.py adds that the releases do not carry: the Leiden ladder and the neighbour-graph
     matrices (the parameters in uns['neighbors_scVI_100'] are kept);
  2. keeps the released .X of the base nuclei (integrate.py recomputes .X as log1p(CP10K) of the counts; the
     recomputation must equal the released values within float32 rounding, X_TOLERANCE), so that the primary
     object's .X is that of 2026_10_29 bit for bit; the supplementary nuclei keep the recomputed values;
  2b. recomputes embedding_outlier on the new embedding over all nuclei with the species' release rule
     (embedding_outliers.py) and writes the flagged nuclei, with their neighbour statistics, to a table;
  3. writes the extended object (all nuclei, obs in_primary) and the primary object (the in_primary nuclei,
     without in_primary, with the category order of the 2026_10_29 object);
  4. checks the primary object against the 2026_10_29 object it was built from: the same nuclei minus the v8
     exclusions, in order; obs values unchanged except embedding_outlier; counts, .X, X_spatial and genes
     unchanged; labels as in the vocabulary; and the extended object restricted to in_primary equals it.
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

snk = snakemake  # noqa: F821 (injected by Snakemake)

CHANGED_OBS = ("embedding_outlier",)
HVG_VAR = ("highly_variable", "highly_variable_rank", "highly_variable_nbatches")
RULE_RETENTION = "Flagged cells are retained in all analyses and omitted from embedding and spatial plots."
X_TOLERANCE = 1e-6        # see rebuild_cells.py: the July mouse .X differs from a recomputation by <= 3.4e-7

LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[finish_rebuild:{snk.params.species}] {msg}", file=LOG, flush=True)


def same_values(a: pd.Series, b: pd.Series) -> bool:
    """Equal element-wise, NaN matching NaN; strings compared as strings."""
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        return bool(np.array_equal(a.to_numpy(float), b.to_numpy(float), equal_nan=True))
    return bool((a.astype(object).where(a.notna(), "<NA>").astype(str).values
                 == b.astype(object).where(b.notna(), "<NA>").astype(str).values).all())


def max_difference(a, b) -> float:
    """Largest absolute difference of two sparse matrices; inf if their non-zero patterns differ."""
    a, b = sp.csr_matrix(a), sp.csr_matrix(b)
    if ((a != 0) != (b != 0)).nnz:
        return float("inf")
    diff = (a.astype(np.float64) - b.astype(np.float64)).tocsr()
    return float(np.abs(diff.data).max()) if diff.nnz else 0.0


def keep_released_x(extended, base):
    """.X of the base nuclei as released; the recomputed .X must agree within X_TOLERANCE."""
    in_primary = extended.obs["in_primary"].to_numpy(bool)
    n_primary = int(in_primary.sum())
    if not in_primary[:n_primary].all():
        raise ValueError("primary nuclei are not first in the re-embedded object")
    rows = base.obs_names.get_indexer(extended.obs_names[:n_primary])
    if (rows < 0).any():
        raise ValueError("primary nuclei missing from the base")
    released = sp.csr_matrix(base.X)[rows]
    recomputed = sp.csr_matrix(extended.X)
    difference = max_difference(recomputed[:n_primary], released)
    if difference > X_TOLERANCE:
        raise ValueError(f"recomputed .X differs from the released .X by {difference:.3g}")
    extended.X = sp.vstack([released, recomputed[n_primary:]], format="csr").astype(recomputed.dtype)
    say(f".X: released values kept for {n_primary:,} base nuclei (recomputation within {difference:.2g}); "
        f"{extended.n_obs - n_primary:,} supplementary nuclei keep the recomputed values")


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


def rule_text(previous: str, n_flagged: int, n_extended_flagged: int, release: str, table_path: str,
              extended: bool) -> str:
    """The release's rule description with this object's count and source."""
    if previous.count(RULE_RETENTION) != 1:
        raise ValueError("uns['embedding_outlier_rule'] does not have the expected wording")
    criterion = previous.split(RULE_RETENTION)[0]
    scope = ("all nuclei of this object" if extended else
             "the nuclei of the extended object (this object plus the supplementary populations)")
    return (f"{criterion}{RULE_RETENTION} {n_flagged} cells flagged in this object ({n_extended_flagged} in the "
            f"extended object), recomputed on the {release} embedding over {scope} by finish_rebuilt_release "
            f"(workflow/scripts/deliverables/embedding_outliers.py); flagged cells and their neighbour "
            f"statistics: {table_path}.")


def check_against_base(primary, base, dropped):
    """The base nuclei minus the v8 exclusions, unchanged except the outlier flags."""
    expected_names = base.obs_names[~base.obs_names.isin(dropped)]
    if not primary.obs_names.equals(expected_names):
        raise ValueError("the primary nuclei are not the base nuclei minus the v8 exclusions, in base order")
    rows = base.obs_names.get_indexer(primary.obs_names)
    base_obs = base.obs.iloc[rows]
    if list(primary.obs.columns) != list(base.obs.columns):
        raise ValueError("primary obs columns differ from the base (names or order)")
    for column in base.obs.columns:
        if column in CHANGED_OBS:
            continue
        if not same_values(primary.obs[column], base_obs[column]):
            raise ValueError(f"obs {column} changed for nuclei of the base")
        if isinstance(base.obs[column].dtype, pd.CategoricalDtype):
            if list(primary.obs[column].cat.categories) != list(base.obs[column].cat.categories):
                raise ValueError(f"obs {column}: categories differ from the base")
    for name, after, before in (("counts", primary.layers["counts"], base.layers["counts"]), (".X", primary.X, base.X)):
        if (sp.csr_matrix(after) != sp.csr_matrix(before)[rows]).nnz:
            raise ValueError(f"{name} changed for nuclei of the base")
    if not np.array_equal(np.asarray(primary.obsm["X_spatial"]), np.asarray(base.obsm["X_spatial"])[rows],
                          equal_nan=True):
        raise ValueError("X_spatial changed for nuclei of the base")
    if not primary.var_names.equals(base.var_names):
        raise ValueError("genes differ from the base")
    for column in base.var.columns:
        if column not in HVG_VAR and not same_values(primary.var[column], base.var[column]):
            raise ValueError(f"var {column} differs from the base")
    extra_obsm = sorted(set(primary.obsm.keys()) - set(base.obsm.keys()))
    if extra_obsm:
        raise ValueError(f"unexpected obsm {extra_obsm}")
    hvg_before = set(base.var_names[base.var["highly_variable"].astype(bool)])
    hvg_after = set(primary.var_names[primary.var["highly_variable"].astype(bool)])
    say(f"HVGs: {len(hvg_after):,} (2026_10_29 fit {len(hvg_before):,}; shared {len(hvg_after & hvg_before):,})")
    say(f"primary = base minus {len(dropped):,} nuclei; unchanged except the recomputed outlier flags")


def check_extended_equals_primary(extended, primary):
    """extended[in_primary] is the primary object (obs except in_primary, X, layers, obsm, var)."""
    sub = extended[extended.obs["in_primary"].to_numpy(bool)]
    if not sub.obs_names.equals(primary.obs_names):
        raise ValueError("extended[in_primary] nuclei differ from the primary object")
    for column in primary.obs.columns:
        if not same_values(sub.obs[column], primary.obs[column]):
            raise ValueError(f"extended[in_primary] obs {column} differs from the primary object")
    for name, a, b in (("counts", sub.layers["counts"], primary.layers["counts"]), (".X", sub.X, primary.X)):
        if (sp.csr_matrix(a) != sp.csr_matrix(b)).nnz:
            raise ValueError(f"extended[in_primary] {name} differs from the primary object")
    if sorted(sub.obsm.keys()) != sorted(primary.obsm.keys()):
        raise ValueError("obsm keys differ between the extended and primary objects")
    for key in primary.obsm.keys():
        if not np.array_equal(np.asarray(sub.obsm[key]), np.asarray(primary.obsm[key]), equal_nan=True):
            raise ValueError(f"extended[in_primary] obsm {key} differs from the primary object")
    if not extended.var.equals(primary.var):
        raise ValueError("var differs between the extended and primary objects")
    say("extended[in_primary] equals the primary object")


def check_counts(obs, vocab, species, key):
    found = {k: int(v) for k, v in obs["cell_type"].astype(str).value_counts().items()}
    expected = {k: int(v) for k, v in vocab[key][species].items()}
    if found != expected or len(obs) != int(vocab[key]["totals"][species]):
        raise ValueError(f"{key}: found {found} ({len(obs)}), vocabulary expects {expected}")
    say(f"  {key}: {found}")


def write(adata, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp.h5ad")
    adata.write_h5ad(tmp_path)
    os.replace(tmp_path, path)
    say(f"wrote {adata.n_obs:,} nuclei -> {path}")


def main():
    species = str(snk.params.species)
    release = str(snk.params.release)
    vocab = yaml.safe_load(open(snk.input.vocab))
    extended = ad.read_h5ad(snk.input.h5ad)
    rebuilt = ad.read_h5ad(snk.input.rebuilt, backed="r")
    base = ad.read_h5ad(snk.input.base)
    say(f"re-embedded {extended.n_obs:,} nuclei; rebuilt {rebuilt.n_obs:,}; base {base.n_obs:,}")
    if not extended.obs_names.equals(rebuilt.obs_names):
        raise ValueError("the re-embedded nuclei are not the rebuilt nuclei in order")
    for column in rebuilt.obs.columns:
        if not same_values(extended.obs[column], rebuilt.obs[column]):
            raise ValueError(f"integrate.py changed obs {column}")
    drop_clustering(extended)
    keep_released_x(extended, base)

    # 2b. Outlier flags over all nuclei of the shared embedding.
    table = outlier_table(extended.obs, extended.obsm, species=species)
    extended.obs["embedding_outlier"] = table["embedding_outlier"].to_numpy(bool)
    flagged = table[table.embedding_outlier].drop(columns="embedding_outlier")
    extra = extended.obs.loc[flagged.index, ["cell_subtype", "sample", "timepoint", "spatial_qc_reason",
                                             "in_primary"]].astype(str)
    umap = pd.DataFrame(np.asarray(extended.obsm["X_umap_scVI"])[extended.obs_names.get_indexer(flagged.index)],
                        index=flagged.index, columns=["umap_1", "umap_2"]).round(3)
    flagged = pd.concat([flagged[["cell_type"]], extra, flagged[["lineage"]], umap,
                         flagged.drop(columns=["cell_type", "lineage"])], axis=1)
    flagged.index.name = "cell"
    Path(snk.output.outliers).parent.mkdir(parents=True, exist_ok=True)
    flagged.round(3).to_csv(snk.output.outliers)
    in_primary = extended.obs["in_primary"].to_numpy(bool)
    n_primary_flagged = int((extended.obs["embedding_outlier"].to_numpy(bool) & in_primary).sum())
    previous = base.obs["embedding_outlier"].astype(bool)
    still = int((previous.reindex(flagged.index).fillna(False)).sum())
    say(f"embedding outliers: {len(flagged)} in the extended object ({n_primary_flagged} primary); "
        f"{still} of the base's {int(previous.sum())} still flagged; by class {flagged['class'].value_counts().to_dict()}; "
        f"by cell type {flagged['cell_type'].value_counts().to_dict()}")

    # 3. The two objects.
    previous_rule = str(base.uns["embedding_outlier_rule"])
    extended.uns["embedding_outlier_rule"] = rule_text(previous_rule, len(flagged), len(flagged), release,
                                                       str(snk.output.outliers), extended=True)
    extended.obs = extended.obs[[*base.obs.columns, "in_primary"]]
    primary = extended[in_primary].copy()
    primary.obs = primary.obs.drop(columns="in_primary")
    for column in base.obs.columns:
        if isinstance(base.obs[column].dtype, pd.CategoricalDtype):
            primary.obs[column] = primary.obs[column].cat.set_categories(base.obs[column].cat.categories)
    primary.uns["embedding_outlier_rule"] = rule_text(previous_rule, n_primary_flagged, len(flagged), release,
                                                      str(snk.output.outliers), extended=False)

    # 4. Checks.
    dropped = base.obs_names.difference(primary.obs_names)
    check_against_base(primary, base, dropped)
    check_extended_equals_primary(extended, primary)
    check_counts(primary.obs, vocab, species, "expected_coarse_counts")
    check_counts(extended.obs, vocab, species, "expected_coarse_counts_extended")
    validate_canonical_annotations(primary, vocab, context=f"{release} {species} all-cells (primary)")
    validate_canonical_annotations(extended, vocab, context=f"{release} {species} all-cells (extended)")
    say(f"uns: {sorted(primary.uns.keys())}")

    write(extended, snk.output.extended)
    write(primary, snk.output.primary)


main()
LOG.close()
