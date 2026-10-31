"""Finish the re-embedded human all-cells object as a release source.

04g finish_reembedded_release. reembed_release ran integrate.py (the stage-02 scVI fit and its
settings) on the output of reinstate_release_cells. This step
  1. removes what integrate.py adds that the release does not carry: the Leiden ladder and the
     neighbour-graph matrices (the release keeps the uns['neighbors_scVI_100'] parameters only);
  2. recomputes embedding_outlier on the new embedding with the release rule
     (embedding_outliers.py, which reproduces the 25 flags of the 2026_10_22 release on the July
     embedding) and writes the flagged nuclei, with their neighbour statistics, to a table;
  3. checks the result against the migrated release it was built from: the release's nuclei
     come first and in order; their obs values are unchanged except cell_subtype of the
     coarse-filled epithelium (now inner/outer) and embedding_outlier; counts, .X, coordinates,
     genes and category orders are unchanged; and the labels match the vocabulary counts.
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

EPITHELIAL = "epithelial_cell"
CHANGED_OBS = ("cell_subtype", "embedding_outlier")
HVG_VAR = ("highly_variable", "highly_variable_rank", "highly_variable_nbatches")
RULE_RETENTION = "Flagged cells are retained in all analyses and omitted from embedding and spatial plots."

LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[finish_reembed] {msg}", file=LOG, flush=True)


def same_values(a: pd.Series, b: pd.Series) -> bool:
    """Equal element-wise, NaN matching NaN; strings compared as strings."""
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        return bool(np.array_equal(a.to_numpy(float), b.to_numpy(float), equal_nan=True))
    return bool((a.astype(object).where(a.notna(), "<NA>").astype(str).values
                 == b.astype(object).where(b.notna(), "<NA>").astype(str).values).all())


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


def rule_text(previous: str, n_flagged: int, release: str, table_path: str) -> str:
    """The release's rule description with this release's count and source."""
    if previous.count(RULE_RETENTION) != 1:
        raise ValueError("uns['embedding_outlier_rule'] does not have the expected wording")
    criterion = previous.split(RULE_RETENTION)[0]
    return (f"{criterion}{RULE_RETENTION} {n_flagged} cells flagged in this object, recomputed on the "
            f"{release} embedding by finish_reembedded_release (workflow/scripts/deliverables/"
            f"embedding_outliers.py); flagged cells and their neighbour statistics: {table_path}.")


def check_against_release(out, base):
    """The release's nuclei are unchanged except the epithelial split and the outlier flags."""
    n = base.n_obs
    if not out.obs_names[:n].equals(base.obs_names):
        raise ValueError("the release's nuclei are not first and in release order")
    for column in base.obs.columns:
        if column not in CHANGED_OBS and not same_values(out.obs[column].iloc[:n], base.obs[column]):
            raise ValueError(f"obs {column} changed for nuclei of the release")
        if isinstance(base.obs[column].dtype, pd.CategoricalDtype):
            before = list(base.obs[column].cat.categories)
            if list(out.obs[column].cat.categories)[:len(before)] != before:
                raise ValueError(f"obs {column}: category order changed")

    subtype_before = base.obs["cell_subtype"].astype(str).values
    subtype_after = out.obs["cell_subtype"].astype(str).values[:n]
    changed = subtype_before != subtype_after
    coarse_epithelial = (base.obs["cell_type"].astype(str).values == EPITHELIAL) & (subtype_before == EPITHELIAL)
    if (changed != coarse_epithelial).any():
        raise ValueError("cell_subtype changed outside the coarse-filled epithelium (or not on all of it)")

    for name, after, before in (("counts", out.layers["counts"], base.layers["counts"]), (".X", out.X, base.X)):
        if (sp.csr_matrix(after)[:n] != sp.csr_matrix(before)).nnz:
            raise ValueError(f"{name} changed for nuclei of the release")
    if not np.array_equal(np.asarray(out.obsm["X_spatial"])[:n], np.asarray(base.obsm["X_spatial"]), equal_nan=True):
        raise ValueError("X_spatial changed for nuclei of the release")

    if not out.var_names.equals(base.var_names):
        raise ValueError("genes differ from the release")
    for column in base.var.columns:
        if column not in HVG_VAR and not same_values(out.var[column], base.var[column]):
            raise ValueError(f"var {column} differs from the release")
    hvg_before = set(base.var_names[base.var["highly_variable"].astype(bool)])
    hvg_after = set(out.var_names[out.var["highly_variable"].astype(bool)])
    # The release var flags predate this fit and are not those of its July embedding (1,691 of 3,000 shared).
    say(f"HVGs: {len(hvg_after):,} (previous release var flags {len(hvg_before):,}; shared {len(hvg_after & hvg_before):,})")
    say(f"release nuclei unchanged except cell_subtype of {int(changed.sum()):,} epithelial nuclei "
        f"and the recomputed outlier flags")


def check_counts_against_vocabulary(obs, vocab):
    tasks = vocab["content_tasks"]
    subtypes = obs["cell_subtype"].astype(str)
    expected_split = tasks["human_epithelial_split"]["expected_counts"]
    split = {label: int((subtypes == label).sum()) for label in expected_split}
    coarse = obs["cell_type"].astype(str).value_counts().to_dict()
    fills = {entry["cell_type"]: entry["expected_coarse_fill_count"]
             for entry in tasks["coarse_subtype_fill"]["applies_to"]["human"]}
    found_fills = {cell_type: int(((obs["cell_type"].astype(str) == cell_type) & (subtypes == cell_type)).sum())
                   for cell_type in fills}
    checks = [
        ("human_epithelial_split", split, expected_split),
        ("expected_coarse_counts", coarse, dict(vocab["expected_coarse_counts"]["human"])),
        ("coarse_subtype_fill", found_fills, fills),
        ("total", len(obs), vocab["expected_coarse_counts"]["totals"]["human"]),
    ]
    for name, found, expected in checks:
        if found != expected:
            raise ValueError(f"{name}: found {found}, vocabulary expects {expected}")
        say(f"  {name}: {found}")


def main():
    vocab = yaml.safe_load(open(snk.input.vocab))
    release = str(snk.params.release)
    out = ad.read_h5ad(snk.input.h5ad)
    base = ad.read_h5ad(snk.input.base)
    say(f"re-embedded {out.n_obs:,} nuclei; release {base.n_obs:,}")
    drop_clustering(out)

    table = outlier_table(out.obs, out.obsm)
    out.obs["embedding_outlier"] = table["embedding_outlier"].to_numpy(bool)
    flagged = table[table.embedding_outlier].drop(columns="embedding_outlier")
    extra = out.obs.loc[flagged.index, ["cell_subtype", "sample", "timepoint", "spatial_qc_reason"]].astype(str)
    umap = pd.DataFrame(np.asarray(out.obsm["X_umap_scVI"])[out.obs_names.get_indexer(flagged.index)],
                        index=flagged.index, columns=["umap_1", "umap_2"]).round(3)
    flagged = pd.concat([flagged[["cell_type"]], extra, flagged[["lineage"]], umap,
                         flagged.drop(columns=["cell_type", "lineage"])], axis=1)
    flagged.index.name = "cell"
    Path(snk.output.outliers).parent.mkdir(parents=True, exist_ok=True)
    flagged.round(3).to_csv(snk.output.outliers)
    previous = base.obs["embedding_outlier"].astype(bool)
    kept_flags = int((previous & out.obs["embedding_outlier"].iloc[:base.n_obs].values).sum())
    say(f"embedding outliers: {len(flagged)} ({flagged['class'].value_counts().to_dict()}); "
        f"{kept_flags} of the release's {int(previous.sum())} still flagged; by cell type "
        f"{flagged['cell_type'].value_counts().to_dict()}")
    out.uns["embedding_outlier_rule"] = rule_text(str(base.uns["embedding_outlier_rule"]), len(flagged),
                                                  release, str(snk.output.outliers))

    check_against_release(out, base)
    check_counts_against_vocabulary(out.obs, vocab)
    validate_canonical_annotations(out, vocab, context=f"{release} human all-cells")

    # Same structure as the release: obs columns in release order, then nothing extra.
    out.obs = out.obs[list(base.obs.columns)]
    extra_obsm = sorted(set(out.obsm.keys()) - set(base.obsm.keys()))
    if extra_obsm:
        raise ValueError(f"unexpected obsm {extra_obsm}")
    say(f"uns added: {sorted(set(out.uns) - set(base.uns))}; uns removed: {sorted(set(base.uns) - set(out.uns))}")

    out_path = Path(snk.output.h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp.h5ad")
    out.write_h5ad(tmp_path)
    os.replace(tmp_path, out_path)
    say(f"wrote {out.n_obs:,} nuclei -> {out_path}")


main()
LOG.close()
