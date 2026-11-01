"""Reuse of the previous release's whole-tissue embedding when a round-2 release keeps exactly its nuclei.

predicted_release_nuclei and same_nuclei are imported at parse time by workflow/rules/04g2_release_round2.smk. They
read obs only (h5py, read-only): the base release's all-cells object, the round-2 lineage fits and the previous
release's all-cells object; and the round-2 decision tables. rebuild_release_round2 keeps the base nuclei outside the
round-2 fits plus the round-2 nuclei whose cluster is not EXCLUDE (apply_annotation.py, fine level), so the prediction
is that set, and the GPU fit is scheduled only for a species whose nuclei change.

reuse_embedding is the copy itself, called by finish_round2.py (04g2 finish_release_round2): it fails unless the
rebuilt nuclei are exactly the previous release's, so a wrong prediction stops the cut.
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "atlas"))
from apply_annotation import load_choices  # noqa: E402


def read_obs_names(path) -> pd.Index:
    """obs_names of an h5ad file, read without loading the object."""
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        return pd.Index(np.asarray(read_elem(obs[obs.attrs["_index"]])).astype(str))


def read_obs_column(path, column: str) -> pd.Series:
    """One obs column of an h5ad file, as strings indexed by obs_names."""
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        if column not in obs:
            raise KeyError(f"{path}: no obs column {column}")
        names = np.asarray(read_elem(obs[obs.attrs["_index"]])).astype(str)
        values = pd.Series(read_elem(obs[column])).astype(str).to_numpy()
    return pd.Series(values, index=pd.Index(names))


def round2_lineage_nuclei(fit_path, table_path, resolution: str) -> tuple:
    """All nuclei of a round-2 lineage fit, and those its decision table keeps (cluster label is not EXCLUDE)."""
    clusters = read_obs_column(fit_path, resolution)
    labels = load_choices(table_path)
    unlabeled = sorted(set(clusters) - set(labels))
    if unlabeled:
        raise ValueError(f"{table_path}: clusters of {fit_path} {resolution} without a row: {unlabeled}")
    is_kept = clusters.map(labels).ne("EXCLUDE").to_numpy()
    return clusters.index, clusters.index[is_kept]


def predicted_release_nuclei(base_path, lineages) -> pd.Index:
    """The nuclei rebuild_release_round2 keeps for one species.

    lineages: one (fit_path, table_path, resolution) per round-2 lineage of the species.
    """
    base = read_obs_names(base_path)
    fit_nuclei, kept_nuclei = [], []
    for fit_path, table_path, resolution in lineages:
        in_fit, kept = round2_lineage_nuclei(fit_path, table_path, resolution)
        fit_nuclei.append(in_fit)
        kept_nuclei.append(kept)
    in_fits = pd.Index(np.concatenate([index.to_numpy() for index in fit_nuclei]))
    kept = pd.Index(np.concatenate([index.to_numpy() for index in kept_nuclei]))
    if not in_fits.is_unique:
        raise ValueError("a nucleus is in two round-2 lineage fits")

    # Kept round-2 nuclei are fit nuclei, so they never repeat a base nucleus outside the fits.
    outside = base[~base.isin(in_fits)]
    return outside.append(kept)


def same_nuclei(predicted: pd.Index, previous_path) -> bool:
    """True when the predicted nuclei are exactly the nuclei of the previous release's object."""
    previous = read_obs_names(previous_path)
    return len(previous) == len(predicted) and set(previous) == set(predicted)


EMBEDDING_OBSM = ("X_scVI", "X_umap_scVI")
EMBEDDING_UNS = ("neighbors_scVI_100", "umap", "hvg", "log1p", "scvi_fit")   # the records rebuild_round2.py removes


def reuse_embedding(rebuilt, previous) -> tuple:
    """The rebuilt object with the previous release's whole-tissue embedding, its records and its outlier flags.

    Returns the new object (rebuilt nuclei and order; obsm, uns and the var columns integrate.py adds copied from
    `previous` by barcode), the copied embedding_outlier flags in that order, and a summary of what was copied.
    """
    missing = rebuilt.obs_names.difference(previous.obs_names)
    extra = previous.obs_names.difference(rebuilt.obs_names)
    if len(missing) or len(extra) or rebuilt.n_obs != previous.n_obs:
        raise ValueError(f"embedding reuse: {len(missing)} rebuilt nuclei are not in the previous release and "
                         f"{len(extra)} of its nuclei are not rebuilt; 04g2 predicted the same nuclei")
    if not rebuilt.var_names.equals(previous.var_names):
        raise ValueError("embedding reuse: genes differ from the previous release")

    adata = rebuilt.copy()
    rows = previous.obs_names.get_indexer(adata.obs_names)
    for key in EMBEDDING_OBSM:
        adata.obsm[key] = np.asarray(previous.obsm[key])[rows]
    copied_uns = [key for key in EMBEDDING_UNS if key in previous.uns]
    for key in copied_uns:
        adata.uns[key] = previous.uns[key]
    added_var = [column for column in previous.var.columns if column not in adata.var.columns]
    for column in added_var:
        adata.var[column] = previous.var[column].to_numpy()
    flags = previous.obs["embedding_outlier"].astype(bool).to_numpy()[rows]
    summary = {"obsm": list(EMBEDDING_OBSM), "uns": copied_uns, "var": added_var, "n_flagged": int(flags.sum())}
    return adata, flags, summary


# The scVI settings that uns['scvi_fit'] records (integrate.py), and where each sits in the fit settings (params.fit).
FIT_RECORD = {
    "seed": ("seed",),
    "n_hvg": ("scvi", "n_hvg"),
    "n_latent": ("scvi", "n_latent"),
    "n_layers": ("scvi", "n_layers"),
}


def fit_settings_differ(record: dict, fit: dict) -> dict:
    """The settings of a recorded scVI fit (uns scvi_fit) that differ from the current fit settings."""
    differ = {}
    for key, path in FIT_RECORD.items():
        current = fit
        for part in path:
            current = current[part]
        if key not in record or int(record[key]) != int(current):
            differ[key] = (record.get(key), current)
    return differ
