"""The embedding-outlier rule of the releases (uns['embedding_outlier_rule']).

A nucleus is an embedding outlier when, among its 30 nearest neighbours (itself excluded) in the
UMAP (obsm X_umap_scVI) or in the scVI latent (obsm X_scVI), >= 70% belong to one other lineage
and <= 20% to its own. Lineages group the coarse cell types, as recorded in each species' objects:
  human: stromal = endometrial_fibroblast + decidual_cell + smooth_muscle_cell + perivascular_cell;
         epithelial, immune and endothelial are their own cell types;
  mouse: stromal = endometrial_fibroblast + decidual_cell; every other cell type is its own lineage
         (the rule of the 2026_10_13 mouse flags).

The flags of the 2026_10_22 release came from remove_barcodes_human_2026_10_14.csv, whose columns
this module reproduces (same_lin_frac_*, top_other_lin_*, top_other_frac_*, misplaced_*, class).
Flagged nuclei stay in every analysis; plots of embeddings and tissue omit them.
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

N_NEIGHBOURS = 30
MAX_OWN_FRACTION = 0.2
MIN_OTHER_FRACTION = 0.7
LINEAGE_OF_CELL_TYPE = {
    "endometrial_fibroblast": "stromal",
    "decidual_cell": "stromal",
    "smooth_muscle_cell": "stromal",
    "perivascular_cell": "stromal",
    "epithelial_cell": "epithelial",
    "immune_cell": "immune",
    "endothelial_cell": "endothelial",
}
MOUSE_STROMAL = ("endometrial_fibroblast", "decidual_cell")
EMBEDDINGS = {"umap": "X_umap_scVI", "scvi": "X_scVI"}


def lineage_of(cell_type: pd.Series, species: str = "human") -> pd.Series:
    """Lineage of each nucleus under the species' rule (see the module docstring)."""
    if species == "human":
        unknown = sorted(set(cell_type) - set(LINEAGE_OF_CELL_TYPE))
        if unknown:
            raise ValueError(f"cell types without a lineage in the outlier rule: {unknown}")
        return cell_type.map(LINEAGE_OF_CELL_TYPE)
    if species == "mouse":
        return cell_type.where(~cell_type.isin(MOUSE_STROMAL), "stromal")
    raise ValueError(f"no outlier rule for species {species!r}")


def neighbours_without_self(coords: np.ndarray, k: int) -> np.ndarray:
    """Indices of each row's k nearest other rows (Euclidean).

    Queries k + 1 and removes the row itself; when an identical point displaces the row from
    its own neighbour list, the farthest of the k + 1 is dropped instead.
    """
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(coords)
    rows = np.arange(len(coords))[:, None]
    is_self = idx == rows
    no_self = ~is_self.any(axis=1)
    is_self[no_self, k] = True
    return idx[~is_self].reshape(len(coords), k)


def lineage_fractions(lineage: pd.Series, coords: np.ndarray, k: int = N_NEIGHBOURS):
    """Own-lineage fraction, top other lineage and its fraction among each nucleus' neighbours."""
    names = np.array(sorted(set(lineage)))
    codes = pd.Categorical(lineage, categories=names).codes
    neighbour_codes = codes[neighbours_without_self(np.asarray(coords, dtype=np.float64), k)]
    counts = np.stack([(neighbour_codes == j).sum(axis=1) for j in range(len(names))], axis=1)

    own = counts[np.arange(len(codes)), codes]
    others = counts.copy()
    others[np.arange(len(codes)), codes] = -1
    top_other = others.argmax(axis=1)
    return own, names[top_other], others.max(axis=1)


def outlier_table(obs: pd.DataFrame, obsm, k: int = N_NEIGHBOURS, species: str = "human") -> pd.DataFrame:
    """Per-nucleus neighbour statistics and flags, columns as in the 2026_10_14 barcode lists."""
    cell_type = obs["cell_type"].astype(str)
    lineage = lineage_of(cell_type, species)

    max_own = int(np.floor(MAX_OWN_FRACTION * k + 1e-9))        # 6 of 30
    min_other = int(np.ceil(MIN_OTHER_FRACTION * k - 1e-9))     # 21 of 30
    table = pd.DataFrame({"cell_type": cell_type.values, "lineage": lineage.values}, index=obs.index)
    for short, key in EMBEDDINGS.items():
        own, top_lineage, top_count = lineage_fractions(lineage, obsm[key], k)
        table[f"same_lin_frac_{short}"] = own / k
        table[f"top_other_lin_{short}"] = top_lineage
        table[f"top_other_frac_{short}"] = top_count / k
        table[f"misplaced_{short}"] = (own <= max_own) & (top_count >= min_other)
    table["class"] = np.select(
        [table.misplaced_umap & table.misplaced_scvi, table.misplaced_umap, table.misplaced_scvi],
        ["both", "umap_only", "scvi_only"], default="")
    table["embedding_outlier"] = table.misplaced_umap | table.misplaced_scvi
    return table
