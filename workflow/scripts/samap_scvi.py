"""External-embedding integration for SAMap (scVI, Harmony, or any latent).

Two insertion points ("slots") into SAMap, both applied to the per-species
SAM objects *before* they are passed to ``SAMAP``. No SAMap internals are
modified — everything goes through state SAMap already consumes.

Slot 1 — :func:`use_embedding_manifold`
    Rebuild a SAM object's within-species manifold (kNN graph + gene
    weights) from an arbitrary cell embedding: a batch-conditioned scVI
    latent, a Harmony-corrected PCA, anything (n_cells x d). SAMap's
    neighborhood coarsening, stitching, and homology refinement consume
    only ``obsp['connectivities']`` / ``var['weights']`` / the Leiden
    clusters derived from them, so this corrects within-species batch
    structure end-to-end.

Slot 2 — :func:`distill_linear_loadings`
    Replace ``varm['PCs_SAMap']`` — the linear gene->latent map SAMap uses
    to build the joint cross-species projection — with a ridge distillation
    of the embedding, fitted in SAMap's own feature space (standardized,
    gene-weighted expression). The joint projection space then approximates
    the batch-corrected latent while remaining a *linear* map, which is the
    property SAMap needs to project homology-translated pseudo-cells
    (``ss_i @ (G_ij . D . L_j)``) into it with sparse matmuls.

    Inspect the returned per-dimension R^2. If it is low, the latent is not
    linearly recoverable from expression and slot 2 will be a poor
    approximation — consider training an LDVAE (``scvi.model.LinearSCVI``,
    ``train_scvi_latent(..., linear=True)``), whose latent is constructed
    to be linearly decodable, or skip slot 2 (slot 1 alone captures most of
    the batch-correction benefit).

Typical use with embeddings you already have::

    sam = SAM(); sam.load_data('species1.h5ad')
    sam.preprocess_data(...); sam.run(...)              # standard SAM
    use_embedding_manifold(sam, my_latent)              # slot 1
    distill_linear_loadings(sam)                        # slot 2 (optional)
    ...
    sm = SAMAP({'s1': sam1, 's2': sam2}, f_maps='maps/')

Or train scVI here from the SAM object's linear-scale matrix::

    latent, model = train_scvi_latent(sam, batch_key='donor')
    use_embedding_manifold(sam, latent)

Notes on ordering
-----------------
* Call ``sam.run()`` before slot 1: SAMap needs ``uns['run_args']`` (for
  the PCA-loadings fallback) and the layers SAM's preprocessing creates.
* Slot 2 must come *after* slot 1 for the same SAM object: it fits in the
  gene-weighted feature space, and slot 1 recomputes the gene weights from
  the corrected graph.
* Slot 1 deletes any stale ``varm['PCs_SAMap']`` so that, if you skip
  slot 2, ``SAMAP`` recomputes the PCA loadings with the corrected weights.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np
import scipy.sparse as sp
from sklearn.preprocessing import StandardScaler

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy.typing import NDArray

    from samap.sam import SAM

__all__ = [
    "use_embedding_manifold",
    "distill_linear_loadings",
    "train_scvi_latent",
    "scvi_manifold",
]


# --------------------------------------------------------------------------- #
# Slot 1: within-species manifold from an external embedding                   #
# --------------------------------------------------------------------------- #


def use_embedding_manifold(
    sam: SAM,
    embedding: NDArray[Any],
    key_added: str = "X_external",
    k: int = 20,
    weight_mode: str = "rms",
    num_norm_avg: int = 50,
    leiden_res: float | None = None,
) -> SAM:
    """Rebuild a SAM object's manifold from an external cell embedding.

    Mirrors exactly what ``SAM.calculate_nnm`` does with its own weighted
    PCA (sam/core.py, ``update_manifold`` branch), but sources the
    coordinates from ``embedding`` instead: cosine kNN graph ->
    ``obsp['distances'/'connectivities'/'nnm']`` -> gene weights
    recomputed on the new graph -> optional Leiden clustering.

    Parameters
    ----------
    sam
        A SAM object that has been through ``preprocess_data`` and
        ``run()`` (the run is needed for SAMap bookkeeping; its manifold
        is overwritten here).
    embedding
        (n_cells x d) array, rows aligned with ``sam.adata.obs_names``.
        E.g. ``scvi_model.get_latent_representation()`` or a
        Harmony-corrected PCA.
    key_added
        ``obsm`` key under which the embedding is stored.
    k
        Neighbors per cell (SAM/SAMap default: 20).
    weight_mode, num_norm_avg
        Passed to ``SAM.dispersion_ranking_NN`` for the gene-weight
        recomputation (SAMAP uses ``weight_mode='rms'``).
    leiden_res
        If given, re-run Leiden clustering on the new graph. ``SAMAP``
        re-runs Leiden itself for the default ``keys``, so this is only
        needed if you use the clusters beforehand.

    Returns
    -------
    The same SAM object, modified in place.
    """
    from samap.sam.core import _csr_setdiag
    from samap.sam.knn import calc_nnm

    emb = np.ascontiguousarray(np.asarray(embedding, dtype=np.float32))
    if emb.ndim != 2 or emb.shape[0] != sam.adata.shape[0]:
        raise ValueError(
            f"embedding must be (n_cells={sam.adata.shape[0]} x d); got {emb.shape}"
        )

    # --- kNN graph, exactly as SAM.calculate_nnm builds it ---------------- #
    edm = calc_nnm(emb, k, "cosine")  # CSR of cosine distances

    edm_dist = edm.copy()
    _csr_setdiag(edm_dist, 0)
    sam.adata.obsp["distances"] = edm_dist

    EDM = edm.copy()
    EDM.data[:] = 1
    _csr_setdiag(EDM, 1)
    sam.adata.obsp["connectivities"] = EDM

    sim = edm.copy()
    sim.data[:] = 1 - sim.data
    _csr_setdiag(sim, 1)
    sim.data[sim.data < 0] = 0.001
    sam.adata.obsp["nnm"] = sim

    # --- gene weights from the corrected graph ---------------------------- #
    W = sam.dispersion_ranking_NN(
        EDM, weight_mode=weight_mode, num_norm_avg=num_norm_avg
    )
    sam.adata.var["weights"] = W

    sam.adata.obsm[key_added] = emb

    # Invalidate stale PCA loadings: if slot 2 is skipped, SAMAP recomputes
    # them (prepare_SAMap_loadings) with the corrected gene weights.
    if "PCs_SAMap" in sam.adata.varm:
        del sam.adata.varm["PCs_SAMap"]

    if leiden_res is not None:
        sam.leiden_clustering(res=leiden_res)

    return sam


# --------------------------------------------------------------------------- #
# Slot 2: linear distillation of the embedding into SAMap loadings            #
# --------------------------------------------------------------------------- #


def distill_linear_loadings(
    sam: SAM,
    embedding: NDArray[Any] | None = None,
    key: str = "X_external",
    ridge: float | str = "auto",
    val_frac: float = 0.2,
    seed: int = 0,
    store: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fit a linear gene->latent map to the embedding; store as SAMap loadings.

    Solves the ridge problem, jointly over latent dimensions::

        L = argmin  || Fc @ L - Zc ||_F^2  +  ridge * ||L||_F^2

    where ``Fc`` is SAMap's own projection feature space — column-
    standardized (``StandardScaler(with_mean=False)``), gene-weighted
    expression, centered — matching ``_projection_precompute`` /
    ``_mapping_window_fast`` in samap/core/projection.py, whose generic
    mean-correction reproduces the centering at projection time. ``Zc`` is
    the centered embedding. Because n_cells << n_genes, the problem is
    solved exactly in dual (kernel) form via an (n x n) system.

    The system is underdetermined (n_genes >> n_cells), so a linear map can
    interpolate *any* target in-sample. The meaningful diagnostic is
    therefore the **held-out** R^2 reported here: ``val_frac`` of cells are
    withheld, the ridge strength is chosen on them (when ``ridge='auto'``),
    and the final map is refitted on all cells. Held-out R^2 near 1 means
    the latent is faithfully linearly recoverable from expression; low
    values mean slot 2 degrades to a rough approximation of the latent
    (consider an LDVAE latent, or slot 1 only).

    The resulting matrix (n_genes x d) replaces ``varm['PCs_SAMap']``, so
    SAMap's joint cross-species projection space becomes (approximately)
    the batch-corrected latent. As with the PCA loadings it replaces, rows
    are fitted on all genes but consumed on the homology-gene subset;
    contributions of unmatched genes are dropped at projection time.

    Parameters
    ----------
    sam
        SAM object; run slot 1 first so ``var['weights']`` reflects the
        corrected manifold.
    embedding
        (n_cells x d) latent. Defaults to ``sam.adata.obsm[key]``.
    ridge
        Ridge strength, or 'auto' to select from a grid (scaled to the
        kernel magnitude) by held-out R^2.
    val_frac, seed
        Size and seed of the held-out split used for the R^2 diagnostic
        and ridge selection.
    store
        If True, write the loadings to ``sam.adata.varm['PCs_SAMap']``.

    Returns
    -------
    dict with ``loadings`` (n_genes x d), ``r2`` (per-dimension held-out
    R^2), ``r2_mean``, and ``ridge`` (the value used).
    """
    from scipy.linalg import cho_factor, cho_solve

    Z = embedding if embedding is not None else sam.adata.obsm[key]
    Z = np.asarray(Z, dtype=np.float64)
    if Z.ndim != 2 or Z.shape[0] != sam.adata.shape[0]:
        raise ValueError(
            f"embedding must be (n_cells={sam.adata.shape[0]} x d); got {Z.shape}"
        )

    X = sam.adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    w = np.asarray(sam.adata.var["weights"].values, dtype=np.float64)

    # SAMap's feature construction, verbatim (projection.py::_projection_precompute)
    F = (
        StandardScaler(with_mean=False)
        .fit_transform(X)
        .multiply(w[None, :])
        .tocsr()
        .astype(np.float64)
    )
    n, G = F.shape
    m = np.asarray(F.mean(0)).flatten()  # column means, for implicit centering

    # Centered Gram matrix Kc = Fc Fc^T without materializing Fc (n x n)
    K = np.asarray((F @ F.T).todense(), dtype=np.float64)
    c = F @ m
    Kc = K - c[:, None] - c[None, :] + float(m @ m)

    Zc = Z - Z.mean(0)

    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_val = max(int(round(val_frac * n)), 1)
    val, tr = perm[:n_val], perm[n_val:]

    K_tt = Kc[np.ix_(tr, tr)]
    K_vt = Kc[np.ix_(val, tr)]
    Z_tr, Z_val = Zc[tr], Zc[val]

    def _val_r2(lam: float) -> NDArray[np.float64]:
        Kr = K_tt + lam * np.eye(tr.size)
        alpha = cho_solve(cho_factor(Kr, lower=True), Z_tr)
        pred = K_vt @ alpha
        ss_res = ((pred - Z_val) ** 2).sum(0)
        ss_tot = (Z_val**2).sum(0)
        return 1.0 - ss_res / np.maximum(ss_tot, 1e-12)

    diag_scale = float(np.mean(np.diag(Kc)))
    if ridge == "auto":
        grid = [diag_scale * f for f in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0)]
        scores = {lam: float(_val_r2(lam).mean()) for lam in grid}
        lam = max(scores, key=scores.get)
    else:
        lam = float(ridge)
    r2 = _val_r2(lam)

    # Refit on all cells at the chosen strength; recover primal loadings:
    #   L = Fc^T alpha = F^T alpha - m * sum(alpha)
    alpha = cho_solve(
        cho_factor(Kc + lam * np.eye(n), lower=True), Zc
    )
    L = np.asarray(F.T @ alpha) - m[:, None] * alpha.sum(0)[None, :]

    if verbose:
        print(
            f"[distill_linear_loadings] held-out R^2: mean={r2.mean():.3f}, "
            f"min={r2.min():.3f}, max={r2.max():.3f} (ridge={lam:.3g})"
        )
    if r2.mean() < 0.5:
        warnings.warn(
            f"Held-out R^2 = {r2.mean():.2f} < 0.5: the embedding is poorly "
            "captured by a linear map of expression. Slot 2 will be a rough "
            "approximation - consider LinearSCVI "
            "(train_scvi_latent(..., linear=True)) or slot 1 only.",
            stacklevel=2,
        )

    if store:
        sam.adata.varm["PCs_SAMap"] = L

    return {"loadings": L, "r2": r2, "r2_mean": float(r2.mean()), "ridge": lam}


# --------------------------------------------------------------------------- #
# Convenience: train scVI / LDVAE from a SAM object                           #
# --------------------------------------------------------------------------- #


def _counts_adata(
    sam: SAM,
    counts: AnnData | NDArray[Any] | sp.spmatrix | None = None,
    layer: str | None = None,
    round_non_integer: bool = True,
) -> AnnData:
    """Assemble the (cells x genes) count matrix scVI trains on.

    Priority: explicit ``counts`` > ``sam.adata.layers[layer]`` >
    ``sam.adata_raw`` (SAM keeps the pre-log, linear-scale matrix there).

    scVI's NB/ZINB likelihoods expect integer counts. SAM itself only
    requires *linear-scale* expression, so ``adata_raw`` may hold
    depth-normalized data (CP10K, estimated counts, ...): with
    ``round_non_integer`` such matrices are rounded, with a warning. Pass
    true counts via ``counts=`` whenever you have them.
    """
    from anndata import AnnData

    obs = sam.adata.obs.copy()

    if counts is not None:
        if hasattr(counts, "X"):  # AnnData
            ad = counts[sam.adata.obs_names].copy()
            for c in obs.columns.difference(ad.obs.columns):
                ad.obs[c] = obs[c].values
        else:
            Xc = counts if sp.issparse(counts) else sp.csr_matrix(np.asarray(counts))
            ad = AnnData(X=Xc.tocsr().astype(np.float32), obs=obs)
    elif layer is not None:
        ad = AnnData(X=sam.adata.layers[layer].copy(), obs=obs)
        ad.var_names = sam.adata.var_names
    else:
        if not hasattr(sam, "adata_raw"):
            raise ValueError("SAM object has no adata_raw; pass counts= explicitly.")
        ad = AnnData(X=sam.adata_raw.X.copy(), obs=obs)
        ad.var_names = sam.adata_raw.var_names

    X = ad.X if sp.issparse(ad.X) else sp.csr_matrix(ad.X)
    X = X.tocsr()
    head = X[: min(200, X.shape[0])].toarray()
    if not np.allclose(head, np.round(head)):
        if not round_non_integer:
            raise ValueError(
                "Count matrix is not integer-valued and round_non_integer=False."
            )
        warnings.warn(
            "Count matrix for scVI is not integer-valued (SAM only stores "
            "linear-scale expression, which may be depth-normalized). Rounding "
            "to integers - fine for exploration; pass counts= with true raw "
            "counts for production runs.",
            stacklevel=2,
        )
        X.data = np.round(X.data)
        X.eliminate_zeros()
    ad.X = X.astype(np.float32)

    # drop all-zero genes (uninformative for scVI; keeps varnames meaningful)
    nz = np.asarray((ad.X > 0).sum(0)).flatten() > 0
    return ad[:, nz].copy()


def train_scvi_latent(
    sam: SAM,
    batch_key: str | None = None,
    linear: bool = False,
    n_latent: int = 30,
    max_epochs: int | None = None,
    counts: AnnData | NDArray[Any] | sp.spmatrix | None = None,
    layer: str | None = None,
    accelerator: str = "cpu",
    train_kwargs: dict[str, Any] | None = None,
    **model_kwargs: Any,
) -> tuple[NDArray[np.floating[Any]], Any]:
    """Train scVI (or LDVAE) on a SAM object's counts; return (latent, model).

    Parameters
    ----------
    batch_key
        Column of ``sam.adata.obs`` holding the within-species batch labels
        the latent should be corrected for.
    linear
        If True train ``scvi.model.LinearSCVI`` (LDVAE) — its latent is
        constructed to be linearly decodable, which makes slot 2's
        distillation more faithful.
    counts, layer
        Where to find integer counts; see :func:`_counts_adata`.
    accelerator
        Lightning accelerator. 'cpu' is the safe default for the small
        per-species datasets SAMap typically sees.
    """
    import scvi

    ad = _counts_adata(sam, counts=counts, layer=layer)
    cls = scvi.model.LinearSCVI if linear else scvi.model.SCVI
    cls.setup_anndata(ad, batch_key=batch_key)
    model = cls(ad, n_latent=n_latent, **model_kwargs)
    model.train(
        max_epochs=max_epochs,
        accelerator=accelerator,
        enable_progress_bar=False,
        **(train_kwargs or {}),
    )
    latent = np.asarray(model.get_latent_representation())
    return latent, model


def scvi_manifold(
    sam: SAM,
    batch_key: str | None = None,
    linear: bool = False,
    n_latent: int = 30,
    max_epochs: int | None = None,
    k: int = 20,
    distill: bool = False,
    key_added: str = "X_scvi",
    **train_kw: Any,
) -> tuple[NDArray[np.floating[Any]], Any]:
    """One-call slot 1 (and optionally slot 2) with a freshly trained scVI.

    Equivalent to ``train_scvi_latent`` + ``use_embedding_manifold``
    (+ ``distill_linear_loadings`` when ``distill=True``).
    Returns (latent, model).
    """
    latent, model = train_scvi_latent(
        sam,
        batch_key=batch_key,
        linear=linear,
        n_latent=n_latent,
        max_epochs=max_epochs,
        **train_kw,
    )
    use_embedding_manifold(sam, latent, key_added=key_added, k=k)
    if distill:
        distill_linear_loadings(sam, key=key_added)
    return latent, model
