"""
Map Ulrich/Hammoud (PNAS 2024) fibroblast subtype centroids onto our human
snRNA-seq fibroblast/decidual nuclei.

Design notes
------------
Primary statistic is a cluster-level centroid comparison, matching the method
Hammoud et al. used to reconcile annotations across studies (Fig. 3A is a 39x39
centroid rank-correlation matrix). Our data are single-nucleus and theirs are
single-cell, so raw correlation carries a shared dataset-level offset; the
double-centered ("contrastive") cosine removes it and asks only about relative
subtype pattern. Both are reported.

Centroids are computed once over every shared gene. Every downstream robustness
check resamples *genes*, and log1p(mean(expm1(x))) is per-gene independent, so
those checks are row-indexing on the cached centroid matrix rather than repeated
passes over the count matrix.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.preprocessing import normalize

# ---------------------------------------------------------------- atlas side

# The six stromal subtypes that are fibroblasts. Stro-1/2/3 are perivascular,
# Stro-4 is a vascular progenitor and Stro-6 is uSMC; our object is
# fibroblast/decidual only, so those four act as negative controls.
FIB_SUBTYPES = [
    "Stro-5_RUNX3+ Fibroblast",
    "Stro-7_C7 Fibroblast",
    "Stro-8_dS1",
    "Stro-9_dS2",
    "Stro-10_dS3",
    "Stro-11_Stromal progenitor",
]

SHORT = {
    "Stro-5_RUNX3+ Fibroblast": "Stro-5 RUNX3+",
    "Stro-7_C7 Fibroblast": "Stro-7 C7",
    "Stro-8_dS1": "Stro-8 dS1",
    "Stro-9_dS2": "Stro-9 dS2",
    "Stro-10_dS3": "Stro-10 dS3",
    "Stro-11_Stromal progenitor": "Stro-11 progenitor",
}


def load_atlas_centroids(path, sheet="5B Subtype centroids"):
    """Parse Dataset S5B. Returns genes x 39 subtypes, log1p(CP10K)."""
    df = pd.ExcelFile(path).parse(sheet, header=2, index_col=0)
    df.columns = (
        pd.Index(df.columns).astype(str)
        .str.replace(r"\s+", " ", regex=True).str.strip()
    )
    df.index = df.index.astype(str).str.strip()
    df = df[~df.index.isin(["nan", ""])]
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df[~df.index.duplicated(keep="first")]

    colsum = np.expm1(df).sum(axis=0).median()
    if not (9e3 < colsum < 1.1e4):
        raise ValueError(f"centroids are not log1p(CP10K): median colsum {colsum:,.0f}")
    return df


def load_atlas_markers(path, sheet="6C Stro markers"):
    """Parse Dataset S6C (the paper's own DE markers per stromal subtype).

    The spreadsheet auto-increments the Stro-3 label down consecutive rows
    (3_PV STEAP4, 3_PV STEAP5, ...); those collapse back to one cluster.
    """
    import re

    d = pd.ExcelFile(path).parse(sheet, header=2)
    d.columns = [str(c).strip() for c in d.columns]
    d = d.dropna(subset=["cluster", "gene"]).copy()
    d["cluster"] = d["cluster"].map(
        lambda x: "3_PV STEAP4"
        if re.fullmatch(r"3_PV STEAP\d+", str(x).strip())
        else str(x).strip()
    )
    d["gene"] = d["gene"].astype(str).str.strip()
    return d


# S6C cluster labels -> S5B centroid column names
S6C_TO_S5B = {
    "5_RUNX3+ Fibroblast": "Stro-5_RUNX3+ Fibroblast",
    "7_C7 Fibroblast": "Stro-7_C7 Fibroblast",
    "8_dS1": "Stro-8_dS1",
    "9_dS2": "Stro-9_dS2",
    "10_dS3": "Stro-10_dS3",
    "11_Stromal progenitor": "Stro-11_Stromal progenitor",
    "1_PV MYH11 FOXC1 hi CTGF lo": "Stro-1_PV MYH11 FOXC1hi CTGFlo",
    "2_PV MYH11 FOXC1 lo CTGF hi": "Stro-2_PV MYH11 FOXC1lo CTGFhi",
    "3_PV STEAP4": "Stro-3_PV STEAP4",
    "4_Tissue resident vascular progenitor": "Stro-4_Tissue resident vascular progenitor",
    "6_uSMC": "Stro-6_uSMC",
}


# ---------------------------------------------------------------- query side

def linear_panel_matrix(adata, genes, layer=None):
    """Expression on `genes` in linear CP10K space. expm1 preserves sparsity."""
    X = adata[:, genes].layers[layer] if layer else adata[:, genes].X
    X = sp.csr_matrix(X, dtype=np.float32).copy()
    X.data = np.expm1(X.data)
    return X


def group_centroids(X_lin, labels, genes, min_n=1):
    """log1p of the within-group mean of linear expression.

    Averaging in linear space then logging is the correct order; averaging
    already-logged values is not the same quantity.
    """
    labels = pd.Series(np.asarray(labels).astype(str))
    cats = sorted(labels.unique())
    codes = pd.Categorical(labels, categories=cats).codes

    G = sp.csr_matrix(
        (np.ones(len(codes), np.float32), (codes, np.arange(len(codes)))),
        shape=(len(cats), len(codes)),
    )
    n = np.asarray(G.sum(axis=1)).ravel()
    means = np.asarray((G @ X_lin).todense()) / n[:, None]

    C = pd.DataFrame(np.log1p(means).T, index=genes, columns=cats)
    ns = pd.Series(n.astype(int), index=cats, name="n")

    keep = ns[ns >= min_n].index
    return C[keep], ns[keep]


# ---------------------------------------------------------------- similarity

def contrastive_cosine(Q, A):
    """Cosine after centering each gene across the columns on each side.

    Removes the shared dataset-level offset (single-nucleus vs single-cell) and
    compares relative subtype pattern. Note the result depends on which atlas
    columns are passed, because the centering is over exactly those columns.
    """
    if Q.shape[1] < 2 or A.shape[1] < 2:
        raise ValueError("contrastive cosine needs >=2 columns on each side")
    Qc = Q.sub(Q.mean(axis=1), axis=0)
    Ac = A.sub(A.mean(axis=1), axis=0)
    Qn = normalize(np.nan_to_num(Qc.T.to_numpy(float)), axis=1)
    An = normalize(np.nan_to_num(Ac.T.to_numpy(float)), axis=1)
    return pd.DataFrame(Qn @ An.T, index=Q.columns, columns=A.columns)


def rank_correlation(Q, A):
    """Spearman between every query column and every atlas column.

    This is the statistic used in the paper's own Fig. 3A.
    """
    Qr = Q.rank(axis=0).to_numpy(float)
    Ar = A.rank(axis=0).to_numpy(float)
    Qz = (Qr - Qr.mean(0)) / Qr.std(0)
    Az = (Ar - Ar.mean(0)) / Ar.std(0)
    return pd.DataFrame(Qz.T @ Az / len(Qr), index=Q.columns, columns=A.columns)


def _margin(row):
    top = row.nlargest(2)
    return top.iloc[0] - top.iloc[1] if len(top) > 1 else np.nan


def correspondence(M, ns=None):
    """Both directions of the mapping, with margins and reciprocity.

    Rows are query clusters, columns are atlas subtypes.
    """
    q_to_a = M.idxmax(axis=1)
    a_to_q = M.idxmax(axis=0)

    by_query = pd.DataFrame({
        "best_atlas": q_to_a,
        "score": M.max(axis=1),
        "margin": M.apply(_margin, axis=1),
        "reciprocal": [a_to_q[q_to_a[q]] == q for q in M.index],
    })
    if ns is not None:
        by_query["n"] = ns.reindex(by_query.index)

    by_atlas = pd.DataFrame({
        "best_query": a_to_q,
        "score": M.max(axis=0),
        "margin": M.apply(_margin, axis=0),
        "reciprocal": [q_to_a[a_to_q[a]] == a for a in M.columns],
    })
    if ns is not None:
        by_atlas["n_in_best_query"] = [ns.get(a_to_q[a], np.nan) for a in M.columns]

    return by_query, by_atlas


# ---------------------------------------------------------------- gene panels

def variance_panel(atlas, cols, genes, n=3000):
    """Genes with the highest variance across the subtypes being discriminated."""
    v = atlas.loc[genes, cols].var(axis=1)
    return v.nlargest(min(n, len(v))).index


def specificity(atlas, subtype, cols):
    """How much higher a gene is in `subtype` than in its best competitor."""
    others = [c for c in cols if c != subtype]
    return (atlas[subtype] - atlas[others].max(axis=1)).sort_values(ascending=False)


def specific_genes(atlas, subtype, cols, genes, n=50, min_expr=0.5):
    """Top specific genes for a subtype, restricted to genes we can detect."""
    spec = specificity(atlas.loc[genes], subtype, cols)
    expressed = atlas.loc[genes, subtype] > min_expr
    return spec[expressed.reindex(spec.index).fillna(False)].head(n).index


# ---------------------------------------------------------------- robustness

def marker_ablation(C, atlas, targets, cols, panel, grid=(0, 3, 5, 10, 20, 50, 100)):
    """Does each subtype keep the same best-matching cluster as its top
    specific genes are progressively removed?"""
    rows = []
    for s in targets:
        spec = specificity(atlas.loc[panel], s, cols)
        for k in grid:
            genes = panel.difference(spec.head(k).index)
            if len(genes) < 50:
                continue
            M = contrastive_cosine(C.loc[genes], atlas.loc[genes, cols])
            best = M[s].idxmax()
            vals = M[s].sort_values(ascending=False)
            rows.append({
                "subtype": s,
                "n_removed": k,
                "best_query": best,
                "score": vals.iloc[0],
                "margin": vals.iloc[0] - vals.iloc[1],
            })
    return pd.DataFrame(rows)


def gene_bootstrap(C, atlas, targets, cols, panel, n_boot=500, frac=0.70, seed=0):
    """How often does each subtype recover the same best cluster across random
    gene subsets?"""
    rng = np.random.default_rng(seed)
    base = panel.to_numpy()
    size = max(50, int(len(base) * frac))

    hits = {s: [] for s in targets}
    scores = {s: [] for s in targets}
    for _ in range(n_boot):
        genes = pd.Index(rng.choice(base, size=size, replace=False))
        M = contrastive_cosine(C.loc[genes], atlas.loc[genes, cols])
        for s in targets:
            hits[s].append(M[s].idxmax())
            scores[s].append(M[s].max())

    rows = []
    for s in targets:
        vc = pd.Series(hits[s]).value_counts(normalize=True)
        rows.append({
            "subtype": s,
            "modal_query": vc.index[0],
            "frequency": vc.iloc[0],
            "runner_up": vc.index[1] if len(vc) > 1 else None,
            "runner_up_freq": vc.iloc[1] if len(vc) > 1 else 0.0,
            "mean_score": np.mean(scores[s]),
            "sd_score": np.std(scores[s]),
        })
    return pd.DataFrame(rows)


def heldout_split(C, atlas, targets, cols, panel, n_rep=500, top_n=200,
                  disc_frac=0.5, seed=0):
    """Pick the best-matching cluster on discovery genes, then test whether it
    survives on genes that played no part in the selection."""
    rng = np.random.default_rng(seed)
    rows = []
    for s in targets:
        informative = specificity(atlas.loc[panel], s, cols).head(top_n).index.to_numpy()
        agree = []
        margins = []
        for _ in range(n_rep):
            perm = rng.permutation(informative)
            cut = max(20, int(len(perm) * disc_frac))
            disc, val = pd.Index(perm[:cut]), pd.Index(perm[cut:])
            if len(val) < 20:
                continue
            Md = contrastive_cosine(C.loc[disc], atlas.loc[disc, cols])
            selected = Md[s].idxmax()
            Mv = contrastive_cosine(C.loc[val], atlas.loc[val, cols])
            vals = Mv[s].sort_values(ascending=False)
            agree.append(Mv[s].idxmax() == selected)
            margins.append(vals.iloc[0] - vals.iloc[1])
        rows.append({
            "subtype": s,
            "validation_agreement": float(np.mean(agree)),
            "median_validation_margin": float(np.median(margins)),
            "n_rep": len(agree),
        })
    return pd.DataFrame(rows)


def per_group_mapping(X_lin, obs, cluster_key, split_key, atlas, cols, panel,
                      genes, min_n=20):
    """Recompute the mapping independently within each sample / phase.

    This is the replication axis: the assignment should not depend on any one
    section, and it needs no smoothing to make that argument.
    """
    rows = []
    for g in sorted(obs[split_key].astype(str).unique()):
        m = (obs[split_key].astype(str) == g).to_numpy()
        if m.sum() < min_n:
            continue
        C_g, ns_g = group_centroids(
            X_lin[m], obs.loc[m, cluster_key], genes, min_n=min_n
        )
        if C_g.shape[1] < 2:
            continue
        M = contrastive_cosine(C_g.loc[panel], atlas.loc[panel, cols])
        for s in cols:
            vals = M[s].sort_values(ascending=False)
            rows.append({
                split_key: g,
                "subtype": s,
                "best_query": vals.index[0],
                "score": vals.iloc[0],
                "margin": vals.iloc[0] - vals.iloc[1],
                "n_in_best_query": int(ns_g[vals.index[0]]),
            })
    return pd.DataFrame(rows)
