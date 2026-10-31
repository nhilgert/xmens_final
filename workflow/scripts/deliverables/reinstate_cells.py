"""Add the reinstated epithelial nuclei to the human all-cells release and label its epithelium.

04g reinstate_release_cells. Inputs:
  base       the human all-cells object relabeled by migrate_release_vocabulary
  lineage    data/interim/latent_human_epithelial.h5ad, the epithelial lineage clustering the
             decisions were written on (pinned by sha256; cluster column from config)
  spatial    data/interim/qc_spatial_human.h5ad: coordinates and spatial QC of every nucleus
  decisions  resources/annotations/human_epithelial.csv (cluster -> label)

Clusters 11, 22 and 44 were excluded before the July whole-tissue scVI fit, so the release has no
row for their nuclei; their decisions now keep them. This step
  1. labels every epithelial nucleus of the release by its cluster's decision
     (epithelial_inner, epithelial_outer, epithelial_ciliated);
  2. appends the nuclei of the reinstated clusters after the release's own: raw counts from the
     lineage object, .X = log1p(CP10K) as in the release, legacy label columns as on the other
     non-ciliated epithelial nuclei, timepoint from their sample, and coordinates, spatial QC
     and tissue_section from the spatial table;
  3. removes the July embedding (X_scVI, X_umap_scVI, graph, HVG flags) and embedding_outlier,
     which reembed_release and finish_reembedded_release recompute for all nuclei.

It fails unless the release and the clustering agree in both directions: every epithelial
nucleus of the release sits in a kept cluster, the clustering's nuclei that are in the release
are all epithelial there, and the clusters absent from the release are exactly the EXCLUDE
clusters plus the reinstated ones. On the nuclei present in both it also checks that the lineage
counts, the spatial table and the tissue_section rule reproduce the release, and at the end that
the labels match the vocabulary counts.
"""
import hashlib
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "atlas"))
from apply_annotation import load_choices, validate_canonical_annotations  # noqa: E402

snk = snakemake  # noqa: F821 (injected by Snakemake)

EPITHELIAL = "epithelial_cell"
CILIATED = "epithelial_ciliated"
EXCLUDE = "EXCLUDE"
# Carried over unchanged from the non-ciliated epithelial nuclei of the release.
LEGACY_COLUMNS = ("coarse_refined", "fine_refined", "origin_lineage", "lineage", "cell_type",
                  "subtype", "annotation_level")
SPATIAL_COLUMNS = ("x_um_dbscan", "y_um_dbscan", "x_um_qc", "y_um_qc", "spatial_qc_reason")
# Recomputed for all nuclei by reembed_release (integrate.py) and finish_reembedded_release.
EMBEDDING_OBSM = ("X_scVI", "X_umap_scVI")
EMBEDDING_UNS = ("neighbors_scVI_100", "umap", "hvg", "log1p")
HVG_VAR = ("highly_variable", "highly_variable_rank", "highly_variable_nbatches")
RECOMPUTED_OBS = ("embedding_outlier",)

LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[reinstate] {msg}", file=LOG, flush=True)


def sha256_of(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def same_values(a: pd.Series, b: pd.Series) -> bool:
    """Equal element-wise, NaN matching NaN; strings compared as strings."""
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        return bool(np.array_equal(a.to_numpy(float), b.to_numpy(float), equal_nan=True))
    return bool((a.astype(object).where(a.notna(), "<NA>").astype(str).values
                 == b.astype(object).where(b.notna(), "<NA>").astype(str).values).all())


def log1p_cp10k(counts):
    """.X of the release: log1p of counts normalised to 10,000 per nucleus (as integrate.py)."""
    tmp = ad.AnnData(X=counts.copy())
    sc.pp.normalize_total(tmp, target_sum=1e4, exclude_highly_expressed=False, inplace=True)
    sc.pp.log1p(tmp)
    return sp.csr_matrix(tmp.X)


def tissue_section(sample: pd.Series, reason: pd.Series) -> pd.Series:
    """Human pucks hold one cross-section: "{sample}__A" on spatially kept nuclei, else missing."""
    section = sample.astype(str) + "__A"
    return section.where(reason.astype(str).values == "kept", other=np.nan)


def extend_column(base_column: pd.Series, added: pd.Series):
    """Base values then added values, in the base column's dtype.

    Categorical columns keep their category order (positional uns colors depend on it) and
    append any new category.
    """
    if isinstance(base_column.dtype, pd.CategoricalDtype):
        values = np.concatenate([base_column.astype(object).values, added.astype(object).values])
        categories = list(base_column.cat.categories)
        categories += sorted(set(added.dropna().astype(str)) - set(categories))
        return pd.Categorical(values, categories=categories)
    return np.concatenate([base_column.to_numpy(), added.to_numpy(dtype=base_column.dtype)])


def check_clusters_against_release(base, cluster_of, decisions, reinstated):
    """The release and the clustering must agree on which nuclei are kept epithelium."""
    unknown = sorted(set(cluster_of) - set(decisions))
    if unknown:
        raise ValueError(f"clusters without a decision in {snk.input.decisions}: {unknown}")
    kept_reinstated = [c for c in reinstated if decisions.get(c, EXCLUDE) != EXCLUDE]
    if kept_reinstated != list(reinstated):
        raise ValueError(f"reinstated clusters must have a kept decision: {reinstated}")

    types = base.obs["cell_type"].astype(str)
    epithelial = base.obs_names[types.values == EPITHELIAL]
    not_clustered = epithelial.difference(cluster_of.index)
    if len(not_clustered):
        raise ValueError(f"{len(not_clustered)} epithelial nuclei of the release are not in the clustering")
    labels = cluster_of.loc[epithelial].map(decisions)
    if labels.eq(EXCLUDE).any():
        raise ValueError(f"{int(labels.eq(EXCLUDE).sum())} epithelial nuclei of the release sit in EXCLUDE clusters")

    in_release = cluster_of.index.isin(base.obs_names)
    other_type = types.loc[cluster_of.index[in_release]].ne(EPITHELIAL)
    if other_type.any():
        raise ValueError(f"{int(other_type.sum())} nuclei of the clustering are in the release as another cell type")

    absent = set(cluster_of[~in_release])
    expected_absent = {c for c, label in decisions.items() if label == EXCLUDE} | set(reinstated)
    if absent != expected_absent:
        raise ValueError(f"clusters absent from the release {sorted(absent)} != EXCLUDE + reinstated "
                         f"{sorted(expected_absent)}")
    partly = sorted(c for c in expected_absent if in_release[cluster_of.values == c].any())
    if partly:
        raise ValueError(f"EXCLUDE or reinstated clusters partly present in the release: {partly}")
    return labels


def check_counts_against_vocabulary(obs, added_labels, vocab):
    """Epithelial labels, reinstated nuclei and coarse counts as recorded in the vocabulary."""
    tasks = vocab["content_tasks"]
    subtypes = obs["cell_subtype"].astype(str)
    split = {label: int((subtypes == label).sum()) for label in tasks["human_epithelial_split"]["expected_counts"]}
    added = {label: int(n) for label, n in added_labels.value_counts().items()}
    coarse = obs["cell_type"].astype(str).value_counts().to_dict()
    expected_coarse = dict(vocab["expected_coarse_counts"]["human"])
    checks = [
        ("human_epithelial_split", split, tasks["human_epithelial_split"]["expected_counts"]),
        ("human_epithelial_reinstatement", added, tasks["human_epithelial_reinstatement"]["expected_counts"]),
        ("expected_coarse_counts", coarse, expected_coarse),
        ("total", len(obs), vocab["expected_coarse_counts"]["totals"]["human"]),
    ]
    for name, found, expected in checks:
        if found != expected:
            raise ValueError(f"{name}: found {found}, vocabulary expects {expected}")
        say(f"  {name}: {found}")


def main():
    vocab = yaml.safe_load(open(snk.input.vocab))
    pinned = str(snk.params.lineage_sha256)
    found_sha = sha256_of(snk.input.lineage)
    if found_sha != pinned:
        raise ValueError(f"{snk.input.lineage}: sha256 {found_sha} != pinned {pinned}")
    reinstated = [str(c) for c in snk.params.reinstated_clusters]

    base = ad.read_h5ad(snk.input.base)
    lineage = ad.read_h5ad(snk.input.lineage)
    spatial = ad.read_h5ad(snk.input.spatial, backed="r").obs[["sample", *SPATIAL_COLUMNS]].copy()
    say(f"release {base.n_obs:,} nuclei; clustering {lineage.n_obs:,}; spatial table {len(spatial):,}")
    for name, names in (("release", base.obs_names), ("clustering", lineage.obs_names),
                        ("spatial table", spatial.index)):
        if not names.is_unique:
            raise ValueError(f"{name} has duplicate nucleus names")
    if not lineage.var_names.equals(base.var_names):
        raise ValueError("clustering genes differ from the release genes (names or order)")

    cluster_of = lineage.obs[snk.params.cluster_column].astype(str)
    decisions = load_choices(snk.input.decisions)
    allowed = set(vocab["allowed"][EPITHELIAL])
    bad = sorted(set(decisions.values()) - allowed - {EXCLUDE})
    if bad:
        raise ValueError(f"decision labels not allowed under {EPITHELIAL}: {bad}")
    labels_in_release = check_clusters_against_release(base, cluster_of, decisions, reinstated)

    # 1. Label the release's epithelium. Ciliated nuclei must keep their label.
    was_ciliated = base.obs["cell_subtype"].astype(str).loc[labels_in_release.index].eq(CILIATED)
    if not (labels_in_release[was_ciliated] == CILIATED).all() or \
            labels_in_release[~was_ciliated].eq(CILIATED).any():
        raise ValueError("ciliated decisions do not match the release's epithelial_ciliated nuclei")
    coarse_epithelial = labels_in_release.index[~was_ciliated.values]
    legacy = {}
    for column in LEGACY_COLUMNS:
        values = set(base.obs.loc[coarse_epithelial, column].astype(str))
        if len(values) != 1:
            raise ValueError(f"{column} is not constant on non-ciliated epithelial nuclei: {sorted(values)}")
        legacy[column] = values.pop()
    say(f"legacy columns of non-ciliated epithelium: {legacy}")

    subtype = base.obs["cell_subtype"].copy()
    subtype = subtype.cat.add_categories(sorted(set(labels_in_release) - set(subtype.cat.categories)))
    subtype.loc[labels_in_release.index] = labels_in_release.values
    base.obs["cell_subtype"] = subtype

    # 2. The reinstated nuclei, checked against the release on the nuclei present in both.
    shared = labels_in_release.index
    added = lineage.obs_names[cluster_of.isin(reinstated).values]
    say(f"reinstating {len(added):,} nuclei from clusters {reinstated}")

    base_counts = sp.csr_matrix(base.layers["counts"])
    lineage_counts = sp.csr_matrix(lineage.layers["counts"])
    rows_base = base.obs_names.get_indexer(shared)
    rows_lineage = lineage.obs_names.get_indexer(shared)
    if (base_counts[rows_base] != lineage_counts[rows_lineage]).nnz:
        raise ValueError("clustering counts differ from release counts on shared nuclei")
    if (log1p_cp10k(base_counts) != sp.csr_matrix(base.X)).nnz:
        raise ValueError("release .X is not log1p(CP10K) of its counts")
    for column in ("sample", "cell_barcode"):
        if not same_values(lineage.obs.loc[shared, column], base.obs.loc[shared, column]):
            raise ValueError(f"{column} differs between clustering and release")

    missing_spatial = base.obs_names.union(added).difference(spatial.index)
    if len(missing_spatial):
        raise ValueError(f"{len(missing_spatial)} nuclei have no row in the spatial table")
    for column in SPATIAL_COLUMNS:
        if not same_values(spatial.loc[base.obs_names, column], base.obs[column]):
            raise ValueError(f"spatial table {column} differs from the release")
    xy = base.obs[["x_um_qc", "y_um_qc"]].to_numpy(float)
    if not np.array_equal(np.asarray(base.obsm["X_spatial"], float), xy, equal_nan=True):
        raise ValueError("release X_spatial is not (x_um_qc, y_um_qc)")
    if not same_values(tissue_section(base.obs["sample"], base.obs["spatial_qc_reason"]),
                       base.obs["tissue_section"]):
        raise ValueError("tissue_section rule does not reproduce the release")

    timepoints = base.obs.groupby("sample", observed=True)["timepoint"].agg(lambda s: set(s.astype(str)))
    if (timepoints.map(len) != 1).any():
        raise ValueError("a sample has more than one timepoint in the release")
    timepoint_of = timepoints.map(lambda s: next(iter(s)))

    counts_added = lineage_counts[lineage.obs_names.get_indexer(added)].astype(base_counts.dtype)
    obs_added = pd.DataFrame(index=added)
    obs_added["sample"] = lineage.obs.loc[added, "sample"].astype(str)
    obs_added["cell_barcode"] = lineage.obs.loc[added, "cell_barcode"].astype(str)
    obs_added["timepoint"] = obs_added["sample"].map(timepoint_of)
    if obs_added["timepoint"].isna().any():
        raise ValueError("reinstated nuclei from a sample not in the release")
    for column, value in legacy.items():
        obs_added[column] = value
    obs_added["cell_subtype"] = cluster_of.loc[added].map(decisions).values
    for column in SPATIAL_COLUMNS:
        obs_added[column] = spatial.loc[added, column].values
    obs_added["tissue_section"] = tissue_section(obs_added["sample"], obs_added["spatial_qc_reason"]).values

    # 3. One object: release nuclei first, in release order, then the reinstated nuclei.
    columns = [c for c in base.obs.columns if c not in RECOMPUTED_OBS]
    missing_columns = sorted(set(columns) - set(obs_added.columns))
    if missing_columns:
        raise ValueError(f"no value for reinstated nuclei in {missing_columns}")
    obs = pd.DataFrame(index=base.obs_names.append(added))
    for column in columns:
        obs[column] = extend_column(base.obs[column], obs_added[column])
    check_counts_against_vocabulary(obs, obs_added["cell_subtype"], vocab)

    uns = {key: value for key, value in base.uns.items() if key not in EMBEDDING_UNS}
    uns["epithelial_reinstatement"] = {
        "clusters": reinstated,
        "n_nuclei": int(len(added)),
        "clustering": str(snk.input.lineage),
        "clustering_sha256": found_sha,
        "cluster_column": str(snk.params.cluster_column),
        "decisions": str(snk.input.decisions),
        "note": ("Epithelial nuclei excluded before the July whole-tissue fit and kept by the "
                 "v7 decisions; appended after the release's nuclei."),
    }
    out = ad.AnnData(
        X=sp.vstack([sp.csr_matrix(base.X), log1p_cp10k(counts_added)], format="csr"),
        obs=obs,
        var=base.var.drop(columns=[c for c in HVG_VAR if c in base.var.columns]),
        uns=uns,
        obsm={"X_spatial": np.vstack([np.asarray(base.obsm["X_spatial"]),
                                      obs_added[["x_um_qc", "y_um_qc"]].to_numpy(base.obsm["X_spatial"].dtype)])},
        layers={"counts": sp.vstack([base_counts, counts_added], format="csr")},
    )
    validate_canonical_annotations(out, vocab, context="reinstated human all-cells")
    say(f"dropped from the release: obsm {list(EMBEDDING_OBSM)}, uns {list(EMBEDDING_UNS)}, "
        f"var {list(HVG_VAR)}, obs {list(RECOMPUTED_OBS)}, obsp {list(base.obsp.keys())}")

    report = (pd.DataFrame({"cluster": cluster_of.values, "in_release": cluster_of.index.isin(base.obs_names)})
              .assign(label=lambda d: d.cluster.map(decisions))
              .groupby(["cluster", "label"], as_index=False)
              .agg(n_nuclei=("in_release", "size"), n_in_release=("in_release", "sum")))
    report["n_added"] = np.where(report.cluster.isin(reinstated), report.n_nuclei, 0)
    report = report.sort_values("cluster", key=lambda s: s.astype(int))
    Path(snk.output.report).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(snk.output.report, sep="\t", index=False)

    out_path = Path(snk.output.h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp.h5ad")
    out.write_h5ad(tmp_path)
    os.replace(tmp_path, out_path)
    say(f"wrote {out.n_obs:,} nuclei x {out.n_vars:,} genes -> {out_path}")


main()
LOG.close()
