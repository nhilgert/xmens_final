"""Build one species' extended all-cells object for release 2026_10_29b (annotation vocabulary v8).

04g rebuild_release_cells. Inputs (config release_rebuild):
  base               the species' all-cells object of release 2026_10_29
  atlas              data/interim/latent_{species}.h5ad: the whole-tissue clustering (pinned by sha256)
  lineage_decisions  resources/annotations/{species}_lineage.csv (whole-tissue cluster -> class or EXCLUDE)
  dropped            class clusterings (pinned by sha256) and decision CSVs of the clusters newly EXCLUDE in v8
  supplementary      resources/annotations/supplementary_populations.csv
  vocab              resources/annotation_vocabulary.yaml

This step
  1. drops the nuclei of the class clusters newly EXCLUDE in v8 (human immune 0, mouse epithelial 28);
  2. appends the nuclei of the supplementary whole-tissue clusters after the base nuclei, labelled from the table:
     raw counts from the atlas, .X = log1p(CP10K) as in the release, timepoint from their sample, legacy label
     columns as on the base nuclei of the same cell type (origin_lineage "whole_tissue"), the raw Slide-tags
     (DBSCAN) coordinates and no spatial QC (spatial_qc_reason "not_assessed", no tissue_section or X_spatial);
     float columns the base carries for other analyses (mouse auc_*) are left missing;
  3. marks every nucleus with obs in_primary (False for the supplementary nuclei);
  4. removes the 2026_10_29 embedding (X_scVI, X_umap_scVI, neighbour/UMAP/HVG metadata, scvi_fit, HVG flags)
     and embedding_outlier, which reembed_release and finish_rebuilt_release recompute for all nuclei.

It fails unless: the dropped clusters are EXCLUDE in their decision tables and all of their nuclei are in the
base, and no nucleus of any other EXCLUDE cluster of those clusterings is; the supplementary clusters are EXCLUDE
at the whole-tissue step, have the recorded sizes, and none of their nuclei is in the base; the base counts equal
the atlas counts; and the labels pass the vocabulary with its primary, extended, exclusion and supplementary
counts.
"""
import csv
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

EXCLUDE = "EXCLUDE"
WHOLE_TISSUE = "whole_tissue"          # origin_lineage of the supplementary nuclei
NOT_ASSESSED = "not_assessed"          # spatial_qc_reason of the supplementary nuclei
# Legacy coarse label for a cell type the base does not contain (mouse pericytes); as in the human objects.
COARSE_REFINED_NEW = {"perivascular_cell": "perivascular cells"}
EMBEDDING_OBSM = ("X_scVI", "X_umap_scVI")
EMBEDDING_UNS = ("neighbors_scVI_100", "umap", "hvg", "log1p", "scvi_fit")
HVG_VAR = ("highly_variable", "highly_variable_rank", "highly_variable_nbatches")
RECOMPUTED_OBS = ("embedding_outlier",)
# .X of the releases is log1p(CP10K) of the counts in float32; the mouse objects were normalised along another
# float path in July and differ from a recomputation by at most 3.4e-7, so .X is compared within this tolerance.
X_TOLERANCE = 1e-6

LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[rebuild:{snk.params.species}] {msg}", file=LOG, flush=True)


def sha256_of(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def check_pinned(path, pinned):
    found = sha256_of(path)
    if found != str(pinned):
        raise ValueError(f"{path}: sha256 {found} != pinned {pinned}")
    return found


def log1p_cp10k(counts):
    """.X of the release: log1p of counts normalised to 10,000 per nucleus (as integrate.py)."""
    tmp = ad.AnnData(X=counts.copy())
    sc.pp.normalize_total(tmp, target_sum=1e4, exclude_highly_expressed=False, inplace=True)
    sc.pp.log1p(tmp)
    return sp.csr_matrix(tmp.X)


def max_difference(a, b) -> float:
    """Largest absolute difference of two sparse matrices; inf if their non-zero patterns differ."""
    a, b = sp.csr_matrix(a), sp.csr_matrix(b)
    if ((a != 0) != (b != 0)).nnz:
        return float("inf")
    diff = (a.astype(np.float64) - b.astype(np.float64)).tocsr()
    return float(np.abs(diff.data).max()) if diff.nnz else 0.0


def extend_column(base_column: pd.Series, added: pd.Series):
    """Base values then added values, in the base column's dtype.

    Categorical columns keep their category order (positional uns colors depend on it) and append any new
    category; float columns take NaN for missing added values.
    """
    if isinstance(base_column.dtype, pd.CategoricalDtype):
        values = np.concatenate([base_column.astype(object).values, added.astype(object).values])
        categories = list(base_column.cat.categories)
        categories += sorted(set(added.dropna().astype(str)) - set(categories))
        return pd.Categorical(values, categories=categories)
    return np.concatenate([base_column.to_numpy(), added.to_numpy(dtype=base_column.dtype)])


def drop_excluded_clusters(base, settings, vocab, species):
    """Names of the base nuclei in the class clusters newly EXCLUDE in v8, and a record per cluster."""
    drop = pd.Index([])
    record = {}
    expected = {key: n for key, n in vocab["content_tasks"]["v8_exclusions"]["expected_counts"].items()
                if key.startswith(f"{species}_")}          # keys are <decision table>_cluster_<id>
    for spec, clusters_path, decisions_path in zip(settings["dropped"], snk.input.dropped_clusters,
                                                   snk.input.dropped_decisions):
        check_pinned(clusters_path, spec["clusters_sha256"])
        cluster_of = ad.read_h5ad(clusters_path, backed="r").obs[spec["cluster_column"]].astype(str)
        decisions = load_choices(decisions_path)
        excluded_now = [str(c) for c in spec["excluded"]]
        for cluster in excluded_now:
            if decisions.get(cluster) != EXCLUDE:
                raise ValueError(f"{decisions_path}: cluster {cluster} is not EXCLUDE")
            names = cluster_of.index[cluster_of.values == cluster]
            missing = names.difference(base.obs_names)
            if len(missing):
                raise ValueError(f"{len(missing)} nuclei of {clusters_path} cluster {cluster} are not in the base")
            key = f"{Path(decisions_path).stem}_cluster_{cluster}"
            if key not in expected or int(expected[key]) != len(names):
                raise ValueError(f"{key}: {len(names)} nuclei, vocabulary v8_exclusions expects {expected.get(key)}")
            record[key] = {"clustering": str(clusters_path), "cluster_column": str(spec["cluster_column"]),
                           "decisions": str(decisions_path), "cluster": cluster, "n_nuclei": int(len(names))}
            say(f"dropping {len(names):,} nuclei of {Path(clusters_path).name} cluster {cluster} "
                f"({base.obs.loc[names, 'cell_subtype'].astype(str).value_counts().to_dict()})")
            drop = drop.append(names)
        other_excluded = {c for c, label in decisions.items() if label == EXCLUDE} - set(excluded_now)
        stale = cluster_of.index[cluster_of.isin(other_excluded).values].intersection(base.obs_names)
        if len(stale):
            raise ValueError(f"{len(stale)} base nuclei sit in other EXCLUDE clusters of {clusters_path}")
    if not drop.is_unique:
        raise ValueError("a nucleus is listed by two dropped clusters")
    if sorted(record) != sorted(expected):
        raise ValueError(f"dropped clusters {sorted(record)} != vocabulary v8_exclusions {sorted(expected)}")
    return drop, record


def supplementary_rows(species):
    rows = [r for r in csv.DictReader(open(snk.input.supplementary)) if r["species"] == species]
    if not rows:
        raise ValueError(f"no supplementary populations for {species} in {snk.input.supplementary}")
    return rows


def legacy_value(base_obs, column, cell_type, fallback=None):
    """The single value the base nuclei of cell_type carry in a legacy column (or the fallback)."""
    values = set(base_obs.loc[base_obs["cell_type"].astype(str) == cell_type, column].astype(str))
    if len(values) == 1:
        return values.pop()
    if not values and fallback is not None:
        return fallback
    raise ValueError(f"legacy column {column} has values {sorted(values)} for {cell_type}; no rule")


def check_counts(obs, vocab, species):
    """Primary, extended, supplementary and dropped counts as recorded in the vocabulary."""
    primary = obs[obs["in_primary"].to_numpy(bool)]
    supplementary = obs[~obs["in_primary"].to_numpy(bool)]
    checks = [
        ("expected_coarse_counts", primary["cell_type"].astype(str).value_counts().to_dict(),
         dict(vocab["expected_coarse_counts"][species])),
        ("total primary", len(primary), vocab["expected_coarse_counts"]["totals"][species]),
        ("expected_coarse_counts_extended", obs["cell_type"].astype(str).value_counts().to_dict(),
         dict(vocab["expected_coarse_counts_extended"][species])),
        ("total extended", len(obs), vocab["expected_coarse_counts_extended"]["totals"][species]),
        ("supplementary_populations", supplementary["cell_subtype"].astype(str).value_counts().to_dict(),
         dict(vocab["content_tasks"]["supplementary_populations"]["expected_counts"][species])),
    ]
    split = {"human": "human_epithelial_split", "mouse": "mouse_epithelial_split"}[species]
    expected_split = vocab["content_tasks"][split]["expected_counts"]
    found_split = {label: int((primary["cell_subtype"].astype(str) == label).sum()) for label in expected_split}
    checks.append((split, found_split, dict(expected_split)))
    for name, found, expected in checks:
        found = {k: int(v) for k, v in found.items()} if isinstance(found, dict) else int(found)
        expected = {k: int(v) for k, v in expected.items()} if isinstance(expected, dict) else int(expected)
        if found != expected:
            raise ValueError(f"{name}: found {found}, vocabulary expects {expected}")
        say(f"  {name}: {found}")


def main():
    species = str(snk.params.species)
    settings = dict(snk.params.settings)
    vocab = yaml.safe_load(open(snk.input.vocab))

    base = ad.read_h5ad(snk.input.base)
    say(f"base {snk.input.base}: {base.n_obs:,} nuclei x {base.n_vars:,} genes")
    if not base.obs_names.is_unique:
        raise ValueError("base has duplicate nucleus names")

    # 1. The class clusters newly EXCLUDE in v8.
    drop, dropped_record = drop_excluded_clusters(base, settings, vocab, species)
    kept = base.obs_names[~base.obs_names.isin(drop)]          # base order
    say(f"keeping {len(kept):,} of {base.n_obs:,} base nuclei")

    # 2. The supplementary whole-tissue clusters.
    atlas_sha = check_pinned(snk.input.atlas, settings["atlas_sha256"])
    atlas = ad.read_h5ad(snk.input.atlas)
    if not atlas.obs_names.is_unique:
        raise ValueError("atlas has duplicate nucleus names")
    if not atlas.var_names.equals(base.var_names):
        raise ValueError("atlas genes differ from the base genes (names or order)")
    atlas_cluster = atlas.obs[settings["atlas_cluster_column"]].astype(str)
    lineage_decisions = load_choices(snk.input.lineage_decisions)

    base_counts = sp.csr_matrix(base.layers["counts"])
    atlas_counts = sp.csr_matrix(atlas.layers["counts"])
    missing = base.obs_names.difference(atlas.obs_names)
    if len(missing):
        raise ValueError(f"{len(missing)} base nuclei are not in the atlas")
    if (base_counts != atlas_counts[atlas.obs_names.get_indexer(base.obs_names)].astype(base_counts.dtype)).nnz:
        raise ValueError("base counts differ from the atlas counts")
    x_difference = max_difference(log1p_cp10k(base_counts), base.X)
    if x_difference > X_TOLERANCE:
        raise ValueError(f"base .X is not log1p(CP10K) of its counts (max difference {x_difference:.3g})")
    say(f"base counts equal the atlas counts; base .X is log1p(CP10K) of them within {x_difference:.2g}")

    timepoints = base.obs.groupby("sample", observed=True)["timepoint"].agg(lambda s: set(s.astype(str)))
    if (timepoints.map(len) != 1).any():
        raise ValueError("a sample has more than one timepoint in the base")
    timepoint_of = timepoints.map(lambda s: next(iter(s)))

    added_frames, added_names, report = [], [], []
    for row in supplementary_rows(species):
        cluster = str(row["atlas_cluster"])
        if lineage_decisions.get(cluster) != EXCLUDE:
            raise ValueError(f"{snk.input.lineage_decisions}: whole-tissue cluster {cluster} is not EXCLUDE")
        names = atlas.obs_names[atlas_cluster.values == cluster]
        if len(names) != int(row["n_nuclei"]):
            raise ValueError(f"whole-tissue cluster {cluster}: {len(names)} nuclei, table says {row['n_nuclei']}")
        if names.isin(base.obs_names).any():
            raise ValueError(f"whole-tissue cluster {cluster}: nuclei already in the base")
        frame = pd.DataFrame(index=names)
        frame["sample"] = atlas.obs.loc[names, "sample"].astype(str).values
        frame["cell_barcode"] = atlas.obs.loc[names, "cell_barcode"].astype(str).values
        frame["timepoint"] = frame["sample"].map(timepoint_of)
        if frame["timepoint"].isna().any():
            raise ValueError(f"whole-tissue cluster {cluster}: nuclei from a sample not in the base")
        frame["cell_type"] = row["cell_type"]
        frame["cell_subtype"] = row["cell_subtype"]
        frame["coarse_refined"] = legacy_value(base.obs, "coarse_refined", row["cell_type"],
                                               COARSE_REFINED_NEW.get(row["cell_type"]))
        frame["fine_refined"] = row["legacy_fine_name"]
        frame["origin_lineage"] = WHOLE_TISSUE
        if "lineage" in base.obs:
            frame["lineage"] = legacy_value(base.obs, "lineage", row["cell_type"])
        if "subtype" in base.obs:
            frame["subtype"] = row["legacy_fine_name"]
        if "annotation_level" in base.obs:
            frame["annotation_level"] = "subtype" if row["legacy_fine_name"] else "cell_type"
        frame["x_um_dbscan"] = atlas.obs.loc[names, "x_um_dbscan"].astype(float).values
        frame["y_um_dbscan"] = atlas.obs.loc[names, "y_um_dbscan"].astype(float).values
        frame["x_um_qc"] = np.nan
        frame["y_um_qc"] = np.nan
        frame["spatial_qc_reason"] = NOT_ASSESSED
        frame["tissue_section"] = np.nan
        added_frames.append(frame)
        added_names.append(names)
        report.append({"kind": "supplementary", "source": f"{Path(snk.input.atlas).name} "
                       f"{settings['atlas_cluster_column']}", "cluster": cluster, "label": row["cell_subtype"],
                       "n_nuclei": len(names)})
        say(f"adding {len(names):,} nuclei of whole-tissue cluster {cluster} as {row['cell_subtype']} "
            f"({int(frame['x_um_dbscan'].notna().sum())} with coordinates)")
    added = pd.concat(added_frames)
    added_index = added.index
    if not added_index.is_unique:
        raise ValueError("a nucleus is in two supplementary clusters")

    # 3. One object: kept base nuclei in base order, then the supplementary nuclei; in_primary marks them.
    columns = [c for c in base.obs.columns if c not in RECOMPUTED_OBS]
    for column in columns:
        if column not in added.columns:
            if pd.api.types.is_float_dtype(base.obs[column]):
                added[column] = np.nan             # values of other analyses (mouse auc_*), not computed here
            else:
                raise ValueError(f"no value for the supplementary nuclei in obs {column}")
    base_kept = base[kept]
    obs = pd.DataFrame(index=kept.append(added_index))
    for column in columns:
        obs[column] = extend_column(base_kept.obs[column], added[column])
    obs["in_primary"] = np.concatenate([np.ones(len(kept), bool), np.zeros(len(added_index), bool)])
    check_counts(obs, vocab, species)

    counts_added = atlas_counts[atlas.obs_names.get_indexer(added_index)].astype(base_counts.dtype)
    rows_kept = base.obs_names.get_indexer(kept)
    spatial_dtype = np.asarray(base.obsm["X_spatial"]).dtype
    uns = {key: value for key, value in base.uns.items() if key not in EMBEDDING_UNS}
    uns["annotation_schema_version"] = vocab["schema_version"]
    uns["display_names"] = vocab["display_names"]
    uns["release_rebuild"] = {
        "base": str(snk.input.base),
        "dropped": dropped_record,
        "n_dropped": int(len(drop)),
        "atlas": str(snk.input.atlas),
        "atlas_sha256": atlas_sha,
        "atlas_cluster_column": str(settings["atlas_cluster_column"]),
        "supplementary_table": str(snk.input.supplementary),
        "supplementary": {f"cluster_{r['atlas_cluster']}": {k: str(v) for k, v in r.items()}
                          for r in supplementary_rows(species)},
        "note": ("The supplementary populations (obs in_primary == False) were excluded before the class "
                 "re-integrations and are not part of any analysis of the paper; they carry raw Slide-tags "
                 "coordinates (x_um_dbscan, y_um_dbscan) but no spatial QC."),
    }
    out = ad.AnnData(
        X=sp.vstack([sp.csr_matrix(base.X)[rows_kept], log1p_cp10k(counts_added)], format="csr"),
        obs=obs,
        var=base.var.drop(columns=[c for c in HVG_VAR if c in base.var.columns]),
        uns=uns,
        obsm={"X_spatial": np.vstack([np.asarray(base.obsm["X_spatial"])[rows_kept],
                                      np.full((len(added_index), 2), np.nan, dtype=spatial_dtype)])},
        layers={"counts": sp.vstack([base_counts[rows_kept], counts_added], format="csr")},
    )
    validate_canonical_annotations(out, vocab, context=f"rebuilt {species} all-cells (extended)")
    say(f"dropped from the base: obsm {list(EMBEDDING_OBSM)}, uns {[k for k in EMBEDDING_UNS if k in base.uns]}, "
        f"var {list(HVG_VAR)}, obs {list(RECOMPUTED_OBS)}, obsp {list(base.obsp.keys())}")

    for entry in dropped_record.values():
        report.append({"kind": "dropped", "source": f"{Path(entry['clustering']).name} {entry['cluster_column']}",
                       "cluster": entry["cluster"], "label": EXCLUDE, "n_nuclei": entry["n_nuclei"]})
    Path(snk.output.report).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(report).to_csv(snk.output.report, sep="\t", index=False)

    out_path = Path(snk.output.h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp.h5ad")
    out.write_h5ad(tmp_path)
    os.replace(tmp_path, out_path)
    say(f"wrote {out.n_obs:,} nuclei ({int(out.obs['in_primary'].sum()):,} primary) x {out.n_vars:,} genes -> {out_path}")


main()
LOG.close()
