"""Build one species' all-cells object for release 2026_10_30 (vocabulary v9): round 2 of the lineage split.

04g2 rebuild_release_round2. Inputs (config release_round2 and annotation_round2):
  base        the species' all-cells object of release 2026_10_29b (its primary nuclei)
  atlas       data/interim/latent_{species}.h5ad, the July whole-tissue fit (pinned by sha256): raw counts, sample,
              barcode and raw coordinates of the nuclei new to the release
  latents     the round-2 lineage fits (pinned by sha256): every nucleus of the round-2 lineages
  lineages    their annotated objects (04c2 reannotate_lineage_round2): the kept nuclei with cell_type, cell_subtype
  spatial_qc  04c2 spatial_qc_round2: the spatial QC of every nucleus the round-2 atlas table keeps
  legacy      resources/annotations/round2/legacy_names.csv: the legacy label columns of each v9 label
  vocab       resources/annotation_vocabulary.yaml (v9)

This step
  1. replaces the base nuclei of the round-2 lineages by the kept nuclei of their round-2 annotation. Base nuclei
     outside these lineages keep their place and every obs value; kept round-2 nuclei already in the base keep their
     place; nuclei new to the release follow, lineage by lineage in config order, in the annotated object's order;
  2. labels every round-2 nucleus from its annotated object (cell_type, cell_subtype) and sets its legacy columns
     (coarse_refined, fine_refined, origin_lineage = its lineage; human also lineage, subtype, annotation_level) from
     the legacy table. Nuclei new to the release take sample, cell_barcode and raw coordinates from the atlas, the
     timepoint their sample has in the base, raw counts from the atlas with .X = log1p(CP10K) of them, and missing
     values in the float columns of other analyses (mouse auc_*);
  3. sets the spatial QC of every nucleus (spatial_qc_reason, x_um_qc, y_um_qc, tissue_section, obsm X_spatial) from
     the round-2 QC and reports the base nuclei whose QC changed;
  4. removes the whole-tissue embedding (X_scVI, X_umap_scVI, neighbour/UMAP/HVG metadata, scvi_fit, HVG flags) and
     embedding_outlier, which reembed_release_round2 and finish_release_round2 recompute for all nuclei.

It fails unless: the atlas and the round-2 fits match their pinned sha256; the base nuclei of the round-2 lineages
(obs origin_lineage) are exactly the base nuclei in the round-2 fits; every annotated nucleus is in its fit and in
the atlas; the base counts equal the atlas counts; the legacy table reproduces the base values of every nucleus whose
label is unchanged; every nucleus has a round-2 spatial QC; and the counts and labels pass the vocabulary
(expected_coarse_counts, content tasks round2_lineages, coarse_subtype_fill and the epithelial splits).
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
from apply_annotation import validate_canonical_annotations  # noqa: E402

snk = snakemake  # noqa: F821 (injected by Snakemake)

EMBEDDING_OBSM = ("X_scVI", "X_umap_scVI")
EMBEDDING_UNS = ("neighbors_scVI_100", "umap", "hvg", "log1p", "scvi_fit")
HVG_VAR = ("highly_variable", "highly_variable_rank", "highly_variable_nbatches")
RECOMPUTED_OBS = ("embedding_outlier",)
SPATIAL_OBS = ("spatial_qc_reason", "x_um_qc", "y_um_qc", "tissue_section")
RETIRED_UNS = ("release_rebuild",)     # the v8 record (extended objects); this release records release_round2
LEGACY_TABLE_COLUMNS = ("coarse_refined", "fine_refined", "subtype", "annotation_level")
X_TOLERANCE = 1e-6                     # as rebuild_cells.py: July mouse .X differs from a recomputation by <= 3.4e-7

LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[rebuild_round2:{snk.params.species}] {msg}", file=LOG, flush=True)


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


def as_column(values, like: pd.Series):
    """values in the dtype of a base column; categorical columns keep their order and append new categories."""
    if isinstance(like.dtype, pd.CategoricalDtype):
        values = pd.Series(values, dtype=object)
        categories = list(like.cat.categories)
        categories += sorted(set(values.dropna().astype(str)) - set(categories))
        return pd.Categorical(values.where(values.notna(), None), categories=categories)
    return np.asarray(values, dtype=like.dtype)


def same(a, b) -> np.ndarray:
    """Element-wise equality, NaN matching NaN, compared as strings unless both are numeric."""
    a, b = pd.Series(a).reset_index(drop=True), pd.Series(b).reset_index(drop=True)
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        x, y = a.to_numpy(float), b.to_numpy(float)
        return (x == y) | (np.isnan(x) & np.isnan(y))
    return (a.astype(object).where(a.notna(), "<NA>").astype(str).values
            == b.astype(object).where(b.notna(), "<NA>").astype(str).values)


def legacy_table(path, species):
    rows = [r for r in csv.DictReader(open(path)) if r["species"] == species]
    table = {r["cell_subtype"]: r for r in rows}
    if len(table) != len(rows):
        raise ValueError(f"{path}: duplicate labels for {species}")
    return table


def round2_nuclei(settings_lineages, species):
    """Per lineage: all nuclei of the round-2 fit, and the annotated (kept) nuclei with their labels."""
    lineages = {}
    for lineage, latent_path, annotated_path in zip(settings_lineages, snk.input.latents, snk.input.lineages):
        check_pinned(latent_path, snk.params.latent_sha256[lineage])
        fit = ad.read_h5ad(latent_path, backed="r")
        all_names = fit.obs_names.copy()
        fit.file.close()
        annotated = ad.read_h5ad(annotated_path, backed="r")
        labels = annotated.obs[["cell_type", "cell_subtype"]].astype(str).copy()
        annotated.file.close()
        if not labels.index.isin(all_names).all():
            raise ValueError(f"{annotated_path}: nuclei that are not in {latent_path}")
        lineages[lineage] = {"all": all_names, "kept": labels}
        say(f"{lineage}: {len(all_names):,} nuclei in the round-2 fit, {len(labels):,} kept "
            f"({len(all_names) - len(labels):,} EXCLUDE)")
    names = pd.Index(np.concatenate([v["all"].values for v in lineages.values()]))
    if not names.is_unique:
        raise ValueError("a nucleus is in two round-2 lineages")
    return lineages


def check_counts(obs, vocab, species, lineages):
    """Counts against the vocabulary: coarse counts, round-2 label counts, coarse fills and epithelial splits."""
    checks = [("expected_coarse_counts", obs["cell_type"].astype(str).value_counts().to_dict(),
               dict(vocab["expected_coarse_counts"][species])),
              ("total", len(obs), vocab["expected_coarse_counts"]["totals"][species])]
    expected_r2 = vocab["content_tasks"]["round2_lineages"]["expected_counts"]
    for lineage, entry in lineages.items():
        found = entry["kept"]["cell_subtype"].value_counts().to_dict()
        found["EXCLUDE"] = len(entry["all"]) - len(entry["kept"])
        checks.append((f"round2_lineages {species}_{lineage}", found, dict(expected_r2[f"{species}_{lineage}"])))
    fills = vocab["content_tasks"]["coarse_subtype_fill"]
    total_fills = 0
    for entry in fills["applies_to"][species]:
        cell_type = entry["cell_type"]
        n = int(((obs["cell_type"].astype(str) == cell_type) & (obs["cell_subtype"].astype(str) == cell_type)).sum())
        total_fills += n
        checks.append((f"coarse fill {cell_type}", n, entry["expected_coarse_fill_count"]))
    checks.append(("coarse fills, total", total_fills, fills["expected_total_fills"][species]))
    split = {"human": "human_epithelial_split", "mouse": "mouse_epithelial_split"}[species]
    for label, n_expected in vocab["content_tasks"][split]["expected_counts"].items():
        checks.append((f"{split} {label}", int((obs["cell_subtype"].astype(str) == label).sum()), n_expected))
    if species == "mouse":
        for label, n_expected in vocab["content_tasks"]["mouse_endothelial_blood_only"]["expected_counts"].items():
            column = "cell_type" if label == "endothelial_cell" else "cell_subtype"
            checks.append((f"mouse_endothelial_blood_only {label}", int((obs[column].astype(str) == label).sum()),
                           n_expected))
    for name, found, expected in checks:
        found = {k: int(v) for k, v in found.items()} if isinstance(found, dict) else int(found)
        expected = {k: int(v) for k, v in expected.items()} if isinstance(expected, dict) else int(expected)
        if found != expected:
            raise ValueError(f"{name}: found {found}, vocabulary expects {expected}")
        say(f"  {name}: {found}")


def main():
    species = str(snk.params.species)
    vocab = yaml.safe_load(open(snk.input.vocab))
    legacy = legacy_table(snk.input.legacy, species)

    base = ad.read_h5ad(snk.input.base)
    say(f"base {snk.input.base}: {base.n_obs:,} nuclei x {base.n_vars:,} genes")
    if not base.obs_names.is_unique:
        raise ValueError("base has duplicate nucleus names")

    # 1. Round-2 lineages: which base nuclei go, which stay, which are new.
    lineages = round2_nuclei(list(snk.params.lineages), species)
    in_round2 = pd.Index(np.concatenate([v["all"].values for v in lineages.values()]))
    base_round2 = base.obs_names[base.obs_names.isin(in_round2)]
    by_origin = base.obs_names[base.obs["origin_lineage"].astype(str).isin(list(lineages)).to_numpy()]
    if set(base_round2) != set(by_origin):
        raise ValueError(f"base nuclei in the round-2 fits ({len(base_round2):,}) are not the base nuclei of the "
                         f"round-2 lineages by origin_lineage ({len(by_origin):,})")
    kept_labels = pd.concat([v["kept"] for v in lineages.values()])
    lineage_of = pd.Series({n: lin for lin, v in lineages.items() for n in v["kept"].index})
    removed = base_round2[~base_round2.isin(kept_labels.index)]
    keep_base = base.obs_names[~base.obs_names.isin(removed)]                        # base order
    new_names = kept_labels.index[~kept_labels.index.isin(base.obs_names)]            # lineage order
    names = keep_base.append(new_names)
    say(f"base: {len(base_round2):,} nuclei in the round-2 lineages, {len(removed):,} of them EXCLUDE in round 2; "
        f"{len(new_names):,} nuclei new to the release; {len(names):,} nuclei in total")

    # 2. Counts and .X.
    atlas_sha = check_pinned(snk.input.atlas, snk.params.atlas_sha256)
    atlas = ad.read_h5ad(snk.input.atlas)
    if not atlas.var_names.equals(base.var_names):
        raise ValueError("atlas genes differ from the base genes (names or order)")
    missing = names.difference(atlas.obs_names)
    if len(missing):
        raise ValueError(f"{len(missing)} nuclei are not in the atlas")
    base_counts = sp.csr_matrix(base.layers["counts"])
    atlas_counts = sp.csr_matrix(atlas.layers["counts"])
    if (base_counts != atlas_counts[atlas.obs_names.get_indexer(base.obs_names)].astype(base_counts.dtype)).nnz:
        raise ValueError("base counts differ from the atlas counts")
    x_difference = max_difference(log1p_cp10k(base_counts), base.X)
    if x_difference > X_TOLERANCE:
        raise ValueError(f"base .X is not log1p(CP10K) of its counts (max difference {x_difference:.3g})")
    rows_base = base.obs_names.get_indexer(keep_base)
    counts_new = atlas_counts[atlas.obs_names.get_indexer(new_names)].astype(base_counts.dtype)
    counts = sp.vstack([base_counts[rows_base], counts_new], format="csr")
    x = sp.vstack([sp.csr_matrix(base.X)[rows_base], log1p_cp10k(counts_new).astype(base.X.dtype)], format="csr")
    say(f"base counts equal the atlas counts; .X of the base within {x_difference:.2g} of a recomputation")

    # 3. obs: base values, round-2 labels and legacy columns, values of the new nuclei, spatial QC.
    columns = [c for c in base.obs.columns if c not in RECOMPUTED_OBS]
    frame = pd.DataFrame(index=names)
    timepoint_of = base.obs.groupby("sample", observed=True)["timepoint"].agg(lambda s: set(s.astype(str)))
    if (timepoint_of.map(len) != 1).any():
        raise ValueError("a sample has more than one timepoint in the base")
    timepoint_of = timepoint_of.map(lambda s: next(iter(s)))
    new_obs = atlas.obs.loc[new_names]
    for column in columns:
        values = pd.Series(index=names, dtype=object)
        values.loc[keep_base] = base.obs.loc[keep_base, column].astype(object).values
        if column in ("sample", "cell_barcode", "x_um_dbscan", "y_um_dbscan"):
            values.loc[new_names] = new_obs[column].values
        elif column == "timepoint":
            values.loc[new_names] = new_obs["sample"].astype(str).map(timepoint_of).values
        elif pd.api.types.is_float_dtype(base.obs[column]) and column not in SPATIAL_OBS:
            values.loc[new_names] = np.nan                  # values of other analyses (mouse auc_*), not computed here
        frame[column] = values
    if frame.loc[new_names, "timepoint"].isna().any():
        raise ValueError("nuclei new to the release come from a sample that is not in the base")

    # labels and legacy columns of every round-2 nucleus
    r2 = kept_labels.index
    before = base.obs.reindex(r2)
    frame.loc[r2, "cell_type"] = kept_labels["cell_type"].values
    frame.loc[r2, "cell_subtype"] = kept_labels["cell_subtype"].values
    unknown = sorted(set(kept_labels["cell_subtype"]) - set(legacy))
    if unknown:
        raise ValueError(f"{snk.input.legacy}: no legacy names for {species} labels {unknown}")
    for column in [c for c in LEGACY_TABLE_COLUMNS if c in columns]:
        frame.loc[r2, column] = kept_labels["cell_subtype"].map({k: v[column] for k, v in legacy.items()}).values
    frame.loc[r2, "origin_lineage"] = lineage_of.reindex(r2).values
    if "lineage" in columns:
        frame.loc[r2, "lineage"] = lineage_of.reindex(r2).values
    unchanged = r2[r2.isin(base.obs_names)
                   & (before["cell_subtype"].astype(str).values == kept_labels["cell_subtype"].values)]
    legacy_columns = [c for c in (*LEGACY_TABLE_COLUMNS, "lineage", "origin_lineage") if c in columns]
    for column in legacy_columns:
        mismatch = ~same(frame.loc[unchanged, column], base.obs.loc[unchanged, column])
        if mismatch.any():
            example = unchanged[np.where(mismatch)[0][0]]
            raise ValueError(f"legacy column {column}: {int(mismatch.sum())} nuclei with an unchanged label get a "
                             f"value other than the base's (e.g. {example}: {frame.at[example, column]!r} vs "
                             f"{base.obs.at[example, column]!r})")
    relabelled = r2[r2.isin(base.obs_names)
                    & (before["cell_subtype"].astype(str).values != kept_labels["cell_subtype"].values)]
    say(f"labels: {len(unchanged):,} round-2 nuclei of the base keep their label (legacy columns reproduce the "
        f"base); {len(relabelled):,} change label; {len(new_names):,} are new")

    # spatial QC of every nucleus
    qc = ad.read_h5ad(snk.input.spatial_qc, backed="r")
    qc_obs = qc.obs[list(SPATIAL_OBS)].copy()
    qc_uns = dict(qc.uns.get("spatial_qc", {}))
    qc.file.close()
    without_qc = names.difference(qc_obs.index)
    if len(without_qc):
        raise ValueError(f"{len(without_qc)} nuclei have no round-2 spatial QC")
    for column in SPATIAL_OBS:
        frame[column] = qc_obs.loc[names, column].astype(object).values
    changed_qc = keep_base[~np.logical_and.reduce(
        [same(frame.loc[keep_base, c], base.obs.loc[keep_base, c]) for c in SPATIAL_OBS])]
    say(f"spatial QC: {len(changed_qc):,} base nuclei change (reasons "
        f"{pd.crosstab(base.obs.loc[changed_qc, 'spatial_qc_reason'].astype(str), frame.loc[changed_qc, 'spatial_qc_reason'].astype(str)).stack().loc[lambda s: s > 0].to_dict()})")

    for column in columns:
        if column in ("tissue_section",) or pd.api.types.is_float_dtype(base.obs[column]):
            continue                                        # missing by design: no section without a kept position
        if frame.loc[new_names, column].isna().any():
            raise ValueError(f"obs {column} has no value for {int(frame.loc[new_names, column].isna().sum())} new nuclei")
    obs = pd.DataFrame(index=names)
    for column in columns:
        obs[column] = as_column(frame[column], base.obs[column])
    check_counts(obs, vocab, species, lineages)

    # 4. The object, without the embedding.
    xy = obs[["x_um_qc", "y_um_qc"]].to_numpy(np.asarray(base.obsm["X_spatial"]).dtype)
    uns = {k: v for k, v in base.uns.items() if k not in EMBEDDING_UNS + RETIRED_UNS}
    for key in [k for k in uns if k.endswith("_colors")]:
        column = key[: -len("_colors")]
        if column in obs and isinstance(obs[column].dtype, pd.CategoricalDtype) \
                and len(obs[column].cat.categories) != len(uns[key]):
            del uns[key]
            say(f"uns {key} dropped: {column} now has {len(obs[column].cat.categories)} categories")
    uns["annotation_schema_version"] = vocab["schema_version"]
    uns["display_names"] = vocab["display_names"]
    uns["release_round2"] = {
        "release": str(snk.params.release),
        "base": str(snk.input.base),
        "atlas": str(snk.input.atlas),
        "atlas_sha256": atlas_sha,
        "lineages": {lin: {"fit": str(p), "fit_sha256": str(snk.params.latent_sha256[lin]), "annotated": str(a),
                           "n_fit": int(len(lineages[lin]["all"])), "n_kept": int(len(lineages[lin]["kept"]))}
                     for lin, p, a in zip(lineages, snk.input.latents, snk.input.lineages)},
        "n_base_removed": int(len(removed)),
        "n_relabelled": int(len(relabelled)),
        "n_new": int(len(new_names)),
        "legacy_names": str(snk.input.legacy),
        "spatial_qc": {"source": str(snk.input.spatial_qc), "n_base_changed": int(len(changed_qc)),
                       **{k: str(v) for k, v in qc_uns.items()}},
        "note": ("Round 2 of the lineage split (vocabulary v9): the nuclei of the round-2 lineages are those kept by "
                 "their round-2 annotation; spatial QC re-run over the nuclei the round-2 atlas tables keep."),
    }
    out = ad.AnnData(X=x, obs=obs, var=base.var.drop(columns=[c for c in HVG_VAR if c in base.var.columns]),
                     uns=uns, obsm={"X_spatial": xy}, layers={"counts": counts})
    validate_canonical_annotations(out, vocab, context=f"rebuilt {species} all-cells (v9)")

    report = []
    for lin, entry in lineages.items():
        gone = removed[removed.isin(entry["all"])]
        report += [{"kind": "removed", "lineage": lin, "from": k, "to": "EXCLUDE", "n_nuclei": int(v)}
                   for k, v in base.obs.loc[gone, "cell_subtype"].astype(str).value_counts().items()]
        added = new_names[new_names.isin(entry["kept"].index)]
        report += [{"kind": "added", "lineage": lin, "from": "not in base", "to": k, "n_nuclei": int(v)}
                   for k, v in obs.loc[added, "cell_subtype"].astype(str).value_counts().items()]
        moved = relabelled[relabelled.isin(entry["kept"].index)]
        pairs = pd.crosstab(base.obs.loc[moved, "cell_subtype"].astype(str), obs.loc[moved, "cell_subtype"].astype(str))
        report += [{"kind": "relabelled", "lineage": lin, "from": a, "to": b, "n_nuclei": int(n)}
                   for (a, b), n in pairs.stack().items() if n > 0]
    pairs = pd.crosstab(base.obs.loc[changed_qc, "spatial_qc_reason"].astype(str),
                        obs.loc[changed_qc, "spatial_qc_reason"].astype(str))
    report += [{"kind": "spatial_qc_changed", "lineage": "", "from": a, "to": b, "n_nuclei": int(n)}
               for (a, b), n in pairs.stack().items() if n > 0]
    Path(snk.output.report).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(report).to_csv(snk.output.report, sep="\t", index=False)

    out_path = Path(snk.output.h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp.h5ad")
    out.write_h5ad(tmp_path)
    os.replace(tmp_path, out_path)
    say(f"wrote {out.n_obs:,} nuclei x {out.n_vars:,} genes -> {out_path}")


main()
LOG.close()
