"""Switch a fibroblast/decidual object to the whole-tissue embedding and spatial QC of the release being cut.

04g2 refresh_fibdec_round2. The fibroblast/decidual objects carry, for their nuclei, the whole-tissue embedding
(human: obsm X_scVI_wholetissue and X_umap_wholetissue; mouse: obsm X_scVI and X_umap_scVI; config release_round2
fibdec_embedding), embedding_outlier and the spatial QC (spatial_qc_reason, x_um_qc, y_um_qc, tissue_section, obsm
X_spatial). This step replaces exactly those with the values of the same nuclei in the finished all-cells object,
rewrites uns['embedding_outlier_rule'], records the sources in uns['embedding_provenance'] and uns['spatial_qc_round2'],
sets uns['display_names'] to the vocabulary's whole map (the v9 fib/dec objects lacked the immune-only entries)
and checks that everything else is unchanged. When the release reuses the previous release's whole-tissue embedding
(params.embedding), the embedding and flags written here must equal those of the previous release's fib/dec object,
nucleus by nucleus. A species with a pseudotime of record (params.pseudotime, config pseudotime_of_record; decision G)
takes that run's pseudotime, diffusion map and record (pseudotime_of_record.py); otherwise the stromal embeddings,
diffusion maps and pseudotimes are not touched.
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
from pseudotime_of_record import (  # noqa: E402
    DIFFMAP_OBSM,
    RANK_OBS,
    RAW_OBS,
    RECORD_UNS,
    apply_pseudotime_of_record,
)

snk = snakemake  # noqa: F821 (injected by Snakemake)

SPATIAL_OBS = ("spatial_qc_reason", "x_um_qc", "y_um_qc", "tissue_section")
RULE_RETENTION = "Flagged cells are retained in all analyses and omitted from embedding and spatial plots."
LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[refresh_fibdec_round2:{snk.params.species}] {msg}", file=LOG, flush=True)


def same_values(a: pd.Series, b: pd.Series) -> np.ndarray:
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        x, y = a.to_numpy(float), b.to_numpy(float)
        return (x == y) | (np.isnan(x) & np.isnan(y))
    return (a.astype(object).where(a.notna(), "<NA>").astype(str).values
            == b.astype(object).where(b.notna(), "<NA>").astype(str).values)


def like(values: pd.Series, column: pd.Series):
    """values in the dtype of an existing column (categorical: its categories, new ones appended)."""
    if isinstance(column.dtype, pd.CategoricalDtype):
        values = values.astype(object)
        categories = list(column.cat.categories)
        categories += sorted(set(values.dropna().astype(str)) - set(categories))
        return pd.Categorical(values.where(values.notna(), None).values, categories=categories)
    return values.to_numpy(dtype=column.dtype)


def check_against_previous(fibdec, mapping: dict, previous_release: str):
    """A reused embedding: the fib/dec embedding and flags equal the previous release's fib/dec object by nucleus."""
    previous = ad.read_h5ad(snk.input.previous[0], backed="r")
    rows = previous.obs_names.get_indexer(fibdec.obs_names)
    if (rows < 0).any() or previous.n_obs != fibdec.n_obs:
        raise ValueError(f"embedding reuse: the fib/dec nuclei differ from those of {previous_release}")
    for key in mapping:
        if not np.array_equal(np.asarray(fibdec.obsm[key]), np.asarray(previous.obsm[key])[rows], equal_nan=True):
            raise ValueError(f"embedding reuse: obsm {key} differs from the {previous_release} fib/dec object")
    previous_flags = previous.obs["embedding_outlier"].astype(bool).to_numpy()[rows]
    if not np.array_equal(fibdec.obs["embedding_outlier"].astype(bool).to_numpy(), previous_flags):
        raise ValueError(f"embedding reuse: embedding_outlier differs from the {previous_release} fib/dec object")
    previous.file.close()
    say(f"embedding reuse: obsm {list(mapping)} and embedding_outlier equal the {previous_release} fib/dec object")


def main():
    release = str(snk.params.release)
    embedding = str(snk.params.embedding)
    previous_release = str(snk.params.previous_release)
    vocab = yaml.safe_load(open(snk.input.vocab))
    mapping = dict(snk.params.obsm_map)          # fib/dec obsm key -> all-cells obsm key
    fibdec = ad.read_h5ad(snk.input.fibdec)
    before = ad.read_h5ad(snk.input.fibdec, backed="r")
    primary = ad.read_h5ad(snk.input.primary, backed="r")
    say(f"fib/dec {fibdec.n_obs:,} nuclei; all-cells {primary.n_obs:,}")

    rows = primary.obs_names.get_indexer(fibdec.obs_names)
    if (rows < 0).any():
        raise ValueError(f"{int((rows < 0).sum())} fib/dec nuclei are not in the all-cells object")
    for column in ("cell_type", "cell_subtype", "sample", "cell_barcode"):
        if column in fibdec.obs and not same_values(fibdec.obs[column], primary.obs[column].iloc[rows]).all():
            raise ValueError(f"obs {column} differs between the fib/dec and all-cells objects")

    # Whole-tissue embedding and outlier flags.
    for key, source in mapping.items():
        if key not in fibdec.obsm:
            raise ValueError(f"fib/dec object has no obsm {key}")
        new = np.asarray(primary.obsm[source])[rows]
        if new.shape != np.asarray(fibdec.obsm[key]).shape:
            raise ValueError(f"obsm {key}: shape {new.shape} != {np.asarray(fibdec.obsm[key]).shape}")
        fibdec.obsm[key] = new.astype(np.asarray(fibdec.obsm[key]).dtype)
        say(f"obsm {key} <- all-cells {source}")
    flags_before = fibdec.obs["embedding_outlier"].astype(bool).to_numpy()
    flags = primary.obs["embedding_outlier"].astype(bool).to_numpy()[rows]
    fibdec.obs["embedding_outlier"] = flags
    say(f"embedding_outlier: {int(flags.sum())} flagged (was {int(flags_before.sum())}; "
        f"{int((flags & flags_before).sum())} in both)")
    if embedding == "reuse":
        check_against_previous(fibdec, mapping, previous_release)

    # Spatial QC.
    spatial = [c for c in SPATIAL_OBS if c in fibdec.obs]
    changed = np.zeros(fibdec.n_obs, bool)
    for column in spatial:
        new = primary.obs[column].iloc[rows].reset_index(drop=True)
        changed |= ~same_values(new, fibdec.obs[column].reset_index(drop=True))
        fibdec.obs[column] = like(new, fibdec.obs[column])
    spatial_dtype = np.asarray(fibdec.obsm["X_spatial"]).dtype
    fibdec.obsm["X_spatial"] = np.asarray(primary.obsm["X_spatial"])[rows].astype(spatial_dtype)
    moved = pd.crosstab(before.obs["spatial_qc_reason"].astype(str).values[changed],
                        fibdec.obs["spatial_qc_reason"].astype(str).values[changed])
    say(f"spatial QC: {int(changed.sum())} nuclei changed ({moved.stack().loc[lambda s: s > 0].to_dict()})")

    previous = str(fibdec.uns["embedding_outlier_rule"])
    if previous.count(RULE_RETENTION) != 1:
        raise ValueError("uns['embedding_outlier_rule'] does not have the expected wording")
    fibdec.uns["embedding_outlier_rule"] = (
        f"{previous.split(RULE_RETENTION)[0]}{RULE_RETENTION} {int(flags.sum())} cells flagged in this object: "
        f"the flags of the same nuclei in the {release} all-cells objects (refresh_fibdec_round2).")
    provenance = dict(fibdec.uns.get("embedding_provenance", {}))
    if embedding == "reuse":
        fit = (f"{previous_release} whole-tissue scVI fit, reused by {release} because its nuclei are the same "
               "(04g2 finish_release_round2)")
    else:
        fit = (f"{release} whole-tissue scVI fit (04g2 reembed_release_round2, integrate.py; all nuclei of the "
               "all-cells objects)")
    for key, source in mapping.items():
        provenance[key] = f"{fit}, obsm {source} row-subset to fib/dec (refresh_fibdec_round2)."
    fibdec.uns["embedding_provenance"] = provenance
    fibdec.uns["spatial_qc_round2"] = {
        "source": f"{release} all-cells object (04c2 spatial_qc_round2: spatial QC over the nuclei the round-2 atlas "
                  f"tables keep)",
        "columns": ", ".join(spatial + ["obsm X_spatial"]),
        "n_changed": str(int(changed.sum())),
    }
    fibdec.uns["annotation_schema_version"] = vocab["schema_version"]
    fibdec.uns["display_names"] = vocab["display_names"]

    # Decision G: the pseudotime of record replaces the release's pseudotime, diffusion map and record.
    pseudotime = dict(snk.params.pseudotime)
    replaced_obs, replaced_obsm = [], []
    if pseudotime:
        summary = apply_pseudotime_of_record(fibdec, pseudotime)
        replaced_obs, replaced_obsm = [RANK_OBS, RAW_OBS], [DIFFMAP_OBSM]
        say(f"pseudotime of record {pseudotime['path']}: obs {replaced_obs}, obsm {replaced_obsm} and uns "
            f"{RECORD_UNS} replaced for {summary['n_cells']:,} nuclei ({summary['n_dcs']} diffusion components; "
            f"root {summary['root']} at the lowest rank; Spearman with the release's rank "
            f"{summary['spearman_with_release']:.4f})")

    # Everything else unchanged.
    if not fibdec.obs_names.equals(before.obs_names) or list(fibdec.obs.columns) != list(before.obs.columns):
        raise ValueError("fib/dec nuclei or obs columns changed")
    for column in before.obs.columns:
        if column != "embedding_outlier" and column not in spatial + replaced_obs \
                and not same_values(fibdec.obs[column], before.obs[column]).all():
            raise ValueError(f"obs {column} changed")
    if (sp.csr_matrix(fibdec.X) != sp.csr_matrix(before.X[:])).nnz:
        raise ValueError(".X changed")
    for layer in before.layers.keys():
        if (sp.csr_matrix(fibdec.layers[layer]) != sp.csr_matrix(before.layers[layer][:])).nnz:
            raise ValueError(f"layer {layer} changed")
    for key in before.obsm.keys():
        if key not in mapping and key not in ["X_spatial", *replaced_obsm] and not np.array_equal(
                np.asarray(fibdec.obsm[key]), np.asarray(before.obsm[key]), equal_nan=True):
            raise ValueError(f"obsm {key} changed")
    if not fibdec.var.equals(before.var.copy()):
        raise ValueError("var changed")
    validate_canonical_annotations(fibdec, vocab, context=f"{release} {snk.params.species} fib/dec")

    out_path = Path(snk.output.h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp.h5ad")
    fibdec.write_h5ad(tmp_path)
    os.replace(tmp_path, out_path)
    say(f"wrote {fibdec.n_obs:,} nuclei -> {out_path}")


main()
LOG.close()
