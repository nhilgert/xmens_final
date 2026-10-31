"""Switch a fibroblast/decidual release object to the 2026_10_29b whole-tissue embedding (vocabulary v8).

04g refresh_fibdec_release. The fibroblast/decidual objects carry the whole-tissue embedding of their nuclei
(human: obsm X_scVI_wholetissue and X_umap_wholetissue; mouse: obsm X_scVI and X_umap_scVI; config
release_rebuild fibdec_embedding) and embedding_outlier. This step replaces exactly those with the values of the
same nuclei in the finished primary all-cells object, rewrites uns['embedding_outlier_rule'] and records the
source in uns['embedding_provenance']; everything else is checked unchanged. The stromal embeddings, diffusion
maps and pseudotimes the object also carries are not touched.
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
from apply_annotation import validate_canonical_annotations  # noqa: E402

snk = snakemake  # noqa: F821 (injected by Snakemake)

RULE_RETENTION = "Flagged cells are retained in all analyses and omitted from embedding and spatial plots."
LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[refresh_fibdec:{snk.params.species}] {msg}", file=LOG, flush=True)


def same_values(a: pd.Series, b: pd.Series) -> bool:
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        return bool(np.array_equal(a.to_numpy(float), b.to_numpy(float), equal_nan=True))
    return bool((a.astype(object).where(a.notna(), "<NA>").astype(str).values
                 == b.astype(object).where(b.notna(), "<NA>").astype(str).values).all())


def main():
    release = str(snk.params.release)
    vocab = yaml.safe_load(open(snk.input.vocab))
    mapping = dict(snk.params.obsm_map)          # fib/dec obsm key -> all-cells obsm key
    fibdec = ad.read_h5ad(snk.input.fibdec)
    before = ad.read_h5ad(snk.input.fibdec, backed="r")
    primary = ad.read_h5ad(snk.input.primary, backed="r")
    say(f"fib/dec {fibdec.n_obs:,} nuclei; primary all-cells {primary.n_obs:,}")

    rows = primary.obs_names.get_indexer(fibdec.obs_names)
    if (rows < 0).any():
        raise ValueError(f"{int((rows < 0).sum())} fib/dec nuclei are not in the primary all-cells object")
    for column in ("cell_type", "cell_subtype", "sample", "cell_barcode"):
        if column in fibdec.obs and not same_values(fibdec.obs[column], primary.obs[column].iloc[rows]):
            raise ValueError(f"obs {column} differs between the fib/dec and all-cells objects")

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

    previous = str(fibdec.uns["embedding_outlier_rule"])
    if previous.count(RULE_RETENTION) != 1:
        raise ValueError("uns['embedding_outlier_rule'] does not have the expected wording")
    fibdec.uns["embedding_outlier_rule"] = (
        f"{previous.split(RULE_RETENTION)[0]}{RULE_RETENTION} {int(flags.sum())} cells flagged in this object: "
        f"the flags of the same nuclei in the {release} all-cells objects (refresh_fibdec_release).")
    provenance = dict(fibdec.uns.get("embedding_provenance", {}))
    for key, source in mapping.items():
        provenance[key] = (f"{release} whole-tissue scVI fit (04g reembed_release, integrate.py; primary and "
                           f"supplementary nuclei of the all-cells objects), obsm {source} row-subset to fib/dec "
                           f"(refresh_fibdec_release).")
    fibdec.uns["embedding_provenance"] = provenance
    fibdec.uns["annotation_schema_version"] = vocab["schema_version"]

    # Everything else unchanged.
    if not fibdec.obs_names.equals(before.obs_names) or list(fibdec.obs.columns) != list(before.obs.columns):
        raise ValueError("fib/dec nuclei or obs columns changed")
    for column in before.obs.columns:
        if column != "embedding_outlier" and not same_values(fibdec.obs[column], before.obs[column]):
            raise ValueError(f"obs {column} changed")
    if (sp.csr_matrix(fibdec.X) != sp.csr_matrix(before.X[:])).nnz:
        raise ValueError(".X changed")
    for layer in before.layers.keys():
        if (sp.csr_matrix(fibdec.layers[layer]) != sp.csr_matrix(before.layers[layer][:])).nnz:
            raise ValueError(f"layer {layer} changed")
    for key in before.obsm.keys():
        if key not in mapping and not np.array_equal(np.asarray(fibdec.obsm[key]), np.asarray(before.obsm[key]),
                                                     equal_nan=True):
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
