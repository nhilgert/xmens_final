"""Extract the per-nucleus plot table behind the human whole-tissue panels (Fig. 5B, 5C, S25, S26).

Called by rule `panels_human_tissue_extract`, once per release (v8, v9).
Reads:  one release object (Cevrim_Human_All_Cells.h5ad) and the release MANIFEST.json; the expression
        matrix is never loaded whole, the counts layer is streamed in row blocks.
Writes: cells.parquet   one row per nucleus: obs, scVI UMAP, total counts, detected genes and
                        log1p(CP10k) of the marker genes
        uns.json        the release palettes and display vocabulary (label_colors, timepoint_colors, ...)
        provenance.json input sha256 (checked against the manifest), counts, script hash
"""

import argparse
import datetime
import hashlib
import json
import logging
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


# Rows per block when streaming the CSR counts layer; ~8k rows keeps one block well under 1 GB.
BLOCK_ROWS = 8192

# Library-size normalization of the S25 dot plot and of the existing marker dot plots
# (xmens_figures/celltype_marker_dotplots): log(1 + counts per 10,000).
TARGET_SUM = 1e4

# uns entries copied to uns.json: palettes, display vocabulary and the release notes the panels cite.
UNS_KEYS = [
    "label_colors",
    "timepoint_colors",
    "annotation_level_colors",
    "display_names",
    "label_relabel",
    "annotation_schema_version",
    "embedding_outlier_rule",
    "spatial_qc_push",
    "epithelial_reinstatement",
    "release_round2",
    "release_rebuild",
]

UMAP_KEY = "X_umap_scVI"

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("extract_human_tissue")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def to_jsonable(value):
    """Convert the nested numpy/pandas values of an AnnData uns entry into plain JSON types."""
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return to_jsonable(value.tolist())
    if isinstance(value, pd.DataFrame):
        return to_jsonable(value.reset_index().to_dict(orient="list"))
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, bytes):
        return value.decode()
    return value


def read_uns(h5: h5py.File) -> dict:
    uns = {}
    for key in UNS_KEYS:
        if key in h5["uns"]:
            uns[key] = to_jsonable(anndata.io.read_elem(h5["uns"][key]))
    return uns


def stream_counts(h5: h5py.File, gene_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Total counts, detected genes and marker-gene counts per nucleus from the CSR counts layer.

    Returns (total_counts [n_obs], n_genes [n_obs], marker_counts [n_obs, n_markers]).
    """
    layer = h5["layers"]["counts"]
    encoding = layer.attrs["encoding-type"]
    if encoding != "csr_matrix":
        raise ValueError(f"layers/counts is {encoding}; this reader streams csr_matrix only")

    n_obs, n_vars = layer.attrs["shape"]
    indptr = layer["indptr"][:]
    total_counts = np.zeros(n_obs, dtype=np.float64)
    n_genes = np.zeros(n_obs, dtype=np.int64)
    marker_counts = np.zeros((n_obs, len(gene_idx)), dtype=np.float32)

    for start in range(0, n_obs, BLOCK_ROWS):
        stop = min(start + BLOCK_ROWS, n_obs)
        lo, hi = indptr[start], indptr[stop]
        block = sp.csr_matrix(
            (layer["data"][lo:hi], layer["indices"][lo:hi], indptr[start : stop + 1] - lo),
            shape=(stop - start, n_vars),
        )

        # Stored zeros would count as detected genes, so they are dropped first.
        block.eliminate_zeros()
        total_counts[start:stop] = np.asarray(block.sum(axis=1)).ravel()
        n_genes[start:stop] = np.diff(block.indptr)
        marker_counts[start:stop] = block[:, gene_idx].toarray()

    return total_counts, n_genes, marker_counts


def build_cells_table(h5: h5py.File, genes: list[str]) -> tuple[pd.DataFrame, list[str]]:
    obs = anndata.io.read_elem(h5["obs"])
    var_names = anndata.io.read_elem(h5["var"]).index.astype(str)

    # Marker genes absent from the object are reported, not silently dropped from the dot plot.
    missing = [gene for gene in genes if gene not in var_names]
    present = [gene for gene in genes if gene in var_names]
    gene_idx = np.array([var_names.get_loc(gene) for gene in present])
    log.info(f"{len(obs):,} nuclei, {len(var_names):,} genes; markers found {len(present)}, missing {missing}")

    umap = anndata.io.read_elem(h5["obsm"][UMAP_KEY])
    total_counts, n_genes, marker_counts = stream_counts(h5, gene_idx)

    cells = obs.copy()
    cells["umap_1"] = umap[:, 0]
    cells["umap_2"] = umap[:, 1]
    cells["total_counts"] = total_counts
    cells["n_genes"] = n_genes

    # log1p(CP10k) per marker gene, from the raw counts layer.
    scale = TARGET_SUM / total_counts
    for j, gene in enumerate(present):
        cells[f"expr_{gene}"] = np.log1p(marker_counts[:, j] * scale).astype(np.float32)

    cells.index.name = "cell"
    return cells, missing


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--release", required=True, help="release tag used in file names, e.g. v9")
    parser.add_argument("--genes", required=True, help="comma-separated marker genes for the S25 dot plot")
    parser.add_argument("--out-cells", required=True, type=Path)
    parser.add_argument("--out-uns", required=True, type=Path)
    parser.add_argument("--out-provenance", required=True, type=Path)
    args = parser.parse_args()
    genes = args.genes.split(",")

    # Input identity: the object must be the one the release manifest describes ----
    manifest = json.loads(args.manifest.read_text())
    expected_sha = manifest["objects"][args.h5ad.stem]["sha256"]
    actual_sha = sha256(args.h5ad)
    if actual_sha != expected_sha:
        raise ValueError(f"{args.h5ad} sha256 {actual_sha} differs from its manifest entry {expected_sha}")
    log.info(f"{args.h5ad.name}: sha256 matches {manifest['release']} manifest")

    with h5py.File(args.h5ad, "r") as h5:
        cells, missing = build_cells_table(h5, genes)
        uns = read_uns(h5)

    args.out_cells.parent.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(args.out_cells)
    args.out_uns.write_text(json.dumps(uns, indent=1, sort_keys=True))
    log.info(f"wrote {args.out_cells} ({cells.shape[0]:,} x {cells.shape[1]}) and {args.out_uns}")

    stat = args.h5ad.stat()
    script = Path(__file__)
    provenance = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "release": args.release,
        "release_manifest": {"path": str(args.manifest.resolve()), "release": manifest["release"]},
        "input": {
            "path": str(args.h5ad.resolve()),
            "bytes": stat.st_size,
            "mtime": datetime.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
            "sha256": actual_sha,
            "sha256_matches_manifest": True,
        },
        "script": {"path": str(script.resolve()), "sha256": sha256(script)},
        "n_obs": int(cells.shape[0]),
        "umap_key": UMAP_KEY,
        "expression": f"log1p(counts / total_counts * {TARGET_SUM:g}) from layers['counts']",
        "marker_genes_requested": genes,
        "marker_genes_missing": missing,
        "uns_keys_exported": sorted(uns),
    }
    args.out_provenance.write_text(json.dumps(provenance, indent=1))


if __name__ == "__main__":
    main()
