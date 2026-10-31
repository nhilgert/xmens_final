"""Per-nucleus plot table from one release object, for the mouse stroma supplementary panels.

Reads a release .h5ad read-only (obs, optional obsm embeddings, and selected gene columns of X,
which holds log-normalized expression) and writes:
  --out   one row per nucleus (index = obs_names): obs columns, embedding coordinates, genes (csv.gz)
  --prov  input path, size, mtime and sha256 (checked against the release MANIFEST), the X
          normalization check, and row counts (json)
Rule: panels_mouse_stroma_table (workflow/rules/91_panels_mouse_stroma.smk).
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd

ROWS_PER_BLOCK = 8192  # CSR rows read per h5py call; keeps memory flat on 2 GB objects
NORM_CHECK_ROWS = 2000  # nuclei used to infer the normalization target of X
NORM_CHECK_SEED = 0


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gene_columns(group: h5py.Group, gene_positions: np.ndarray, n_obs: int) -> np.ndarray:
    """Dense (n_obs x n_genes) array of the given var positions from a CSR matrix group."""
    lookup = {int(pos): j for j, pos in enumerate(gene_positions)}
    wanted = np.zeros(int(group.attrs["shape"][1]), dtype=bool)
    wanted[gene_positions] = True
    column_of = np.full(wanted.size, -1, dtype=np.int64)
    column_of[list(lookup)] = list(lookup.values())

    dense = np.zeros((n_obs, len(gene_positions)), dtype=np.float32)
    indptr = group["indptr"][:]
    for start in range(0, n_obs, ROWS_PER_BLOCK):
        stop = min(start + ROWS_PER_BLOCK, n_obs)
        lo, hi = indptr[start], indptr[stop]
        indices = group["indices"][lo:hi]
        values = group["data"][lo:hi]
        rows = np.repeat(np.arange(start, stop), np.diff(indptr[start:stop + 1]))
        keep = wanted[indices]
        dense[rows[keep], column_of[indices[keep]]] = values[keep]
    return dense


def row_sums(group: h5py.Group, start: int, stop: int, transform=None) -> np.ndarray:
    indptr = group["indptr"][start:stop + 1]
    values = group["data"][indptr[0]:indptr[-1]].astype(np.float64)
    if transform is not None:
        values = transform(values)
    rows = np.repeat(np.arange(stop - start), np.diff(indptr))
    return np.bincount(rows, weights=values, minlength=stop - start)


def check_normalization(handle: h5py.File, n_obs: int) -> dict:
    """Infer the per-nucleus target sum of X = log1p(counts / total * target) on a block of nuclei."""
    rng = np.random.default_rng(NORM_CHECK_SEED)
    start = int(rng.integers(0, max(1, n_obs - NORM_CHECK_ROWS)))
    stop = min(start + NORM_CHECK_ROWS, n_obs)
    expm1_sums = row_sums(handle["X"], start, stop, transform=np.expm1)
    count_sums = row_sums(handle["layers/counts"], start, stop)
    placed = count_sums > 0
    return {
        "rows_checked": [start, stop],
        "median_expm1_rowsum": round(float(np.median(expm1_sums[placed])), 3),
        "min_expm1_rowsum": round(float(expm1_sums[placed].min()), 3),
        "max_expm1_rowsum": round(float(expm1_sums[placed].max()), 3),
        "median_counts_rowsum": float(np.median(count_sums[placed])),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--object-name", required=True)
    parser.add_argument("--obs-cols", required=True)
    parser.add_argument("--obsm", default="")
    parser.add_argument("--genes", default="")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--prov", required=True, type=Path)
    args = parser.parse_args()

    obs_cols = [c for c in args.obs_cols.split(",") if c]
    obsm_keys = [k for k in args.obsm.split(",") if k]
    genes = [g for g in args.genes.split(",") if g]

    # The release MANIFEST records each object's sha256; a mismatch means the file on disk is not the release.
    manifest = json.loads(args.manifest.read_text())
    expected_sha = manifest["objects"][args.object_name]["sha256"]
    observed_sha = sha256_of(args.h5ad)
    if observed_sha != expected_sha:
        raise ValueError(f"{args.h5ad}: sha256 {observed_sha} does not match MANIFEST {expected_sha}")

    with h5py.File(args.h5ad, "r") as handle:
        obs = anndata.io.read_elem(handle["obs"])
        var_names = pd.Index(anndata.io.read_elem(handle["var"]).index.astype(str))
        n_obs = len(obs)

        missing_cols = [c for c in obs_cols if c not in obs.columns]
        if missing_cols:
            raise KeyError(f"obs columns not in {args.h5ad.name}: {missing_cols}")
        table = obs[obs_cols].copy()

        for key in obsm_keys:
            coords = np.asarray(anndata.io.read_elem(handle["obsm"][key]))
            table[f"{key}_1"] = coords[:, 0]
            table[f"{key}_2"] = coords[:, 1]

        missing_genes = [g for g in genes if g not in var_names]
        if missing_genes:
            raise KeyError(f"genes not in var of {args.h5ad.name}: {missing_genes}")
        if genes:
            positions = var_names.get_indexer(genes)
            expression = read_gene_columns(handle["X"], positions, n_obs)
            for j, gene in enumerate(genes):
                table[gene] = expression[:, j]

        normalization = check_normalization(handle, n_obs)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.index.name = "cell"
    table.to_csv(args.out, compression="gzip")

    provenance = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "input": {
            "path": str(args.h5ad.resolve()),
            "bytes": args.h5ad.stat().st_size,
            "mtime": datetime.datetime.fromtimestamp(args.h5ad.stat().st_mtime).isoformat(timespec="seconds"),
            "sha256": observed_sha,
            "sha256_matches_manifest": True,
            "manifest": str(args.manifest.resolve()),
            "release": manifest.get("release"),
        },
        "n_obs": n_obs,
        "obs_columns": obs_cols,
        "obsm_keys": obsm_keys,
        "genes": genes,
        "expression_source": "X (log-normalized; see x_normalization)",
        "x_normalization": normalization,
        "output": str(args.out),
    }
    args.prov.write_text(json.dumps(provenance, indent=1))
    print(f"wrote {args.out} ({n_obs} nuclei, {len(genes)} genes, {len(obsm_keys)} embeddings)")


if __name__ == "__main__":
    main()
