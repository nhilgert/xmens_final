"""Per-nucleus plot table of one fibroblast/decidual release object, as read by the SAMap and Fig. 6 panels.

Called by rule samap_fibdec_plot_table (workflow/rules/10c_samap_figures.smk), once per species and release.
Reads:  data/deliverables/<release>/Cevrim_{XMens,Human}_Only_Fib_Dec_Cells.h5ad (h5py, read-only) and the
        release MANIFEST.json, whose sha256 for the object must match.
Writes: <species>_fibdec_plot_table.parquet, one row per nucleus in object order: every obs column as stored,
        two columns per embedding (<key>__0, <key>__1), the diffusion map of the harmonized pseudotime as
        dc0..dc<n-1>, and per gene <gene>__counts, <gene>__X and <gene>__log1p_cp10k plus total_counts_layer;
        and <table>.provenance.json with input checksums and counts.
The column layout is that of the 2026_10_22_flagged_harmonized plot tables, which an analysis session saved
before this rule existed.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem


# Counts are streamed in blocks of rows so a 2 GB object never sits in memory as a dense matrix.
ROW_BLOCK = 8192

# CP10k: counts per 10,000 counts of the nucleus, the normalization of the release's X.
CP10K = 1e4


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def check_manifest(h5ad: Path, manifest: Path) -> str:
    """sha256 of the object, which must equal the release MANIFEST entry for it."""
    # MANIFEST objects are keyed by object name, the file name without .h5ad.
    expected = json.loads(manifest.read_text())["objects"][h5ad.stem]["sha256"]
    observed = sha256_of(h5ad)
    if observed != expected:
        raise ValueError(f"{h5ad}: sha256 {observed} differs from {manifest} ({expected})")
    return observed


def gene_columns(group: h5py.Group, gene_index: int, n_obs: int) -> tuple[np.ndarray, np.ndarray]:
    """One gene's values and every nucleus's total, from a CSR matrix group (data, indices, indptr)."""
    if group.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError(f"{group.name} is {group.attrs.get('encoding-type')!r}, expected csr_matrix")
    indptr = group["indptr"][:]
    values = np.zeros(n_obs)
    totals = np.zeros(n_obs)

    for start in range(0, n_obs, ROW_BLOCK):
        stop = min(start + ROW_BLOCK, n_obs)
        lo, hi = indptr[start], indptr[stop]
        data = group["data"][lo:hi].astype(np.float64)
        indices = group["indices"][lo:hi]
        rows = np.repeat(np.arange(start, stop), np.diff(indptr[start:stop + 1]))
        totals[start:stop] = np.bincount(rows - start, weights=data, minlength=stop - start)
        hit = indices == gene_index
        values[rows[hit]] = data[hit]
    return values, totals


def build_table(h5ad: Path, embeddings: list[str], diffmap_key: str, genes: list[str]) -> pd.DataFrame:
    with h5py.File(h5ad, "r") as handle:
        obs = read_elem(handle["obs"])
        var_names = pd.Index(read_elem(handle["var"]).index.astype(str))
        table = obs.copy()

        # Embeddings as two columns each, as the panels read them.
        for key in embeddings:
            coords = np.asarray(read_elem(handle["obsm"][key]))
            table[f"{key}__0"] = coords[:, 0]
            table[f"{key}__1"] = coords[:, 1]

        # Diffusion components of the harmonized pseudotime, named dc0..dc<n-1> as in the 22 Oct tables.
        diffmap = np.asarray(read_elem(handle["obsm"][diffmap_key]))
        for component in range(diffmap.shape[1]):
            table[f"dc{component}"] = diffmap[:, component]

        # Expression of the named genes: raw counts, the release's X and log1p CP10k from the counts layer.
        missing = [gene for gene in genes if gene not in var_names]
        if missing:
            raise ValueError(f"genes not in var: {missing}")
        for gene in genes:
            gene_index = int(var_names.get_loc(gene))
            counts, totals = gene_columns(handle["layers"]["counts"], gene_index, len(table))
            normalized, _ = gene_columns(handle["X"], gene_index, len(table))
            table[f"{gene}__counts"] = counts
            table[f"{gene}__X"] = normalized
            table[f"{gene}__log1p_cp10k"] = np.log1p(counts / np.where(totals > 0, totals, 1) * CP10K)
            table["total_counts_layer"] = totals
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--embeddings", required=True, help="comma-separated obsm keys, e.g. X_spatial,X_umap_scVI")
    parser.add_argument("--diffmap-key", default="X_diffmap_harmonized")
    parser.add_argument("--genes", default="", help="comma-separated gene names (may be empty)")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()

    object_sha256 = check_manifest(args.h5ad, args.manifest)
    embeddings = [key for key in args.embeddings.split(",") if key]
    genes = [gene for gene in args.genes.split(",") if gene]
    table = build_table(args.h5ad, embeddings, args.diffmap_key, genes)

    # Index named as in the 2026_10_22 plot tables: export_samap_tables.py writes it to CSV and samap_edge_share.py /
    # ptnull_permutation.py read index_col="cell" (bundle a66e05b wrote it unnamed and the S30 chain failed).
    table.index.name = "cell"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(args.out)

    flags = table["embedding_outlier"].astype(bool)
    record = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        # Paths as the rule passes them, relative to the repository, so the record reads the same in any checkout.
        "script": os.path.relpath(__file__),
        "script_sha256": sha256_of(Path(__file__)),
        "release": json.loads(args.manifest.read_text())["release"],
        "inputs": {
            "h5ad": {"path": str(args.h5ad), "sha256": object_sha256},
            "manifest": {"path": str(args.manifest), "sha256": sha256_of(args.manifest)},
        },
        "output": {"path": str(args.out), "sha256": sha256_of(args.out)},
        "n_nuclei": int(len(table)),
        "n_columns": int(table.shape[1]),
        "embeddings": embeddings,
        "diffmap_key": args.diffmap_key,
        "genes": genes,
        "n_embedding_outlier": int(flags.sum()),
        "cell_subtype_counts": {str(k): int(v) for k, v in table["cell_subtype"].value_counts().sort_index().items()},
        "timepoint_counts": {str(k): int(v) for k, v in table["timepoint"].value_counts().sort_index().items()},
    }
    args.provenance.write_text(json.dumps(record, indent=1))
    print(f"{args.out}: {len(table)} nuclei, {table.shape[1]} columns, {int(flags.sum())} embedding outliers")


if __name__ == "__main__":
    main()
