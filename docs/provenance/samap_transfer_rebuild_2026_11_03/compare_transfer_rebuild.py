"""Compare the SAMap pseudotime transfer rebuilt by its Snakemake rules with the transfer of record.

The rules samap_raw_pt_transfer and samap_raw_attach_pt (workflow/rules/21_samap_raw.smk, variant _harmonized) were
run into a scratch root. This script compares their outputs with the files of record made by hand on 24 Oct 2026:
  pt_transfer_harmonized/pt_transferred_hs.csv   per human nucleus, every transfer variant and the edge weight
  pt_transfer_harmonized/summary.json            join tiers, coverage and ranges
  samap_transferred_pt_harmonized.h5ad           the joint object with the transfer attached as obs columns
Writes one JSON report (--out): file checksums, per-column equality and the largest absolute difference, and every
HDF5 dataset or attribute outside obs that differs between the two objects.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem


# Obs columns written by samap_attach_pt.py; the figures draw the _ranknorm one.
TRANSFER_COLUMNS = [
    "transferred_mm_dpt_pseudotime",
    "transferred_mm_dpt_pseudotime_ranknorm",
    "transferred_mm_dpt_pseudotime_ranknorm_insample",
    "transferred_mm_dpt_pseudotime_weight",
]

# Rows read at a time when hashing a dataset, so a large matrix is never held whole.
HASH_ROWS = 1 << 20


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def compare_columns(original: pd.DataFrame, rebuilt: pd.DataFrame) -> dict:
    """Per shared column: whether the values are identical (NaN equal to NaN) and, for numbers, the largest gap."""
    report = {
        "same_index_and_order": bool(original.index.equals(rebuilt.index)),
        "columns_only_in_original": sorted(set(original.columns) - set(rebuilt.columns)),
        "columns_only_in_rebuilt": sorted(set(rebuilt.columns) - set(original.columns)),
        "columns": {},
    }
    rebuilt = rebuilt.reindex(original.index)

    for column in [c for c in original.columns if c in rebuilt.columns]:
        a, b = original[column], rebuilt[column]
        if a.dtype.kind in "fiu" and b.dtype.kind in "fiu":
            x, y = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
            same_nan = np.isnan(x) == np.isnan(y)
            both = ~np.isnan(x) & ~np.isnan(y)
            gap = float(np.max(np.abs(x[both] - y[both]))) if both.any() else 0.0
            report["columns"][column] = {
                "identical": bool(same_nan.all() and np.array_equal(x[both], y[both])),
                "nan_pattern_equal": bool(same_nan.all()),
                "max_abs_diff": gap,
                "n_values": int(both.sum()),
            }
        else:
            differ = int((a.astype(str).to_numpy() != b.astype(str).to_numpy()).sum())
            report["columns"][column] = {"identical": differ == 0, "n_differ": differ}
    return report


def dataset_digest(dataset: h5py.Dataset) -> str:
    digest = hashlib.sha256()
    if dataset.shape == ():
        digest.update(repr(dataset[()]).encode())
        return digest.hexdigest()
    for start in range(0, dataset.shape[0], HASH_ROWS):
        block = dataset[start:start + HASH_ROWS]
        if dataset.dtype.kind == "O":
            digest.update("\x00".join(map(repr, np.ravel(block))).encode())
        else:
            digest.update(np.ascontiguousarray(block).tobytes())
    return digest.hexdigest()


def inventory(path: Path) -> dict:
    """Every group and dataset outside obs, keyed by HDF5 name: attributes, and dtype, shape and content digest."""
    items = {}
    with h5py.File(path, "r") as handle:
        items["/"] = {"attrs": sorted((k, repr(v)) for k, v in handle.attrs.items())}

        def visit(name, obj):
            if name == "obs" or name.startswith("obs/"):
                return
            entry = {"attrs": sorted((k, repr(v)) for k, v in obj.attrs.items())}
            if isinstance(obj, h5py.Dataset):
                entry.update(dtype=str(obj.dtype), shape=list(obj.shape), sha256=dataset_digest(obj))
            items[name] = entry

        handle.visititems(visit)
    return items


def compare_objects(original: Path, rebuilt: Path) -> dict:
    with h5py.File(original, "r") as handle:
        obs_original = read_elem(handle["obs"])
    with h5py.File(rebuilt, "r") as handle:
        obs_rebuilt = read_elem(handle["obs"])

    before, after = inventory(original), inventory(rebuilt)
    differing = sorted(name for name in set(before) & set(after) if before[name] != after[name])
    return {
        "obs": compare_columns(obs_original, obs_rebuilt),
        "non_obs_items": len(before),
        "non_obs_only_in_original": sorted(set(before) - set(after)),
        "non_obs_only_in_rebuilt": sorted(set(after) - set(before)),
        "non_obs_differing": differing,
        "non_obs_differing_detail": {name: {"original": before[name], "rebuilt": after[name]} for name in differing},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--original-dir", type=Path, required=True, help="runs_raw/vanilla_raw of record")
    parser.add_argument("--rebuilt-dir", type=Path, required=True, help="runs_raw/vanilla_raw of the scratch root")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    csv_name = "pt_transfer_harmonized/pt_transferred_hs.csv"
    summary_name = "pt_transfer_harmonized/summary.json"
    object_name = "samap_transferred_pt_harmonized.h5ad"

    report = {"original_dir": str(args.original_dir), "rebuilt_dir": str(args.rebuilt_dir), "files": {}}
    for name in (csv_name, summary_name, object_name):
        original, rebuilt = args.original_dir / name, args.rebuilt_dir / name
        report["files"][name] = {
            "original_sha256": sha256_of(original),
            "rebuilt_sha256": sha256_of(rebuilt),
            "original_bytes": original.stat().st_size,
            "rebuilt_bytes": rebuilt.stat().st_size,
        }
        report["files"][name]["byte_identical"] = (
            report["files"][name]["original_sha256"] == report["files"][name]["rebuilt_sha256"]
        )

    report["csv"] = compare_columns(
        pd.read_csv(args.original_dir / csv_name, index_col="cell"),
        pd.read_csv(args.rebuilt_dir / csv_name, index_col="cell"),
    )

    summary_original = json.loads((args.original_dir / summary_name).read_text())
    summary_rebuilt = json.loads((args.rebuilt_dir / summary_name).read_text())
    report["summary_json_differing_keys"] = sorted(
        key for key in set(summary_original) | set(summary_rebuilt)
        if summary_original.get(key) != summary_rebuilt.get(key)
    )
    report["summary_json_differing"] = {
        key: {"original": summary_original.get(key), "rebuilt": summary_rebuilt.get(key)}
        for key in report["summary_json_differing_keys"]
    }

    report["object"] = compare_objects(args.original_dir / object_name, args.rebuilt_dir / object_name)

    transfer = {column: report["object"]["obs"]["columns"].get(column) for column in TRANSFER_COLUMNS}
    report["transfer_columns"] = transfer
    args.out.write_text(json.dumps(report, indent=1) + "\n")

    print(json.dumps({
        "byte_identical": {name: entry["byte_identical"] for name, entry in report["files"].items()},
        "csv_columns_identical": {c: v["identical"] for c, v in report["csv"]["columns"].items()},
        "transfer_columns": transfer,
        "obs_columns_not_identical": [
            c for c, v in report["object"]["obs"]["columns"].items() if not v["identical"]
        ],
        "non_obs_differing": report["object"]["non_obs_differing"],
        "summary_json_differing_keys": report["summary_json_differing_keys"],
    }, indent=1))


if __name__ == "__main__":
    main()
