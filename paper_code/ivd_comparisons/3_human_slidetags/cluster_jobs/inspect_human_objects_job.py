"""Diagnostic: describe the four human fib/dec objects and compare their labels with the current deliverable's.

For each object: size, obs and var columns, the matrices it holds (type, and whether the first 200,000 stored values
are integers), and value counts of label-like obs columns. Then, for each other object: how many nuclei it shares
with the 2026_09_15 deliverable (matched by obs name) and how often each shared label column agrees.

Only obs, var and a slice of each matrix are read; the objects are not modified. Writes human_objects_report.json
into the working directory.
"""
import json

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

DATA = "/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data"
OBJECTS = {
    "deliv_2026_09_15": f"{DATA}/deliverables/2026_09_15/Cevrim_Human_Only_Fib_Dec_Cells.h5ad",
    "deliv_2026_08_28": f"{DATA}/deliverables/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad",
    "final_human": f"{DATA}/interim/final_human.h5ad",
    "old_20251008": ("/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/"
                     "20251008_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"),
}
REFERENCE = "deliv_2026_09_15"

# An obs column counts as a label column if its name contains one of these and it has at most 60 distinct values.
LABEL_NAME_PARTS = ("type", "annot", "label", "coarse", "fine", "subtype", "refined", "timepoint", "donor", "sample",
                    "phase", "section", "puck", "library", "batch")
MAX_LABEL_VALUES = 60
INTEGER_CHECK_VALUES = 200_000


def looks_like_label(column, values):
    return any(part in column.lower() for part in LABEL_NAME_PARTS) and values.nunique() <= MAX_LABEL_VALUES


def matrix_info(node):
    """Type and shape of a stored matrix; for sparse matrices, also whether the first stored values are integers."""
    if isinstance(node, h5py.Dataset):
        return {"kind": "dense", "dtype": str(node.dtype), "shape": list(node.shape)}
    data = node["data"]
    first_values = data[: min(INTEGER_CHECK_VALUES, data.shape[0])]
    return {"kind": str(node.attrs.get("encoding-type", "")),
            "dtype": str(data.dtype),
            "shape": [int(x) for x in list(node.attrs.get("shape", []))],
            "nnz": int(data.shape[0]),
            "first_200k_integer": bool(np.allclose(first_values, np.round(first_values))),
            "first_200k_max": float(first_values.max()) if first_values.size else None}


def describe(path):
    """Everything the report says about one object, and its obs table."""
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        layers = list(f["layers"].keys()) if "layers" in f else []
        info = {"path": path,
                "n_obs": int(obs.shape[0]),
                "n_var": int(var.shape[0]),
                "obs_columns": list(obs.columns),
                "var_columns": list(var.columns),
                "var_index_head": [str(x) for x in var.index[:6]],
                "X": matrix_info(f["X"]),
                "layers": {layer: matrix_info(f["layers"][layer]) for layer in layers},
                "raw": ({"X": matrix_info(f["raw/X"]), "n_var": int(read_elem(f["raw/var"]).shape[0])}
                        if "raw" in f else None)}
    info["categorical_summaries"] = {
        column: {str(value): int(n) for value, n in obs[column].astype(str).value_counts().items()}
        for column in obs.columns if looks_like_label(column, obs[column])}
    return info, obs


def compare_labels(reference_obs, other_obs):
    """Shared nuclei, and the fraction of shared nuclei on which each shared label column agrees."""
    shared = reference_obs.index.astype(str).intersection(other_obs.index.astype(str))
    result = {"shared_barcodes": int(len(shared)),
              "only_in_new": int(len(reference_obs) - len(shared)),
              "only_in_other": int(len(other_obs) - len(shared))}
    for column in sorted(set(reference_obs.columns) & set(other_obs.columns)):
        if not looks_like_label(column, reference_obs[column]):
            continue
        new_labels = reference_obs.loc[shared, column].astype(str)
        other_labels = other_obs.reindex(shared)[column].astype(str)
        agree = new_labels.to_numpy() == other_labels.to_numpy()
        result[f"agree::{column}"] = round(float(agree.mean()), 6)
        if not agree.all():
            changes = pd.crosstab(other_labels, new_labels).stack()
            result[f"changes::{column}"] = {f"{old} -> {new}": int(n) for (old, new), n in changes.items()
                                            if n and old != new}
    return result


def main():
    report, obs_tables = {}, {}
    for name, path in OBJECTS.items():
        try:
            report[name], obs_tables[name] = describe(path)
            print(f"read {name}: {obs_tables[name].shape}", flush=True)
        except Exception as error:
            report[name] = {"path": path, "error": repr(error)}
            print(f"failed {name}: {error!r}", flush=True)

    comparisons = {}
    if REFERENCE in obs_tables:
        for name, obs in obs_tables.items():
            if name != REFERENCE:
                comparisons[name] = compare_labels(obs_tables[REFERENCE], obs)
    report[f"comparisons_vs_{REFERENCE}"] = comparisons

    with open("human_objects_report.json", "w") as handle:
        json.dump(report, handle, indent=1, default=str)
    print("done", flush=True)


if __name__ == "__main__":
    main()
