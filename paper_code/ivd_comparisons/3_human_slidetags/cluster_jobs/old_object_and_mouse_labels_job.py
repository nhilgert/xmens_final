"""Step 2 (only needed for the regression check), plus a check of the mouse object's label columns.

1. Pseudobulks of the October 2025 human object, built the way our_human_invivo_DE.ipynb built them: layer "raw",
   one pseudobulk per phase (timepoint) and arm, with the arm taken from cell_type_annotations ("decidual cells" or
   "endometrial fibroblasts"; other labels are left out). Writes old_pb_meta.csv and old_pb_counts.csv.gz.
2. The current mouse fib/dec deliverable: do its label columns agree with coarse_refined?

Both parts are summarized in old_object_and_mouse_labels.json. The objects are only read; everything is written into
the working directory.
"""
import glob
import json

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

OLD_OBJECT = ("/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/"
              "20251008_Cevrim_Human_Only_Fib_Dec_Cells.h5ad")
DELIVERABLES = "/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/deliverables"
MOUSE_OBJECT_PATTERNS = [f"{DELIVERABLES}/2026_09_15/*XMens*Fib_Dec*.h5ad",
                         f"{DELIVERABLES}/2026_09_02_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"]
MOUSE_LABEL_COLUMNS = ("coarse_refined", "cell_type", "cell_subtype", "deliverable_cell_type")


def nonzero_crosstab(rows, columns, key_format):
    """Counts of each combination of two labels, as {"<row> <sep> <column>": count}, omitting zeros."""
    table = pd.crosstab(rows, columns)
    return {key_format.format(row, column): int(n) for (row, column), n in table.stack().items() if n}


def old_object_pseudobulks():
    """Sum the old object's raw counts per phase and arm, and write the pseudobulks."""
    with h5py.File(OLD_OBJECT, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        counts = sp.csr_matrix(read_elem(f["layers/raw"]))
    assert counts.shape == (len(obs), len(var))
    assert np.array_equal(counts.data, np.round(counts.data)), "layers/raw holds non-integer values"
    counts = counts.astype(np.int64)

    label = obs["cell_type_annotations"].astype(str).to_numpy()
    phase = obs["timepoint"].astype(str).to_numpy()
    arm = np.where(label == "decidual cells", "decidual",
                   np.where(label == "endometrial fibroblasts", "fibroblast", "other"))

    descriptions, summed_counts = [], []
    for phase_name in sorted(set(phase)):
        for arm_name in ("decidual", "fibroblast"):
            rows = np.flatnonzero((phase == phase_name) & (arm == arm_name))
            descriptions.append(dict(id=f"{phase_name}__{arm_name}", timepoint=phase_name, group=arm_name,
                                     n_cells=len(rows)))
            summed_counts.append(np.asarray(counts[rows].sum(axis=0)).ravel())

    pd.DataFrame(descriptions).to_csv("old_pb_meta.csv", index=False)
    pd.DataFrame(np.vstack(summed_counts), index=[d["id"] for d in descriptions],
                 columns=var.index.astype(str)).to_csv("old_pb_counts.csv.gz")
    return {"old_groups": descriptions, "old_other": int((arm == "other").sum())}


def mouse_label_agreement():
    """Label counts in each mouse fib/dec object, and how each label column cross-tabulates with coarse_refined."""
    paths = []
    for pattern in MOUSE_OBJECT_PATTERNS:
        paths += sorted(glob.glob(pattern))
    report = {"mouse_candidates": paths}
    for path in paths:
        with h5py.File(path, "r") as f:
            obs = read_elem(f["obs"])
        entry = {"n_obs": int(len(obs)), "columns": list(obs.columns)}
        label_columns = [column for column in MOUSE_LABEL_COLUMNS if column in obs.columns]
        for column in label_columns:
            entry[f"counts::{column}"] = {str(k): int(v) for k, v in obs[column].astype(str).value_counts().items()}
        if "coarse_refined" in obs.columns:
            for column in label_columns:
                if column != "coarse_refined":
                    entry[f"coarse_refined x {column}"] = nonzero_crosstab(obs.coarse_refined.astype(str),
                                                                           obs[column].astype(str), "{} | {}")
        report[path] = entry
    return report


def main():
    report = old_object_pseudobulks()
    report.update(mouse_label_agreement())
    with open("old_object_and_mouse_labels.json", "w") as handle:
        json.dump(report, handle, indent=1, default=str)
    print("done", flush=True)


if __name__ == "__main__":
    main()
