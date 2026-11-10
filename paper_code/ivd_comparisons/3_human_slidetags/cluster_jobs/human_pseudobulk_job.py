"""Step 1: pseudobulks of the human fib/dec object, for the decidual vs fibroblast DE rerun.

Reads, without modifying:
  NEW_OBJECT  the current object, layers/counts
  OLD_OBJECT  the October 2025 object that the old notebook used, layers/raw (only for the old vs new comparison)

Writes into the working directory:
  human_pb_meta.csv       one row per pseudobulk: id, label column, unit (donor or section), donor, phase, arm,
                          number of nuclei, total UMIs
  human_pb_counts.csv.gz  UMIs summed over each pseudobulk's nuclei (pseudobulks x genes)
  human_pb_detect.csv.gz  number of each pseudobulk's nuclei with at least one UMI of the gene (pseudobulks x genes)
  human_var.csv           Ensembl gene ID and symbol, in column order
  human_old_vs_new.json   comparison with the old object, matching nuclei by sample and barcode

Pseudobulks are built for both label columns (coarse_refined and cell_type), per donor and per Slide-tags section,
and per arm. A few more hold, for each donor, the nuclei on which the two label columns disagree.
"""
import json
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

NEW_OBJECT = ("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/deliverables/2026_09_15/"
              "Cevrim_Human_Only_Fib_Dec_Cells.h5ad")
OLD_OBJECT = ("/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data/"
              "20251008_Cevrim_Human_Only_Fib_Dec_Cells.h5ad")

LABEL_COLUMNS = ("coarse_refined", "cell_type")

# The label values used in either object, mapped to the two arms of the comparison
ARM = {
    "decidual cells": "decidual",
    "decidual_cell": "decidual",
    "endometrial fibroblasts": "fibroblast",
    "endometrial_fibroblast": "fibroblast",
}

# Compared nucleus by nucleus between the two objects
MARKER_GENES = ("PRL", "IGFBP1", "LEFTY2")

START = time.time()


def log(message):
    print(f"[{time.time() - START:7.1f}s] {message}", flush=True)


def read_object(path, layer):
    """obs, var and one count layer of an h5ad file, as an integer sparse matrix (nuclei x genes)."""
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        counts = sp.csr_matrix(read_elem(f[f"layers/{layer}"]))
    assert counts.shape == (len(obs), len(var))
    assert np.array_equal(counts.data, np.round(counts.data)), f"layers/{layer} holds non-integer values"
    return obs, var, counts.astype(np.int64)


def single_value(values):
    """The one value a column takes within a group of nuclei; fails if there is more than one."""
    unique = values.astype(str).unique()
    assert len(unique) == 1, unique
    return unique[0]


def nonzero_crosstab(rows, columns, key_format):
    """Counts of each combination of two labels, as {"<row> <sep> <column>": count}, omitting zeros."""
    table = pd.crosstab(rows, columns)
    return {key_format.format(row, column): int(n) for (row, column), n in table.stack().items() if n}


def main():
    # ---- the new object -----------------------------------------------------------------------------------------
    obs, var, counts = read_object(NEW_OBJECT, "counts")
    gene_ids = var["gene_ids"].astype(str).to_numpy()
    symbols = var.index.astype(str).to_numpy()
    assert len(set(gene_ids)) == len(gene_ids)
    log(f"new object {counts.shape}, nnz {counts.nnz}, stored zeros {int((counts.data == 0).sum())}")

    is_detected = (counts > 0).astype(np.int64)                  # 1 where a nucleus has at least one UMI of a gene
    umis_per_nucleus = np.asarray(counts.sum(axis=1)).ravel()

    descriptions, summed_counts, detecting_nuclei = [], [], []

    def add_pseudobulk(in_pseudobulk, **description):
        """Sum the selected nuclei into one pseudobulk and record what it contains."""
        rows = np.flatnonzero(in_pseudobulk)
        descriptions.append(dict(description, n_cells=int(len(rows)), total_counts=int(umis_per_nucleus[rows].sum())))
        summed_counts.append(np.asarray(counts[rows].sum(axis=0)).ravel())
        detecting_nuclei.append(np.asarray(is_detected[rows].sum(axis=0)).ravel())

    arm_by_label = {label: obs[label].astype(str).map(ARM) for label in LABEL_COLUMNS}
    for label, arm in arm_by_label.items():
        assert arm.notna().all(), (label, obs[label].unique())

    # Donor x arm and section x arm pseudobulks, for both label columns
    for label, arm in arm_by_label.items():
        for unit in ("donor", "sample"):
            for unit_value in sorted(obs[unit].astype(str).unique()):
                in_unit = (obs[unit].astype(str) == unit_value).to_numpy()
                donor = single_value(obs.loc[in_unit, "donor"])
                phase = single_value(obs.loc[in_unit, "timepoint"])
                for arm_name in ("fibroblast", "decidual"):
                    add_pseudobulk(in_unit & (arm == arm_name).to_numpy(),
                                   id=f"{label}|{unit}={unit_value}|{arm_name}", label_source=label,
                                   unit_type=unit, unit=unit_value, donor=donor, timepoint=phase, arm=arm_name)

    # Nuclei on which the two label columns disagree, per donor
    coarse_arm = arm_by_label["coarse_refined"].to_numpy()
    labels_disagree = coarse_arm != arm_by_label["cell_type"].to_numpy()
    for donor in sorted(obs.donor.astype(str).unique()):
        of_donor = (obs.donor.astype(str) == donor).to_numpy()
        for kind, arm_in_coarse_refined in (("coarseFib_celltypeDec", "fibroblast"),
                                            ("coarseDec_celltypeFib", "decidual")):
            selected = labels_disagree & (coarse_arm == arm_in_coarse_refined) & of_donor
            if selected.sum():
                add_pseudobulk(selected, id=f"discordant_{kind}|donor={donor}", label_source=f"discordant_{kind}",
                               unit_type="donor", unit=donor, donor=donor,
                               timepoint=single_value(obs.loc[selected, "timepoint"]), arm=kind)

    pseudobulks = pd.DataFrame(descriptions)
    pseudobulks.to_csv("human_pb_meta.csv", index=False)
    pd.DataFrame(np.vstack(summed_counts), index=pseudobulks.id, columns=gene_ids).to_csv("human_pb_counts.csv.gz")
    pd.DataFrame(np.vstack(detecting_nuclei), index=pseudobulks.id, columns=gene_ids).to_csv("human_pb_detect.csv.gz")
    pd.DataFrame({"gene_id": gene_ids, "gene_symbol": symbols}).to_csv("human_var.csv", index=False)
    log(f"wrote {len(pseudobulks)} pseudobulks")

    # ---- old vs new, matching nuclei by sample and cell barcode -------------------------------------------------
    report = {"coarse_refined_x_cell_type": nonzero_crosstab(obs.coarse_refined.astype(str),
                                                             obs.cell_type.astype(str), "{} | {}")}
    old_obs, old_var, old_counts = read_object(OLD_OBJECT, "raw")
    new_keys = (obs["sample"].astype(str) + "|" + obs["cell_barcode"].astype(str)).to_numpy()
    old_keys = (old_obs["sample"].astype(str) + "|" + old_obs["cell_barcode"].astype(str)).to_numpy()
    shared_keys, new_rows, old_rows = np.intersect1d(new_keys, old_keys, return_indices=True)
    report.update(n_new=int(len(new_keys)), n_old=int(len(old_keys)), shared=int(len(shared_keys)),
                  new_only=int(len(new_keys) - len(shared_keys)), old_only=int(len(old_keys) - len(shared_keys)),
                  same_feature_order=bool(list(old_var.index.astype(str)) == list(symbols)))

    umis_new = umis_per_nucleus[new_rows]
    umis_old = np.asarray(old_counts.sum(axis=1)).ravel()[old_rows]
    report["total_counts_shared"] = {
        "identical_fraction": float((umis_new == umis_old).mean()),
        "median_ratio_new_over_old": float(np.median(umis_new / np.maximum(umis_old, 1))),
        "pearson": float(np.corrcoef(umis_new, umis_old)[0, 1]),
        "new_median": float(np.median(umis_new)),
        "old_median": float(np.median(umis_old)),
    }

    old_arm = old_obs["cell_type_annotations"].astype(str).map(ARM).to_numpy()[old_rows]
    for label, arm in arm_by_label.items():
        report[f"old_label_vs_new_{label}"] = nonzero_crosstab(old_arm, arm.to_numpy()[new_rows], "old {} -> new {}")

    for gene in MARKER_GENES:
        new_column = int(np.flatnonzero(symbols == gene)[0])
        old_column = int(np.flatnonzero(old_var.index.astype(str) == gene)[0])
        new_values = np.asarray(counts[new_rows][:, new_column].todense()).ravel()
        old_values = np.asarray(old_counts[old_rows][:, old_column].todense()).ravel()
        report[f"{gene}_shared_nuclei"] = {
            "detected_new": int((new_values > 0).sum()),
            "detected_old": int((old_values > 0).sum()),
            "umi_new": int(new_values.sum()),
            "umi_old": int(old_values.sum()),
            "identical_counts_fraction": float((new_values == old_values).mean()),
        }

    new_only_rows = np.setdiff1d(np.arange(len(new_keys)), new_rows)
    report["new_only_by_sample_and_coarse"] = nonzero_crosstab(obs["sample"].astype(str).to_numpy()[new_only_rows],
                                                               coarse_arm[new_only_rows], "{} | {}")
    with open("human_old_vs_new.json", "w") as handle:
        json.dump(report, handle, indent=1, default=str)
    log("done")


if __name__ == "__main__":
    main()
