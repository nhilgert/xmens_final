"""Diagnostic: which human label columns did the 2 September "11/24 relabel" change?

Compares the obs table of each pre-relabel backup (*.pre_relabel_11_24.bak) with the current object, for the human
fib/dec deliverable and for final_human.h5ad. For the current objects it also counts the nuclei on which
coarse_refined and cell_type disagree, and where they sit (Leiden clusters, samples). Finally it cross-tabulates the
same two columns in the mouse object.

Only obs tables are read; the objects are not modified. Writes relabel_report.json into the working directory.
"""
import json

import h5py
import pandas as pd
from anndata.io import read_elem

DATA = "/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data"
FIBDEC = f"{DATA}/deliverables/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"
BEFORE_AND_AFTER = {
    "deliverable_fibdec": (f"{FIBDEC}.pre_relabel_11_24.bak", FIBDEC),
    "final_human": (f"{DATA}/interim/final_human.h5ad.pre_relabel_11_24.bak", f"{DATA}/interim/final_human.h5ad"),
}
MOUSE_OBJECT = f"{DATA}/interim/final_mouse.h5ad"
LABEL_COLUMNS = ("coarse_refined", "cell_type", "cell_subtype", "fine_refined", "subtype")


def read_obs(path):
    with h5py.File(path, "r") as f:
        return read_elem(f["obs"])


def nonzero_crosstab(rows, columns, key_format, skip_diagonal=False):
    """Counts of each combination of two labels, as {"<row> <sep> <column>": count}, omitting zeros (and, if asked,
    combinations where both labels are the same)."""
    table = pd.crosstab(rows, columns)
    return {key_format.format(row, column): int(n) for (row, column), n in table.stack().items()
            if n and not (skip_diagonal and row == column)}


def same_style(labels):
    """Labels in one naming style: "decidual cells" and "decidual_cell" both become "decidual_cell", so that
    coarse_refined and cell_type, which use different styles, can be compared."""
    return labels.astype(str).str.replace(" ", "_").str.rstrip("s")


def compare_before_after(before, after):
    shared = before.index.intersection(after.index)
    result = {"n_pre": len(before), "n_cur": len(after), "shared": len(shared)}

    # Which label columns changed, and how
    for column in LABEL_COLUMNS:
        if column in before.columns and column in after.columns:
            old_labels = before.loc[shared, column].astype(str)
            new_labels = after.loc[shared, column].astype(str)
            changed = old_labels != new_labels
            result[f"{column}: changed"] = int(changed.sum())
            if changed.any():
                result[f"{column}: transitions"] = nonzero_crosstab(old_labels, new_labels, "{} -> {}",
                                                                    skip_diagonal=True)

    # coarse_refined vs cell_type, before and after
    for when, obs in (("pre", before), ("cur", after)):
        if {"coarse_refined", "cell_type"} <= set(obs.columns):
            result[f"{when}: coarse_refined x cell_type"] = nonzero_crosstab(obs.coarse_refined.astype(str),
                                                                             obs.cell_type.astype(str), "{} | {}")

    # Nuclei on which the two columns now disagree, and where they are
    discordant = after.index[same_style(after.coarse_refined) != same_style(after.cell_type)]
    result["current discordant nuclei"] = int(len(discordant))
    leiden_columns = [column for column in after.columns if column.startswith("leiden")][:12]
    result["discordant by leiden"] = {}
    for column in leiden_columns:
        top_clusters = after.loc[discordant, column].astype(str).value_counts().head(5)
        result["discordant by leiden"][column] = {str(cluster): int(n) for cluster, n in top_clusters.items()}
    result["discordant by sample"] = (
        {str(sample): int(n) for sample, n in after.loc[discordant, "sample"].astype(str).value_counts().items()}
        if "sample" in after.columns else None)
    return result


def main():
    report = {}
    for name, (before_path, after_path) in BEFORE_AND_AFTER.items():
        report[name] = compare_before_after(read_obs(before_path), read_obs(after_path))
        print(name, "done", flush=True)

    mouse = read_obs(MOUSE_OBJECT)
    if "cell_type" in mouse.columns:
        crosstab = nonzero_crosstab(mouse.coarse_refined.astype(str), mouse.cell_type.astype(str), "{} | {}")
        report["mouse_final"] = {"coarse_refined x cell_type": crosstab}
    else:
        report["mouse_final"] = {"columns": list(mouse.columns)}

    with open("relabel_report.json", "w") as handle:
        json.dump(report, handle, indent=1, default=str)
    print("done", flush=True)


if __name__ == "__main__":
    main()
