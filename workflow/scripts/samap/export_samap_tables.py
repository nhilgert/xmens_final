"""Per-nucleus tables behind the SAMap pseudotime figures and the permutation control.

Called by rule samap_figure_tables (workflow/rules/10c_samap_figures.smk). Reads obs and the UMAP of the joint
SAMap object that carries the harmonized pseudotime transfer of record (no expression, no graph), the release plot
tables of the two fibroblast/decidual objects, and optionally a human pseudotime computed separately (--human-pt).
Writes, in --out-dir:
  joint.parquet    one row per joint nucleus, in the object's order: species, sample, timepoint,
                   transferred_mm_dpt_pseudotime*, X_umap__0, X_umap__1
  mouse_cells.csv  one row per mouse nucleus: cell_type, cell_subtype, pt_harmonized_rank, embedding_outlier
  human_cells.csv  one row per human nucleus: sample, timepoint, cell_type, x_um_qc, y_um_qc, pt_harmonized_rank,
                   embedding_outlier, pt_transferred (the transfer of record), pt_human (the human pseudotime the
                   figures and the permutation control compare against: --human-pt's column if given, else the
                   release's pt_harmonized_rank)
"""
import argparse
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem


TRANSFERRED = [
    "transferred_mm_dpt_pseudotime",
    "transferred_mm_dpt_pseudotime_ranknorm",
    "transferred_mm_dpt_pseudotime_ranknorm_insample",
    "transferred_mm_dpt_pseudotime_weight",
]

# Mouse pt_harmonized_rank averaged over each human nucleus's SAMap edges (samap_transfer_pt.py, variant _harmonized).
TRANSFER_OF_RECORD = "transferred_mm_dpt_pseudotime_ranknorm"

MOUSE_COLUMNS = ["cell_type", "cell_subtype", "pt_harmonized_rank", "embedding_outlier"]
HUMAN_COLUMNS = ["sample", "timepoint", "cell_type", "x_um_qc", "y_um_qc", "pt_harmonized_rank", "embedding_outlier"]


def read_joint(path: Path) -> pd.DataFrame:
    """obs columns and obsm['X_umap'] of the joint object, one row per nucleus, in the object's order."""
    with h5py.File(path, "r") as f:
        index = pd.Index(read_elem(f["obs/_index"]), name="cell")
        columns = {key: read_elem(f[f"obs/{key}"]) for key in ["species", "sample", "timepoint", *TRANSFERRED]}
        umap = f["obsm/X_umap"][:]
    joint = pd.DataFrame(columns, index=index)
    joint["X_umap__0"] = umap[:, 0]
    joint["X_umap__1"] = umap[:, 1]
    return joint


def human_pseudotime(human: pd.DataFrame, table_path, column: str):
    """The human pseudotime to compare against, per nucleus of `human`, and a note of its source.

    A separately computed pseudotime must cover every human nucleus and, where it carries labels, agree with the
    release's cell types; otherwise it describes another cell set.
    """
    if table_path is None:
        return human["pt_harmonized_rank"].to_numpy(), "release pt_harmonized_rank"
    table = pd.read_parquet(table_path)
    missing = human.index.difference(table.index)
    if len(missing):
        raise ValueError(f"{table_path} lacks {len(missing)} human nuclei, e.g. {list(missing[:3])}")
    if "label" in table.columns:
        disagree = int((table.loc[human.index, "label"].astype(str) != human["cell_type"].astype(str)).sum())
        if disagree:
            raise ValueError(f"{table_path}: {disagree} nuclei carry another label than the release's cell_type")
    return table.loc[human.index, column].to_numpy(), f"{table_path}:{column}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--joint", type=Path, required=True, help="samap_transferred_pt_harmonized.h5ad")
    parser.add_argument("--mouse-table", type=Path, required=True, help="release mouse_fibdec_plot_table.parquet")
    parser.add_argument("--human-table", type=Path, required=True, help="release human_fibdec_plot_table.parquet")
    parser.add_argument("--human-pt", type=Path, help="table of a separately computed human pseudotime, indexed by nucleus")
    parser.add_argument("--human-pt-col", default="pt_rank", help="its pseudotime column (default pt_rank)")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    joint = read_joint(args.joint)
    mouse = pd.read_parquet(args.mouse_table, columns=MOUSE_COLUMNS)
    human = pd.read_parquet(args.human_table, columns=HUMAN_COLUMNS)

    # The release tables and the joint object must hold the same nuclei, species by species; a mismatch would
    # mean the transfer and the labels come from different cell sets.
    for species, table in (("mm", mouse), ("hs", human)):
        in_joint = joint.index[joint["species"].astype(str) == species]
        if len(in_joint) != len(table) or not in_joint.isin(table.index).all():
            raise ValueError(f"{species}: the release table ({len(table)}) and the joint object ({len(in_joint)}) "
                             "hold different nuclei")

    human["pt_transferred"] = joint.loc[human.index, TRANSFER_OF_RECORD].to_numpy()
    human["pt_human"], human_pt_source = human_pseudotime(human, args.human_pt, args.human_pt_col)

    joint.to_parquet(args.out_dir / "joint.parquet")
    mouse.to_csv(args.out_dir / "mouse_cells.csv")
    human.to_csv(args.out_dir / "human_cells.csv")

    print(f"joint: {len(joint)} nuclei; mouse: {len(mouse)} ({int(mouse.pt_harmonized_rank.notna().sum())} with "
          f"pseudotime, {int(mouse.embedding_outlier.sum())} embedding outliers); human: {len(human)} "
          f"({int(human.pt_transferred.notna().sum())} with a transferred value, {int(human.x_um_qc.notna().sum())} "
          f"positioned, {int(human.embedding_outlier.sum())} embedding outliers); human pseudotime: {human_pt_source}")


if __name__ == "__main__":
    main()
