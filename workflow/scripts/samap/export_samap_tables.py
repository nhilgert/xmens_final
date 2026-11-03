"""Per-nucleus tables behind the SAMap pseudotime figures and the permutation control.

Called by rule samap_figure_tables (workflow/rules/10c_samap_figures.smk). Reads obs and the UMAP of the joint
SAMap object that carries the harmonized pseudotime transfer of record (no expression, no graph), the release plot
tables of the two fibroblast/decidual objects (rule samap_fibdec_plot_table) and the human pseudotime of record
(--human-pt, sha256 checked against --human-pt-sha256).
Writes, in --out-dir:
  joint.parquet    one row per joint nucleus, in the object's order: species, sample, timepoint,
                   transferred_mm_dpt_pseudotime_ranknorm (the transfer of record, read from --transfer-col),
                   X_umap__0, X_umap__1
  mouse_cells.csv  one row per mouse nucleus: cell_type, cell_subtype, pt_harmonized_rank, embedding_outlier
  human_cells.csv  one row per human nucleus: sample, timepoint, cell_type, x_um_qc, y_um_qc, embedding_outlier,
                   pt_transferred (the transfer of record), pt_human (the human pseudotime of record, which the
                   figures and the permutation control compare against)
  provenance.json  input paths and sha256, the release of the plot tables, the human pseudotime source and counts
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem


# Mouse pt_harmonized_rank averaged over each human nucleus's SAMap edges (samap_transfer_pt.py, variant _harmonized),
# the one transfer the figures draw. joint.parquet carries it under this name whatever the joint object calls it
# (Samap_Final.h5ad: transferred_mouse_pseudotime); the object's other transfer variants have no reader.
TRANSFER_OF_RECORD = "transferred_mm_dpt_pseudotime_ranknorm"

MOUSE_COLUMNS = ["cell_type", "cell_subtype", "pt_harmonized_rank", "embedding_outlier"]
HUMAN_COLUMNS = ["sample", "timepoint", "cell_type", "x_um_qc", "y_um_qc", "embedding_outlier"]


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def read_joint(path: Path, transfer_col: str) -> pd.DataFrame:
    """species, sample, timepoint, the transfer of record and obsm['X_umap'] of the joint object, in its order."""
    with h5py.File(path, "r") as f:
        index = pd.Index(read_elem(f["obs/_index"]), name="cell")
        columns = {key: read_elem(f[f"obs/{key}"]) for key in ["species", "sample", "timepoint"]}
        columns[TRANSFER_OF_RECORD] = read_elem(f[f"obs/{transfer_col}"])
        umap = f["obsm/X_umap"][:]
    joint = pd.DataFrame(columns, index=index)
    joint["X_umap__0"] = umap[:, 0]
    joint["X_umap__1"] = umap[:, 1]
    return joint


def human_pseudotime(human: pd.DataFrame, table_path: Path, column: str):
    """The human pseudotime of record per nucleus of `human`, and a note of its source.

    The table must cover every human nucleus and, where it carries labels, agree with the release's cell types;
    otherwise it describes another cell set.
    """
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
    parser.add_argument("--joint-sha256", required=True, help="sha256 the joint object must have")
    parser.add_argument("--mouse-table", type=Path, required=True, help="release mouse_fibdec_plot_table.parquet")
    parser.add_argument("--human-table", type=Path, required=True, help="release human_fibdec_plot_table.parquet")
    parser.add_argument("--human-pt", type=Path, required=True, help="human pseudotime of record, indexed by nucleus")
    parser.add_argument("--human-pt-col", required=True, help="its pseudotime column (pt_rank)")
    parser.add_argument("--human-pt-sha256", required=True, help="sha256 the human pseudotime table must have")
    parser.add_argument("--release", required=True, help="release folder of the plot tables (provenance)")
    parser.add_argument("--mouse-pt-col", default="pt_harmonized_rank",
                        help="mouse table column of the pseudotime transferred; written as pt_harmonized_rank")
    parser.add_argument("--human-pt-release-col", default="",
                        help="take pt_human from this human table column (must equal --human-pt); empty: --human-pt")
    parser.add_argument("--transfer-col", default=TRANSFER_OF_RECORD,
                        help=f"joint obs column holding the transfer of record; written as {TRANSFER_OF_RECORD}")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    joint_sha256 = sha256_of(args.joint)
    if joint_sha256 != args.joint_sha256:
        raise ValueError(f"{args.joint}: sha256 {joint_sha256}, expected {args.joint_sha256}")

    human_pt_sha256 = sha256_of(args.human_pt)
    if human_pt_sha256 != args.human_pt_sha256:
        raise ValueError(f"{args.human_pt}: sha256 {human_pt_sha256}, expected {args.human_pt_sha256}")

    joint = read_joint(args.joint, args.transfer_col)
    # The mouse pseudotime transferred is written as pt_harmonized_rank whatever the release calls it.
    mouse_columns = [args.mouse_pt_col if column == "pt_harmonized_rank" else column for column in MOUSE_COLUMNS]
    mouse = pd.read_parquet(args.mouse_table, columns=mouse_columns)
    mouse = mouse.rename(columns={args.mouse_pt_col: "pt_harmonized_rank"})
    release_pt = [args.human_pt_release_col] if args.human_pt_release_col else []
    human = pd.read_parquet(args.human_table, columns=HUMAN_COLUMNS + release_pt)

    # The release tables and the joint object must hold the same nuclei, species by species; a mismatch would
    # mean the transfer and the labels come from different cell sets.
    for species, table in (("mm", mouse), ("hs", human)):
        in_joint = joint.index[joint["species"].astype(str) == species]
        if len(in_joint) != len(table) or not in_joint.isin(table.index).all():
            raise ValueError(f"{species}: the release table ({len(table)}) and the joint object ({len(in_joint)}) "
                             "hold different nuclei")

    human["pt_transferred"] = joint.loc[human.index, TRANSFER_OF_RECORD].to_numpy()
    human["pt_human"], human_pt_source = human_pseudotime(human, args.human_pt, args.human_pt_col)
    if args.human_pt_release_col:
        # The release object's own copy is drawn; it must equal the pseudotime of record nucleus by nucleus.
        from_release = human.pop(args.human_pt_release_col).to_numpy(dtype=float)
        of_record = human["pt_human"].to_numpy(dtype=float)
        differ = ~((from_release == of_record) | (np.isnan(from_release) & np.isnan(of_record)))
        if differ.any():
            raise ValueError(f"{args.human_table}:{args.human_pt_release_col} differs from {human_pt_source} "
                             f"on {int(differ.sum())} nuclei")
        human["pt_human"] = from_release
        human_pt_source = f"{args.human_table}:{args.human_pt_release_col} (equal to {human_pt_source})"

    joint.to_parquet(args.out_dir / "joint.parquet")
    mouse.to_csv(args.out_dir / "mouse_cells.csv")
    human.to_csv(args.out_dir / "human_cells.csv")

    inputs = {"joint": args.joint, "mouse_table": args.mouse_table, "human_table": args.human_table,
              "human_pt": args.human_pt}
    outputs = [args.out_dir / name for name in ("joint.parquet", "mouse_cells.csv", "human_cells.csv")]
    record = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "script": str(Path(__file__)),
        "script_sha256": sha256_of(Path(__file__)),
        "release": args.release,
        "inputs": {name: {"path": str(path), "sha256": sha256_of(path)} for name, path in inputs.items()},
        "human_pt_source": human_pt_source,
        "transfer_of_record": TRANSFER_OF_RECORD,
        "transfer_column_in_joint": args.transfer_col,
        "outputs": {path.name: sha256_of(path) for path in outputs},
        "n_nuclei": {"joint": int(len(joint)), "mouse": int(len(mouse)), "human": int(len(human))},
        "n_embedding_outlier": {
            "mouse": int(mouse.embedding_outlier.sum()),
            "human": int(human.embedding_outlier.sum()),
        },
    }
    (args.out_dir / "provenance.json").write_text(json.dumps(record, indent=1) + "\n")

    print(f"joint: {len(joint)} nuclei; mouse: {len(mouse)} ({int(mouse.pt_harmonized_rank.notna().sum())} with "
          f"pseudotime, {int(mouse.embedding_outlier.sum())} embedding outliers); human: {len(human)} "
          f"({int(human.pt_transferred.notna().sum())} with a transferred value, {int(human.x_um_qc.notna().sum())} "
          f"positioned, {int(human.embedding_outlier.sum())} embedding outliers); human pseudotime: {human_pt_source}")


if __name__ == "__main__":
    main()
