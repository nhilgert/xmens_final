"""Share of cross-species SAMap edge weight that runs to each compartment, both directions (SAMap supplement, panel d).

Called by rule samap_edge_share (workflow/rules/10c_samap_figures.smk). For each source species and source compartment
(fibroblast, decidual), the weights of graph edges to nuclei of the other species are summed and split by the target
nucleus's compartment. Computed for every graph in obsp, with the compartment labels stored in the SAMap object (when
it has them; --samap-labels-col) and with those of the release, as shares of the summed weight and as per-nucleus mean
shares.
Writes, in --out-dir:
  samap_compartment_edge_share_all.csv  every variant above
  samap_compartment_edge_share.csv      the variant the figure shows: connectivities, release labels, summed weight
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


COMPARTMENTS = ["endometrial_fibroblast", "decidual_cell"]
FIGURE_VARIANT = {"graph": "connectivities", "labels": "release"}


def edge_shares(graph: sp.csr_matrix, species: np.ndarray, compartment: np.ndarray) -> list[dict]:
    """One row per direction and source compartment: summed and per-nucleus shares of cross-species edge weight."""
    rows = []
    for source, target, direction in (("hs", "mm", "human_to_mouse"), ("mm", "hs", "mouse_to_human")):
        source_idx = np.flatnonzero(species == source)
        target_idx = np.flatnonzero(species == target)
        cross = graph[source_idx][:, target_idx]  # source nuclei x target nuclei
        target_compartment = compartment[target_idx]
        to_compartment = np.column_stack([
            np.asarray(cross[:, target_compartment == c].sum(axis=1)).ravel() for c in COMPARTMENTS
        ])

        for c in COMPARTMENTS:
            in_compartment = compartment[source_idx] == c
            summed = to_compartment[in_compartment].sum(axis=0)
            totals = to_compartment[in_compartment].sum(axis=1)

            # Nuclei without a cross-species edge have no share; they count only in n_source_cells.
            with_edges = totals > 0
            per_cell = (to_compartment[in_compartment][with_edges] / totals[with_edges, None]).mean(axis=0)
            rows.append({
                "direction": direction,
                "source_compartment": c,
                "n_source_cells": int(in_compartment.sum()),
                "n_with_cross_edges": int(with_edges.sum()),
                **{f"weight_share_{t}": float(s) for t, s in zip(COMPARTMENTS, summed / summed.sum())},
                **{f"percell_share_{t}": float(s) for t, s in zip(COMPARTMENTS, per_cell)},
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--joint", required=True, help="joint SAMap h5ad (obs species, cell_type_annotations; obsp)")
    parser.add_argument("--human-cells", required=True, help="human_cells.csv (release cell_type)")
    parser.add_argument("--mouse-cells", required=True, help="mouse_cells.csv (release cell_type)")
    parser.add_argument("--samap-labels-col", default="cell_type_annotations",
                        help="joint obs column of SAMap's own compartment labels; empty for an object without them")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.joint, backed="r")
    obs = adata.obs.copy()
    graphs = {key: sp.csr_matrix(adata.obsp[key]) for key in adata.obsp.keys()}
    print("graphs:", {key: (g.shape, g.nnz) for key, g in graphs.items()})
    species = obs["species"].astype(str).to_numpy()

    # Release labels, joined by nucleus barcode. Mouse decidual subtypes (apoptotic included) are the decidual
    # compartment, as in the SAMap object's own labels.
    human = pd.read_csv(args.human_cells, index_col="cell", usecols=["cell", "cell_type"])["cell_type"]
    mouse = pd.read_csv(args.mouse_cells, index_col="cell", usecols=["cell", "cell_type"])["cell_type"]
    release = pd.concat([human, mouse]).reindex(obs.index)
    if release.isna().any():
        raise ValueError(f"{int(release.isna().sum())} joint nuclei have no release label")
    # The deposited Samap_Final.h5ad dropped SAMap's own labels (its cell_type is the final objects'), so its variants
    # table has no samap_object rows; the figure's variant uses the release labels either way.
    label_sets = {"release": release.astype(str)}
    if args.samap_labels_col:
        label_sets = {"samap_object": obs[args.samap_labels_col].astype(str), **label_sets}
    else:
        print("joint object without SAMap compartment labels: release labels only")

    rows = []
    for graph_name, graph in graphs.items():
        for label_name, labels in label_sets.items():
            for row in edge_shares(graph, species, labels.to_numpy()):
                rows.append({"graph": graph_name, "labels": label_name, **row})
    table = pd.DataFrame(rows)
    table.to_csv(args.out_dir / "samap_compartment_edge_share_all.csv", index=False)
    print(table.to_string())

    # The figure's table: rows = direction x source compartment, columns = target compartment.
    chosen = table[(table["graph"] == FIGURE_VARIANT["graph"]) & (table["labels"] == FIGURE_VARIANT["labels"])]
    figure_table = chosen.set_index(["direction", "source_compartment"])[[f"weight_share_{c}" for c in COMPARTMENTS]]
    figure_table.columns = COMPARTMENTS
    figure_table.to_csv(args.out_dir / "samap_compartment_edge_share.csv")


if __name__ == "__main__":
    main()
