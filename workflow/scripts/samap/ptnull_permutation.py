"""Permutation control for the SAMap transfer of mouse pseudotime onto human nuclei.

Called by rule samap_ptnull (workflow/rules/10c_samap_figures.smk). Reads the human x mouse block of the joint SAMap
graph (C_hs_mm.npz with its row and column names) and the per-nucleus tables of export_samap_tables.py. Transfers mouse
pt_harmonized_rank as samap_transfer_pt.py does, checks the result against the stored transfer of record, then repeats
the transfer with mouse pseudotime permuted (a) across all mouse nuclei with a value and (b) only among nuclei of the
same mouse compartment, with the graph held fixed. Per replicate: Spearman rho with the human pseudotime pt_human
(all nuclei, decidual, fibroblast) and Moran's I of the transferred values (over all nuclei, and within each
compartment).
Writes null_replicates.csv, example_transfers.csv, null_summary.json and ptnull_statistics_table.csv in --out-dir.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata
from sklearn.neighbors import NearestNeighbors


MOUSE_PT = "pt_harmonized_rank"  # transferred
HUMAN_PT = "pt_human"  # compared with: export_samap_tables.py's human pseudotime (see its --human-pt)
SPATIAL_K = 10  # spatial neighbors per nucleus for Moran's I, within one 10x lane
CHUNK = 250  # replicates transferred per sparse-dense product (memory ~ 47k x CHUNK doubles)

# Statistics reported in the figure and the methods, in display order.
REPORTED = [
    "rho_all", "rho_decidual", "rho_fibroblast",
    "moran_all", "moran_within_decidual", "moran_within_fibroblast",
]


def transfer(C: sp.csr_matrix, values: np.ndarray, has_value: np.ndarray) -> np.ndarray:
    """Kernel-weighted average of mouse values onto human rows; NaN where no valued neighbor.

    values is (n_mouse, R); rows without a value must already be 0 so they drop out of the numerator, and has_value
    removes them from the denominator, as in samap_transfer_pt.py.
    """
    numerator = C @ values
    denominator = C @ has_value.astype(float)
    out = np.full(numerator.shape, np.nan)
    ok = denominator > 1e-10
    out[ok] = numerator[ok] / denominator[ok][:, None]
    return out


def rank_unit(values: np.ndarray, has_value: np.ndarray) -> np.ndarray:
    """Average rank / N over nuclei with a value; 0 elsewhere, so they drop out of the numerator."""
    out = np.zeros(len(values))
    out[has_value] = rankdata(values[has_value]) / has_value.sum()
    return out


def spearman_columns(transferred: np.ndarray, human_pt: np.ndarray) -> np.ndarray:
    """Spearman rho between human_pt and each column of transferred (rows already filtered)."""
    rank_human = rankdata(human_pt)
    rank_human = (rank_human - rank_human.mean()) / rank_human.std()
    ranks = rankdata(transferred, axis=0)
    ranks = (ranks - ranks.mean(0)) / ranks.std(0)
    return (ranks * rank_human[:, None]).mean(0)


def lane_weights(coords: np.ndarray) -> sp.csr_matrix:
    """Row-standardized k-nearest-neighbor weights on spatial coordinates (microns)."""
    nn = NearestNeighbors(n_neighbors=SPATIAL_K + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(len(coords)), SPATIAL_K)
    cols = idx[:, 1:].ravel()  # column 0 is the nucleus itself
    data = np.full(rows.size, 1.0 / SPATIAL_K)
    return sp.csr_matrix((data, (rows, cols)), shape=(len(coords), len(coords)))


def morans_i(W: sp.csr_matrix, values: np.ndarray) -> np.ndarray:
    """Moran's I per column; with row-standardized W this is z'Wz / z'z."""
    z = values - values.mean(0)
    return (z * (W @ z)).sum(0) / (z * z).sum(0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--mouse-cells", type=Path, required=True)
    parser.add_argument("--human-cells", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-perm", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # Load inputs, in the graph's row and column order ----
    C = sp.load_npz(args.graph_dir / "C_hs_mm.npz").tocsr()
    hs_names = (args.graph_dir / "hs_names.txt").read_text().split()
    mm_names = (args.graph_dir / "mm_names.txt").read_text().split()

    mouse = pd.read_csv(args.mouse_cells, index_col="cell").loc[mm_names]
    human = pd.read_csv(args.human_cells, index_col="cell").loc[hs_names]
    print(f"graph {C.shape}, nnz {C.nnz}; mouse {len(mouse)}; human {len(human)}")

    # Observed transfer: global rank / N over mouse nuclei with a value. Apoptotic nuclei have no value and stay out
    # of both numerator and denominator in every replicate.
    has_value = mouse[MOUSE_PT].notna().to_numpy()
    pt_rank = rank_unit(mouse[MOUSE_PT].to_numpy(), has_value)
    observed = transfer(C, pt_rank[:, None], has_value)[:, 0]
    print(f"{MOUSE_PT}: {has_value.sum()} of {len(mouse)} mouse nuclei carry a value")

    # The operator must reproduce the stored transfer of record before any null is run.
    stored = human["pt_transferred"].to_numpy()
    same_nan = np.array_equal(np.isnan(observed), np.isnan(stored))
    max_diff = float(np.nanmax(np.abs(observed - stored)))
    print(f"operator check: NaN pattern identical {same_nan}, max |diff| {max_diff:.2e}")
    if not same_nan or max_diff > 1e-9:
        raise SystemExit("the recomputed transfer does not match the stored transfer of record (pt_transferred)")

    # Human nuclei entering each statistic ----
    human_pt = human[HUMAN_PT].to_numpy()
    scored = np.isfinite(observed) & np.isfinite(human_pt)
    subsets = {
        "all": scored,
        "decidual": scored & (human["cell_type"] == "decidual_cell").to_numpy(),
        "fibroblast": scored & (human["cell_type"] == "endometrial_fibroblast").to_numpy(),
    }

    # Spatial neighbors are found within each 10x lane. The three lanes of a phase share one tissue section and one
    # coordinate frame; Moran's I is computed per lane and averaged over the six lanes.
    positioned = scored & human["x_um_qc"].notna().to_numpy()
    lanes = {
        lane: np.flatnonzero(positioned & (human["sample"] == lane).to_numpy())
        for lane in sorted(human["sample"].unique())
    }
    weights = {lane: lane_weights(human.iloc[ix][["x_um_qc", "y_um_qc"]].to_numpy()) for lane, ix in lanes.items()}
    print({k: int(v.sum()) for k, v in subsets.items()}, {lane: len(ix) for lane, ix in lanes.items()})

    # Spatial autocorrelation within one compartment: neighbors are drawn only among that compartment's nuclei in a
    # lane, so the fibroblast/decidual boundary cannot contribute.
    compartment_cells = {
        (k, lane): ix[subsets[k][ix]] for k in ("decidual", "fibroblast") for lane, ix in lanes.items()
    }
    compartment_weights = {
        key: lane_weights(human.iloc[ix][["x_um_qc", "y_um_qc"]].to_numpy()) for key, ix in compartment_cells.items()
    }
    print({f"{k}/{lane}": len(ix) for (k, lane), ix in compartment_cells.items()})

    def statistics(transferred: np.ndarray) -> pd.DataFrame:
        stats = {f"rho_{k}": spearman_columns(transferred[m], human_pt[m]) for k, m in subsets.items()}
        for lane, ix in lanes.items():
            stats[f"moran_{lane}"] = morans_i(weights[lane], transferred[ix])
        stats["moran_all"] = np.mean([stats[f"moran_{lane}"] for lane in lanes], axis=0)
        for k in ("decidual", "fibroblast"):
            per_lane = [
                morans_i(compartment_weights[(k, lane)], transferred[compartment_cells[(k, lane)]]) for lane in lanes
            ]
            stats[f"moran_within_{k}"] = np.mean(per_lane, axis=0)
        stats["sd_raw"] = np.nanstd(transferred[scored], axis=0)
        return pd.DataFrame(stats)

    # Permutation designs ----
    # uniform: values shuffled across all valued mouse nuclei (destroys all structure).
    # within_compartment: shuffled only among valued nuclei of the same mouse cell type, which keeps the
    # fibroblast/decidual offset and destroys the order within compartments.
    groups = {
        "uniform": [np.flatnonzero(has_value)],
        "within_compartment": [
            np.flatnonzero(has_value & (mouse["cell_type"] == ct).to_numpy()) for ct in sorted(mouse["cell_type"].unique())
        ],
    }
    tables = [statistics(observed[:, None]).assign(design="observed", replicate=0)]
    examples = {}
    for design, index_groups in groups.items():
        for start in range(0, args.n_perm, CHUNK):
            n = min(CHUNK, args.n_perm - start)
            permuted = np.zeros((len(mouse), n))
            for r in range(n):
                for ix in index_groups:
                    permuted[ix, r] = pt_rank[rng.permutation(ix)]
            transferred = transfer(C, permuted, has_value)
            if start == 0:
                examples[design] = transferred[:, 0]
            tables.append(statistics(transferred).assign(design=design, replicate=np.arange(start, start + n) + 1))
        print(f"{design}: {args.n_perm} replicates done")

    human_stats = statistics(human_pt[:, None])  # spatial autocorrelation of human pseudotime itself
    replicates = pd.concat(tables, ignore_index=True)
    replicates.to_csv(args.out_dir / "null_replicates.csv", index=False)

    # Example transfers for the spatial panels: observed and the first replicate of each design.
    pd.DataFrame(
        {"cell": hs_names, "observed": observed, **{f"{d}_example": v for d, v in examples.items()}}
    ).to_csv(args.out_dir / "example_transfers.csv", index=False)

    metric_cols = [c for c in replicates.columns if c.startswith(("rho_", "moran_", "sd_"))]
    summary = {
        "n_perm": args.n_perm,
        "seed": args.seed,
        "spatial_k": SPATIAL_K,
        "mouse_pt_col": MOUSE_PT,
        "human_pt_col": HUMAN_PT,
        "reproduction_max_abs_diff": max_diff,
        "n_scored": {k: int(v.sum()) for k, v in subsets.items()},
        "human_pt_moran": human_stats.filter(like="moran_").iloc[0].round(4).to_dict(),
        "by_design": {
            d: {
                c: {
                    "mean": round(float(g[c].mean()), 4),
                    "sd": round(float(g[c].std()), 4),
                    "max": round(float(g[c].max()), 4),
                }
                for c in metric_cols
            }
            for d, g in replicates.groupby("design")
        },
    }
    (args.out_dir / "null_summary.json").write_text(json.dumps(summary, indent=1))

    # One row per reported statistic: the observed value against both nulls.
    observed_row = replicates[replicates["design"] == "observed"].iloc[0]
    rows = []
    for stat in REPORTED:
        row = {"statistic": stat, "observed": observed_row[stat]}
        for design in groups:
            vals = replicates.loc[replicates["design"] == design, stat]
            row |= {
                f"{design}_mean": vals.mean(),
                f"{design}_max": vals.max(),
                f"{design}_n_ge_obs": int((vals >= observed_row[stat]).sum()),
            }
        rows.append(row)
    pd.DataFrame(rows).round(3).to_csv(args.out_dir / "ptnull_statistics_table.csv", index=False)
    print("wrote null_replicates.csv, example_transfers.csv, null_summary.json, ptnull_statistics_table.csv")


if __name__ == "__main__":
    main()
