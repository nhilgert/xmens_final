"""Select the one-to-one orthologs of the cross-species pseudotime heatmap, Fig. 5F (rule mi_heatmap_selection).

Vendored from results/mi_dpt/overlap_harmonized_nall/crossspecies_heatmap_selection/crossspecies_gene_selection.py
(sha256 94fc7ffe...), with the pseudotime read from the MI arm's input object (obs --*-pt-key, the nuclei and
pseudotime the arm scored) instead of a pseudotime run directory. Procedure (Methods, pseudotime expression heatmaps):
  1. Pool: one-to-one pairs measured in both MI arms; each species' top fraction by MI over that shared universe;
     pool = pairs in the top fraction of both.
  2. Per species: log1p(CP10K) from raw counts, mean in 50 equal-count bins of the pseudotime rank
     (pd.qcut on rank(method="first")), z-scored per gene (population SD); r = Pearson(human, mouse z-profile).
  3. Centroid per species = sum_b max(z_b, 0) * b / (nbin - 1) / sum_b max(z_b, 0); gene centroid = mean of species.
  4. Keep r > --rmin; take the --n-early early and --n-late late genes with the highest r; order rows by centroid.
Reads the two mi_all.csv, the pairs table and the two objects (read only). Writes to --outdir the per-nucleus count
tables of the pool (input to the Monocle smoother), the pool concordance table, the selected genes and a JSON summary
that compares the selection with --panel-genes (the genes drawn in the panel).
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

N_BINS_SELECTION = 50   # bins of the concordance r used for the selection
N_BINS_CHECK = 20       # coarser bins, reported as profile_r_20bins
CP10K = 1e4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mi-dir", required=True, type=Path, help="directory holding one folder per MI arm")
    parser.add_argument("--mouse-arm", required=True)
    parser.add_argument("--human-arm", required=True)
    parser.add_argument("--pairs", required=True, type=Path, help="one-to-one pairs (rule mi_overlap)")
    parser.add_argument("--mouse-h5ad", required=True, type=Path, help="object the mouse arm scored")
    parser.add_argument("--mouse-pt-key", required=True, help="obs column the mouse arm scored")
    parser.add_argument("--human-h5ad", required=True, type=Path, help="object the human arm scored")
    parser.add_argument("--human-pt-key", required=True, help="obs column the human arm scored")
    parser.add_argument("--panel-genes", required=True, type=Path, help="genes drawn in the panel (csv, human)")
    parser.add_argument("--top-fraction", type=float, required=True)
    parser.add_argument("--rmin", type=float, required=True)
    parser.add_argument("--n-early", type=int, required=True)
    parser.add_argument("--n-late", type=int, required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def conserved_pool(mi_dir, arms, pairs, top_fraction):
    """Pairs measured in both arms, and the pool of pairs in the top fraction of MI in both species."""
    mi = {
        species: pd.read_csv(mi_dir / arm / "mi_all.csv").drop_duplicates("gene").set_index("gene")["mi"]
        for species, arm in arms.items()
    }
    pairs = pairs[pairs["human"].isin(mi["human"].index) & pairs["mouse"].isin(mi["mouse"].index)]
    pairs = pairs.reset_index(drop=True)
    mi_mouse = mi["mouse"].reindex(pairs["mouse"]).to_numpy()
    mi_human = mi["human"].reindex(pairs["human"]).to_numpy()
    is_finite = np.isfinite(mi_mouse) & np.isfinite(mi_human)
    pairs = pairs[is_finite].reset_index(drop=True)
    mi_mouse, mi_human = mi_mouse[is_finite], mi_human[is_finite]

    k_top = int(round(top_fraction * len(pairs)))
    in_top = (
        (pd.Series(mi_mouse).rank(ascending=False).to_numpy() <= k_top)
        & (pd.Series(mi_human).rank(ascending=False).to_numpy() <= k_top)
    )
    pool = pairs[in_top].assign(mi_mouse=mi_mouse[in_top], mi_human=mi_human[in_top]).reset_index(drop=True)
    return pairs, k_top, pool


def per_nucleus_counts(h5ad: Path, pt_key: str, genes: list[str]) -> pd.DataFrame:
    """Total counts, pseudotime and the raw counts of the pool genes, one row per nucleus with a pseudotime.

    Nuclei without a pseudotime (mouse decidual_apoptotic) are dropped, as in the MI arm; rows keep object order.
    """
    adata = ad.read_h5ad(h5ad, backed="r")
    pt = adata.obs[pt_key].to_numpy(dtype=float)
    has_pt = np.isfinite(pt)
    columns = pd.Index(adata.var_names.astype(str)).get_indexer(pd.Index(genes))
    if (columns < 0).any():
        raise ValueError(f"{h5ad}: {int((columns < 0).sum())} pool genes missing from var")

    counts = sp.csr_matrix(adata.layers["counts"])[np.flatnonzero(has_pt)]
    total = np.asarray(counts.sum(axis=1)).ravel()
    table = pd.DataFrame(
        counts[:, columns].toarray(),
        columns=[f"cnt_{gene}" for gene in genes],
        index=pd.Index(adata.obs_names[has_pt].astype(str), name="cell"),
    )
    adata.file.close()
    table.insert(0, "total_counts", total)
    table.insert(1, "pt_rank_norm", pt[has_pt])
    print(f"{h5ad}: {len(table):,} nuclei with {pt_key} ({int((~has_pt).sum()):,} without)", flush=True)
    return table


def profiles(per_nucleus: pd.DataFrame, genes: list[str], nbin: int):
    """Mean log1p(CP10K) per equal-count pseudotime bin (genes x bins) and its per-gene z-score."""
    counts = per_nucleus[[f"cnt_{gene}" for gene in genes]].to_numpy(float)
    expr = np.log1p(counts / per_nucleus["total_counts"].to_numpy(float)[:, None] * CP10K)
    bins = pd.qcut(per_nucleus["pt_rank_norm"].rank(method="first"), nbin, labels=False).to_numpy()
    means = np.vstack([expr[bins == b].mean(axis=0) for b in range(nbin)]).T
    sd = means.std(axis=1, ddof=0)
    sd[sd == 0] = 1
    return means, (means - means.mean(axis=1, keepdims=True)) / sd[:, None]


def centroid(z: np.ndarray) -> np.ndarray:
    """Pseudotime centroid on a 0-1 scale, weighted by the positive part of the z-scored profile."""
    weights = np.clip(z, 0, None)
    position = np.arange(z.shape[1]) / (z.shape[1] - 1)
    return (weights * position).sum(1) / weights.sum(1)


def rowwise_pearson(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.array([np.corrcoef(a[i], b[i])[0, 1] for i in range(len(a))])


def concordance_table(pool, per_nucleus_mouse, per_nucleus_human):
    """Profile concordance, centroid, peak, range and phase of every pool pair."""
    human_genes, mouse_genes = list(pool["human"]), list(pool["mouse"])
    Mh50, Zh50 = profiles(per_nucleus_human, human_genes, N_BINS_SELECTION)
    Mm50, Zm50 = profiles(per_nucleus_mouse, mouse_genes, N_BINS_SELECTION)
    _, Zh20 = profiles(per_nucleus_human, human_genes, N_BINS_CHECK)
    _, Zm20 = profiles(per_nucleus_mouse, mouse_genes, N_BINS_CHECK)

    table = pool.copy()
    table["r50"] = rowwise_pearson(Zh50, Zm50)
    table["r20"] = rowwise_pearson(Zh20, Zm20)
    table["centroid_human"] = centroid(Zh50)
    table["centroid_mouse"] = centroid(Zm50)
    table["pt_centroid"] = (table["centroid_human"] + table["centroid_mouse"]) / 2
    table["peak_h"] = Zh50.argmax(1)
    table["peak_m"] = Zm50.argmax(1)
    table["range_h"] = Mh50.max(1) - Mh50.min(1)
    table["range_m"] = Mm50.max(1) - Mm50.min(1)

    # Phase cutoffs as written in the Methods (< 0.4 early, > 0.6 late). The literals stay in this statement:
    # the methods-number registry reads them from the code (registry_methods.py, xs_heatmap_early/late).
    table["phase"] = np.where(table.pt_centroid < 0.4, "early", np.where(table.pt_centroid > 0.6, "late", "mid"))
    return table


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    arms = {"mouse": args.mouse_arm, "human": args.human_arm}
    pairs = pd.read_csv(args.pairs)[["human", "mouse"]]
    universe, k_top, pool = conserved_pool(args.mi_dir, arms, pairs, args.top_fraction)
    print(f"universe {len(universe):,} one-to-one pairs; top {args.top_fraction:.0%} = {k_top:,} per species; "
          f"pool {len(pool):,}", flush=True)

    per_nucleus_mouse = per_nucleus_counts(args.mouse_h5ad, args.mouse_pt_key, list(pool["mouse"]))
    per_nucleus_human = per_nucleus_counts(args.human_h5ad, args.human_pt_key, list(pool["human"]))
    per_nucleus_mouse.to_csv(args.outdir / "percell_mouse_pool.tsv.gz", sep="\t")
    per_nucleus_human.to_csv(args.outdir / "percell_human_pool.tsv.gz", sep="\t")

    table = concordance_table(pool, per_nucleus_mouse, per_nucleus_human)
    concordance_path = args.outdir / "conserved_profiles_pool_concordance.csv"
    table.sort_values("r50", ascending=False).to_csv(concordance_path, index=False)

    # Selection: the most concordant early and late genes above rmin, drawn in centroid order.
    concordant = table[table["r50"] > args.rmin]
    picked = pd.concat([
        concordant[concordant["phase"] == "early"].nlargest(args.n_early, "r50"),
        concordant[concordant["phase"] == "late"].nlargest(args.n_late, "r50"),
    ]).sort_values("pt_centroid")
    selected = picked.rename(columns={"r50": "profile_r_50bins", "r20": "profile_r_20bins"})[
        ["human", "mouse", "profile_r_50bins", "profile_r_20bins", "phase", "pt_centroid", "peak_h", "peak_m"]
    ]
    selected.to_csv(args.outdir / "conserved_profiles_selected_genes.csv", index=False)

    panel = pd.read_csv(args.panel_genes)
    summary = {
        "mi_arms": arms,
        "universe": int(len(universe)),
        "k_top": k_top,
        "pool": int(len(pool)),
        "n_r_gt": int(len(concordant)),
        "rmin": args.rmin,
        "phase_counts_r_gt": concordant["phase"].value_counts().to_dict(),
        "median_r50_pool": round(float(table["r50"].median()), 3),
        "n_selected": int(len(selected)),
        "r50_range_selected": [
            round(float(selected["profile_r_50bins"].min()), 3),
            round(float(selected["profile_r_50bins"].max()), 3),
        ],
        "centroid_definition": "per species: sum_b max(z_b,0)*b/(nbin-1) / sum_b max(z_b,0) on 50-bin z-profiles; "
                               "gene = mean of species",
        "panel_genes": str(args.panel_genes),
        "overlap_with_previous_28": int(len(set(selected["human"]) & set(panel["human"]))),
        "kept": sorted(set(selected["human"]) & set(panel["human"])),
        "dropped": sorted(set(panel["human"]) - set(selected["human"])),
        "added": sorted(set(selected["human"]) - set(panel["human"])),
        "previous_genes_in_new_pool": int(panel["human"].isin(pool["human"]).sum()),
    }
    (args.outdir / "selection_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1), flush=True)


if __name__ == "__main__":
    main()
