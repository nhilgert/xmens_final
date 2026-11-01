"""Cross-species overlap of the pseudotime MI arms (rule mi_overlap, 11_mi_pseudotime.smk).

Extracted from notebooks/mi_dpt/mi_dpt_overlap_harmonized_nall.ipynb (sections 1-7 and 9; the figures, the
trajectory clustering of section 8 and the regulator tables of section 10 are not ported). Reads the mi_all.csv,
manifest.json and build_meta.json of the mouse and human main arms and of the mouse control arms, and the Ensembl
Compara component table. Writes to --outdir:
  ortholog_pairs_one2one.csv       strict one-to-one protein-coding pairs (human, mouse)
  overlap_summary.csv              top-decile overlap, hypergeometric p and Spearman rho per comparison
  overlap_permutation_null.csv     overlap under shuffled ortholog pairing (main comparisons)
  overlap_cutoff_sweep.csv         the same overlap at top fractions 1-30%
  shared_top_decile_genes.json     first 40 genes in both top deciles, per comparison
  specificity_shared_universe.csv  mouse arms against human on the pairs measurable in every arm
  methods_overlap.json             Methods overlap: N, K, n, k, expected, fold, hypergeometric p (log scale)
  overlap_percentile_sweep.csv     Methods overlap at every percentile cutoff 40-99 (Fig. S16 C-D)
  overlap_report.json              run records, ortholog source, overlap and null records
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

CUTOFF_FRACTIONS = (0.01, 0.02, 0.05, 0.10, 0.20, 0.30)   # top fractions of the cutoff sweep (notebook section 5)
SHARED_TOP_N = (100, 200, 500, 1000, 2000)                 # fixed cutoffs on the shared universe (section 9)
SHARED_FRACTIONS = (0.02, 0.05, 0.10, 0.20)
PERCENTILES = range(40, 100)                               # Fig. S16 C-D: every percentile cutoff, 40th to 99th
N_SHARED_GENES = 40                                        # genes listed per comparison in the JSON


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mi-dir", required=True, type=Path, help="directory holding one folder per MI arm")
    parser.add_argument("--mouse-arm", required=True, help="mouse main arm (global pseudotime rank)")
    parser.add_argument("--human-arm", required=True, help="human main arm")
    parser.add_argument("--control-arms", nargs="+", required=True, help="mouse arms: uninduced, induced_matched")
    parser.add_argument("--orthologs", required=True, type=Path, help="Ensembl Compara component table (tsv)")
    parser.add_argument("--top-fraction", type=float, required=True, help="top fraction for the overlap (0.10)")
    parser.add_argument("--n-perm", type=int, required=True, help="shuffled pairings for the permutation null")
    parser.add_argument("--seed", type=int, required=True, help="seed of the permutation null")
    parser.add_argument("--q-thresh", type=float, required=True, help="MI significance (Methods: q < 0.05)")
    parser.add_argument("--top-percentile", type=float, required=True, help="Methods top set (90th percentile)")
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def load_arm(arm_dir: Path) -> tuple[pd.DataFrame, dict]:
    """MI table indexed by gene, and the arm's build and prep records."""
    table = pd.read_csv(arm_dir / "mi_all.csv")
    if not {"gene", "mi", "p", "q"} <= set(table.columns):
        raise ValueError(f"{arm_dir}/mi_all.csv: unexpected columns {list(table.columns)}")
    table = table.drop_duplicates("gene").set_index("gene")

    build_path = arm_dir / "build_meta.json"
    build = json.loads(build_path.read_text()) if build_path.exists() else {}
    prep = json.loads((arm_dir / "manifest.json").read_text())
    return table, {"build": build, "prep": prep}


def one_to_one_pairs(components: pd.DataFrame) -> pd.DataFrame:
    """Strict one-to-one protein-coding (human, mouse) pairs, pivoted on component_id."""
    needed = {"component_id", "species", "gene_symbol", "homology_type", "biotype"}
    if not needed <= set(components.columns):
        raise ValueError(f"ortholog table lacks {sorted(needed - set(components.columns))}")
    if pd.api.types.is_numeric_dtype(components["gene_symbol"]):
        raise ValueError("gene_symbol is numeric: the table's columns are not the expected ones")

    one2one = components[
        (components["homology_type"] == "ortholog_one2one") & (components["biotype"] == "protein_coding")
    ]
    genes_per_component = one2one.groupby(["component_id", "species"]).size()
    if not (genes_per_component == 1).all():
        raise ValueError("a one-to-one component has more than one gene in a species")

    wide = one2one.pivot_table(index="component_id", columns="species", values="gene_symbol", aggfunc="first")
    pairs = wide.dropna()[["human", "mouse"]].drop_duplicates()
    is_unique = ~pairs["human"].duplicated(keep=False) & ~pairs["mouse"].duplicated(keep=False)
    return pairs[is_unique]


def aligned_mi(pairs, tables, species, arm_a, arm_b):
    """MI of the two arms over the genes (same species) or pairs (cross-species) measured in both, finite only."""
    table_a, table_b = tables[arm_a], tables[arm_b]
    if species[arm_a] == species[arm_b]:
        genes_a = genes_b = pd.Index(table_a.index.intersection(table_b.index))
    else:
        human_table = table_b if species[arm_b] == "human" else table_a
        mouse_table = table_a if species[arm_a] == "mouse" else table_b
        measured = pairs[pairs["human"].isin(human_table.index)]
        measured = measured[measured["mouse"].isin(mouse_table.index)]
        genes_a = pd.Index(measured["mouse"] if species[arm_a] == "mouse" else measured["human"])
        genes_b = pd.Index(measured["human"] if species[arm_b] == "human" else measured["mouse"])

    mi_a = table_a["mi"].reindex(genes_a).to_numpy()
    mi_b = table_b["mi"].reindex(genes_b).to_numpy()
    is_finite = np.isfinite(mi_a) & np.isfinite(mi_b)
    return mi_a[is_finite], mi_b[is_finite], genes_a[is_finite]


def compare(pairs, tables, species, arm_a, arm_b, top_fraction):
    """Top-fraction overlap of two arms against the hypergeometric expectation (notebook section 3)."""
    mi_a, mi_b, genes_a = aligned_mi(pairs, tables, species, arm_a, arm_b)
    n = len(mi_a)
    if n < 100:
        return None

    k = int(round(top_fraction * n))
    top_a = set(np.argsort(-mi_a)[:k])
    top_b = set(np.argsort(-mi_b)[:k])
    overlap = len(top_a & top_b)
    expected = k * k / n
    return {
        "a": arm_a,
        "b": arm_b,
        "universe": n,
        "k": k,
        "overlap": overlap,
        "expected": round(expected, 1),
        "fold": round(overlap / expected, 3) if expected else np.nan,
        "p": float(hypergeom.sf(overlap - 1, n, k, k)),
        "rho": round(float(spearmanr(mi_a, mi_b).statistic), 4),
        "shared_genes_a": [str(genes_a[i]) for i in sorted(top_a & top_b)][:N_SHARED_GENES],
    }


def permutation_null(pairs, tables, species, arm_a, arm_b, top_fraction, n_perm, seed):
    """Overlap when the ortholog pairing is shuffled: keeps each species' MI distribution (section 4)."""
    mi_a, mi_b, _ = aligned_mi(pairs, tables, species, arm_a, arm_b)
    n = len(mi_a)
    k = int(round(top_fraction * n))
    top_a = set(np.argsort(-mi_a)[:k])
    observed = len(top_a & set(np.argsort(-mi_b)[:k]))

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, int)
    for i in range(n_perm):
        null[i] = len(top_a & set(np.argsort(-rng.permutation(mi_b))[:k]))

    return {
        "a": arm_a,
        "b": arm_b,
        "observed": observed,
        "null_mean": round(float(null.mean()), 1),
        "null_sd": round(float(null.std()), 1),
        "null_max": int(null.max()),
        "z": round(float((observed - null.mean()) / max(null.std(), 1e-9)), 2),
        "p_emp": float((null >= observed).mean()),
        "n_perm": n_perm,
    }


def shared_universe_specificity(pairs, tables, human_arm, mouse_arms):
    """Each mouse arm against human on the pairs measurable in every arm, fixed and fractional cutoffs (sec. 9)."""
    universe = pairs[pairs["human"].isin(tables[human_arm].index)]
    for arm in mouse_arms:
        universe = universe[universe["mouse"].isin(tables[arm].index)]
    mi_human = tables[human_arm]["mi"].reindex(pd.Index(universe["human"])).to_numpy()
    mi_mouse = {arm: tables[arm]["mi"].reindex(pd.Index(universe["mouse"])).to_numpy() for arm in mouse_arms}
    is_finite = np.isfinite(mi_human)
    for arm in mouse_arms:
        is_finite &= np.isfinite(mi_mouse[arm])
    mi_human = mi_human[is_finite]
    mi_mouse = {arm: values[is_finite] for arm, values in mi_mouse.items()}
    n_universe = int(is_finite.sum())

    fixed = [k for k in SHARED_TOP_N if k < n_universe // 2]
    fractional = [int(round(f * n_universe)) for f in SHARED_FRACTIONS]
    rows = []
    for k in sorted(set(fixed + fractional)):
        top_human = set(np.argsort(-mi_human)[:k])
        for arm in mouse_arms:
            overlap = len(top_human & set(np.argsort(-mi_mouse[arm])[:k]))
            expected = k * k / n_universe
            rows.append({
                "arm": arm,
                "k": k,
                "N": n_universe,
                "overlap": overlap,
                "expected": round(expected, 1),
                "fold": round(overlap / expected, 3) if expected else np.nan,
                "p": float(hypergeom.sf(overlap - 1, n_universe, k, k)),
                "rho": round(float(spearmanr(mi_human, mi_mouse[arm]).statistic), 4),
            })
    table = pd.DataFrame(rows)
    table["cutoff"] = np.where(
        table["k"].isin(fixed),
        "top-" + table["k"].astype(str),
        (100 * table["k"] / n_universe).round(0).astype(int).astype(str) + "%",
    )
    return table


def methods_overlap(pairs, mi_mouse, mi_human, q_thresh, percentile):
    """Methods overlap: pairs with both genes significant (N); top sets above the percentile of significant MI.

    K and n are the mouse and human top genes among the N pairs, k the pairs in both; p is the upper tail
    P(X >= k) of the hypergeometric distribution, computed on the log scale (hypergeom.logsf) so it cannot underflow.
    """
    sets = {}
    for name, table in (("mouse", mi_mouse), ("human", mi_human)):
        significant = table[table["q"] < q_thresh]
        threshold = float(np.percentile(significant["mi"], percentile))
        sets[name] = {
            "tested": len(table),
            "significant": len(significant),
            "threshold_nats": threshold,
            "top": set(significant.index[significant["mi"] > threshold]),
            "significant_genes": set(significant.index),
        }

    is_significant = (
        pairs["mouse"].isin(sets["mouse"]["significant_genes"])
        & pairs["human"].isin(sets["human"]["significant_genes"])
    )
    background = pairs[is_significant]
    is_top_mouse = background["mouse"].isin(sets["mouse"]["top"])
    is_top_human = background["human"].isin(sets["human"]["top"])
    n_pairs = len(background)
    top_mouse = int(is_top_mouse.sum())
    top_human = int(is_top_human.sum())
    in_both = int((is_top_mouse & is_top_human).sum())
    expected = top_mouse * top_human / n_pairs
    log_p = float(hypergeom.logsf(in_both - 1, n_pairs, top_mouse, top_human))

    record = {
        "q_thresh": q_thresh,
        "top_percentile": percentile,
        "N": n_pairs,
        "K": top_mouse,
        "n": top_human,
        "k": in_both,
        "expected": expected,
        "fold": in_both / expected,
        "log10_p": log_p / np.log(10),
    }
    for name, values in sets.items():
        record[name] = {
            "tested": values["tested"],
            "significant": values["significant"],
            "threshold_nats": values["threshold_nats"],
            "threshold_bits": values["threshold_nats"] / np.log(2),
            "top": len(values["top"]),
        }
    return record


def percentile_sweep(pairs, mi_mouse, mi_human, q_thresh):
    """Methods overlap at every percentile cutoff, with the hypergeometric 2.5-97.5% interval of the fold."""
    sig_mouse = mi_mouse[mi_mouse["q"] < q_thresh]
    sig_human = mi_human[mi_human["q"] < q_thresh]
    background = pairs[pairs["mouse"].isin(sig_mouse.index) & pairs["human"].isin(sig_human.index)]
    n_pairs = len(background)

    rows = []
    for percentile in PERCENTILES:
        thr_mouse = float(np.percentile(sig_mouse["mi"], percentile))
        thr_human = float(np.percentile(sig_human["mi"], percentile))
        top_mouse = set(sig_mouse.index[sig_mouse["mi"] > thr_mouse])
        top_human = set(sig_human.index[sig_human["mi"] > thr_human])
        is_top_mouse = background["mouse"].isin(top_mouse)
        is_top_human = background["human"].isin(top_human)
        K = int(is_top_mouse.sum())
        n = int(is_top_human.sum())
        k = int((is_top_mouse & is_top_human).sum())
        expected = K * n / n_pairs
        rows.append({
            "percentile": percentile,
            "thr_mouse": thr_mouse,
            "thr_human": thr_human,
            "top_mouse_genomewide": len(top_mouse),
            "top_human_genomewide": len(top_human),
            "N": n_pairs,
            "K_mouse": K,
            "n_human": n,
            "overlap": k,
            "expected": expected,
            "fold": k / expected,
            "neglog10_p": -float(hypergeom.logsf(k - 1, n_pairs, K, n)) / np.log(10),
            "null_fold_lo": hypergeom.ppf(0.025, n_pairs, K, n) / expected,
            "null_fold_hi": hypergeom.ppf(0.975, n_pairs, K, n) / expected,
            "pct_of_human_top": 100 * k / n,
        })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Arms and their species; the comparisons are those of the notebook, main ones first.
    arms = [*args.control_arms, args.mouse_arm, args.human_arm]
    species = {arm: "mouse" for arm in [*args.control_arms, args.mouse_arm]} | {args.human_arm: "human"}
    unind_arm, indmatch_arm = args.control_arms
    comparisons = [
        (unind_arm, args.human_arm),
        (indmatch_arm, args.human_arm),
        (args.mouse_arm, args.human_arm),
        (unind_arm, indmatch_arm),
        (unind_arm, args.mouse_arm),
    ]

    tables, records = {}, {}
    for arm in arms:
        tables[arm], records[arm] = load_arm(args.mi_dir / arm)
        print(f"{arm}: {len(tables[arm]):,} genes", flush=True)

    components = pd.read_csv(args.orthologs, sep="\t")
    pairs = one_to_one_pairs(components)
    pairs.to_csv(args.outdir / "ortholog_pairs_one2one.csv", index=False)
    print(f"{len(pairs):,} one-to-one pairs from {args.orthologs}", flush=True)

    # Top-decile overlap per comparison and the list of shared genes.
    overlap_rows = [compare(pairs, tables, species, a, b, args.top_fraction) for a, b in comparisons]
    overlap = pd.DataFrame([row for row in overlap_rows if row])
    overlap.drop(columns=["shared_genes_a"]).to_csv(args.outdir / "overlap_summary.csv", index=False)
    shared = {f"{row['a']} vs {row['b']}": row["shared_genes_a"] for row in overlap_rows if row}
    (args.outdir / "shared_top_decile_genes.json").write_text(json.dumps(shared, indent=1))

    # Shuffled-pairing null for the three comparisons against human; each starts from the same seed.
    null_rows = [
        permutation_null(pairs, tables, species, a, b, args.top_fraction, args.n_perm, args.seed)
        for a, b in comparisons[:3]
    ]
    null = pd.DataFrame(null_rows)
    null.to_csv(args.outdir / "overlap_permutation_null.csv", index=False)
    if not (null["observed"] > null["null_max"]).all():
        raise SystemExit(f"observed overlap within the shuffled-pairing null:\n{null.to_string()}")

    sweep_rows = []
    for fraction in CUTOFF_FRACTIONS:
        for a, b in comparisons:
            row = compare(pairs, tables, species, a, b, fraction)
            if row:
                sweep_rows.append({
                    "frac": fraction,
                    "pair": f"{a} vs {b}",
                    "fold": row["fold"],
                    "overlap": row["overlap"],
                    "k": row["k"],
                    "p": row["p"],
                })
    pd.DataFrame(sweep_rows).to_csv(args.outdir / "overlap_cutoff_sweep.csv", index=False)

    mouse_arms = [*args.control_arms, args.mouse_arm]
    specificity = shared_universe_specificity(pairs, tables, args.human_arm, mouse_arms)
    specificity.to_csv(args.outdir / "specificity_shared_universe.csv", index=False)

    # Methods numbers: overlap of the significant top sets, and the same at every percentile (Fig. S16 C-D).
    methods = methods_overlap(
        pairs, tables[args.mouse_arm], tables[args.human_arm], args.q_thresh, args.top_percentile
    )
    methods["arms"] = {"mouse": args.mouse_arm, "human": args.human_arm}
    (args.outdir / "methods_overlap.json").write_text(json.dumps(methods, indent=1) + "\n")
    sweep = percentile_sweep(pairs, tables[args.mouse_arm], tables[args.human_arm], args.q_thresh)
    sweep.to_csv(args.outdir / "overlap_percentile_sweep.csv", index=False)
    print(
        f"Methods overlap: N {methods['N']:,}, K {methods['K']:,}, n {methods['n']:,}, k {methods['k']:,}, "
        f"fold {methods['fold']:.3f}, log10 p {methods['log10_p']:.4f}",
        flush=True,
    )

    runs = {
        arm: {
            **records[arm]["build"],
            "prep_layer": records[arm]["prep"].get("layer"),
            "prep_n_genes": records[arm]["prep"].get("n_genes"),
            "prep_n_chunks": records[arm]["prep"].get("n_chunks"),
        }
        for arm in arms
    }
    release = sorted(components["source_release"].unique()) if "source_release" in components else None
    report = {
        "runs": runs,
        "ortholog_source": str(args.orthologs),
        "ortholog_sha256": hashlib.sha256(args.orthologs.read_bytes()).hexdigest(),
        "ortholog_release": release,
        "n_one2one_pairs": int(len(pairs)),
        "topfrac": args.top_fraction,
        "overlap": overlap.drop(columns=["shared_genes_a"]).to_dict("records"),
        "permutation_null": null.to_dict("records"),
        "methods_overlap": {key: value for key, value in methods.items() if key not in ("mouse", "human")},
    }
    (args.outdir / "overlap_report.json").write_text(json.dumps(report, indent=1, default=str))
    print(f"wrote {args.outdir}", flush=True)


if __name__ == "__main__":
    main()
