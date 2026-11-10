#!/usr/bin/env python3
"""Step 3: compare the X-Mens decidual vs fibroblast calls with the reference gene lists.

Input: xmens_broad_all_genes.csv from step 2, and the reference lists:
  consensus      the 123-gene human in vitro consensus and its mouse orthology table (corrected list of 22 October 2026)
  Wang et al.    the published direction-only lists of Wang et al. 2020: oil-induced decidualization, natural pregnancy
                 (day 8) and in vitro decidualization

A gene counts as called in X-Mens when it passed the expression filter, its fold change is estimable, its adjusted
p-value is below 0.05 and its absolute unshrunk log2 fold change is above 1. Each comparison counts the reference genes
that X-Mens calls and how many of those change in the same direction.

How genes are matched, as in the release of 19 October 2026:
  consensus         through strict one-to-one Ensembl 116 orthologs (114 of the 123 genes), by mouse Ensembl gene ID;
                    built by the release's own function (release_comparison_functions.build_consensus_ortholog_support)
  Wang in vitro     by current symbol or source feature ID, choosing among several candidates as the release did
                    (release_comparison_functions.build_wang_support)
  Wang oil and      by symbol, among X-Mens genes that passed the filter, are estimable and have a unique symbol
  pregnancy

Expected by chance: the agreement two unrelated lists would show, given the share of up genes on each side among the
genes compared: p_up_reference * p_up_xmens + (1 - p_up_reference) * (1 - p_up_xmens).

Writes into --out:
  consensus_ortholog_xmens_support_rows.csv   one row per consensus gene with a one-to-one ortholog (114)
  wang_ivd_support_rows.csv                   one row per gene of Wang's in vitro list (1,501)
  genes_xmens_vs_<reference>.csv              the genes compared, for consensus, oil, pregnancy and in vitro
  SUMMARY_xmens_comparisons.csv               genes compared, same direction, % same direction, % expected by chance

Usage
  python compare_xmens_with_references.py --xmens xmens_broad_all_genes.csv --consensus CONSENSUS.csv
         --orthology ORTHOLOGY.csv --wang-dir WANG_DIR --out OUT_DIR

  WANG_DIR  wang_oil.csv, wang_pregnancy.csv, wang_ivd.csv (columns gene, direction)
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# release_comparison_functions.py sits next to this script; import it from there whatever the working directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import release_comparison_functions as release  # noqa: E402

CALL_RULE = {"padj_lt": 0.05, "abs_unshrunk_lfc_gt": 1.0}
N_CONSENSUS_ORTHOLOGS = 114
LABELS = {
    "consensus": "X-Mens vs human consensus in vitro (123 genes, 114 with mouse ortholog)",
    "oil": "X-Mens vs Wang oil-induced decidualization",
    "pregnancy": "X-Mens vs Wang natural pregnancy (day 8)",
    "in_vitro": "X-Mens vs Wang in vitro decidualization",
}


def read_xmens(path):
    """The all-genes table, with floats read back exactly as written."""
    return pd.read_csv(path, keep_default_na=False, na_values=[""], float_precision="round_trip")


def percent_expected_by_chance(reference_up, xmens_up):
    p, q = np.mean(reference_up), np.mean(xmens_up)
    return 100 * (p * q + (1 - p) * (1 - q))


def summary_row(label, reference_direction, xmens_direction):
    reference_direction, xmens_direction = np.asarray(reference_direction), np.asarray(xmens_direction)
    same = int((reference_direction == xmens_direction).sum())
    return {"comparison": label, "genes compared": len(reference_direction), "same direction": same,
            "% same direction": round(100 * same / len(reference_direction), 1),
            "% expected by chance": round(percent_expected_by_chance(reference_direction == "up", xmens_direction == "up"), 1)}


def compare_consensus(xmens, consensus_path, orthology_path):
    """Support rows for the consensus orthologs, and the called genes."""
    rows = release.build_consensus_ortholog_support(Path(consensus_path), Path(orthology_path), xmens, CALL_RULE,
                                                    N_CONSENSUS_ORTHOLOGS)
    rows["source_id"] = xmens.source_id.iloc[0]         # the release's function writes its own design's ID
    called = rows[rows.support_status.isin(["supported", "reversed"])].sort_values("human_gene")
    genes = pd.DataFrame({"gene": called.human_gene, "mouse gene": called.mouse_gene,
                          "in vitro (human consensus)": called.consensus_direction, "X-Mens": called.observed_direction,
                          "X-Mens log2 fold change": called.log2FoldChange_unshrunk.round(3),
                          "X-Mens adjusted p": called.padj_for_call.map(lambda p: float(f"{p:.3g}")),
                          "same direction": called.support_status.eq("supported")})
    return rows, genes, summary_row(LABELS["consensus"], called.consensus_direction, called.observed_direction)


def compare_wang_by_symbol(xmens, wang_list, name):
    """Wang oil or pregnancy list against the X-Mens genes with a unique symbol."""
    testable = xmens[xmens.production_eligible.astype(bool) & xmens.estimable.astype(bool)]
    testable = testable.drop_duplicates("gene_symbol", keep=False).set_index("gene_symbol")
    called = testable[testable.significant_fdr05_lfc1.astype(bool)]
    reference = wang_list.drop_duplicates("gene", keep=False).set_index("gene").direction
    shared = sorted(set(reference.index) & set(called.index))
    lfc = called.log2FoldChange_unshrunk.reindex(shared)
    genes = pd.DataFrame({"gene": shared, f"{name} (Wang)": reference.reindex(shared).to_numpy(),
                          "X-Mens": np.where(lfc > 0, "up", "down"), "X-Mens log2 fold change": lfc.round(3).to_numpy(),
                          "X-Mens adjusted p": called.padj_for_call.reindex(shared).map(lambda p: float(f"{p:.3g}")).to_numpy()})
    genes["same direction"] = genes[f"{name} (Wang)"] == genes["X-Mens"]
    return genes, summary_row(LABELS[name], genes[f"{name} (Wang)"], genes["X-Mens"])


def compare_wang_in_vitro(xmens, wang_dir):
    """Wang's in vitro list, matched as in the release; the oil and pregnancy directions ride along."""
    rows = release.build_wang_support(Path(f"{wang_dir}/wang_ivd.csv"), Path(f"{wang_dir}/wang_oil.csv"),
                                      Path(f"{wang_dir}/wang_pregnancy.csv"), xmens, CALL_RULE)
    called = rows[rows.xmens_called_strict_lfc1.astype(bool)].sort_values("wang_ivd_gene")
    matched = called.xmens_current_gene_symbol.where(called.xmens_current_gene_symbol != called.wang_ivd_gene, "")
    # The release's function returns object columns (it mixes in pd.NA for unmatched genes); make them numeric.
    lfc = pd.to_numeric(called.xmens_effect_log2_fold_change_unshrunk)
    padj = pd.to_numeric(called.xmens_adjusted_p_value_for_call)
    genes = pd.DataFrame({"gene": called.wang_ivd_gene, "in vitro (Wang)": called.wang_ivd_direction,
                          "X-Mens gene": matched, "X-Mens": called.xmens_direction,
                          "X-Mens log2 fold change": lfc.round(3),
                          "X-Mens adjusted p": padj.map(lambda p: float(f"{p:.3g}"))})
    genes["same direction"] = genes["in vitro (Wang)"] == genes["X-Mens"]
    return rows, genes, summary_row(LABELS["in_vitro"], genes["in vitro (Wang)"], genes["X-Mens"])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xmens", required=True, help="xmens_broad_all_genes.csv from step 2")
    parser.add_argument("--consensus", required=True, help="consensus123_stadtmauer_authors_criteria.csv")
    parser.add_argument("--orthology", required=True, help="consensus123_stadtmauer_authors_criteria_mouse_orthology.csv")
    parser.add_argument("--wang-dir", required=True, help="folder with wang_oil.csv, wang_pregnancy.csv, wang_ivd.csv")
    parser.add_argument("--out", required=True, help="output folder")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    xmens = read_xmens(args.xmens)
    consensus_rows, consensus_genes, consensus_summary = compare_consensus(xmens, args.consensus, args.orthology)
    oil_genes, oil_summary = compare_wang_by_symbol(xmens, pd.read_csv(f"{args.wang_dir}/wang_oil.csv"), "oil")
    pregnancy_genes, pregnancy_summary = compare_wang_by_symbol(xmens, pd.read_csv(f"{args.wang_dir}/wang_pregnancy.csv"),
                                                                "pregnancy")
    in_vitro_rows, in_vitro_genes, in_vitro_summary = compare_wang_in_vitro(xmens, args.wang_dir)

    consensus_rows.to_csv(f"{args.out}/consensus_ortholog_xmens_support_rows.csv", index=False)
    in_vitro_rows.to_csv(f"{args.out}/wang_ivd_support_rows.csv", index=False)
    for name, genes in (("consensus", consensus_genes), ("oil", oil_genes), ("pregnancy", pregnancy_genes),
                        ("in_vitro", in_vitro_genes)):
        genes.to_csv(f"{args.out}/genes_xmens_vs_{name}.csv", index=False)
    summary = pd.DataFrame([consensus_summary, oil_summary, pregnancy_summary, in_vitro_summary])
    summary.to_csv(f"{args.out}/SUMMARY_xmens_comparisons.csv", index=False)
    print(summary.to_string(index=False))
    print("consensus support:", consensus_rows.support_status.value_counts().to_dict())


if __name__ == "__main__":
    main()
