#!/usr/bin/env python3
"""Step 4: write the human decidual vs fibroblast DE tables, in the column format of the X-Mens release's tables.

There is one row for every gene in the object. Genes that did not pass the expression filter have
production_eligible = False and empty statistics. Besides the fit's statistics, each row has:

  log2_cpm_ratio_<donor>_<phase>      decidual vs fibroblast within one donor:
                                      log2((CPM decidual + 0.5) / (CPM fibroblast + 0.5)).
                                      Filled in only for genes that passed the filter.
  pct_detected_<donor>_<phase>_<arm>  percent of nuclei with at least one UMI of the gene

Outputs
  human_decidual_vs_fibroblast_all_genes.csv                        fit B, `cell_type` labels (primary)
  human_decidual_vs_fibroblast_all_genes_coarse_refined_labels.csv  fit C, `coarse_refined` labels (sensitivity)

Usage
  python build_human_de_tables.py --pb-dir PB_DIR --fits FIT_DIR --out OUT_DIR

  PB_DIR   human_pb_meta.csv, human_pb_counts.csv.gz, human_pb_detect.csv.gz, human_var.csv (step 1)
  FIT_DIR  DE_B_new_celltype_release.csv.gz, DE_C_new_coarse_release.csv.gz, fit_diagnostics.json (step 3)
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

OBJECT = "xmens_final/data/deliverables/2026_09_15/Cevrim_Human_Only_Fib_Dec_Cells.h5ad (layers/counts)"
PSEUDOCOUNT_CPM = 0.5

TABLES = [  # fit, label column, output file
    ("B_new_celltype_release", "cell_type", "human_decidual_vs_fibroblast_all_genes.csv"),
    ("C_new_coarse_release", "coarse_refined", "human_decidual_vs_fibroblast_all_genes_coarse_refined_labels.csv"),
]


def donor_pseudobulks(pseudobulks, label_column):
    """The four donor x arm pseudobulks built from one label column."""
    return pseudobulks[(pseudobulks.label_source == label_column) & (pseudobulks.unit_type == "donor")]


def per_donor_log2_ratios(pseudobulks, counts, label_column):
    """Decidual vs fibroblast log2 CPM ratio within each donor, one column per donor."""
    selected = donor_pseudobulks(pseudobulks, label_column)
    selected_counts = counts.loc[selected.index]
    cpm = selected_counts.div(selected_counts.sum(axis=1), axis=0) * 1e6
    ratios = {}
    for donor in sorted(selected.donor.unique()):
        of_donor = selected.donor == donor
        decidual = selected.index[of_donor & (selected.arm == "decidual")][0]
        fibroblast = selected.index[of_donor & (selected.arm == "fibroblast")][0]
        phase = selected.loc[decidual, "timepoint"]
        ratios[f"log2_cpm_ratio_{donor}_{phase}"] = np.log2((cpm.loc[decidual] + PSEUDOCOUNT_CPM)
                                                            / (cpm.loc[fibroblast] + PSEUDOCOUNT_CPM))
    return pd.DataFrame(ratios)


def build_table(fit_name, label_column, pseudobulks, counts, detected, symbols, diagnostics, fits_dir):
    fitted = pd.read_csv(f"{fits_dir}/DE_{fit_name}.csv.gz", index_col=0, float_precision="round_trip")

    # Per-donor ratios, for the genes that were fitted
    ratios = per_donor_log2_ratios(pseudobulks, counts, label_column)
    for column in ratios.columns:
        fitted[column] = ratios[column].reindex(fitted.index).round(3)

    # One row per gene in the object; genes that were not fitted get empty statistics
    was_fitted = pd.Series(False, index=symbols.index)
    was_fitted[fitted.index] = True
    table = pd.DataFrame({"ensembl_gene_id": symbols.index,
                          "gene_symbol": symbols.values,
                          "production_eligible": was_fitted.values})
    statistics = fitted.drop(columns=["gene_symbol"]).rename(columns={"log2FoldChange": "log2FoldChange_unshrunk",
                                                                      "call": "significant_fdr05_lfc1"})
    table = table.merge(statistics, left_on="ensembl_gene_id", right_index=True, how="left")
    table["significant_fdr05_lfc1"] = table.significant_fdr05_lfc1.eq(True)    # empty for unfitted genes: False
    table["estimable"] = table.estimable.eq(True)
    table["padj_for_call"] = table.padj.where(np.isfinite(table.padj), 1.0)

    # Detection in each donor x arm pseudobulk
    for pseudobulk_id, pseudobulk in donor_pseudobulks(pseudobulks, label_column).iterrows():
        percent = 100 * detected.loc[pseudobulk_id] / pseudobulk.n_cells
        column = f"pct_detected_{pseudobulk.donor}_{pseudobulk.timepoint}_{pseudobulk.arm}"
        table[column] = percent.reindex(table.ensembl_gene_id).round(3).to_numpy()

    # Where the table came from, in the release's columns
    table["source_id"] = f"human_decidual_vs_fibroblast_{label_column}_by_donor"
    table["label_column"] = label_column
    table["positive_lfc_meaning"] = "decidual_minus_fibroblast"
    table["analysis_method"] = "pydeseq2_wald_fixed_block"
    table["object"] = OBJECT

    assert int(table.significant_fdr05_lfc1.sum()) == diagnostics["calls"]
    assert int(table.production_eligible.sum()) == diagnostics["eligible"]
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pb-dir", required=True, help="folder with step 1's outputs")
    parser.add_argument("--fits", required=True, help="folder with step 3's outputs")
    parser.add_argument("--out", required=True, help="output folder")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    pseudobulks = pd.read_csv(f"{args.pb_dir}/human_pb_meta.csv").set_index("id")
    counts = pd.read_csv(f"{args.pb_dir}/human_pb_counts.csv.gz", index_col=0)
    detected = pd.read_csv(f"{args.pb_dir}/human_pb_detect.csv.gz", index_col=0)
    symbols = pd.read_csv(f"{args.pb_dir}/human_var.csv").set_index("gene_id").gene_symbol
    with open(f"{args.fits}/fit_diagnostics.json") as handle:
        diagnostics = json.load(handle)

    for fit_name, label_column, filename in TABLES:
        table = build_table(fit_name, label_column, pseudobulks, counts, detected, symbols,
                            diagnostics[fit_name], args.fits)
        table.to_csv(f"{args.out}/{filename}", index=False)
        print(filename, table.shape, "| eligible", int(table.production_eligible.sum()),
              "| called", int(table.significant_fdr05_lfc1.sum()))


if __name__ == "__main__":
    main()
