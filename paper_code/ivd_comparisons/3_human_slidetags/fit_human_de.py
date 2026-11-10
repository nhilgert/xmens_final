#!/usr/bin/env python3
"""Step 3: pseudobulk differential expression, human decidual vs fibroblast nuclei (PyDESeq2).

Input: the pseudobulks that cluster_jobs/human_pseudobulk_job.py built from
xmens_final/data/deliverables/2026_09_15/Cevrim_Human_Only_Fib_Dec_Cells.h5ad (layers/counts).

Every fit is a Wald test of decidual vs fibroblast. A positive log2 fold change means higher in decidual nuclei.
A gene is called when its adjusted p-value is below 0.05 and its absolute log2 fold change is above 1.

The fits
  B_new_celltype_release             PRIMARY. `cell_type` labels, one pseudobulk per donor and arm (4 in total),
                                     design ~ donor + group, release settings.
  C_new_coarse_release               As B, with `coarse_refined` labels (sensitivity analysis).
  D_new_celltype_release_by_section  As B, with one pseudobulk per Slide-tags section and arm (12), ~ sample + group.
                                     Sections from one donor are not independent replicates; for reference only.
  E_new_coarse_release_by_section    As D, with `coarse_refined` labels.
  A_new_celltype_old_settings        As B, with the old notebook's settings. A vs the old table shows the effect of
                                     the new object; A vs B shows the effect of the settings.
  R0_old_object_old_settings         Optional regression check: the October 2025 object with the old notebook's
                                     settings should reproduce the old notebook's table (our_human_invivo_DE.csv).

The two groups of settings
  release       Shared by the X-Mens, human and Hammoud fits (from ulrich_ivd_pipeline.py). A pseudobulk needs at least
                20 nuclei (25 until 29 October 2026; no human pseudobulk has fewer than 1,347). A gene is
                kept if its CPM is at least 1 in at least half (rounded up) of either arm's pseudobulks: 1 of the 2
                donors, or 3 of the 6 sections. (The release's docstring calls this a majority, which it is for the
                release's 5 pseudobulks per arm.) Cook's outlier refit on; no independent filtering.
  old_notebook  From our_human_invivo_DE.ipynb. A pseudobulk needs at least 10 nuclei. A gene is kept if it has
                counts in at least 3 pseudobulks. Cook's outlier refit on, Cook's filtering off, PyDESeq2's default
                independent filtering.

With two donors the donor design leaves 1 residual degree of freedom. Fold changes and per-donor directions are
informative; the p-values are not evidence that a difference replicates across women.

Usage
  python fit_human_de.py --pb-dir PB_DIR --out OUT_DIR [--old-pb-dir OLD_DIR --old-table our_human_invivo_DE.csv]

  PB_DIR   human_pb_meta.csv, human_pb_counts.csv.gz, human_pb_detect.csv.gz, human_var.csv (step 1)
  OLD_DIR  old_pb_meta.csv, old_pb_counts.csv.gz (step 2; only needed for the regression check)

Writes DE_<fit>.csv.gz for each fit, detection_pct_by_donor_and_arm.csv.gz and fit_diagnostics.json.
"""
import argparse
import json
import os
import warnings
from typing import NamedTuple

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# A gene is called differentially expressed when both thresholds are met.
MAX_PADJ = 0.05
MIN_ABS_LOG2FC = 1

# Minimum number of nuclei in a pseudobulk, for each group of settings. No human pseudobulk comes close to either.
# 20 is the cell floor of the Hammoud (Ulrich) release; the X-Mens code uses the same value (both were 25 until 29 October 2026).
MIN_NUCLEI = {"release": 20, "old_notebook": 10}


class Fit(NamedTuple):
    name: str          # the output file is DE_<name>.csv.gz
    label_column: str  # the obs column that assigns each nucleus to decidual or fibroblast
    unit: str          # what one pseudobulk is: "donor" (one per donor and arm) or "sample" (one per section and arm)
    block: str         # the term the design adjusts for, alongside group
    settings: str      # "release" or "old_notebook"


FITS = [
    # Fit A uses ~ timepoint + group, as the old notebook did. Each phase is one donor, so it is the same model as
    # ~ donor + group.
    Fit("A_new_celltype_old_settings", "cell_type", "donor", "timepoint", "old_notebook"),
    Fit("B_new_celltype_release", "cell_type", "donor", "donor", "release"),
    Fit("C_new_coarse_release", "coarse_refined", "donor", "donor", "release"),
    Fit("D_new_celltype_release_by_section", "cell_type", "sample", "sample", "release"),
    Fit("E_new_coarse_release_by_section", "coarse_refined", "sample", "sample", "release"),
]


# ---------------------------------------------------------------------------------------------------------------
# Inputs

def load_pseudobulks(pb_dir):
    """Step 1's outputs: one row per pseudobulk, its summed counts, its detection counts, and the gene symbols."""
    pseudobulks = pd.read_csv(f"{pb_dir}/human_pb_meta.csv").set_index("id")
    # Both matrices are pseudobulks x genes: summed UMIs, and the number of nuclei with at least one UMI.
    counts = pd.read_csv(f"{pb_dir}/human_pb_counts.csv.gz", index_col=0)
    detected = pd.read_csv(f"{pb_dir}/human_pb_detect.csv.gz", index_col=0)
    symbols = pd.read_csv(f"{pb_dir}/human_var.csv").set_index("gene_id").gene_symbol
    assert list(counts.index) == list(pseudobulks.index) and list(counts.columns) == list(symbols.index)
    assert list(detected.index) == list(counts.index) and list(detected.columns) == list(counts.columns)
    return pseudobulks, counts, detected, symbols


def select_pseudobulks(pseudobulks, label_column, unit):
    """The pseudobulks built from one label column and one unit, with the "group" and "sample" columns the
    designs use."""
    selected = pseudobulks[(pseudobulks.label_source == label_column) & (pseudobulks.unit_type == unit)]
    return selected.assign(group=selected.arm, sample=selected.unit)


# ---------------------------------------------------------------------------------------------------------------
# Gene filters

def genes_passing_release_filter(counts, is_decidual, min_cpm=1.0):
    """Genes with CPM >= min_cpm in at least half (rounded up) of the decidual pseudobulks, or of the fibroblast
    pseudobulks.

    Same rule as arm_aware_eligible() in the release's ulrich_ivd_pipeline.py. Each arm is checked on its own, so a
    gene that is switched on in only one arm is kept.
    """
    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    needed_in_decidual = (int(is_decidual.sum()) + 1) // 2
    needed_in_fibroblast = (int((~is_decidual).sum()) + 1) // 2
    passes_in_decidual = cpm.loc[is_decidual].ge(min_cpm).sum(axis=0).ge(needed_in_decidual)
    passes_in_fibroblast = cpm.loc[~is_decidual].ge(min_cpm).sum(axis=0).ge(needed_in_fibroblast)
    return passes_in_decidual | passes_in_fibroblast


def genes_passing_old_notebook_filter(counts):
    """Genes with at least one count in at least 3 pseudobulks."""
    return (counts > 0).sum(axis=0) >= 3


# ---------------------------------------------------------------------------------------------------------------
# The fit

def fit_decidual_vs_fibroblast(counts, pseudobulks, block, settings, inference):
    """Fit ~ block + group with PyDESeq2 and test decidual vs fibroblast.

    counts       pseudobulks x genes, integer counts
    pseudobulks  one row per pseudobulk, in the same order as counts, with the columns `block` and "group"
    Returns the per-gene results and a dictionary of diagnostics.
    """
    design = f"~ {block} + group"
    metadata = pd.DataFrame({block: pseudobulks[block].astype(str).to_numpy(),
                             "group": pseudobulks["group"].astype(str).to_numpy()},
                            index=counts.index.astype(str))
    is_decidual = (metadata.group == "decidual").to_numpy()

    if settings == "release":
        kept_genes = genes_passing_release_filter(counts, is_decidual)
        test_options = dict(independent_filter=False)
    elif settings == "old_notebook":
        kept_genes = genes_passing_old_notebook_filter(counts)
        # The old notebook set fibroblast as the reference level. The contrast below fixes the direction anyway.
        metadata["group"] = pd.Categorical(metadata.group, categories=["fibroblast", "decidual"])
        test_options = dict(cooks_filter=False)
    else:
        raise ValueError(f"unknown settings: {settings}")

    kept_counts = counts.loc[:, kept_genes.to_numpy()].astype(int).set_axis(metadata.index)
    dds = DeseqDataSet(counts=kept_counts, metadata=metadata, design=design,
                       refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()

    # The design must have full rank, i.e. group must not be confounded with the block term.
    design_matrix = dds.obsm["design_matrix"]
    rank = int(np.linalg.matrix_rank(design_matrix.to_numpy(dtype=float)))
    assert rank == design_matrix.shape[1], (rank, list(design_matrix.columns))

    stats = DeseqStats(dds, contrast=["group", "decidual", "fibroblast"], inference=inference, quiet=True,
                       **test_options)
    stats.summary()

    results = stats.results_df[["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]].copy()
    results["estimable"] = np.isfinite(results.log2FoldChange) & np.isfinite(results.pvalue)
    # A missing padj (for example, a gene removed by independent filtering) counts as not significant.
    padj_or_1 = results.padj.where(np.isfinite(results.padj), 1.0)
    results["call"] = (results.estimable
                       & padj_or_1.lt(MAX_PADJ)
                       & results.log2FoldChange.abs().gt(MIN_ABS_LOG2FC))
    results["direction"] = np.where(results.log2FoldChange > 0, "up",
                                    np.where(results.log2FoldChange < 0, "down", ""))

    called = results.call
    diagnostics = dict(
        settings=settings,
        design=design,
        n_pseudobulks=len(metadata),
        residual_df=int(len(metadata) - rank),
        eligible=int(kept_genes.sum()),
        calls=int(called.sum()),
        calls_up=int((called & (results.log2FoldChange > 0)).sum()),
        calls_down=int((called & (results.log2FoldChange < 0)).sum()),
        padj_na=int(results.padj.isna().sum()),
    )
    return results, diagnostics


# ---------------------------------------------------------------------------------------------------------------
# Regression check and detection

def log10_padj(padj):
    """log10 of adjusted p-values; values of 0 are raised to 1e-300 first."""
    return np.log10(padj.clip(lower=1e-300))


def regression_check(old_pb_dir, old_table_path, inference):
    """Refit the October 2025 object's pseudobulks with the old notebook's settings and compare the result with
    the old notebook's table."""
    pseudobulks = pd.read_csv(f"{old_pb_dir}/old_pb_meta.csv").set_index("id")
    counts = pd.read_csv(f"{old_pb_dir}/old_pb_counts.csv.gz", index_col=0)
    pseudobulks = pseudobulks[pseudobulks.n_cells >= MIN_NUCLEI["old_notebook"]]
    results, diagnostics = fit_decidual_vs_fibroblast(counts.loc[pseudobulks.index], pseudobulks,
                                                      "timepoint", "old_notebook", inference)
    results.index = results.index.astype(str)

    old = pd.read_csv(old_table_path, index_col=0)
    shared = results.index.intersection(old.index)
    old_calls = set(old.index[(old.padj < MAX_PADJ) & (old.log2FoldChange.abs() > MIN_ABS_LOG2FC)])
    lfc_difference = results.log2FoldChange[shared] - old.log2FoldChange[shared]
    log10_padj_difference = log10_padj(results.padj[shared]) - log10_padj(old.padj[shared])
    diagnostics.update(
        regression_genes_old=len(old),
        regression_genes_refit=len(results),
        regression_shared=len(shared),
        regression_max_abs_dlfc=float(lfc_difference.abs().max()),
        regression_max_abs_dlog10padj=float(log10_padj_difference.abs().max()),
        regression_calls_old=len(old_calls),
        regression_calls_refit=int(results.call.sum()),
        regression_calls_identical=bool(set(results.index[results.call]) == old_calls),
    )
    return results, diagnostics


def detection_percentages(pseudobulks, detected, symbols):
    """Percent of nuclei with at least one UMI of each gene, in every donor x arm pseudobulk, for both label
    columns."""
    columns = {}
    for label_column in ("cell_type", "coarse_refined"):
        for pseudobulk_id, pseudobulk in select_pseudobulks(pseudobulks, label_column, "donor").iterrows():
            name = f"{label_column}|{pseudobulk.donor}|{pseudobulk.timepoint}|{pseudobulk.arm}"
            columns[name] = 100 * detected.loc[pseudobulk_id] / pseudobulk.n_cells
    table = pd.DataFrame(columns)
    table.insert(0, "gene_symbol", symbols.reindex(table.index).to_numpy())
    return table


# ---------------------------------------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pb-dir", required=True, help="folder with step 1's outputs")
    parser.add_argument("--out", required=True, help="output folder")
    parser.add_argument("--old-pb-dir",
                        help="folder with step 2's outputs; with --old-table, adds the regression check")
    parser.add_argument("--old-table", help="the old notebook's table, our_human_invivo_DE.csv")
    parser.add_argument("--n-cpus", type=int, default=1, help="1 and 8 CPUs give identical results (checked 29 October 2026)")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # PyDESeq2 warns on every fit that the residual degrees of freedom are few; the docstring above covers this.
    warnings.filterwarnings("ignore", message="As the residual degrees of freedom")
    inference = DefaultInference(n_cpus=args.n_cpus)

    pseudobulks, counts, detected, symbols = load_pseudobulks(args.pb_dir)

    results, diagnostics = {}, {}
    if args.old_pb_dir and args.old_table:
        name = "R0_old_object_old_settings"
        results[name], diagnostics[name] = regression_check(args.old_pb_dir, args.old_table, inference)

    for fit in FITS:
        selected = select_pseudobulks(pseudobulks, fit.label_column, fit.unit)
        selected = selected[selected.n_cells >= MIN_NUCLEI[fit.settings]]
        fit_results, fit_diagnostics = fit_decidual_vs_fibroblast(counts.loc[selected.index], selected,
                                                                  fit.block, fit.settings, inference)
        fit_results.insert(0, "gene_symbol", symbols.reindex(fit_results.index).to_numpy())
        fit_diagnostics.update(label_source=fit.label_column, unit=fit.unit)
        results[fit.name], diagnostics[fit.name] = fit_results, fit_diagnostics

    for name, table in results.items():
        table.to_csv(f"{args.out}/DE_{name}.csv.gz")
    detection_percentages(pseudobulks, detected, symbols).to_csv(f"{args.out}/detection_pct_by_donor_and_arm.csv.gz")
    with open(f"{args.out}/fit_diagnostics.json", "w") as handle:
        json.dump(diagnostics, handle, indent=1, default=str)

    for name, fit_diagnostics in diagnostics.items():
        shown = ("design", "n_pseudobulks", "residual_df", "eligible", "calls")
        print(name, {key: fit_diagnostics[key] for key in shown})


if __name__ == "__main__":
    main()
