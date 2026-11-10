#!/usr/bin/env python3
"""Step 2: X-Mens decidual vs fibroblast differential expression (PyDESeq2), and the all-genes results table.

Input: the pseudobulks that cluster_jobs/mouse_pseudobulk_job.py (step 1) built. With its default design these are
11: one decidual and one fibroblast pseudobulk for each of the five induced timepoints, and one fibroblast pseudobulk
from the uninduced sample. (Its induced_only design gives the 10 pseudobulks of the release of 19 October 2026.)

The fit
  design     ~ timepoint + group. Timepoint is a fixed block, so decidual is compared with fibroblast within
             timepoints. Each induced timepoint has one pseudobulk of each arm; the uninduced timepoint has a
             fibroblast pseudobulk only. Its timepoint coefficient fits it exactly, so it adds nothing to the
             decidual vs fibroblast estimate directly; it enters the size factors, the dispersion estimates and the
             expression filter. 11 pseudobulks and 7 coefficients leave 4 residual degrees of freedom (10 and 6
             without the uninduced pseudobulk).
  test       Wald test, decidual vs fibroblast. A positive log2 fold change means higher in decidual nuclei.
  genes      a gene is fitted if its CPM is at least 1 in at least half (rounded up) of either arm's pseudobulks:
             3 of the 5 decidual or 3 of the 6 fibroblast pseudobulks (3 of 5 without the uninduced one). Each arm
             is checked on its own, so genes switched on in only one arm are kept.
  settings   Cook's outlier refit on and Cook's filtering on (PyDESeq2's defaults); no independent filtering, so
             adjusted p-values are Benjamini-Hochberg over exactly the fitted genes.
  shrinkage  log2FoldChange_shrunk is PyDESeq2's apeglm shrinkage of the same effect. It is for display only;
             calls use the unshrunk fold change.
  calls      adjusted p < 0.05 and |unshrunk log2 fold change| > 1.

The settings are those of the release of 19 October 2026 (release_20261019_v6_dS_vs_C7fib_luteal, made by
ulrich_ivd_pipeline.py). On step 1's induced_only pseudobulks, this script reproduces that release's
xmens_broad_all_genes.csv byte for byte; see README_mouse_DE_code.md.

Usage
  python fit_mouse_de.py --pb-dir PB_DIR --out OUT_DIR [--n-cpus N]

  PB_DIR  xmens_broad_carrier_counts.tsv, xmens_broad_metadata.tsv and xmens_broad_genes.tsv, from step 1

Writes xmens_broad_all_genes.csv (one row per gene in the object) and fit_diagnostics.json.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# A gene is called differentially expressed when both thresholds are met.
MAX_PADJ = 0.05
MIN_ABS_LOG2FC = 1.0

# Expression filter: CPM of at least MIN_CPM in at least half (rounded up) of either arm's pseudobulks.
MIN_CPM = 1.0

DESIGN = "~ timepoint + group"

# PyDESeq2 orders factor levels alphabetically, so "decidual" is the reference level and the design's group
# coefficient is fibroblast vs decidual. The Wald test below is set up as decidual vs fibroblast directly; only the
# shrinkage, which works on a coefficient, needs this coefficient's sign flipped.
GROUP_COEFFICIENT = "group[T.fibroblast]"

# Fixed columns of the results table, in the release's format. The source ID names the design.
SOURCE_ID = {"induced_only": "mouse_all_decidual_vs_fibroblast_timepoint_blocked",
             "with_uninduced": "mouse_all_decidual_vs_fibroblast_timepoint_blocked_with_uninduced_fibroblasts"}
UNINDUCED = "uninduced"
POSITIVE_LFC_MEANING = "decidual_minus_fibroblast"
ANALYSIS_METHOD = "pydeseq2_wald_fixed_block"

STATISTICS = ["log2FoldChange_unshrunk", "log2FoldChange_shrunk", "lfcSE", "stat", "pvalue", "padj"]


def load_pseudobulks(pb_dir):
    """The pseudobulk counts (pseudobulks x genes), their descriptions, and the genes, all in matching order."""
    counts = pd.read_csv(f"{pb_dir}/xmens_broad_carrier_counts.tsv", sep="\t", index_col="pseudobulk_id")
    pseudobulks = pd.read_csv(f"{pb_dir}/xmens_broad_metadata.tsv", sep="\t")
    # Empty fields are missing values; everything else is read literally, so no gene symbol is mistaken for "NA".
    genes = pd.read_csv(f"{pb_dir}/xmens_broad_genes.tsv", sep="\t", keep_default_na=False, na_values=[""])

    assert list(counts.index) == list(pseudobulks.pseudobulk_id), "counts and descriptions list different pseudobulks"
    assert list(counts.columns) == list(genes.ensembl_gene_id), "counts and gene table list different genes"
    assert set(pseudobulks.group) == {"decidual", "fibroblast"}
    # Every induced timepoint has one pseudobulk of each arm; the uninduced timepoint, if present, a fibroblast one.
    arms = pseudobulks.groupby("timepoint").group.apply(lambda g: tuple(sorted(g)))
    for timepoint, timepoint_arms in arms.items():
        expected = ("fibroblast",) if timepoint == UNINDUCED else ("decidual", "fibroblast")
        assert timepoint_arms == expected, f"{timepoint}: pseudobulks {timepoint_arms}, expected {expected}"
    return counts, pseudobulks, genes


def design_name(pseudobulks):
    return "with_uninduced" if (pseudobulks.timepoint == UNINDUCED).any() else "induced_only"


def genes_passing_filter(counts, is_decidual):
    """True for genes with CPM >= MIN_CPM in at least half (rounded up) of the decidual pseudobulks or of the
    fibroblast pseudobulks. CPM is computed over all genes in the object, not only the ones that pass."""
    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    needed_decidual = (int(is_decidual.sum()) + 1) // 2
    needed_fibroblast = (int((~is_decidual).sum()) + 1) // 2
    passes_in_decidual = cpm.loc[is_decidual].ge(MIN_CPM).sum(axis=0).ge(needed_decidual)
    passes_in_fibroblast = cpm.loc[~is_decidual].ge(MIN_CPM).sum(axis=0).ge(needed_fibroblast)
    return passes_in_decidual | passes_in_fibroblast


def fit_decidual_vs_fibroblast(counts, pseudobulks, n_cpus):
    """PyDESeq2 fit of ~ timepoint + group and Wald test of decidual vs fibroblast.

    Returns one row per fitted gene (unshrunk and shrunk log2 fold changes, standard error, Wald statistic, p-value,
    adjusted p-value, whether the gene was estimable) and a dictionary of fit diagnostics."""
    design_data = pd.DataFrame({"timepoint": pseudobulks.timepoint.astype(str).to_numpy(),
                                "group": pseudobulks.group.astype(str).to_numpy()},
                               index=pseudobulks.pseudobulk_id.astype(str).to_numpy())
    counts = counts.copy()
    counts.index = design_data.index

    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(counts=counts, metadata=design_data, design=DESIGN, refit_cooks=True,
                       inference=inference, quiet=True)
    dds.deseq2()

    design_matrix = dds.obsm["design_matrix"]
    rank = int(np.linalg.matrix_rank(design_matrix.to_numpy(dtype=float)))
    assert rank == design_matrix.shape[1], f"the design {DESIGN} is rank-deficient"
    group_coefficients = [column for column in design_matrix.columns if column.startswith("group[T.")]
    assert group_coefficients == [GROUP_COEFFICIENT], group_coefficients

    # Wald test, decidual vs fibroblast; the expression filter above is the only gene filter.
    stats = DeseqStats(dds, contrast=["group", "decidual", "fibroblast"], independent_filter=False,
                       inference=inference, quiet=True)
    stats.summary()
    results = stats.results_df.copy()

    # apeglm shrinkage of the same effect. lfc_shrink overwrites stats.results_df, which is why the unshrunk
    # results were copied above. The coefficient is fibroblast vs decidual, so flip its sign.
    stats.lfc_shrink(coeff=GROUP_COEFFICIENT)
    shrunk = stats.results_df["log2FoldChange"] * -1.0

    # The shrunk and unshrunk fold changes must agree in sign wherever both are clearly away from zero;
    # if they did not, the coefficient would be oriented the wrong way round.
    clear = results.log2FoldChange.abs().gt(1) & shrunk.abs().gt(0.5)
    disagree = np.sign(results.loc[clear, "log2FoldChange"]) != np.sign(shrunk[clear])
    assert not (clear.any() and disagree.mean() > 0.01), "shrunk and unshrunk fold changes disagree in sign"

    fitted = pd.DataFrame({
        "ensembl_gene_id": results.index.astype(str),
        "baseMean": results.baseMean.to_numpy(),
        "log2FoldChange_unshrunk": results.log2FoldChange.to_numpy(),
        "log2FoldChange_shrunk": shrunk.reindex(results.index).to_numpy(),
        "lfcSE": results.lfcSE.to_numpy(),
        "stat": results.stat.to_numpy(),
        "pvalue": results.pvalue.to_numpy(),
        "padj": results.padj.to_numpy(),
    })
    fitted["estimable"] = np.isfinite(fitted.log2FoldChange_unshrunk) & np.isfinite(fitted.pvalue)

    size_factors = dds.obs["size_factors"]
    diagnostics = {
        "design": DESIGN,
        "design_columns": list(design_matrix.columns),
        "design_rank": rank,
        "residual_df": int(design_matrix.shape[0] - rank),
        "coefficient": GROUP_COEFFICIENT,
        "n_pseudobulks": int(len(design_data)),
        "n_genes_fitted": int(len(fitted)),
        "size_factor_min": float(size_factors.min()),
        "size_factor_median": float(size_factors.median()),
        "size_factor_max": float(size_factors.max()),
        "n_padj_na": int(fitted.padj.isna().sum()),
        "n_not_estimable": int((~fitted.estimable).sum()),
    }
    return fitted, diagnostics


def build_results_table(genes, counts, passes_filter, fitted, pseudobulks):
    """One row per gene in the object, in the release's column format. Genes that were not fitted keep their
    annotation and counts, with empty statistics and production_eligible = False."""
    table = genes.reset_index(drop=True).copy()
    table["carrier"] = True                                   # every gene in the object is carried into the table
    table["production_eligible"] = passes_filter.to_numpy()   # passed the expression filter, so was fitted
    table["total_count"] = counts.sum(axis=0).to_numpy(dtype=np.int64)
    table["nonzero_pseudobulks"] = counts.gt(0).sum(axis=0).to_numpy(dtype=np.int64)

    fitted_ids = table.ensembl_gene_id[table.production_eligible].astype(str).tolist()
    assert fitted.ensembl_gene_id.tolist() == fitted_ids, "fitted genes are not the genes that passed the filter"
    table = table.merge(fitted, on="ensembl_gene_id", how="left", validate="one_to_one")

    # Genes that were not fitted have no estimability flag; they count as not estimable.
    table["estimable"] = table.production_eligible & table.estimable.astype("boolean").fillna(False).astype(bool)
    table.loc[~table.estimable, STATISTICS] = np.nan
    # A missing adjusted p-value counts as 1 (not significant) when calling genes.
    table["padj_for_call"] = table.padj.where(np.isfinite(table.padj), 1.0)
    table["significant_fdr05_lfc1"] = (table.estimable
                                       & table.padj_for_call.lt(MAX_PADJ)
                                       & table.log2FoldChange_unshrunk.abs().gt(MIN_ABS_LOG2FC))

    table["source_id"] = SOURCE_ID[design_name(pseudobulks)]
    table["positive_lfc_meaning"] = POSITIVE_LFC_MEANING
    table["analysis_method"] = ANALYSIS_METHOD
    is_decidual = pseudobulks.group.astype(str) == "decidual"
    table["n_pseudobulks_numerator"] = int(is_decidual.sum())
    table["n_pseudobulks_denominator"] = int((~is_decidual).sum())
    return table


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pb-dir", required=True, help="folder with step 1's outputs")
    parser.add_argument("--out", required=True, help="output folder")
    parser.add_argument("--n-cpus", type=int, default=8, help="the release used 8")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    counts, pseudobulks, genes = load_pseudobulks(args.pb_dir)
    is_decidual = (pseudobulks.group == "decidual").to_numpy()
    passes_filter = genes_passing_filter(counts, is_decidual)

    fitted, diagnostics = fit_decidual_vs_fibroblast(counts.loc[:, passes_filter.to_numpy()], pseudobulks,
                                                     args.n_cpus)
    table = build_results_table(genes, counts, passes_filter, fitted, pseudobulks)

    diagnostics["pseudobulk_design"] = design_name(pseudobulks)
    diagnostics["pseudobulks_per_arm"] = {arm: int((pseudobulks.group == arm).sum()) for arm in ("decidual", "fibroblast")}
    diagnostics["nuclei_per_arm"] = {arm: int(pseudobulks.n_cells[pseudobulks.group == arm].sum())
                                     for arm in ("decidual", "fibroblast")}
    diagnostics["eligible_genes"] = int(passes_filter.sum())
    diagnostics["strict_calls"] = int(table.significant_fdr05_lfc1.sum())
    diagnostics["calls_up"] = int((table.significant_fdr05_lfc1 & (table.log2FoldChange_unshrunk > 0)).sum())
    diagnostics["calls_down"] = int((table.significant_fdr05_lfc1 & (table.log2FoldChange_unshrunk < 0)).sum())

    # Same number format as the release: 15 significant digits, empty fields for missing values.
    table.to_csv(f"{args.out}/xmens_broad_all_genes.csv", index=False, na_rep="", float_format="%.15g",
                 lineterminator="\n")
    with open(f"{args.out}/fit_diagnostics.json", "w") as handle:
        json.dump(diagnostics, handle, indent=1)
    shown = ("pseudobulk_design", "design", "n_pseudobulks", "residual_df", "eligible_genes", "strict_calls",
             "calls_up", "calls_down")
    print({key: diagnostics[key] for key in shown})


if __name__ == "__main__":
    main()
