#!/usr/bin/env python3
"""Step 5: compare the old and new human decidual vs fibroblast DE tables, and redo every downstream comparison with
each human table against the same reference lists, so that any difference comes from the human rerun alone.

Reference lists
  - X-Mens decidual vs fibroblast, genes called there: the primary fit of 27 October 2026, which includes the
    uninduced fibroblasts (xmens_broad_all_genes_with_uninduced.csv), or the release of 19 October 2026
    (xmens_broad_all_genes.csv); the summary row names whichever table is given
  - Wang et al. 2020 mouse lists: oil-induced decidualization, natural pregnancy (day 8), in vitro decidualization
  - the corrected human in vitro consensus (123 genes; Stadtmauer genes called by the authors' criteria)

Mouse genes are matched to human genes through ortholog_pairs_one2one.csv (xmens_final/results/mi_dpt/overlap/):
Ensembl 116 Compara one-to-one, protein-coding pairs, keyed by symbol, built from
xmens_final/resources/orthology/ensembl_human_mouse_components.tsv. Pairs whose mouse symbol belongs to two genes
(Tsc22d4) are left out.

A human gene is called at adjusted p < 0.05 and |log2FC| > 1. "Same direction" counts the genes called on both
sides. "Chance" is the agreement expected if the directions were unrelated:
    p_up(human) * p_up(reference) + p_down(human) * p_down(reference), over the genes compared.

Usage
  python compare_downstream.py --fits FIT_DIR --var human_var.csv --old-table our_human_invivo_DE.csv \\
      --xmens xmens_broad_all_genes_with_uninduced.csv --orthology ortholog_pairs_one2one.csv \\
      --consensus consensus123_stadtmauer_authors_criteria.csv --wang-dir WANG_DIR --out OUT_DIR

  FIT_DIR   DE_<fit>.csv.gz from step 3
  WANG_DIR  wang_oil.csv, wang_pregnancy.csv, wang_ivd.csv (columns gene, direction)

Writes SUMMARY_human_tables_downstream.csv, SUMMARY_human_table_agreement.csv, TABLE_marker_genes_by_version.csv
and TABLE_consensus123_vs_human_by_version.csv.
"""
import argparse
import os

import numpy as np
import pandas as pd

MAX_PADJ = 0.05
MIN_ABS_LOG2FC = 1

# The old notebook's table first, then the fits from step 3
HUMAN_TABLES = [
    "old notebook table",
    "A_new_celltype_old_settings",
    "B_new_celltype_release",
    "C_new_coarse_release",
    "D_new_celltype_release_by_section",
    "E_new_coarse_release_by_section",
]

# Pairs of human tables compared with each other
TABLE_PAIRS = [
    ("old notebook table", "A_new_celltype_old_settings"),        # effect of the new object
    ("A_new_celltype_old_settings", "B_new_celltype_release"),    # effect of the settings
    ("old notebook table", "B_new_celltype_release"),             # both together
    ("B_new_celltype_release", "C_new_coarse_release"),           # effect of the label column
    ("B_new_celltype_release", "D_new_celltype_release_by_section"),
    ("C_new_coarse_release", "E_new_coarse_release_by_section"),
]

WANG_LISTS = [
    ("Wang oil", "wang_oil.csv"),
    ("Wang pregnancy", "wang_pregnancy.csv"),
    ("Wang in vitro", "wang_ivd.csv"),
]

# Consensus genes whose published symbol differs from the object's symbol (renamed since publication)
RENAMED_CONSENSUS_GENES = {
    "C10ORF10": ["DEPP1", "C10orf10"],
    "FAM19A2": ["TAFA2"],
    "MUM1L1": ["PWWP3B"],
    "PPAP2B": ["PLPP3"],
    "ZCCHC6": ["TUT7"],
}

MARKERS = ["PRL", "IGFBP1", "LEFTY2", "FOXO1", "HAND2", "ZBTB16", "WNT4", "IL15", "PRLR", "BMP2", "CEBPB", "DKK1",
           "SST", "GADD45A"]


# ---------------------------------------------------------------------------------------------------------------
# Inputs

def load_human_tables(fits_dir, symbols, old):
    """Every human table as a frame indexed by gene symbol, with columns lfc, padj, called, tested."""
    old_called = (old.padj < MAX_PADJ) & (old.log2FoldChange.abs() > MIN_ABS_LOG2FC)
    tables = {"old notebook table": pd.DataFrame({"lfc": old.log2FoldChange,
                                                  "padj": old.padj,
                                                  "called": old_called.fillna(False),
                                                  "tested": old.log2FoldChange.notna()},
                                                 index=old.index.astype(str))}
    for name in HUMAN_TABLES[1:]:
        fitted = pd.read_csv(f"{fits_dir}/DE_{name}.csv.gz", index_col=0)
        tables[name] = pd.DataFrame({"lfc": fitted.log2FoldChange.to_numpy(),
                                     "padj": fitted.padj.to_numpy(),
                                     "called": fitted.call.to_numpy(),
                                     "tested": fitted.estimable.to_numpy()},
                                    index=symbols.reindex(fitted.index).astype(str).to_numpy())
    return tables


def without_missing_or_repeated_symbols(directions):
    """Drop genes without a human ortholog, and human symbols that occur more than once."""
    directions = directions[directions.index.notna()]
    return directions[~directions.index.duplicated(keep=False)]


def reference_lists(xmens, orthology, wang_dir):
    """Each reference list as +1 (up) or -1 (down), indexed by the human ortholog's symbol."""
    mouse_to_human = orthology.set_index("mouse").human

    # X-Mens: genes called in the given fit. Symbols that occur more than once are dropped.
    xmens_called = xmens[xmens.significant_fdr05_lfc1.astype(bool)].drop_duplicates("gene_symbol", keep=False)
    xmens_directions = pd.Series(np.sign(xmens_called.log2FoldChange_unshrunk.to_numpy()),
                                 index=xmens_called.gene_symbol.map(mouse_to_human).to_numpy())
    with_uninduced = str(xmens.source_id.iloc[0]).endswith("with_uninduced_fibroblasts")
    label = "X-Mens decidual vs fibroblast" + (" (with uninduced fibroblasts)" if with_uninduced else " (release)")
    lists = {label: without_missing_or_repeated_symbols(xmens_directions)}

    # Wang et al.: genes as published, with their direction. Genes listed more than once are dropped.
    for name, filename in WANG_LISTS:
        wang = pd.read_csv(f"{wang_dir}/{filename}").drop_duplicates("gene", keep=False)
        wang_directions = pd.Series(wang.direction.str.lower().map({"up": 1, "down": -1}).to_numpy(),
                                    index=wang.gene.map(mouse_to_human).to_numpy())
        lists[name] = without_missing_or_repeated_symbols(wang_directions)
    return lists


# ---------------------------------------------------------------------------------------------------------------
# Comparisons

def direction_agreement(human, reference_direction):
    """How many of the human calls that are also in the reference list point the same way, and the share
    expected by chance."""
    compared = reference_direction.index.intersection(human.index[human.called])
    human_sign = np.sign(human.lfc[compared])
    reference_sign = reference_direction[compared]
    n_same = int((human_sign == reference_sign).sum())
    n_compared = len(compared)
    if not n_compared:
        return dict(same=0, compared=0, pct_same=np.nan, pct_chance=np.nan)
    reference_up = (reference_sign > 0).mean()
    human_up = (human_sign > 0).mean()
    chance = reference_up * human_up + (1 - reference_up) * (1 - human_up)
    return dict(same=n_same, compared=n_compared,
                pct_same=round(100 * n_same / n_compared, 1), pct_chance=round(100 * chance, 1))


def find_symbol(gene, symbols):
    """The consensus gene's symbol in the human table, trying the current names of renamed genes."""
    if gene in symbols:
        return gene
    return next((alias for alias in RENAMED_CONSENSUS_GENES.get(gene, []) if alias in symbols), None)


def consensus_status(human, consensus):
    """For each consensus gene: same direction, opposite direction, tested but not called, or not testable."""
    rows = []
    for gene, consensus_direction in zip(consensus.gene, consensus.consensus_direction):
        symbol = find_symbol(gene, human.index)
        if symbol is None or not bool(human.tested.get(symbol, False)):
            status = "not testable"
        elif not bool(human.called[symbol]):
            status = "tested, not called"
        else:
            expected_sign = 1 if consensus_direction == "up" else -1
            status = "same direction" if np.sign(human.lfc[symbol]) == expected_sign else "opposite direction"
        lfc = human.lfc.get(symbol, np.nan) if symbol else np.nan
        rows.append(dict(gene=gene, matched_symbol=symbol, consensus_direction=consensus_direction,
                         status=status, lfc=lfc))
    return pd.DataFrame(rows)


def compare_two_tables(tables, first, second):
    """Agreement between two human tables, over the genes tested in both."""
    a, b = tables[first], tables[second]
    tested_in_both = a.index[a.tested].intersection(b.index[b.tested])
    calls_a = set(a.index[a.called]) & set(tested_in_both)
    calls_b = set(b.index[b.called]) & set(tested_in_both)
    called_in_both = sorted(calls_a & calls_b)
    same_direction = int((np.sign(a.lfc[called_in_both]) == np.sign(b.lfc[called_in_both])).sum())
    return {
        "shared tested genes": len(tested_in_both),
        "log2FC Pearson r": round(float(np.corrcoef(a.lfc[tested_in_both], b.lfc[tested_in_both])[0, 1]), 4),
        "calls (first / second)": f"{len(calls_a):,} / {len(calls_b):,}",
        "called in both": len(called_in_both),
        "Jaccard": round(len(called_in_both) / max(1, len(calls_a | calls_b)), 3),
        "called in both, same direction": same_direction,
    }


def summary_row(name, human, references, consensus_table):
    """One column of SUMMARY_human_tables_downstream.csv."""
    n = consensus_table.status.value_counts()
    row = {
        "human table": name,
        "genes tested": int(human.tested.sum()),
        "calls": int(human.called.sum()),
        "calls up": int((human.called & (human.lfc > 0)).sum()),
        "calls down": int((human.called & (human.lfc < 0)).sum()),
        "consensus 123 (corrected): same / opposite / tested not called / not testable": (
            f"{n.get('same direction', 0)} / {n.get('opposite direction', 0)} / "
            f"{n.get('tested, not called', 0)} / {n.get('not testable', 0)}"),
    }
    for reference_name, reference_direction in references.items():
        agreement = direction_agreement(human, reference_direction)
        row[reference_name] = (f"{agreement['same']:,} of {agreement['compared']:,} "
                               f"({agreement['pct_same']}%; chance {agreement['pct_chance']}%)")
    return row


# ---------------------------------------------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for option in ("--fits", "--var", "--old-table", "--xmens", "--orthology", "--consensus", "--wang-dir", "--out"):
        parser.add_argument(option, required=True)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    symbols = pd.read_csv(args.var).set_index("gene_id").gene_symbol
    assert symbols.is_unique
    human_tables = load_human_tables(args.fits, symbols, pd.read_csv(args.old_table, index_col=0))
    references = reference_lists(pd.read_csv(args.xmens), pd.read_csv(args.orthology), args.wang_dir)
    consensus = pd.read_csv(args.consensus)

    consensus_tables = {name: consensus_status(human, consensus) for name, human in human_tables.items()}
    summary = pd.DataFrame([summary_row(name, human, references, consensus_tables[name])
                            for name, human in human_tables.items()]).set_index("human table").T
    agreement = pd.DataFrame({f"{first} vs {second}": compare_two_tables(human_tables, first, second)
                              for first, second in TABLE_PAIRS})
    markers = pd.concat({name: human_tables[name].reindex(MARKERS)[["lfc", "padj", "called"]]
                         for name in HUMAN_TABLES}, axis=1)
    consensus_by_table = pd.concat({name: table.set_index("gene") for name, table in consensus_tables.items()},
                                   axis=1)

    summary.to_csv(f"{args.out}/SUMMARY_human_tables_downstream.csv")
    agreement.to_csv(f"{args.out}/SUMMARY_human_table_agreement.csv")
    markers.to_csv(f"{args.out}/TABLE_marker_genes_by_version.csv")
    consensus_by_table.to_csv(f"{args.out}/TABLE_consensus123_vs_human_by_version.csv")
    print(summary.loc[:, ["old notebook table", "B_new_celltype_release"]].to_string())


if __name__ == "__main__":
    main()
