#!/usr/bin/env python
"""Build IVD_comparisons_with_genes.xlsx, the gene-level workbook of the in vitro / in vivo comparisons, from tables alone.

Until 2 Nov 2026 this workbook was built in layers: the Wang sheets and the layout by session code (19 Oct), the consensus
sheets by build_outputs_stadtmauer_fix.py (22 Oct), the X-Mens sheets by update_ivd_workbook.py (27 Oct) and the two
human sheets by add_human_comparisons_to_ivd_workbook.py (29 Oct), each starting from the previous file. This script
builds the same 12 sheets, Summary formulas and Notes in one pass from the pipeline's tables, with the same functions
where they exist (build_house_style_workbooks.py for the layout; add_human_comparisons_to_ivd_workbook.py for the two
human comparisons and their notes; build_outputs_stadtmauer_fix.py for the study labels of the consensus sheet).

Sheets, in order: Summary, Notes, consensus 123 genes, X-Mens vs consensus in vitro, X-Mens vs oil, X-Mens vs pregnancy,
X-Mens vs in vitro, oil vs pregnancy, in vitro vs oil, in vitro vs pregnancy, Slide-tags vs consensus, Hammoud vs consensus.
Summary: genes compared = COUNTA over the gene column, same direction = COUNTIF(...,TRUE) over the "same direction"
column, % same direction = ROUND(100*C/B,1), % expected by chance = p_up(A)*p_up(B) + p_down(A)*p_down(B) over the genes
compared (a fixed value), and the sheet name.

Inputs
  --consensus-dir   output folder of build_consensus123.py (consensus list, orthology)
  --consensus-tables  folder with TABLE_consensus123_per_gene_calls.csv (make_consensus_tables.py)
  --xmens-compare   output folder of compare_xmens_with_references.py (genes_xmens_vs_*.csv, SUMMARY_xmens_comparisons.csv)
  --wang-dir        wang_oil.csv, wang_pregnancy.csv, wang_ivd.csv
  --slidetags       human_decidual_vs_fibroblast_all_genes.csv
  --hammoud-release, --hammoud-release-name, --hammoud-config, --hammoud-donors, --hammoud-contrast, --hammoud-label
                    as in add_human_comparisons_to_ivd_workbook.py
  --out             the workbook (must not exist)
"""
from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
READ = dict(keep_default_na=False, na_values=[""], float_precision="round_trip")

XMENS_SHEETS = [  # sheet, step-3 file key, Summary label as the comparison script writes it
    ("X-Mens vs consensus in vitro", "consensus"),
    ("X-Mens vs oil", "oil"),
    ("X-Mens vs pregnancy", "pregnancy"),
    ("X-Mens vs in vitro", "in_vitro"),
]
WANG_SHEETS = [  # sheet, Summary label, (list A, column A), (list B, column B)
    ("oil vs pregnancy", "Wang oil-induced vs Wang natural pregnancy", ("oil", "oil (Wang)"), ("pregnancy", "pregnancy (Wang)")),
    ("in vitro vs oil", "Wang in vitro vs Wang oil-induced", ("ivd", "in vitro (Wang)"), ("oil", "oil (Wang)")),
    ("in vitro vs pregnancy", "Wang in vitro vs Wang natural pregnancy", ("ivd", "in vitro (Wang)"), ("pregnancy", "pregnancy (Wang)")),
]
# Notes 4-11 of the workbook of record. Rows 10 and 11 record the two revisions the workbook went through before it was
# built by this script; they are kept as history.
FIXED_NOTES = [
    "X-Mens rows: genes in the Wang list that X-Mens calls (adjusted p < 0.05, fold change > 2x, decidual cells vs endometrial "
    "fibroblasts); % same direction = share changed the same way as in Wang.",
    "Wang vs Wang rows: genes in both Wang lists; % same direction = share listed the same way in both.",
    "% expected by chance: agreement if the two lists were unrelated, given how many genes go up vs down in each.",
    "Each comparison has its own sheet listing the genes compared.",
    "Wang lists: Wang et al. 2020, Genes 11:935 (published gene lists, direction only).",
    'Summary columns "genes compared", "same direction" and "% same direction" are formulas over the gene sheets (COUNTA / '
    'COUNTIF on the "same direction" column). "% expected by chance" is a fixed value.',
    "Consensus list revised 2026-10-22: Stadtmauer genes with tabulated values are counted only where a treatment meets the "
    "authors' criteria (FDR < 1e-6, fold change > 1.5, > 3 TPM vs base medium). This removed KLF4 and TSHZ1 and added ADAMTS5 "
    "and CLU; X-Mens changed from 19 of 49 to 20 of 50 genes in the same direction.",
    "X-Mens fit revised 2026-10-27: the uninduced sample's endometrial fibroblasts (9,270 nuclei) are now included as a sixth "
    "timepoint, giving 11 pseudobulks; decidual nuclei are still compared with fibroblasts within each timepoint. X-Mens rows "
    "changed from 1,429 of 1,578 to 1,439 of 1,588 (Wang oil), from 1,092 of 1,505 to 1,105 of 1,518 (Wang pregnancy) and "
    "from 459 of 654 to 461 of 657 (Wang in vitro); the consensus row stays at 20 of 50.",
]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def chance(a: pd.Series, b: pd.Series) -> float:
    """Agreement expected if the two direction lists were unrelated, given their up/down proportions."""
    p, q = (a == "up").mean(), (b == "up").mean()
    return round(100 * (p * q + (1 - p) * (1 - q)), 1)


def wang_pair(lists: dict, a: tuple, b: tuple) -> pd.DataFrame:
    """Genes in both Wang lists, with both directions and whether they agree; sorted by gene."""
    A = lists[a[0]].set_index("gene").direction
    B = lists[b[0]].set_index("gene").direction
    shared = sorted(set(A.index) & set(B.index))
    d = pd.DataFrame({"gene": shared, a[1]: A[shared].to_numpy(), b[1]: B[shared].to_numpy()})
    d["same direction"] = d[a[1]] == d[b[1]]
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("consensus-dir", "consensus-tables", "xmens-compare", "wang-dir", "slidetags", "hammoud-release",
                 "hammoud-release-name", "hammoud-config", "hammoud-donors", "out"):
        ap.add_argument(f"--{name}", required=True)
    ap.add_argument("--hammoud-contrast", default=None)
    ap.add_argument("--hammoud-label", default=None)
    a = ap.parse_args()
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; refusing to overwrite")
    hs = load_module("house_style", HERE / "build_house_style_workbooks.py")
    AH = load_module("add_human", HERE / "add_human_comparisons_to_ivd_workbook.py")
    BO = load_module("build_outputs", HERE.parent / "consensus" / "build_outputs_stadtmauer_fix.py")

    # ---- consensus list and orthology
    best = pd.read_csv(f"{a.consensus_dir}/consensus123_stadtmauer_authors_criteria.csv")
    orth = pd.read_csv(f"{a.consensus_dir}/consensus123_stadtmauer_authors_criteria_mouse_orthology.csv", dtype=str, keep_default_na=False)
    n123, n_orth = len(best), int((orth.strict_one_to_one == "True").sum())
    per_gene = pd.read_csv(f"{a.consensus_tables}/TABLE_consensus123_per_gene_calls.csv", **READ)
    cons_sheet = per_gene[["gene", "consensus direction", "studies supporting"]].rename(
        columns={"consensus direction": "direction in vitro", "studies supporting": "studies agreeing (of 3)"})
    mature = best.set_index("gene").mature_papers
    cons_sheet["studies"] = cons_sheet.gene.map(lambda g: BO.studies_label(mature.loc[g]))
    cons_sheet["mouse ortholog (1:1)"] = per_gene["mouse ortholog (1:1)"].fillna("").values
    cons_sheet = cons_sheet.sort_values(["studies agreeing (of 3)", "gene"], ascending=[False, True]).reset_index(drop=True)
    cons_sheet["mouse ortholog (1:1)"] = cons_sheet["mouse ortholog (1:1)"].replace("", np.nan)

    # ---- X-Mens comparison sheets (step 3 of the X-Mens DE code, read as written)
    xmens_summary = pd.read_csv(f"{a.xmens_compare}/SUMMARY_xmens_comparisons.csv").set_index("comparison")
    xmens_sheets, xmens_rows = [], []
    for sheet_name, key in XMENS_SHEETS:
        table = pd.read_csv(f"{a.xmens_compare}/genes_xmens_vs_{key}.csv", **READ)
        label = [c for c in xmens_summary.index if c.startswith("X-Mens vs human consensus")][0] if key == "consensus" else \
                {"oil": "X-Mens vs Wang oil-induced decidualization", "pregnancy": "X-Mens vs Wang natural pregnancy (day 8)",
                 "in_vitro": "X-Mens vs Wang in vitro decidualization"}[key]
        row = xmens_summary.loc[label]
        if int(row["genes compared"]) != len(table) or int(row["same direction"]) != int(table["same direction"].astype(bool).sum()):
            raise SystemExit(f"{sheet_name}: SUMMARY_xmens_comparisons.csv does not match the gene table")
        xmens_sheets.append((sheet_name, table))
        xmens_rows.append((label, sheet_name, table, float(row["% expected by chance"])))
    n_called = len(xmens_sheets[0][1])

    # ---- Wang vs Wang sheets, from the published lists
    lists = {k: pd.read_csv(f"{a.wang_dir}/wang_{k}.csv") for k in ("oil", "pregnancy", "ivd")}
    for k, v in lists.items():
        if v.gene.duplicated().any() or not set(v.direction) <= {"up", "down"}:
            raise SystemExit(f"wang_{k}.csv: duplicated genes or unexpected directions")
    wang_sheets, wang_rows = [], []
    for sheet_name, label, A, B in WANG_SHEETS:
        table = wang_pair(lists, A, B)
        wang_sheets.append((sheet_name, table))
        wang_rows.append((label, sheet_name, table, chance(table[A[1]], table[B[1]])))

    # ---- the two human comparisons, with the functions and note templates of add_human_comparisons_to_ivd_workbook.py
    consensus = best.merge(orth[["human_gene", "human_ensembl_gene_id"]], left_on="gene", right_on="human_gene", how="left", validate="one_to_one")
    if len(consensus) != n123 or consensus.human_ensembl_gene_id.isna().any():
        raise SystemExit("expected every consensus gene to carry a human Ensembl ID")
    s_cmp, values = AH.slidetags_values(consensus, pd.read_csv(a.slidetags, float_precision="round_trip"))
    h_args = SimpleNamespace(hammoud_release=a.hammoud_release, hammoud_release_name=a.hammoud_release_name,
                             hammoud_config=a.hammoud_config, hammoud_donors=a.hammoud_donors,
                             hammoud_contrast=a.hammoud_contrast, hammoud_label=a.hammoud_label)
    hammoud, u_values, old_list_rows = AH.hammoud_inputs(h_args)
    u_cmp = AH.compare(consensus, hammoud)
    values.update(u_values, u_testable=int(u_cmp.tested.sum()), u_untested=AH.join_words(u_cmp.gene[~u_cmp.tested]),
                  u_called=int(u_cmp.called.sum()))
    s_sheet, u_sheet = AH.gene_sheet(s_cmp, AH.SLIDETAGS_COL), AH.gene_sheet(u_cmp, AH.HAMMOUD_COL)
    values["u_same"] = int(u_sheet["same direction"].sum())
    old = old_list_rows.set_index("human_gene").support_status
    if old.reindex(list(AH.REMOVED_BY_CORRECTION)).isin(["supported", "reversed"]).any():
        raise SystemExit("a gene removed by the 22 October correction is DE in the Hammoud fit; revise the Hammoud notes")
    if u_cmp[u_cmp.gene.isin(AH.ADDED_BY_CORRECTION)].called.any():
        raise SystemExit("a gene added by the 22 October correction is DE in the Hammoud fit; revise the Hammoud notes")
    human_rows = [
        (f"Human Slide-tags DE vs human consensus in vitro ({n123} genes, {values['s_testable']} testable)", AH.SLIDETAGS_SHEET, s_sheet, AH.chance(s_sheet, AH.SLIDETAGS_COL)),
        (f"Hammoud DE vs human consensus in vitro ({n123} genes, {values['u_testable']} testable)", AH.HAMMOUD_SHEET, u_sheet, AH.chance(u_sheet, AH.HAMMOUD_COL)),
    ]

    # ---- Summary and Notes
    gene_sheets = [("consensus 123 genes", cons_sheet)] + xmens_sheets + wang_sheets + [(AH.SLIDETAGS_SHEET, s_sheet), (AH.HAMMOUD_SHEET, u_sheet)]
    rows = []
    for r, (label, sheet_name, table, chance_value) in enumerate(xmens_rows + wang_rows + human_rows, start=2):
        rows.append({"comparison": label,
                     "genes compared": f"=COUNTA({hs.column_range(sheet_name, table, 'gene')})",
                     "same direction": f"=COUNTIF({hs.column_range(sheet_name, table, 'same direction')},TRUE)",
                     "% same direction": f"=ROUND(100*C{r}/B{r},1)",
                     "% expected by chance": chance_value,
                     "genes listed on": f"gene sheet: {sheet_name}"})
    summary = pd.DataFrame(rows)
    notes = [
        f"Human consensus in vitro list: {n123} genes called by at least 2 of 3 human in vitro decidualization studies (Aghajanova, "
        f"Rytkönen, Stadtmauer), with every call in the same direction; {n_orth} have a one-to-one mouse ortholog. Full list on the "
        '"consensus 123 genes" sheet.',
        f"X-Mens vs human consensus: of the {n_orth} orthologs, {n_called} are called by X-Mens (adjusted p < 0.05, fold change > 2x); "
        "% same direction = share changed the same way as in vitro.",
        *FIXED_NOTES,
        *[line.format(**values) for line in AH.SLIDETAGS_NOTES + AH.HAMMOUD_NOTES],
    ]
    # the workbook of record fixes the width of the Summary's "genes listed on" column at 40 (session code of 20 Oct 2026)
    hs.write_workbook(a.out, summary, notes, gene_sheets, summary_widths={"F": 40})

    # ---- check the saved formulas against the sheets, and report
    saved = load_workbook(a.out)
    for r, (label, sheet_name, table, chance_value) in enumerate(xmens_rows + wang_rows + human_rows, start=2):
        n = hs.evaluate(saved, saved["Summary"][f"B{r}"].value)
        same = hs.evaluate(saved, saved["Summary"][f"C{r}"].value)
        if (n, same) != (len(table), int(table["same direction"].astype(bool).sum())):
            raise SystemExit(f"Summary row {r}: formulas give {(n, same)}, the table gives {(len(table), int(table['same direction'].sum()))}")
        print(f"{label}: {same} of {n} same direction ({100 * same / n:.1f}%; {chance_value}% expected by chance)")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
