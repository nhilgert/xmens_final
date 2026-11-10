#!/usr/bin/env python3
"""Step 4: write a new IVD_comparisons_with_genes.xlsx from step 3's comparison tables.

Starts from the previous workbook and replaces only what depends on the X-Mens fit:
  - the four X-Mens gene sheets (vs consensus in vitro, vs oil, vs pregnancy, vs in vitro), from step 3's
    genes_xmens_vs_<reference>.csv
  - Summary rows 2 to 5 (X-Mens vs consensus, oil, pregnancy and in vitro): the formula ranges and the fixed
    "% expected by chance"
  - Notes: one revision line is appended
Every other sheet, row and note is kept as it was; the script checks this, and checks that each Summary formula
evaluates to step 3's counts, before it reports success. The previous workbook is never modified.

Usage
  python update_ivd_workbook.py --previous IVD_comparisons_with_genes.xlsx --comparisons STEP3_DIR
         --note "revision note" --out NEW_WORKBOOK.xlsx
"""
import argparse
import re

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

FONT, BOLD = Font(name="Calibri", size=11), Font(name="Calibri", size=11, bold=True)

# gene sheet, step 3 file, Summary row, and the column of the sheet that holds "same direction"
SHEETS = [("X-Mens vs consensus in vitro", "consensus", 2, "G"),
          ("X-Mens vs oil", "oil", 3, "F"),
          ("X-Mens vs pregnancy", "pregnancy", 4, "F"),
          ("X-Mens vs in vitro", "in_vitro", 5, "G")]


def write_sheet(workbook, name, table):
    """Replace a sheet in place, in the workbook's style: Calibri 11, bold header, fitted widths, frozen header."""
    index = workbook.sheetnames.index(name)
    del workbook[name]
    sheet = workbook.create_sheet(name, index)
    sheet.append(list(table.columns))
    for record in table.itertuples(index=False):
        sheet.append([None if (isinstance(v, float) and np.isnan(v)) else (bool(v) if isinstance(v, np.bool_) else v)
                      for v in record])
    for row in sheet.iter_rows():
        for cell in row:
            cell.font = BOLD if cell.row == 1 else FONT
    for j, column in enumerate(table.columns, 1):
        width = max([len(str(column))] + [len(str(v)) for v in table[column].head(500)])
        sheet.column_dimensions[get_column_letter(j)].width = min(width + 2, 110)
    sheet.freeze_panes = "A2"


def evaluate(workbook, formula):
    """Evaluate the Summary's COUNTA / COUNTIF(...,TRUE) formulas against the saved sheets."""
    match = re.match(r"=(COUNTA|COUNTIF)\('([^']+)'!\$([A-Z])\$2:\$[A-Z]\$(\d+)(?:,TRUE)?\)", formula)
    function, sheet, column, last = match.group(1), match.group(2), match.group(3), int(match.group(4))
    values = [workbook[sheet][f"{column}{k}"].value for k in range(2, last + 1)]
    return sum(v is not None for v in values) if function == "COUNTA" else sum(v is True for v in values)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--previous", required=True, help="the previous IVD_comparisons_with_genes.xlsx")
    parser.add_argument("--comparisons", required=True, help="step 3's output folder")
    parser.add_argument("--note", required=True, help="revision line appended to the Notes sheet")
    parser.add_argument("--out", required=True, help="the new workbook")
    args = parser.parse_args()

    summary = pd.read_csv(f"{args.comparisons}/SUMMARY_xmens_comparisons.csv").set_index("comparison")
    workbook = load_workbook(args.previous)
    previous = load_workbook(args.previous)
    order = workbook.sheetnames
    expected = {}
    for sheet_name, file_key, summary_row, same_column in SHEETS:
        table = pd.read_csv(f"{args.comparisons}/genes_xmens_vs_{file_key}.csv", keep_default_na=False, na_values=[""],
                            float_precision="round_trip")
        write_sheet(workbook, sheet_name, table)
        label = workbook["Summary"][f"A{summary_row}"].value
        assert label in summary.index, f"Summary row {summary_row} label not in step 3's summary: {label}"
        last = len(table) + 1
        workbook["Summary"][f"B{summary_row}"] = f"=COUNTA('{sheet_name}'!$A$2:$A${last})"
        workbook["Summary"][f"C{summary_row}"] = f"=COUNTIF('{sheet_name}'!${same_column}$2:${same_column}${last},TRUE)"
        workbook["Summary"][f"E{summary_row}"] = float(summary.loc[label, "% expected by chance"])
        expected[summary_row] = (int(summary.loc[label, "genes compared"]), int(summary.loc[label, "same direction"]))
    notes = workbook["Notes"]
    notes.cell(row=notes.max_row + 1, column=1, value=args.note).font = FONT
    assert workbook.sheetnames == order
    workbook.save(args.out)

    # Checks on the saved file
    saved = load_workbook(args.out)
    for row, (compared, same) in expected.items():
        got = (evaluate(saved, saved["Summary"][f"B{row}"].value), evaluate(saved, saved["Summary"][f"C{row}"].value))
        assert got == (compared, same), f"Summary row {row}: formulas give {got}, step 3 gives {(compared, same)}"
    replaced = {name for name, *_ in SHEETS} | {"Summary", "Notes"}
    for name in saved.sheetnames:
        if name not in replaced:
            old = [[c.value for c in r] for r in previous[name].iter_rows()]
            new = [[c.value for c in r] for r in saved[name].iter_rows()]
            assert old == new, f"sheet {name} changed"
    for row in range(1, saved["Summary"].max_row + 1):
        if row not in expected:
            assert all(previous["Summary"].cell(row, col).value == saved["Summary"].cell(row, col).value for col in range(1, 7)), row
    old_notes = [r[0] for r in previous["Notes"].iter_rows(values_only=True)]
    assert [r[0] for r in saved["Notes"].iter_rows(values_only=True)] == old_notes + [args.note]
    for row, (compared, same) in expected.items():
        print(saved["Summary"][f"A{row}"].value, "|", compared, same, round(100 * same / compared, 1),
              "| chance", saved["Summary"][f"E{row}"].value)


if __name__ == "__main__":
    main()
