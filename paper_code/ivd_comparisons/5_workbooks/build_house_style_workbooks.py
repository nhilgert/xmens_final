"""Write workbooks in the layout of IVD_comparisons_with_genes.xlsx.

Layout: a Summary sheet first (one row per comparison, counts as Excel formulas over the gene sheets), then a Notes
sheet (header "note", one note per row, column width 110), then the gene sheets. Every cell is Calibri 11, the header
row is bold and frozen, and there are no filters, fills, number formats or wrapped text. Each column is as wide as its
longest entry plus 2 (formulas are not counted), which reproduces the widths of IVD_comparisons_with_genes.xlsx.
"""
import re

import numpy as np
import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

FONT, BOLD = Font(name="Calibri", size=11), Font(name="Calibri", size=11, bold=True)
NOTES_WIDTH = 110


def cell_value(value):
    """A pandas/numpy value as a plain Python value openpyxl writes as-is (missing values become empty cells)."""
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if np.isnan(value) else float(value)
    return value


def add_sheet(workbook, name, table, widths=None):
    """widths: optional {column letter: width} that replaces the rule for those columns."""
    sheet = workbook.create_sheet(name)
    sheet.append(list(table.columns))
    for row in table.itertuples(index=False):
        sheet.append([cell_value(v) for v in row])
    for row in sheet.iter_rows():
        for cell in row:
            cell.font = BOLD if cell.row == 1 else FONT
    for j, column in enumerate(table.columns, 1):
        entries = [str(column)] + [str(cell_value(v)) for v in table[column]
                                   if cell_value(v) is not None and not (isinstance(v, str) and v.startswith("="))]
        sheet.column_dimensions[get_column_letter(j)].width = max(len(e) for e in entries) + 2
    for letter, width in (widths or {}).items():
        sheet.column_dimensions[letter].width = width
    sheet.freeze_panes = "A2"
    return sheet


def add_notes(workbook, lines):
    sheet = workbook.create_sheet("Notes")
    sheet.append(["note"])
    for line in lines:
        sheet.append([line])
    for row in sheet.iter_rows():
        for cell in row:
            cell.font = BOLD if cell.row == 1 else FONT
    sheet.column_dimensions["A"].width = NOTES_WIDTH
    sheet.freeze_panes = "A2"


def column_range(sheet_name, table, column_name):
    """Absolute range of one column of a gene sheet, data rows only, for use in a Summary formula."""
    letter = get_column_letter(list(table.columns).index(column_name) + 1)
    return f"'{sheet_name}'!${letter}$2:${letter}${len(table) + 1}"


def write_workbook(path, summary, notes, gene_sheets, summary_widths=None):
    """summary: DataFrame for the Summary sheet; notes: list of strings; gene_sheets: list of (name, DataFrame);
    summary_widths: optional fixed widths for Summary columns (IVD_comparisons_with_genes.xlsx fixes B-F)."""
    for name, _ in gene_sheets:
        assert len(name) <= 31, f"sheet name longer than Excel allows: {name}"
    workbook = Workbook()
    workbook.remove(workbook.active)
    add_sheet(workbook, "Summary", summary, summary_widths)
    add_notes(workbook, notes)
    for name, table in gene_sheets:
        add_sheet(workbook, name, table)
    workbook.save(path)


def evaluate(workbook, formula):
    """Evaluate the COUNTA / COUNTIF / ROUND formulas these Summary sheets use, to check them against the data."""
    m = re.fullmatch(r"=(COUNTA|COUNTIF)\('([^']+)'!\$([A-Z]+)\$(\d+):\$[A-Z]+\$(\d+)(?:,(TRUE|\"[^\"]*\"))?\)", formula)
    if m:
        kind, sheet, column, first, last, criterion = m.groups()
        values = [workbook[sheet][f"{column}{r}"].value for r in range(int(first), int(last) + 1)]
        if kind == "COUNTA":
            return sum(v is not None for v in values)
        return sum(v is True for v in values) if criterion == "TRUE" else sum(v == criterion.strip('"') for v in values)
    m = re.fullmatch(r"=ROUND\(100\*([A-Z]+\d+)/([A-Z]+\d+),1\)", formula)
    if m:
        return m.groups()
    raise ValueError(formula)


def format_fingerprint(path):
    """Fonts, fills, number formats, wrapping, filters and frozen panes across all sheets."""
    workbook = load_workbook(path)
    cells = [c for s in workbook for r in s.iter_rows() for c in r if c.value is not None]
    return dict(fonts=sorted({(c.font.name, c.font.size, bool(c.font.bold)) for c in cells}),
                fills=sum(bool(c.fill.fill_type) for c in cells),
                number_formats=sorted({c.number_format for c in cells} - {"General"}),
                wrapped=sum(bool(c.alignment.wrap_text) for c in cells),
                filters=[s.title for s in workbook if s.auto_filter.ref], tables=[s.title for s in workbook if s.tables],
                frozen={s.freeze_panes for s in workbook}, bold_outside_header=sum(bool(c.font.bold) for c in cells if c.row > 1))
