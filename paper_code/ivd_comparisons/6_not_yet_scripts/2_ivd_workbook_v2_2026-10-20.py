# Code as run inside the analysis session that produced IVD_comparisons_with_genes.xlsx, version 2 (artifact version cd577b71-954c-4171-8e8a-0630541cf624).
# What it did: second version, rebuilt from version 1 and the Wang comparison workbook.
# Extracted from the session record for reading. It is not a standalone script: paths point at files in the
# Claude Science artifact store, and it will not run elsewhere as is.

import pandas as pd
import numpy as np
from openpyxl import load_workbook, Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

OUT = 'wang_oil_summary'
import os
os.makedirs(OUT, exist_ok=True)

src = '/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/018eec7a-3f6a-45b7-a16e-09648bcdb252/v8a05c1e7_Wang_comparisons_simple_with_genes.xlsx'
W = {ws.title: pd.DataFrame(list(ws.iter_rows(values_only=True))[1:], columns=list(next(ws.iter_rows(values_only=True)))) for ws in load_workbook(src, read_only=True)}
print({k: v.shape for k, v in W.items()})

# Load the consensus vs XMens data from the simple workbook
wb_src = load_workbook(src, read_only=True)

# We need the consensus ortholog vs XMens data - load from the Wang simple workbook
# The OX and C123 data came from an earlier workbook - reconstruct from available data
# Load from the simple workbook's X-Mens vs consensus in vitro sheet if it exists,
# otherwise we need to rebuild from the release tables

# Load the existing IVD comparisons workbook to get the data
src2 = '/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/c0e13e9f-9222-48cd-9a8a-e96f51d5ab5b/v5b859e95_IVD_comparisons_with_genes.xlsx'
wb2 = load_workbook(src2, read_only=True)
ws_cx = wb2['X-Mens vs consensus in vitro']
rows_cx = list(ws_cx.iter_rows(values_only=True))
cxs = pd.DataFrame(rows_cx[1:], columns=list(rows_cx[0]))

ws_sum = wb2['Summary']
hdr_sum = [c.value for c in ws_sum[1]]
sum_rows = list(ws_sum.iter_rows(values_only=True))[1:]
SUMM_orig = pd.DataFrame(sum_rows, columns=hdr_sum)

ws_notes = wb2['Notes']
notes_rows = list(ws_notes.iter_rows(values_only=True))
NOTES_orig = pd.DataFrame(notes_rows[1:], columns=list(notes_rows[0]))

ws_list123 = wb2['consensus 123 genes']
rows_list123 = list(ws_list123.iter_rows(values_only=True))
list123 = pd.DataFrame(rows_list123[1:], columns=list(rows_list123[0]))

# Verify counts
assert len(cxs) == 49 and int(cxs['same direction'].astype(bool).sum()) == 19

SHEET = {'X-Mens vs human consensus in vitro (123 genes, 114 with mouse ortholog)': 'X-Mens vs consensus in vitro',
         'X-Mens vs Wang oil-induced decidualization': 'X-Mens vs oil', 'X-Mens vs Wang natural pregnancy (day 8)': 'X-Mens vs pregnancy',
         'X-Mens vs Wang in vitro decidualization': 'X-Mens vs in vitro', 'Wang oil-induced vs Wang natural pregnancy': 'oil vs pregnancy',
         'Wang in vitro vs Wang oil-induced': 'in vitro vs oil', 'Wang in vitro vs Wang natural pregnancy': 'in vitro vs pregnancy'}

def write(wb, name, df):
    ws = wb.create_sheet(name)
    ws.append(list(df.columns))
    for r in df.itertuples(index=False):
        ws.append([None if (isinstance(v, float) and np.isnan(v)) else (bool(v) if isinstance(v, np.bool_) else v) for v in r])
    for row in ws.iter_rows():
        for c in row:
            c.font = Font(name='Calibri', size=11, bold=(c.row == 1))
    for j, col in enumerate(df.columns, 1):
        ws.column_dimensions[get_column_letter(j)].width = min(max([len(str(col))] + [len(str(v)) for v in df[col].head(500)]) + 2, 110)
    ws.freeze_panes = 'A2'

wb_out = Workbook()
wb_out.remove(wb_out.active)
write(wb_out, 'Summary', SUMM_orig)
write(wb_out, 'Notes', NOTES_orig)
write(wb_out, 'consensus 123 genes', list123)
write(wb_out, 'X-Mens vs consensus in vitro', cxs)
for name in ('X-Mens vs oil', 'X-Mens vs pregnancy', 'X-Mens vs in vitro', 'oil vs pregnancy', 'in vitro vs oil', 'in vitro vs pregnancy'):
    write(wb_out, name, W[name])

fn = f'{OUT}/IVD_comparisons_with_genes.xlsx'
wb_out.save(fn)

# Now add formulas to the Summary sheet
wb = load_workbook(fn)
ws = wb['Summary']
hdr = [c.value for c in ws[1]]
print(hdr)

before = {}
for r in range(2, ws.max_row + 1):
    comp = ws.cell(r, 1).value
    sh = SHEET[comp]
    g = wb[sh]
    col = [c.value for c in g[1]].index('same direction') + 1
    L = get_column_letter(col)
    n = g.max_row
    before[comp] = (ws.cell(r, 2).value, ws.cell(r, 3).value, ws.cell(r, 4).value)
    ws.cell(r, 2).value = f"=COUNTA('{sh}'!$A$2:$A${n})"
    ws.cell(r, 3).value = f"=COUNTIF('{sh}'!${L}$2:${L}${n},TRUE)"
    ws.cell(r, 4).value = f"=ROUND(100*C{r}/B{r},1)"
    ws.cell(r, 6).value = f"gene sheet: {sh}"
    for c in (2, 3, 4, 6):
        ws.cell(r, c).font = Font(name='Calibri', size=11)

ws.cell(1, 6).value = 'genes listed on'
ws.cell(1, 6).font = Font(name='Calibri', size=11, bold=True)
ws.column_dimensions['F'].width = 40

notes = wb['Notes']
notes.append(['Summary columns "genes compared", "same direction" and "% same direction" are formulas over the gene sheets (COUNTA / COUNTIF on the "same direction" column). "% expected by chance" is a fixed value.'])
notes.cell(notes.max_row, 1).font = Font(name='Calibri', size=11)

wb.save(fn)

# verify: evaluate the formulas independently from the saved gene sheets
chk = load_workbook(fn, read_only=True)
for comp, sh in SHEET.items():
    d = pd.DataFrame(list(chk[sh].iter_rows(values_only=True))[1:], columns=list(next(chk[sh].iter_rows(values_only=True))))
    n, same = len(d), int(d['same direction'].astype(bool).sum())
    assert (n, same, round(100 * same / n, 1)) == before[comp], (comp, n, same, before[comp])
print('formulas reproduce the stored values for all', len(SHEET), 'rows')
print([ws.cell(2, c).value for c in range(1, 7)])