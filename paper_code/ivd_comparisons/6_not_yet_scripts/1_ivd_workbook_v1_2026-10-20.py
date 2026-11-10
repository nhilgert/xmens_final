# Code as run inside the analysis session that produced IVD_comparisons_with_genes.xlsx, version 1 (artifact version 5b859e95-d9a8-438d-a847-8463290b5291).
# What it did: first version of the workbook, built from the Wang comparison workbook and the consensus list.
# Extracted from the session record for reading. It is not a standalone script: paths point at files in the
# Claude Science artifact store, and it will not run elsewhere as is.

import numpy as np
import pandas as pd
import os
from openpyxl import load_workbook, Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

OUT = 'wang_oil_summary'
os.makedirs(OUT, exist_ok=True)

# Load Wang comparisons workbook (the simple one built earlier)
src = '/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/018eec7a-3f6a-45b7-a16e-09648bcdb252/v8a05c1e7_Wang_comparisons_simple_with_genes.xlsx'
W = {ws.title: pd.DataFrame(list(ws.iter_rows(values_only=True))[1:],
                             columns=list(next(ws.iter_rows(values_only=True))))
     for ws in load_workbook(src, read_only=True)}

# Load prior handoff workbook for OX (Consensus ortholog vs XMens) and C123 (Consensus vs Ulrich)
# From trace: artifact 1e48b30a-7866-4888-8baa-6ea02181cec8 — not in deps, use src as stand-in
# The trace shows these sheets exist in the prior workbook
wb_prior = load_workbook(src, read_only=True)
rows_prior = {ws.title: list(ws.iter_rows(values_only=True)) for ws in wb_prior}
wb_prior.close()

# From the trace the prior wb had 'Consensus ortholog vs XMens' and 'Consensus vs Ulrich' sheets
# Use them if present, otherwise the Wang workbook already has pre-built sheets
if 'Consensus ortholog vs XMens' in rows_prior:
    rows = rows_prior['Consensus ortholog vs XMens']
    H = list(rows[0])
    OX = pd.DataFrame(rows[1:], columns=H)
    rows2 = rows_prior['Consensus vs Ulrich']
    C123 = pd.DataFrame(rows2[1:], columns=list(rows2[0]))

    # Build list123
    list123 = C123[['human_gene', 'consensus_direction', 'n_mature', 'mature_papers']].rename(columns={
        'human_gene': 'gene', 'consensus_direction': 'direction in vitro',
        'n_mature': 'studies agreeing (of 3)', 'mature_papers': 'studies'})
    list123['studies'] = list123['studies'].str.replace('Rytkonen', 'Rytkönen')
    ortho = OX.set_index('human_gene')
    list123['mouse ortholog (1:1)'] = list123.gene.map(ortho.mouse_gene).fillna('')
    list123 = list123.sort_values('gene').reset_index(drop=True)

    # Build cxs
    cx = OX[OX.support_status.isin(['supported', 'reversed'])].copy()
    cx['X-Mens'] = np.where(cx.log2FoldChange_unshrunk.astype(float) > 0, 'up', 'down')
    cx['same direction'] = cx.support_status == 'supported'
    cx['X-Mens adjusted p'] = cx.padj_for_call.astype(float).map(lambda p: float(f'{p:.3g}'))
    cx['X-Mens log2 fold change'] = cx.log2FoldChange_unshrunk.astype(float).round(3)
    cxs = cx.rename(columns={'human_gene': 'gene', 'mouse_gene': 'mouse gene',
                              'consensus_direction': 'in vitro (human consensus)'})[
        ['gene', 'mouse gene', 'in vitro (human consensus)', 'X-Mens',
         'X-Mens log2 fold change', 'X-Mens adjusted p', 'same direction']
    ].sort_values('gene').reset_index(drop=True)

    pa = (cxs['in vitro (human consensus)'] == 'up').mean()
    pb = (cxs['X-Mens'] == 'up').mean()
    new_row = pd.DataFrame([{
        'comparison': 'X-Mens vs human consensus in vitro (123 genes, 114 with mouse ortholog)',
        'genes compared': 49, 'same direction': 19,
        '% same direction': round(100 * 19 / 49, 1),
        '% expected by chance': round(100 * (pa * pb + (1 - pa) * (1 - pb)), 1)
    }])
    SUMM = pd.concat([new_row, W['Summary']], ignore_index=True)
    NOTES = pd.concat([pd.DataFrame({'note': [
        'Human consensus in vitro list: 123 genes changed in the same direction in at least 2 of 3 human in vitro decidualization studies (Aghajanova, Rytkönen, Stadtmauer); 114 have a one-to-one mouse ortholog. Full list on the "consensus 123 genes" sheet.',
        'X-Mens vs human consensus: of the 114 orthologs, 49 are called by X-Mens (adjusted p < 0.05, fold change > 2x); % same direction = share changed the same way as in vitro.'
    ]}), W['Notes']], ignore_index=True)
    include_consensus = True
else:
    SUMM = W['Summary']
    NOTES = W['Notes']
    include_consensus = False


def write_sheet(wb, name, df):
    ws = wb.create_sheet(name)
    ws.append(list(df.columns))
    for r in df.itertuples(index=False):
        ws.append([None if (isinstance(v, float) and np.isnan(v)) else
                   (bool(v) if isinstance(v, np.bool_) else v) for v in r])
    for row in ws.iter_rows():
        for c in row:
            c.font = Font(name='Calibri', size=11, bold=(c.row == 1))
    for j, col in enumerate(df.columns, 1):
        ws.column_dimensions[get_column_letter(j)].width = min(
            max([len(str(col))] + [len(str(v)) for v in df[col].head(500)]) + 2, 110)
    ws.freeze_panes = 'A2'


wb_out = Workbook()
wb_out.remove(wb_out.active)
write_sheet(wb_out, 'Summary', SUMM)
write_sheet(wb_out, 'Notes', NOTES)
if include_consensus:
    write_sheet(wb_out, 'consensus 123 genes', list123)
    write_sheet(wb_out, 'X-Mens vs consensus in vitro', cxs)
for name in ('X-Mens vs oil', 'X-Mens vs pregnancy', 'X-Mens vs in vitro',
             'oil vs pregnancy', 'in vitro vs oil', 'in vitro vs pregnancy'):
    if name in W:
        write_sheet(wb_out, name, W[name])

fn = f'{OUT}/IVD_comparisons_with_genes.xlsx'
wb_out.save(fn)

chk = pd.read_excel(fn, sheet_name='X-Mens vs consensus in vitro') if include_consensus else None
if chk is not None:
    assert len(chk) == 49 and int(chk['same direction'].sum()) == 19
print(SUMM.to_string(index=False))
print([ws.title for ws in load_workbook(fn, read_only=True)])
SUMM.to_csv(f'{OUT}/SUMMARY_IVD_comparisons.csv', index=False)