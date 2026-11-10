# Code as run inside the analysis session that produced Wang_comparisons_simple_with_genes.xlsx (artifact version 8a05c1e7-1948-426e-9a2e-2e56c1077abb).
# What it did: the workbook the IVD workbook started from (X-Mens vs Wang and Wang vs Wang sheets).
# Extracted from the session record for reading. It is not a standalone script: paths point at files in the
# Claude Science artifact store, and it will not run elsewhere as is.

import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

xm = pd.read_csv('/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/35be2137-4052-45b6-a4b9-c154e97cb5d6/v0fbf3296_xmens_broad_all_genes.csv')
xm = xm[xm.production_eligible & xm.estimable].drop_duplicates('gene_symbol', keep=False).set_index('gene_symbol')

lists = {k: pd.read_csv(f'{{artifact:{aid}}}').drop_duplicates('gene', keep=False).set_index('gene').direction
         for k, aid in (('oil', '13eb8264-00c6-4e54-8a8b-5c4978a374d4'),
                        ('pregnancy', '0de469d2-778b-4d1f-ae23-800bb37e1be1'),
                        ('ivd', '595656b6-9de1-4846-8c01-66e2ef38663d'))}
oil, preg = lists['oil'], lists['pregnancy']
ivd_rows = pd.read_csv('/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/e231587d-bf90-4ff9-ba49-27e51d7b42e1/ve08e9bb3_wang_ivd_support_rows.csv')

def oil_pregnancy_class(gene):
    in_o, in_p = gene in oil.index, gene in preg.index
    if in_o and in_p: return 'both, same direction' if oil[gene] == preg[gene] else 'both, opposite directions'
    return 'oil only' if in_o else 'pregnancy only' if in_p else 'neither'

def xmens_sheet(ref, name):
    d = pd.DataFrame({'gene': ref.index, f'wang_{name}_direction': ref.values})
    d['oil_pregnancy_class'] = d.gene.map(oil_pregnancy_class)
    d['tested_in_xmens'] = d.gene.isin(xm.index)
    d['xmens_log2fc'] = d.gene.map(xm.log2FoldChange_unshrunk).round(3)
    d['xmens_adjusted_p'] = d.gene.map(xm.padj_for_call)
    d['xmens_called'] = d.gene.map(xm.significant_fdr05_lfc1).fillna(False).astype(bool)
    d['xmens_direction'] = np.where(d.xmens_log2fc > 0, 'up', np.where(d.xmens_log2fc < 0, 'down', ''))
    d['same_direction'] = np.where(d.xmens_called, d.xmens_direction == d[f'wang_{name}_direction'], None)
    d['counted_in_percentage'] = d.xmens_called
    return d.sort_values(['counted_in_percentage', 'gene'], ascending=[False, True]).reset_index(drop=True)

S = {}
S['XMens vs Wang oil'] = xmens_sheet(oil, 'oil')
S['XMens vs Wang pregnancy'] = xmens_sheet(preg, 'pregnancy')

iv = ivd_rows.rename(columns={'wang_ivd_gene': 'gene', 'wang_ivd_direction': 'wang_in_vitro_direction', 'xmens_current_gene_symbol': 'xmens_gene_matched',
                              'xmens_effect_log2_fold_change_unshrunk': 'xmens_log2fc', 'xmens_adjusted_p_value_for_call': 'xmens_adjusted_p',
                              'xmens_called_strict_lfc1': 'xmens_called'})
iv = iv[['gene', 'wang_in_vitro_direction', 'xmens_match_route', 'xmens_gene_matched', 'xmens_log2fc', 'xmens_adjusted_p', 'xmens_called', 'xmens_direction']].copy()
iv['xmens_called'] = iv.xmens_called.fillna(False).astype(bool); iv['xmens_log2fc'] = iv.xmens_log2fc.round(3)
iv['same_direction'] = np.where(iv.xmens_called, iv.xmens_direction == iv.wang_in_vitro_direction, None); iv['counted_in_percentage'] = iv.xmens_called
S['XMens vs Wang in vitro'] = iv.sort_values(['counted_in_percentage', 'gene'], ascending=[False, True]).reset_index(drop=True)

def wang_pair(a, a_name, b, b_name):
    d = pd.DataFrame({'gene': a.index, f'wang_{a_name}_direction': a.values})
    d[f'wang_{b_name}_direction'] = d.gene.map(b).fillna('')
    d[f'in_{b_name}_list'] = d.gene.isin(b.index)
    d['same_direction'] = np.where(d[f'in_{b_name}_list'], d[f'wang_{a_name}_direction'] == d[f'wang_{b_name}_direction'], None)
    d['counted_in_percentage'] = d[f'in_{b_name}_list']
    return d.sort_values(['counted_in_percentage', 'gene'], ascending=[False, True]).reset_index(drop=True)

wv = ivd_rows.set_index('wang_ivd_gene')
iv_oil = pd.DataFrame({'gene': wv.index, 'wang_in_vitro_direction': wv.wang_ivd_direction.values, 'wang_oil_direction': wv.wang_oil_direction.fillna('').values,
                       'in_oil_list': wv.wang_oil_called.astype(bool).values})
iv_oil['same_direction'] = np.where(iv_oil.in_oil_list, iv_oil.wang_in_vitro_direction == iv_oil.wang_oil_direction, None); iv_oil['counted_in_percentage'] = iv_oil.in_oil_list
iv_preg = pd.DataFrame({'gene': wv.index, 'wang_in_vitro_direction': wv.wang_ivd_direction.values, 'wang_pregnancy_direction': wv.wang_pregnancy_direction.fillna('').values,
                        'in_pregnancy_list': wv.wang_pregnancy_called.astype(bool).values})
iv_preg['same_direction'] = np.where(iv_preg.in_pregnancy_list, iv_preg.wang_in_vitro_direction == iv_preg.wang_pregnancy_direction, None); iv_preg['counted_in_percentage'] = iv_preg.in_pregnancy_list
S['Wang in vitro vs oil'] = iv_oil.sort_values(['counted_in_percentage', 'gene'], ascending=[False, True]).reset_index(drop=True)
S['Wang in vitro vs pregnancy'] = iv_preg.sort_values(['counted_in_percentage', 'gene'], ascending=[False, True]).reset_index(drop=True)
S['Wang oil vs pregnancy'] = wang_pair(oil, 'oil', preg, 'pregnancy')

RENAME = {'gene': 'gene', 'wang_oil_direction': 'oil (Wang)', 'wang_pregnancy_direction': 'pregnancy (Wang)', 'wang_in_vitro_direction': 'in vitro (Wang)',
          'xmens_gene_matched': 'X-Mens gene', 'xmens_direction': 'X-Mens', 'xmens_log2fc': 'X-Mens log2 fold change',
          'xmens_adjusted_p': 'X-Mens adjusted p', 'same_direction': 'same direction'}

PLAN = [('X-Mens vs Wang oil-induced decidualization', 'X-Mens vs oil', 'XMens vs Wang oil', 'wang_oil_direction', 'xmens_direction'),
        ('X-Mens vs Wang natural pregnancy (day 8)', 'X-Mens vs pregnancy', 'XMens vs Wang pregnancy', 'wang_pregnancy_direction', 'xmens_direction'),
        ('X-Mens vs Wang in vitro decidualization', 'X-Mens vs in vitro', 'XMens vs Wang in vitro', 'wang_in_vitro_direction', 'xmens_direction'),
        ('Wang oil-induced vs Wang natural pregnancy', 'oil vs pregnancy', 'Wang oil vs pregnancy', 'wang_oil_direction', 'wang_pregnancy_direction'),
        ('Wang in vitro vs Wang oil-induced', 'in vitro vs oil', 'Wang in vitro vs oil', 'wang_in_vitro_direction', 'wang_oil_direction'),
        ('Wang in vitro vs Wang natural pregnancy', 'in vitro vs pregnancy', 'Wang in vitro vs pregnancy', 'wang_in_vitro_direction', 'wang_pregnancy_direction')]

sheets, rows = {}, []
for label, name, src, a_col, b_col in PLAN:
    d = S[src][S[src].counted_in_percentage.astype(bool)].copy()
    keep = ['gene', a_col] + (['xmens_gene_matched'] if 'xmens_gene_matched' in d else []) + \
           (['xmens_direction', 'xmens_log2fc', 'xmens_adjusted_p'] if b_col == 'xmens_direction' else [b_col]) + ['same_direction']
    d = d[keep]; d['same_direction'] = d.same_direction.astype(bool)
    if 'xmens_gene_matched' in d:
        d['xmens_gene_matched'] = np.where(d.xmens_gene_matched == d.gene, '', d.xmens_gene_matched)
    if 'xmens_adjusted_p' in d:
        d['xmens_adjusted_p'] = d['xmens_adjusted_p'].map(lambda p: float(f'{p:.3g}'))
    d = d.rename(columns=RENAME).sort_values('gene').reset_index(drop=True)
    sheets[name] = d
    pa, pb = (S[src].loc[S[src].counted_in_percentage.astype(bool), a_col] == 'up').mean(), (S[src].loc[S[src].counted_in_percentage.astype(bool), b_col] == 'up').mean()
    rows.append({'comparison': label, 'genes compared': len(d), 'same direction': int(d['same direction'].sum()),
                 '% same direction': round(100 * d['same direction'].mean(), 1), '% expected by chance': round(100 * (pa * pb + (1 - pa) * (1 - pb)), 1)})
SUMMARY_MIN = pd.DataFrame(rows)
NOTES = pd.DataFrame({'note': [
    'X-Mens rows: genes in the Wang list that X-Mens calls (adjusted p < 0.05, fold change > 2x, decidual cells vs endometrial fibroblasts); % same direction = share changed the same way as in Wang.',
    'Wang vs Wang rows: genes in both Wang lists; % same direction = share listed the same way in both.',
    '% expected by chance: agreement if the two lists were unrelated, given how many genes go up vs down in each.',
    'Each comparison has its own sheet listing the genes compared.',
    'Wang lists: Wang et al. 2020, Genes 11:935 (published gene lists, direction only).']})
SUMMARY_MIN.to_csv('SUMMARY_wang_comparisons_simple.csv', index=False)

wb = Workbook(); wb.remove(wb.active)
def write(name, df):
    ws = wb.create_sheet(name); ws.append(list(df.columns))
    for r in df.itertuples(index=False): ws.append([None if (isinstance(v, float) and np.isnan(v)) else (bool(v) if isinstance(v, np.bool_) else v) for v in r])
    for row in ws.iter_rows():
        for c in row: c.font = Font(name='Calibri', size=11, bold=(c.row == 1))
    for j, col in enumerate(df.columns, 1):
        w = max([len(str(col))] + [len(str(v)) for v in df[col].head(500)])
        ws.column_dimensions[get_column_letter(j)].width = min(w + 2, 110)
    ws.freeze_panes = 'A2'
write('Summary', SUMMARY_MIN); write('Notes', NOTES)
for name, d in sheets.items(): write(name, d)
wb.save('Wang_comparisons_simple_with_genes.xlsx')