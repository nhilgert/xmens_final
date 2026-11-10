# Code as run inside the analysis session that produced IVD_comparisons_with_genes.xlsx, version 4 (artifact version cac1e20c-1a14-4f2c-8286-39f204aedfba).
# What it did: rewording of the consensus-list note.
# Extracted from the session record for reading. It is not a standalone script: paths point at files in the
# Claude Science artifact store, and it will not run elsewhere as is.

import os
import re
from openpyxl import load_workbook

NEW = '/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/c0e13e9f-9222-48cd-9a8a-e96f51d5ab5b/v5cda8f85_IVD_comparisons_with_genes.xlsx'

os.makedirs('ivd_workbook_v4', exist_ok=True)
OUTW = 'IVD_comparisons_with_genes.xlsx'

wb4 = load_workbook(NEW)
hits_ = [c for r in wb4['Notes'].iter_rows() for c in r if isinstance(c.value, str) and c.value.startswith('Human consensus in vitro list:')]
assert len(hits_) == 1
cell = hits_[0]
cell.value = ('Human consensus in vitro list: 123 genes called by at least 2 of 3 human in vitro decidualization studies (Aghajanova, Rytkönen, Stadtmauer), '
              'with every call in the same direction; 114 have a one-to-one mouse ortholog. Full list on the "consensus 123 genes" sheet.')
wb4.save(OUTW)