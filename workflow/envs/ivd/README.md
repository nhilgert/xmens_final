# ivd environment

The fits of the in vitro / in vivo decidualization comparison (X-Mens, our endometrium, the Ulrich et al. 2024 atlas)
ran in `/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111`: Python 3.11.11, NumPy 2.2.0, pandas 2.3.3,
SciPy 1.17.1, anndata 0.12.11, PyDESeq2 0.5.4, scikit-learn 1.8.0, formulaic 1.2.2 (`requirements-frozen.txt`, frozen
2 Nov 2026). `config/config.yaml` names the interpreter as `ivd.python_fit`.

Why it is pinned: PyDESeq2 0.5.4 gives the same calls but different standard errors and p-values for a few genes under
NumPy 2.4.6 (X-Mens: 13 and 6 genes with standard errors differing by more than 1% in the two designs; our endometrium:
3 genes moved across the adjusted-p 0.05 line in a local run), so the tables of record reproduce exactly only in this
environment (paper_code/ivd_comparisons/2_xmens/README_mouse_DE_code.md and 3_human_slidetags/README_human_DE_code.md).

Recreate it with uv:

```
uv venv --python 3.11.11 /path/to/ivd_env
uv pip install --python /path/to/ivd_env/bin/python -r workflow/envs/ivd/requirements-frozen.txt
```

The table and workbook steps (consensus, comparisons, S9/S10, workbooks, verification) need only pandas, numpy, openpyxl
and PyYAML and run in the pipeline's rsc_uv environment (`ivd.python`).
