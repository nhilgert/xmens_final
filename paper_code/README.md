# paper_code: Methods code that no rule runs yet

Everything here is cited in the Methods but is not yet a step of the pipeline. Each folder is ported into
workflow/ (see docs/PORTING.md) and then leaves; when this folder is empty, the pipeline builds everything in the
paper. Paths inside these files still point to their old locations (shown in brackets).

- ulrich_mapping/: comparison of the fibroblast and decidual clusters with the Ulrich et al. 2024 atlas, run by
  hand with its submit*.sh scripts (analysis/hammoud/fib_mapping).
- epithelial_notebooks/: generators of the two epithelial annotation notebooks (tools/); the notebooks are in
  notebooks/public/ and step 28 renders them.
- mi_overlap/: the cross-species MI overlap notebook and its executed copy (notebooks/mi_dpt).
- notebooks/: eight notebooks (analysis/notebooks, notebooks/, notebooks/spatial_qc).
