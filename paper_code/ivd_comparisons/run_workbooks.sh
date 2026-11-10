#!/bin/bash
# Step 5: build the two workbooks in xmens_figures/DESeq2_DE_and_IVD_comparisons/ from the three fits' tables.
#   pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx   (build_de_workbook.py)
#   IVD_comparisons_with_genes.xlsx                             (add_human_comparisons_to_ivd_workbook.py)
# These are the commands that built the 29 October 2026 files, with paths made relative.
#
# Usage, from this folder:  bash run_workbooks.sh INPUT_DIR OUTPUT_DIR
# INPUT_DIR must contain (README.md, "Inputs for step 5", says where each file comes from):
#   xmens_broad_all_genes_with_uninduced.csv
#   human_decidual_vs_fibroblast_all_genes.csv
#   hammoud_release/                       the Hammoud (Ulrich) release tables, including
#                                          config_ulrich_dS_vs_C7fib_luteal_follicular.yaml and
#                                          intermediate_pseudobulk_tables/ulrich_dS_vs_C7fib_all_donors_metadata.tsv
#   IVD_comparisons_with_genes_without_human.xlsx
#   consensus123_stadtmauer_authors_criteria.csv
#   consensus123_stadtmauer_authors_criteria_mouse_orthology.csv
# Needs Python with pandas, numpy, openpyxl and PyYAML (set PYTHON to choose the interpreter; default: python).
# Neither script overwrites an existing output file.
set -euo pipefail
IN=${1:?input folder}
OUT=${2:?output folder}
PY=${PYTHON:-python}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"

HAMMOUD=(--hammoud-release "$IN/hammoud_release"
         --hammoud-release-name release_20261028_v7_dS_vs_C7fib_luteal_follicular
         --hammoud-config "$IN/hammoud_release/config_ulrich_dS_vs_C7fib_luteal_follicular.yaml"
         --hammoud-donors "$IN/hammoud_release/intermediate_pseudobulk_tables/ulrich_dS_vs_C7fib_all_donors_metadata.tsv"
         --hammoud-contrast dS_vs_C7fib_all_donors
         --hammoud-label "all decidual stroma (dS1-3) vs C7-Fibroblast, within donor, no phase filter"
         --house-style "$HERE/5_workbooks/build_house_style_workbooks.py")

"$PY" "$HERE/5_workbooks/build_de_workbook.py" \
  --xmens "$IN/xmens_broad_all_genes_with_uninduced.csv" \
  --slidetags "$IN/human_decidual_vs_fibroblast_all_genes.csv" \
  "${HAMMOUD[@]}" \
  --out "$OUT/pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx"

"$PY" "$HERE/5_workbooks/add_human_comparisons_to_ivd_workbook.py" \
  --workbook "$IN/IVD_comparisons_with_genes_without_human.xlsx" \
  --consensus "$IN/consensus123_stadtmauer_authors_criteria.csv" \
  --orthology "$IN/consensus123_stadtmauer_authors_criteria_mouse_orthology.csv" \
  --slidetags "$IN/human_decidual_vs_fibroblast_all_genes.csv" \
  "${HAMMOUD[@]}" \
  --out "$OUT/IVD_comparisons_with_genes.xlsx"
