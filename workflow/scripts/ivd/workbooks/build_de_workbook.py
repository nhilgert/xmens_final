#!/usr/bin/env python3
"""Build the pseudobulk DE workbook: decidual vs fibroblast in X-Mens, human Slide-tags and Hammoud (Ulrich et al. 2024).

Sheets, in the layout of IVD_comparisons_with_genes.xlsx (helpers in build_house_style_workbooks.py):
  Summary              one row per dataset: genes tested, DE genes, up and down, as formulas over the gene sheets
  Notes                one note per row: what each fit compares, the DE rule, sources
  X-Mens DE            every gene that passed the expression filter, alphabetically
  human Slide-tags DE  the same, plus per-donor log2 ratios and % detected
  Hammoud DE           the same, for one fit of the Ulrich DESeq2 release
Gene sheet columns: gene, Ensembl gene ID, mean normalized count, <dataset> log2 fold change (unshrunk), <dataset>
standard error, <dataset> p, <dataset> adjusted p, DE, sign (up or down, DE genes only). The two p-value columns are
shown in scientific notation.

DE: adjusted p < 0.05 and |unshrunk log2 fold change| > 1; checked against each table's significant_fdr05_lfc1 flag.

Inputs
  --xmens                  xmens_broad_all_genes_with_uninduced.csv (xmens_with_uninduced_20261027.zip)
  --slidetags              human_decidual_vs_fibroblast_all_genes.csv (human_DE_rerun_2026_09_15_object.zip, 29 October 2026 version)
  --hammoud-release        Ulrich release directory (ulrich_contrasts.csv, ulrich_<contrast>_all_genes.csv, run_metadata.json)
  --hammoud-release-name   the release's name, quoted in the notes
  --hammoud-config         the release's config (YAML)
  --hammoud-donors         the contrast's pseudobulk table (intermediate/ulrich_<contrast>_metadata.tsv of the release)
  --hammoud-contrast       contrast id (default: the release's primary contrast)
  --hammoud-label          description of the fit for the notes (default: the contrast's label in the release)
  --house-style            build_house_style_workbooks.py
  --out                    output workbook path (must not exist)
"""
import argparse
import importlib.util
import json
import os
from collections import Counter

import numpy as np
import pandas as pd
import yaml
from openpyxl import Workbook

SIGN = "sign"
P_FORMAT = "0.00E+00"
READ = dict(keep_default_na=False, na_values=[""], float_precision="round_trip")
# The X-Mens and Slide-tags notes below describe these two fits; stop if the input tables come from other fits.
XMENS_FIT = dict(tested=16784, de=4703)
SLIDETAGS_FIT = dict(tested=20211, de=3198)   # fitted in the cluster release environment, 29 October 2026 (3,195 in the 26 October local fit)
SLIDETAGS_DONOR_COLUMNS = {
    "mck5 (menstrual) log2 ratio": "log2_cpm_ratio_mck5_menstrual",
    "mck6 (secretory) log2 ratio": "log2_cpm_ratio_mck6_secretory",
    "mck5 decidual % detected": "pct_detected_mck5_menstrual_decidual",
    "mck5 fibroblast % detected": "pct_detected_mck5_menstrual_fibroblast",
    "mck6 decidual % detected": "pct_detected_mck6_secretory_decidual",
    "mck6 fibroblast % detected": "pct_detected_mck6_secretory_fibroblast",
}
PHASE_WORDS = {"luteal": "luteal", "follicular": "follicular", "not applicable": "with no recorded phase"}

XMENS_NOTE = (
    "X-Mens: decidual nuclei (all decidual subtypes) and endometrial fibroblast nuclei (all fibroblast subtypes), summed "
    "per timepoint over the five induced timepoints, plus the uninduced sample's fibroblasts (included since 27 October "
    "2026): 11 pseudobulks, 43,133 nuclei. PyDESeq2 with ~ timepoint + group, so decidual nuclei are compared with "
    "fibroblasts of the same timepoint (4 residual degrees of freedom).")
SLIDETAGS_NOTE = (
    "Human Slide-tags: decidual and fibroblast nuclei of the two donors (mck5, menstrual; mck6, secretory), summed per "
    "donor (cell_type labels): 4 pseudobulks, 47,533 nuclei. PyDESeq2 with ~ donor + group, which leaves 1 residual degree "
    "of freedom, so fold changes and the per-donor columns are the informative part; the p-values are not evidence that a "
    "difference replicates across women.")
HAMMOUD_NOTE = (
    "Hammoud (Ulrich et al. 2024, PNAS 121:e2404775121; human uterus single-cell atlas): {label}: {donors} donors "
    "({phases}), each with at least {floor} cells of both cell types, summed per donor and cell type: {pseudobulks} "
    "pseudobulks, {cells:,} cells. PyDESeq2 with {design}, which leaves {residual_df} residual degrees of freedom. \"up\" = "
    "higher in decidual stroma cells.")


def load_house_style(path):
    spec = importlib.util.spec_from_file_location("house_style", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sig3(p):
    return float(f"{p:.3g}") if pd.notna(p) else None


def tested_genes(table, name):
    """The genes a fit tested, after checking the DE flag against the rule."""
    t = table[table.production_eligible.astype(bool)].copy()
    rule = (t.padj_for_call < 0.05) & (t.log2FoldChange_unshrunk.abs() > 1)
    if not (rule == t.significant_fdr05_lfc1.astype(bool)).all():
        raise ValueError(f"{name}: significant_fdr05_lfc1 does not follow the DE rule")
    return t


def de_table(t, who, extra):
    d = pd.DataFrame({"gene": t.gene_symbol.fillna("").astype(str).values,
                      "Ensembl gene ID": t.ensembl_gene_id.values,
                      "mean normalized count": t.baseMean.round(2).values,
                      f"{who} log2 fold change": t.log2FoldChange_unshrunk.round(3).values,
                      f"{who} standard error": t.lfcSE.round(3).values,
                      f"{who} p": t.pvalue.map(sig3).values,
                      f"{who} adjusted p": t.padj_for_call.map(sig3).values,
                      "DE": t.significant_fdr05_lfc1.astype(bool).values})
    d[SIGN] = np.where(~d.DE, None, np.where(t.log2FoldChange_unshrunk.values > 0, "up", "down"))
    for column, source in extra.items():
        d[column] = t[source].round(3).values
    d = d.assign(_blank=d.gene.eq("")).sort_values(["_blank", "gene", "Ensembl gene ID"], kind="stable")
    return d.drop(columns="_blank").reset_index(drop=True)


def check_fit(t, expected, name):
    got = dict(tested=len(t), de=int(t.significant_fdr05_lfc1.astype(bool).sum()))
    if got != expected:
        raise ValueError(f"{name}: table has {got}, the notes describe {expected}; update the notes for the new fit")


def hammoud_inputs(a):
    contrasts = pd.read_csv(os.path.join(a.hammoud_release, "ulrich_contrasts.csv")).set_index("contrast_id")
    primary = contrasts.index[contrasts.primary.astype(str).str.lower() == "true"]
    if len(primary) != 1:
        raise ValueError("the release must have exactly one primary contrast")
    cid = a.hammoud_contrast or primary[0]
    info = contrasts.loc[cid]
    table = pd.read_csv(os.path.join(a.hammoud_release, f"ulrich_{cid}_all_genes.csv"), **READ)
    tested = tested_genes(table, "Hammoud")
    if len(tested) != int(info.eligible_genes) or int(tested.significant_fdr05_lfc1.astype(bool).sum()) != int(info.strict_calls):
        raise ValueError(f"{cid}: all-genes table does not match its row in ulrich_contrasts.csv")
    if not tested.estimable.astype(bool).all():
        raise ValueError(f"{cid}: some tested genes are not estimable")
    config = yaml.safe_load(open(a.hammoud_config))
    spec = [c for c in config["human"]["contrasts"] if c["id"] == cid]
    if len(spec) != 1:
        raise ValueError(f"{cid} is not defined in {a.hammoud_config}")
    floor = int(spec[0].get("cell_floor", config["human"]["cell_floor"]))
    if float(config["human"]["eligibility_cpm_threshold"]) != 1.0:
        raise ValueError("the notes say CPM >= 1; the release used another threshold")
    run = json.load(open(os.path.join(a.hammoud_release, "run_metadata.json")))
    if run["packages"].get("pydeseq2") != "0.5.4":
        raise ValueError("the notes say PyDESeq2 0.5.4; the release used another version")
    donors = pd.read_csv(a.hammoud_donors, sep="\t")
    if donors.donor.nunique() != int(info.donors) or len(donors) != int(info.pseudobulks) or int(donors.n_cells.min()) < floor:
        raise ValueError(f"{a.hammoud_donors} does not match contrast {cid}")
    phases = Counter(donors.drop_duplicates("donor").menstrual_phase_at_collection.astype(str))
    label = a.hammoud_label or str(info.label)
    note = dict(label=label[0].upper() + label[1:], donors=int(info.donors),
                phases=", ".join(f"{n} {PHASE_WORDS.get(p, p)}" for p, n in sorted(phases.items(), key=lambda kv: (-kv[1], kv[0]))),
                floor=floor, pseudobulks=int(info.pseudobulks), cells=int(donors.n_cells.sum()), design=str(info.design),
                residual_df=int(info.residual_df), table=f"ulrich_{cid}_all_genes.csv", release=a.hammoud_release_name)
    return tested, note


def de_row(hs, label, sheet_name, table):
    return {"comparison": label,
            "genes tested": f"=COUNTA({hs.column_range(sheet_name, table, 'Ensembl gene ID')})",
            "DE genes": f"=COUNTIF({hs.column_range(sheet_name, table, 'DE')},TRUE)",
            "higher in decidual (up)": f"=COUNTIF({hs.column_range(sheet_name, table, SIGN)},\"up\")",
            "higher in fibroblasts (down)": f"=COUNTIF({hs.column_range(sheet_name, table, SIGN)},\"down\")",
            "genes listed on": f"gene sheet: {sheet_name}"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("xmens", "slidetags", "hammoud-release", "hammoud-release-name", "hammoud-config", "hammoud-donors",
                 "house-style", "out"):
        ap.add_argument(f"--{name}", required=True)
    ap.add_argument("--hammoud-contrast", default=None)
    ap.add_argument("--hammoud-label", default=None)
    a = ap.parse_args()
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; refusing to overwrite")
    hs = load_house_style(a.house_style)

    xmens = tested_genes(pd.read_csv(a.xmens, **READ), "X-Mens")
    check_fit(xmens, XMENS_FIT, "X-Mens")
    slidetags = tested_genes(pd.read_csv(a.slidetags, **READ), "human Slide-tags")
    check_fit(slidetags, SLIDETAGS_FIT, "human Slide-tags")
    hammoud, h = hammoud_inputs(a)

    sheets = [("X-Mens DE", "X-Mens", de_table(xmens, "X-Mens", {}),
               "X-Mens decidual vs fibroblast nuclei (5 induced timepoints, plus uninduced fibroblasts)"),
              ("human Slide-tags DE", "human", de_table(slidetags, "human", SLIDETAGS_DONOR_COLUMNS),
               "Human Slide-tags decidual vs fibroblast nuclei (2 donors)"),
              ("Hammoud DE", "Hammoud", de_table(hammoud, "Hammoud", {}),
               f"Hammoud decidual stroma vs C7-Fibroblast cells ({h['donors']} donors)")]
    summary = pd.DataFrame([de_row(hs, label, name, table) for name, _, table, label in sheets])

    zero = [f"{dataset}: {', '.join(t.gene[t[f'{who} adjusted p'] == 0])}"
            for (_, who, t, _), dataset in zip(sheets, ("X-Mens", "human Slide-tags", "Hammoud")) if (t[f"{who} adjusted p"] == 0).any()]
    notes = [
        "Pseudobulk differential expression of decidual vs fibroblast cells in three datasets: X-Mens (mouse, single-nucleus), "
        "human Slide-tags (2 donors, single-nucleus) and the Hammoud-lab human uterus atlas (Ulrich et al. 2024, single-cell). "
        "Each gene sheet lists every gene that passed the expression filter, alphabetically.",
        'DE: adjusted p < 0.05 (Benjamini-Hochberg) and fold change > 2x (|log2 fold change| > 1, unshrunk). A positive log2 '
        'fold change and "up" in the "sign" column mean higher in decidual cells or nuclei; "down" means higher in '
        'fibroblasts. "sign" is filled in for DE genes only.',
        XMENS_NOTE,
        SLIDETAGS_NOTE,
        HAMMOUD_NOTE.format(**h),
        "All three fits (PyDESeq2, the Python implementation of DESeq2): genes with CPM >= 1 in at least half (rounded up) of "
        "either cell type's pseudobulks were tested; Cook's outlier refit on, no independent filtering; PyDESeq2 0.5.4.",
        "Human Slide-tags per-donor columns: log2 ratio = log2((CPM decidual + 0.5) / (CPM fibroblast + 0.5)) within one "
        "donor; % detected = percent of nuclei with at least one UMI of the gene.",
    ]
    if zero:
        notes.append(f"An adjusted p of 0 ({'; '.join(zero)}) is below the smallest number the software represents.")
    notes += [
        "The sheets are not matched to each other. The two human sheets can be joined on Ensembl gene ID; for cross-species "
        "comparisons, match genes through one-to-one orthologs (Ensembl 116); many decidual genes, PRL among them, have none "
        "(PRL has 26 mouse prolactin-family co-orthologs).",
        'Summary columns "genes tested", "DE genes", "higher in decidual (up)" and "higher in fibroblasts (down)" are formulas '
        'over the gene sheets (COUNTA / COUNTIF on the "Ensembl gene ID", "DE" and "sign" columns).',
        "Sources: X-Mens, xmens_broad_all_genes_with_uninduced.csv (xmens_with_uninduced_20261027.zip); human Slide-tags, "
        "human_decidual_vs_fibroblast_all_genes.csv (human_DE_rerun_2026_09_15_object.zip); "
        f"Hammoud, {h['table']} ({h['release']}).",
    ]

    workbook = Workbook()
    workbook.remove(workbook.active)
    hs.add_sheet(workbook, "Summary", summary)
    hs.add_notes(workbook, notes)
    for name, who, table, _ in sheets:
        if len(name) > 31:
            raise ValueError(f"sheet name longer than Excel allows: {name}")
        sheet = hs.add_sheet(workbook, name, table)
        for column in (f"{who} p", f"{who} adjusted p"):
            j = list(table.columns).index(column) + 1
            for (cell,) in sheet.iter_rows(min_row=2, max_row=len(table) + 1, min_col=j, max_col=j):
                cell.number_format = P_FORMAT
    workbook.save(a.out)
    for name, who, table, _ in sheets:
        print(f"{name}: {len(table)} genes tested, {int(table.DE.sum())} DE "
              f"({int((table[SIGN] == 'up').sum())} up, {int((table[SIGN] == 'down').sum())} down)")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
