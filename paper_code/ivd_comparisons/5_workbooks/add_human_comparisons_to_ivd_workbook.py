#!/usr/bin/env python3
"""Add the two human in vivo comparisons to the IVD comparisons workbook.

Adds, after the existing sheets of IVD_comparisons_with_genes.xlsx:
  Summary  two rows: human Slide-tags DE vs the human consensus in vitro list, and Hammoud (Ulrich et al. 2024) DE vs
           the same list. Same columns and formulas as the existing rows.
  Notes    lines describing the two fits and where the values come from; with --previous-hammoud, one more line on
           how the Hammoud row changed from an earlier fit.
  Sheets   "Slide-tags vs consensus" and "Hammoud vs consensus": the consensus genes each fit calls, with direction,
           log2 fold change, adjusted p and whether the direction matches in vitro.
Existing cells are not changed.

Consensus genes are matched to each DE table by human Ensembl gene ID (from the consensus orthology table). A gene is
compared when the fit tests it (production_eligible and estimable) and calls it: adjusted p < 0.05 and |unshrunk
log2 fold change| > 1, the same rule as the X-Mens rows.

The Hammoud fit is read from a release of the Ulrich DESeq2 pipeline (ulrich_ivd_pipeline.py): the contrast's
all-genes table and its row in ulrich_contrasts.csv, the contrast's pseudobulk table, the release's support rows for
the consensus list before the 22 October correction, and the release config (cell floor). The fit description in the
notes is built from these files.

Inputs
  --workbook               IVD_comparisons_with_genes.xlsx without the human rows (artifact version fb72892e, 27 Oct 2026)
  --consensus              consensus123_stadtmauer_authors_criteria.csv (consensus123_stadtmauer_fix_bundle.zip)
  --orthology              consensus123_stadtmauer_authors_criteria_mouse_orthology.csv (same bundle; human Ensembl IDs)
  --slidetags              human_decidual_vs_fibroblast_all_genes.csv (human_DE_rerun_2026_09_15_object.zip)
  --hammoud-release        Ulrich release directory (ulrich_contrasts.csv, ulrich_<contrast>_all_genes.csv, support rows)
  --hammoud-release-name   the release's name, quoted in the notes
  --hammoud-config         the release's config (YAML)
  --hammoud-donors         the contrast's pseudobulk table (intermediate/ulrich_<contrast>_metadata.tsv of the release)
  --hammoud-contrast       contrast id (default: the release's primary contrast)
  --hammoud-label          description of the fit for the notes (default: the contrast's label in the release)
  --previous-hammoud       optional: all-genes table of the Hammoud fit an earlier version of this workbook showed
  --previous-hammoud-change what changed since that fit, as a sentence fragment for the notes (required with
                           --previous-hammoud); the donor counts are added from the tables
  --house-style            build_house_style_workbooks.py (shared workbook formatting helpers)
  --out                    output workbook path (must not exist)
"""
import argparse
import importlib.util
import os
from collections import Counter

import numpy as np
import pandas as pd
import yaml
from openpyxl import load_workbook

SLIDETAGS_SHEET, HAMMOUD_SHEET = "Slide-tags vs consensus", "Hammoud vs consensus"
SLIDETAGS_COL, HAMMOUD_COL = "human Slide-tags", "Hammoud"
# genes the 22 October consensus correction removed (KLF4, TSHZ1) and added (ADAMTS5, CLU)
REMOVED_BY_CORRECTION, ADDED_BY_CORRECTION = ("KLF4", "TSHZ1"), ("ADAMTS5", "CLU")
PHASE_WORDS = {"luteal": "luteal", "follicular": "follicular", "not applicable": "with no recorded phase"}
FILTER_RULE = "at least 1 CPM in at least half of the pseudobulks of one group"

SLIDETAGS_NOTES = [
    'Human Slide-tags DE vs human consensus (row and sheet re-added 2026-10-28): decidual vs fibroblast nuclei of the '
    'two Slide-tags donors (mck5, menstrual; mck6, secretory), one pseudobulk per donor and cell type (cell_type labels), '
    'PyDESeq2 with ~ donor + group. Of the 123 consensus genes, {s_testable} pass the expression filter ({s_untested} '
    'do not) and {s_called} are DE (adjusted p < 0.05, fold change > 2x); % same direction = share changed the same way '
    'as in vitro. "human Slide-tags" up = higher in decidual nuclei.',
    'With two donors the Slide-tags fit has 1 residual degree of freedom, so its adjusted p-values are not evidence of '
    'replication across women; {s_both} of the {s_called} DE genes change in the same direction in both donors '
    '({s_not_both} does not). Values from human_decidual_vs_fibroblast_all_genes.csv '
    '(human_DE_rerun_2026_09_15_object.zip).',
]
HAMMOUD_NOTES = [
    'Hammoud DE vs human consensus (row and sheet added 2026-10-28): Ulrich et al. 2024, PNAS 121:e2404775121 (Hammoud '
    'lab), human uterus single-cell atlas. {u_label}: {u_donors} donors ({u_phases}), one pseudobulk per donor and cell '
    'type (at least {u_floor} cells), PyDESeq2 with {u_design}: {u_pseudobulks} pseudobulks, {u_residual_df} residual '
    'degrees of freedom; {u_tested:,} genes tested, {u_all_calls} DE. Of the 123 consensus genes, {u_testable} pass the '
    'expression filter ({u_untested} do not) and {u_called} are DE (adjusted p < 0.05, fold change > 2x). "Hammoud" up '
    '= higher in decidual stroma.',
    'Hammoud values from {u_table} ({u_release}). That release compared the consensus list from before the 22 October '
    'correction; none of the four genes the correction changed (KLF4, TSHZ1, ADAMTS5, CLU) is DE in this fit, so the '
    'counts are the same with either list.',
]
PREVIOUS_NOTE = (
    'Hammoud fit revised 2026-10-28: {p_change} ({p_donors} donors before, {u_donors} now). The Hammoud row changed from '
    '{p_same} of {p_called} to {u_same} of {u_called} genes in the same direction; {p_kept_text} consensus genes DE in the '
    'earlier fit are DE in the same direction here.{p_lost_text}')


def load_house_style(path):
    spec = importlib.util.spec_from_file_location("house_style", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def join_words(words):
    words = list(words)
    if not words:
        return "none"
    return words[0] if len(words) == 1 else ", ".join(words[:-1]) + " and " + words[-1]


def compare(consensus, table):
    """One row per consensus gene: tested or not, called or not, direction."""
    t = table.drop_duplicates("ensembl_gene_id", keep=False).set_index("ensembl_gene_id")
    rows = []
    for r in consensus.itertuples():
        row = dict(gene=r.gene, ensembl_gene_id=r.human_ensembl_gene_id, consensus=r.consensus_direction, tested=False, called=False)
        if r.human_ensembl_gene_id in t.index:
            x = t.loc[r.human_ensembl_gene_id]
            lfc, padj = float(x.log2FoldChange_unshrunk), float(x.padj_for_call)
            tested = bool(x.production_eligible) and bool(x.estimable) and np.isfinite(lfc) and np.isfinite(padj)
            called = tested and padj < 0.05 and abs(lfc) > 1
            if called != bool(x.significant_fdr05_lfc1):
                raise ValueError(f"call rule disagrees with the table's own flag for {r.gene}")
            row.update(tested=tested, called=called, observed="up" if lfc > 0 else "down", lfc=lfc, padj=padj)
        rows.append(row)
    return pd.DataFrame(rows)


def gene_sheet(compared, column):
    d = compared[compared.called].sort_values("gene")
    return pd.DataFrame({"gene": d.gene.values, "in vitro (human consensus)": d.consensus.values, column: d.observed.values,
                         f"{column} log2 fold change": d.lfc.round(3).values,
                         f"{column} adjusted p": [float(f"{p:.3g}") for p in d.padj],
                         "same direction": (d.observed == d.consensus).values})


def chance(sheet, column):
    """Agreement expected if the two direction lists were unrelated, given their up/down proportions (as in the other rows)."""
    p = (sheet["in vitro (human consensus)"] == "up").mean()
    q = (sheet[column] == "up").mean()
    return round(100 * (p * q + (1 - p) * (1 - q)), 1)


def slidetags_values(consensus, slidetags):
    cmp_ = compare(consensus, slidetags)
    # per-donor consistency of the called genes (two donors, so this is the only replication available)
    ratio_cols = sorted(c for c in slidetags.columns if c.startswith("log2_cpm_ratio_"))
    if len(ratio_cols) != 2:
        raise ValueError(f"expected two per-donor ratio columns, got {ratio_cols}")
    st = slidetags.drop_duplicates("ensembl_gene_id", keep=False).set_index("ensembl_gene_id")
    called_ids = cmp_[cmp_.called].set_index("gene").ensembl_gene_id
    d1, d2 = np.sign(st.loc[called_ids, ratio_cols[0]].values), np.sign(st.loc[called_ids, ratio_cols[1]].values)
    fit = np.sign(st.loc[called_ids, "log2FoldChange_unshrunk"].values)
    consistent = (d1 == d2) & (d1 == fit)
    values = dict(s_testable=int(cmp_.tested.sum()), s_untested=join_words(cmp_.gene[~cmp_.tested]), s_called=int(cmp_.called.sum()),
                  s_both=int(consistent.sum()), s_not_both=join_words(list(called_ids.index[~consistent])))
    return cmp_, values


def hammoud_inputs(a):
    """The contrast's table, its description and the release's old-list support rows."""
    contrasts = pd.read_csv(os.path.join(a.hammoud_release, "ulrich_contrasts.csv")).set_index("contrast_id")
    primary = contrasts.index[contrasts.primary.astype(str).str.lower() == "true"]
    if len(primary) != 1:
        raise ValueError("the release must have exactly one primary contrast")
    cid = a.hammoud_contrast or primary[0]
    info = contrasts.loc[cid]
    table = pd.read_csv(os.path.join(a.hammoud_release, f"ulrich_{cid}_all_genes.csv"), float_precision="round_trip")
    eligible = table.production_eligible.astype(bool) & table.estimable.astype(bool)
    calls = table.significant_fdr05_lfc1.astype(bool)
    if int(eligible.sum()) != int(info.eligible_genes) or int(calls.sum()) != int(info.strict_calls):
        raise ValueError(f"{cid}: all-genes table does not match its row in ulrich_contrasts.csv")
    rows_name = "consensus123_ulrich_support_rows.csv" if cid == primary[0] else f"consensus123_ulrich_{cid}_support_rows.csv"
    old_list_rows = pd.read_csv(os.path.join(a.hammoud_release, rows_name))
    config = yaml.safe_load(open(a.hammoud_config))
    spec = [c for c in config["human"]["contrasts"] if c["id"] == cid]
    if len(spec) != 1:
        raise ValueError(f"{cid} is not defined in {a.hammoud_config}")
    floor = int(spec[0].get("cell_floor", config["human"]["cell_floor"]))
    donors = pd.read_csv(a.hammoud_donors, sep="\t")
    if donors.donor.nunique() != int(info.donors) or len(donors) != int(info.pseudobulks) or int(donors.n_cells.min()) < floor:
        raise ValueError(f"{a.hammoud_donors} does not match contrast {cid}")
    phases = Counter(donors.drop_duplicates("donor").menstrual_phase_at_collection.astype(str))
    phase_text = ", ".join(f"{n} {PHASE_WORDS.get(p, p)}" for p, n in sorted(phases.items(), key=lambda kv: (-kv[1], kv[0])))
    label = a.hammoud_label or str(info.label)
    values = dict(u_label=label[0].upper() + label[1:], u_label_lower=label, u_donors=int(info.donors), u_phases=phase_text,
                  u_floor=floor, u_design=str(info.design), u_pseudobulks=int(info.pseudobulks), u_residual_df=int(info.residual_df),
                  u_tested=int(eligible.sum()), u_all_calls=int(calls.sum()), u_table=f"ulrich_{cid}_all_genes.csv",
                  u_release=a.hammoud_release_name)
    return table, values, old_list_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("workbook", "consensus", "orthology", "slidetags", "hammoud-release", "hammoud-release-name",
                 "hammoud-config", "hammoud-donors", "house-style", "out"):
        ap.add_argument(f"--{name}", required=True)
    ap.add_argument("--hammoud-contrast", default=None)
    ap.add_argument("--hammoud-label", default=None)
    ap.add_argument("--previous-hammoud", default=None)
    ap.add_argument("--previous-hammoud-change", default=None)
    a = ap.parse_args()
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists; refusing to overwrite")
    if bool(a.previous_hammoud) != bool(a.previous_hammoud_change):
        raise SystemExit("--previous-hammoud and --previous-hammoud-change go together")
    hs = load_house_style(a.house_style)

    cons = pd.read_csv(a.consensus)
    orth = pd.read_csv(a.orthology)[["human_gene", "human_ensembl_gene_id"]]
    consensus = cons.merge(orth, left_on="gene", right_on="human_gene", how="left", validate="one_to_one")
    if len(consensus) != 123 or consensus.human_ensembl_gene_id.isna().any():
        raise ValueError("expected 123 consensus genes, each with a human Ensembl ID")

    s_cmp, values = slidetags_values(consensus, pd.read_csv(a.slidetags, float_precision="round_trip"))
    hammoud, u_values, old_list_rows = hammoud_inputs(a)
    u_cmp = compare(consensus, hammoud)
    values.update(u_values, u_testable=int(u_cmp.tested.sum()), u_untested=join_words(u_cmp.gene[~u_cmp.tested]),
                  u_called=int(u_cmp.called.sum()))
    s_sheet, u_sheet = gene_sheet(s_cmp, SLIDETAGS_COL), gene_sheet(u_cmp, HAMMOUD_COL)
    values["u_same"] = int(u_sheet["same direction"].sum())

    # the second Hammoud note says the 22 October correction does not change the counts; check it
    old = old_list_rows.set_index("human_gene").support_status
    if old.reindex(list(REMOVED_BY_CORRECTION)).isin(["supported", "reversed"]).any():
        raise ValueError("a gene removed by the 22 October correction is DE in this fit; revise the Hammoud notes")
    if u_cmp[u_cmp.gene.isin(ADDED_BY_CORRECTION)].called.any():
        raise ValueError("a gene added by the 22 October correction is DE in this fit; revise the Hammoud notes")

    notes = [line.format(**values) for line in SLIDETAGS_NOTES + HAMMOUD_NOTES]
    if a.previous_hammoud:
        p_table = pd.read_csv(a.previous_hammoud, float_precision="round_trip")
        p_donors = p_table.n_pseudobulks_numerator.dropna().unique()  # within-donor fit: one decidual pseudobulk per donor
        if len(p_donors) != 1:
            raise ValueError("cannot read the earlier fit's donor count")
        p_cmp = compare(consensus, p_table).set_index("gene")
        u_idx = u_cmp.set_index("gene")
        p_called = p_cmp.index[p_cmp.called]
        kept = [g for g in p_called if u_idx.loc[g, "called"] and u_idx.loc[g, "observed"] == p_cmp.loc[g, "observed"]]
        lost = [g for g in p_cmp.index[p_cmp.tested] if not u_idx.loc[g, "tested"]]
        lost_text = (f" {join_words(lost)} {'was' if len(lost) == 1 else 'were'} tested in the earlier fit but "
                     f"{'does' if len(lost) == 1 else 'do'} not pass the expression filter here ({FILTER_RULE})." if lost else "")
        p_sheet = gene_sheet(p_cmp.reset_index(), HAMMOUD_COL)
        notes.append(PREVIOUS_NOTE.format(**values, p_change=a.previous_hammoud_change, p_donors=int(p_donors[0]),
                                          p_same=int(p_sheet["same direction"].sum()), p_called=len(p_called),
                                          p_kept_text=(f"all {len(p_called)}" if len(kept) == len(p_called) else f"{len(kept)} of the {len(p_called)}"),
                                          p_lost_text=lost_text))

    wb = load_workbook(a.workbook)
    for name in (SLIDETAGS_SHEET, HAMMOUD_SHEET):
        if name in wb.sheetnames:
            raise ValueError(f"workbook already has a sheet named {name}")
    hs.add_sheet(wb, SLIDETAGS_SHEET, s_sheet)
    hs.add_sheet(wb, HAMMOUD_SHEET, u_sheet)

    summary = wb["Summary"]
    rows = (
        (f"Human Slide-tags DE vs human consensus in vitro (123 genes, {values['s_testable']} testable)", SLIDETAGS_SHEET, s_sheet, SLIDETAGS_COL),
        (f"Hammoud DE vs human consensus in vitro (123 genes, {values['u_testable']} testable)", HAMMOUD_SHEET, u_sheet, HAMMOUD_COL),
    )
    for label, sheet_name, sheet, column in rows:
        r = summary.max_row + 1
        cells = [label,
                 f"=COUNTA({hs.column_range(sheet_name, sheet, 'gene')})",
                 f"=COUNTIF({hs.column_range(sheet_name, sheet, 'same direction')},TRUE)",
                 f"=ROUND(100*C{r}/B{r},1)",
                 chance(sheet, column),
                 f"gene sheet: {sheet_name}"]
        for j, value in enumerate(cells, 1):
            summary.cell(r, j, value).font = hs.FONT
    summary.column_dimensions["A"].width = max(len(str(summary.cell(r, 1).value)) for r in range(1, summary.max_row + 1)) + 2

    notes_sheet = wb["Notes"]
    for line in notes:
        notes_sheet.append([line])
        notes_sheet.cell(notes_sheet.max_row, 1).font = hs.FONT
    wb.save(a.out)

    for label, sheet_name, sheet, _ in rows:
        print(f"{label}: {int(sheet['same direction'].sum())} of {len(sheet)} same direction")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
