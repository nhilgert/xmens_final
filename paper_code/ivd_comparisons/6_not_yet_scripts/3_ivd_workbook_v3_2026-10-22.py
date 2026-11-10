# Code as run inside the analysis session that produced IVD_comparisons_with_genes.xlsx, version 3 (artifact version 5cda8f85-f3d5-4fa0-9d75-5a4eaaf6a94d).
# What it did: consensus sheet and X-Mens vs consensus sheet rebuilt for the corrected consensus list (Stadtmauer criteria).
# Extracted from the session record for reading. It is not a standalone script: paths point at files in the
# Claude Science artifact store, and it will not run elsewhere as is.

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

ARMS = ("MPA+PGE2", "MPA+cAMP")
TP_LABEL = {"6h": "6 h", "48h": "48 h", "14d": "14 d", "3d": "3 d", "8d": "8 d"}
XLAB = {"supported": "same direction", "reversed": "opposite direction",
        "tested_not_called": "tested, not called", "ineligible_or_unavailable": "not testable"}
STUDY_NAME = {"Aghajanova": "Aghajanova", "Rytkonen": "Rytkönen", "Stadtmauer": "Stadtmauer"}
COL_AGH = "Aghajanova 2011 (E2+P4 vs E2)"
COL_RYT = "Rytkönen 2019 (cAMP+MPA vs untreated)"
COL_STAD = "Stadtmauer 2022 (MPA+PGE2 or MPA+cAMP vs base)"
FONT = Font(name="Calibri", size=11)
BOLD = Font(name="Calibri", size=11, bold=True)

best = pd.read_csv('/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/0300d4de-4525-4bbb-a115-e38fac12d56e/v68d85693_consensus123_stadtmauer_authors_criteria.csv')
sup = pd.read_csv('/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/d6a060b9-9e7e-43a5-b952-9df775742ced/v340acabc_consensus_ortholog_xmens_broad_support_rows_stadtmauer_authors_criteria.csv').set_index("human_gene")
agh = pd.read_csv('/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/c81e24b2-4be7-4ad3-9187-21a3b7af8718/v219bf010_aghajanova2011_invitro_decidualization_nonendo.csv')
ryt = pd.read_csv('/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/bf156b1d-b8ad-4933-9fdf-fc3e66502323/v428064de_rytkonen2019_complete_signature.csv')

def timepoint_call(sig, gene, order):
    d = sig[sig.g == gene]
    if d.empty:
        return ""
    dirs = sorted(set(d.direction))
    return f"{dirs[0] if len(dirs) == 1 else 'conflicting'} ({', '.join(TP_LABEL[t] for t in order if t in set(d.timepoint))})"

def studies_label(mature_papers):
    return "; ".join(STUDY_NAME[s.strip()] for s in mature_papers.split(";"))

agh["g"] = agh.gene
ryt["g"] = ryt.gene

# Build per-gene table
rows = []
for r in best.itertuples(index=False):
    g = r.gene
    stad = ""
    strict = g in sup.index and str(sup.loc[g, "strict_one_to_one"]) == "True"
    rows.append({
        "gene": g,
        "consensus direction": r.consensus_direction,
        "studies supporting": int(r.n_mature),
        COL_AGH: timepoint_call(agh, g, ["6h", "48h", "14d"]),
        COL_RYT: timepoint_call(ryt, g, ["3d", "8d"]),
        COL_STAD: stad,
        "mouse ortholog (1:1)": sup.loc[g, "mouse_ortholog"] if strict else "",
        "X-Mens": XLAB[sup.loc[g, "support_status"]] if g in sup.index else "no 1:1 mouse ortholog",
        "X-Mens log2 fold change": round(float(sup.loc[g, "log2FoldChange_unshrunk"]), 2)
        if g in sup.index and pd.notna(sup.loc[g, "log2FoldChange_unshrunk"]) else None
    })

per_gene = pd.DataFrame(rows).sort_values(["studies supporting", "gene"], ascending=[False, True]).reset_index(drop=True)

called = sup[sup.support_status.isin(["supported", "reversed"])].reset_index()
called = called.sort_values("human_gene")

cons_sheet = per_gene[["gene", "consensus direction", "studies supporting"]].rename(
    columns={"consensus direction": "direction in vitro", "studies supporting": "studies agreeing (of 3)"})
best_idx = best.set_index("gene")
cons_sheet["studies"] = cons_sheet.gene.map(lambda g: studies_label(best_idx.loc[g, "mature_papers"]))
cons_sheet["mouse ortholog (1:1)"] = per_gene["mouse ortholog (1:1)"].values
cons_sheet = cons_sheet.sort_values(["studies agreeing (of 3)", "gene"], ascending=[False, True])

xs = pd.DataFrame({
    "gene": called.human_gene,
    "mouse gene": called.mouse_gene,
    "in vitro (human consensus)": called.consensus_direction,
    "X-Mens": called.observed_direction,
    "X-Mens log2 fold change": called.log2FoldChange_unshrunk.round(3),
    "X-Mens adjusted p": called.padj_for_call.map(lambda p: float(f"{p:.3g}")),
    "same direction": called.support_status.eq("supported").astype(bool)
})

def write_sheet(wb, name, df, index):
    if name in wb.sheetnames:
        del wb[name]
    ws = wb.create_sheet(name, index)
    ws.append(list(df.columns))
    for rec in df.itertuples(index=False):
        ws.append([None if (isinstance(v, float) and pd.isna(v)) else v for v in rec])
    for row in ws.iter_rows():
        for c in row:
            c.font = BOLD if c.row == 1 else FONT
    for j, col in enumerate(df.columns, 1):
        w = max([len(str(col))] + [len(str(v)) for v in df[col].head(500)])
        ws.column_dimensions[get_column_letter(j)].width = min(w + 2, 110)
    ws.freeze_panes = "A2"

wb = load_workbook('/Users/nickhilgert/.claude-science/orgs/06f88a99-9105-4941-a8cf-87e6b02aff65/artifacts/proj_3111dcb6f5e8/c0e13e9f-9222-48cd-9a8a-e96f51d5ab5b/vcd577b71_IVD_comparisons_with_genes.xlsx')
order = wb.sheetnames

write_sheet(wb, "consensus 123 genes", cons_sheet, order.index("consensus 123 genes"))
write_sheet(wb, "X-Mens vs consensus in vitro", xs, order.index("X-Mens vs consensus in vitro"))

ws = wb["Summary"]
last = len(xs) + 1
n_orth = int((sup.strict_one_to_one.astype(str) == "True").sum())
label = f"X-Mens vs human consensus in vitro ({len(best)} genes, {n_orth} with mouse ortholog)"
ws["A2"] = label
ws["B2"] = f"=COUNTA('X-Mens vs consensus in vitro'!$A$2:$A${last})"
ws["C2"] = f"=COUNTIF('X-Mens vs consensus in vitro'!$G$2:$G${last},TRUE)"
called_sup = sup[sup.support_status.isin(["supported", "reversed"])]
p_cons = (called_sup.consensus_direction == "up").mean()
p_x = (called_sup.observed_direction == "up").mean()
chance = round(100 * (p_cons * p_x + (1 - p_cons) * (1 - p_x)), 1)
ws["E2"] = chance
for cell in [ws["A2"], ws["B2"], ws["C2"], ws["E2"]]:
    cell.font = FONT

notes = wb["Notes"]
n_called = len(called_sup)
n_same = int((called_sup.support_status == "supported").sum())
notes["A2"] = (f"Human consensus in vitro list: {len(best)} genes changed in the same direction in at least 2 of 3 human in vitro "
               f"decidualization studies (Aghajanova, Rytkönen, Stadtmauer); {n_orth} have a one-to-one mouse ortholog. "
               "Full list on the \"consensus 123 genes\" sheet.")
notes["A3"] = (f"X-Mens vs human consensus: of the {n_orth} orthologs, {n_called} are called by X-Mens (adjusted p < 0.05, "
               "fold change > 2x); % same direction = share changed the same way as in vitro.")
revision_note = ("Consensus list revised 2026-10-22: Stadtmauer genes with tabulated values are counted only where a treatment meets the authors' criteria "
                 "(FDR < 1e-6, fold change > 1.5, > 3 TPM vs base medium). This removed KLF4 and TSHZ1 and added ADAMTS5 and CLU; X-Mens changed from "
                 f"19 of 49 to {n_same} of {n_called} genes in the same direction.")
notes.cell(row=notes.max_row + 1, column=1, value=revision_note).font = FONT

wb.save('IVD_comparisons_with_genes.xlsx')