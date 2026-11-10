#!/usr/bin/env python3
"""Build the user-facing outputs of the Stadtmauer correction from the rebuild job's outputs.

Reads the folder written by rebuild_consensus123_stadtmauer_thresholds.py (--job-out) and
writes, into --outdir:
  TABLE_consensus123_per_gene_calls.csv        every consensus gene: call per study, ortholog, X-Mens result
  TABLE_consensus123_sources_by_study.csv      what each study contributed
  TABLE_consensus123_changes_vs_released.csv   genes removed and added by the correction
  TABLE_consensus123_sensitivity.csv           X-Mens agreement for the released, corrected and sensitivity lists
  IVD_comparisons_with_genes.xlsx              the comparison workbook with the consensus sheets replaced
  SUMMARY_IVD_comparisons.csv                  its summary table as plain values
  METHODS_IVD_comparison_draft.md              the methods draft with the consensus numbers updated
The previous workbook, summary CSV and methods draft are inputs and are not modified.

Usage:
  python build_outputs_stadtmauer_fix.py --job-out out --build-consensus build_consensus.py \\
      --prev-workbook IVD_comparisons_with_genes.xlsx --prev-summary SUMMARY_IVD_comparisons.csv \\
      --prev-methods METHODS_IVD_comparison_draft.md --outdir outputs
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

ARMS = ("MPA+PGE2", "MPA+cAMP")
TP_LABEL = {"6h": "6 h", "48h": "48 h", "14d": "14 d", "3d": "3 d", "8d": "8 d"}
XLAB = {"supported": "same direction", "reversed": "opposite direction",
        "tested_not_called": "tested, not called", "ineligible_or_unavailable": "not testable"}
STUDY_NAME = {"Aghajanova": "Aghajanova", "Rytkonen": "Rytkönen", "Stadtmauer": "Stadtmauer"}
COL_AGH, COL_RYT, COL_STAD = ("Aghajanova 2011 (E2+P4 vs E2)", "Rytkönen 2019 (cAMP+MPA vs untreated)",
                              "Stadtmauer 2022 (MPA+PGE2 or MPA+cAMP vs base)")
FONT, BOLD = Font(name="Calibri", size=11), Font(name="Calibri", size=11, bold=True)


def load_norm(path: Path):
    spec = importlib.util.spec_from_file_location("build_consensus_upstream", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.norm


def timepoint_call(sig: pd.DataFrame, gene: str, order: list[str]) -> str:
    d = sig[sig.g == gene]
    if d.empty:
        return ""
    dirs = sorted(set(d.direction))
    return f"{dirs[0] if len(dirs) == 1 else 'conflicting'} ({', '.join(TP_LABEL[t] for t in order if t in set(d.timepoint))})"


def studies_label(mature_papers: str) -> str:
    return "; ".join(STUDY_NAME[s.strip()] for s in mature_papers.split(";"))


def per_gene_table(job: Path, norm) -> pd.DataFrame:
    best = pd.read_csv(job / "consensus123_stadtmauer_authors_criteria.csv")
    orth = pd.read_csv(job / "consensus123_stadtmauer_authors_criteria_mouse_orthology.csv").set_index("human_gene")
    sup = pd.read_csv(job / "consensus_ortholog_xmens_broad_support_rows_stadtmauer_authors_criteria.csv").set_index("human_gene")
    tab = pd.read_csv(job / "stadtmauer_tables_S2_S3_parsed_and_called.csv", keep_default_na=False).set_index("gene")
    hits = pd.read_csv(job / "stadtmauer2022_hits_authors_criteria.csv")
    text_genes = set(hits[hits.source.str.contains("Results text")].gene)
    vote_dir = job / "work" / "vote_corrected"
    agh = pd.read_csv(vote_dir / "aghajanova2011_invitro_decidualization_nonendo.csv")
    ryt = pd.read_csv(vote_dir / "rytkonen2019_complete_signature.csv")
    agh["g"], ryt["g"] = agh.gene.map(norm), ryt.gene.map(norm)
    rows = []
    for r in best.itertuples(index=False):
        g = r.gene
        if g in tab.index and tab.loc[g, "new_study_call"] != "no call":
            arms = [a for a in ARMS if tab.loc[g, f"paper_call_{a}"]]
            stad = f"{tab.loc[g, 'new_study_call']} (table: {'both treatments' if len(arms) == 2 else arms[0] + ' only'})"
        elif g in text_genes:
            stad = f"{sorted(set(hits[hits.gene == g].direction))[0]} (Results text)"
        else:
            stad = ""
        strict = g in orth.index and str(orth.loc[g, "strict_one_to_one"]) == "True"
        rows.append({"gene": g, "consensus direction": r.consensus_direction, "studies supporting": int(r.n_mature),
                     COL_AGH: timepoint_call(agh, g, ["6h", "48h", "14d"]), COL_RYT: timepoint_call(ryt, g, ["3d", "8d"]), COL_STAD: stad,
                     "mouse ortholog (1:1)": orth.loc[g, "mouse_ortholog"] if strict else "",
                     "X-Mens": XLAB[sup.loc[g, "support_status"]] if g in sup.index else "no 1:1 mouse ortholog",
                     "X-Mens log2 fold change": round(float(sup.loc[g, "log2FoldChange_unshrunk"]), 2)
                     if g in sup.index and pd.notna(sup.loc[g, "log2FoldChange_unshrunk"]) else None})
    out = pd.DataFrame(rows).sort_values(["studies supporting", "gene"], ascending=[False, True]).reset_index(drop=True)
    assert len(out) == len(best) and (out[COL_RYT] != "").all() and (out[[COL_AGH, COL_RYT, COL_STAD]] != "").sum(axis=1).ge(2).all()
    return out


def per_study_table(manifest: dict) -> pd.DataFrame:
    t = manifest["tallies"]["corrected"]
    st = manifest["stadtmauer"]
    no_call = sorted(st.get("table_genes_without_a_call", []))
    base = [
        dict(study="Aghajanova et al. 2011, Biol Reprod 84:801", cells="primary ESF, 4 women without endometriosis", comparison="E2 + P4 vs E2 alone",
             timepoints="6 h, 48 h, 14 d", platform_and_published_threshold="Affymetrix Human Gene 1.0 ST; two-way ANOVA, BH-adjusted P < 0.05, fold change > 1.5",
             lists_used="Tables 4 and 5, women without endometriosis",
             not_used="Table 3 (endometriosis vs control, untreated); endometriosis arms of Tables 4 and 5", key="Aghajanova"),
        dict(study="Rytkönen et al. 2019, Reprod Sci", cells="immortalized ESF (T HESC, ATCC CRL-4003)", comparison="8-Br-cAMP + MPA vs untreated",
             timepoints="3 d, 8 d", platform_and_published_threshold="RNA-seq; edgeR, FDR < 0.05, |fold change| > 2, TPM > 2",
             lists_used="Supplemental Table S1A-D (3 d and 8 d vs untreated, up and down)", not_used="3 d vs 8 d comparisons (S1F-S1H)", key="Rytkonen"),
        dict(study="Stadtmauer and Wagner 2022, Biol Reprod 106:155", cells="immortalized ESF (ATCC CRL-4003)",
             comparison="MPA + PGE2 and MPA + 8-Br-cAMP vs base medium", timepoints="6 d",
             platform_and_published_threshold="RNA-seq; edgeR, FDR < 1e-6, fold change > 1.5, > 3 TPM",
             lists_used=(f"Supplementary Tables S2-S3: {st['table_genes']} genes with TPM and FDR, a treatment counted only if it meets the authors' "
                         f"criteria vs base medium ({st['table_genes_with_a_call']} genes); directions for {st['text_genes']} genes transcribed from the Results text"),
             not_used=(f"{len(no_call)} tabulated genes meeting the criteria in neither treatment ({', '.join(no_call)}); "
                       "genes shown only in single-cell figures without a stated direction"), key="Stadtmauer"),
    ]
    rows = []
    for b in base:
        k = b.pop("key")
        rows.append({**b, "genes_with_a_call": t[k]["genes_listed"], "consistent_up": t[k]["consistent_up"], "consistent_down": t[k]["consistent_down"],
                     "conflicting_excluded": t[k]["conflicting_excluded"], "genes_in_consensus": t[k]["genes_in_list"]})
    return pd.DataFrame(rows)


def write_sheet(wb, name: str, df: pd.DataFrame, index: int) -> None:
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


def update_workbook(prev: Path, dst: Path, per_gene: pd.DataFrame, job: Path, summary_row: dict, revision_note: str) -> dict:
    wb = load_workbook(prev)
    order = wb.sheetnames
    sup = pd.read_csv(job / "consensus_ortholog_xmens_broad_support_rows_stadtmauer_authors_criteria.csv")
    called = sup[sup.support_status.isin(["supported", "reversed"])].sort_values("human_gene")
    cons_sheet = per_gene[["gene", "consensus direction", "studies supporting"]].rename(
        columns={"consensus direction": "direction in vitro", "studies supporting": "studies agreeing (of 3)"})
    best = pd.read_csv(job / "consensus123_stadtmauer_authors_criteria.csv").set_index("gene")
    cons_sheet["studies"] = cons_sheet.gene.map(lambda g: studies_label(best.loc[g, "mature_papers"]))
    cons_sheet["mouse ortholog (1:1)"] = per_gene["mouse ortholog (1:1)"].values
    cons_sheet = cons_sheet.sort_values(["studies agreeing (of 3)", "gene"], ascending=[False, True])
    xs = pd.DataFrame({"gene": called.human_gene, "mouse gene": called.mouse_gene, "in vitro (human consensus)": called.consensus_direction,
                       "X-Mens": called.observed_direction, "X-Mens log2 fold change": called.log2FoldChange_unshrunk.round(3),
                       "X-Mens adjusted p": called.padj_for_call.map(lambda p: float(f"{p:.3g}")),
                       "same direction": called.support_status.eq("supported").astype(bool)})
    write_sheet(wb, "consensus 123 genes", cons_sheet, order.index("consensus 123 genes"))
    write_sheet(wb, "X-Mens vs consensus in vitro", xs, order.index("X-Mens vs consensus in vitro"))
    ws = wb["Summary"]
    assert str(ws["A2"].value).startswith("X-Mens vs human consensus"), ws["A2"].value
    last = len(xs) + 1
    ws["A2"] = summary_row["label"]
    ws["B2"] = f"=COUNTA('X-Mens vs consensus in vitro'!$A$2:$A${last})"
    ws["C2"] = f"=COUNTIF('X-Mens vs consensus in vitro'!$G$2:$G${last},TRUE)"
    ws["E2"] = summary_row["chance"]
    notes = wb["Notes"]
    notes["A2"] = summary_row["note_consensus"]
    notes["A3"] = summary_row["note_xmens"]
    notes.cell(row=notes.max_row + 1, column=1, value=revision_note).font = FONT
    assert wb.sheetnames == order, (wb.sheetnames, order)
    wb.save(dst)
    # re-read and evaluate the row-2 formulas against the saved sheet
    chk = load_workbook(dst)
    rows = list(chk["X-Mens vs consensus in vitro"].iter_rows(min_row=2, max_row=last, values_only=True))
    n, same = sum(r[0] is not None for r in rows), sum(r[6] is True for r in rows)
    other_formulas_same = all(load_workbook(prev)["Summary"].cell(r, c).value == chk["Summary"].cell(r, c).value
                              for r in range(3, chk["Summary"].max_row + 1) for c in range(1, 7))
    return {"genes_compared": n, "same_direction": same, "pct_same": round(100 * same / n, 1), "other_summary_rows_unchanged": other_formulas_same}


def methods_text(prev: str, m: dict) -> str:
    """Targeted replacements in the methods draft; each old passage must occur exactly once."""
    reps = m["replacements"]
    for old, new in reps:
        assert prev.count(old) == 1, f"methods passage not found exactly once: {old[:90]!r}"
        prev = prev.replace(old, new)
    return prev


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    for a in ("--job-out", "--build-consensus", "--prev-workbook", "--prev-summary", "--prev-methods", "--outdir"):
        ap.add_argument(a, type=Path, required=True)
    args = ap.parse_args(argv)
    job, outdir = args.job_out, args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = json.load(open(job / "manifest.json"))
    assert all(c["pass"] for c in manifest["checks"].values()), manifest["checks"]
    norm = load_norm(args.build_consensus)
    tab = pd.read_csv(job / "stadtmauer_tables_S2_S3_parsed_and_called.csv", keep_default_na=False)
    manifest["stadtmauer"]["table_genes_without_a_call"] = tab.loc[tab.new_study_call == "no call", "gene"].tolist()

    per_gene = per_gene_table(job, norm)
    per_gene.to_csv(outdir / "TABLE_consensus123_per_gene_calls.csv", index=False)
    per_study = per_study_table(manifest)
    per_study.to_csv(outdir / "TABLE_consensus123_sources_by_study.csv", index=False)
    changes = pd.read_csv(job / "consensus123_changes_vs_released.csv")
    changes.to_csv(outdir / "TABLE_consensus123_changes_vs_released.csv", index=False)
    S = pd.read_csv(job / "summary_released_corrected_and_sensitivity.csv").set_index("variant")
    S.reset_index().to_csv(outdir / "TABLE_consensus123_sensitivity.csv", index=False)
    rel, cor = S.loc["released"], S.loc["corrected"]

    # sensitivity facts for the methods text
    lists = {v: pd.read_csv(job / "sensitivity" / f"consensus_{v}.csv").set_index("gene") for v in
             ("corrected_aghajanova_14d_only", "corrected_stadtmauer_table_genes_only", "corrected_TPM_filter_read_as_mean", "corrected_plus_SOX4_from_text")}
    best_c = pd.read_csv(job / "consensus123_stadtmauer_authors_criteria.csv").set_index("gene")
    a14 = lists["corrected_aghajanova_14d_only"]
    a14_removed = sorted(set(best_c.index) - set(a14.index))
    a14_fewer = sorted(g for g in set(best_c.index) & set(a14.index) if a14.loc[g, "n_mature"] < best_c.loc[g, "n_mature"])
    tpm = lists["corrected_TPM_filter_read_as_mean"]
    tpm_fewer = sorted(g for g in set(best_c.index) & set(tpm.index) if tpm.loc[g, "n_mature"] < best_c.loc[g, "n_mature"])
    tpm_same_list = set(tpm.index) == set(best_c.index)
    tpm_arm_genes = manifest["stadtmauer"]["genes_whose_arm_call_changes_under_mean_TPM_reading"]
    hits = pd.read_csv(job / "stadtmauer2022_hits_authors_criteria.csv")
    text_in = sorted(set(hits[hits.source.str.contains("Results text")].gene) & set(best_c.index))
    assert all("Rytkonen" in best_c.loc[g, "mature_papers"] and best_c.loc[g, "n_mature"] == 2 for g in text_in)
    sox = pd.read_csv(job / "sensitivity" / "xmens_support_rows_corrected_plus_SOX4_from_text.csv").set_index("human_gene")
    sox_status = XLAB.get(sox.loc["SOX4", "support_status"], "not scored") if "SOX4" in sox.index else "no 1:1 mouse ortholog"
    sox_lfc = float(sox.loc["SOX4", "log2FoldChange_unshrunk"]) if "SOX4" in sox.index else float("nan")
    sox_fdr = float(tab.set_index("gene").loc["SOX4", "fdr_MPA+cAMP"])
    T = manifest["tallies"]["corrected"]
    combos = T["vote"]["by_study_combination"]
    n_ag_ry, n_ry_st, n_all = combos.get("Aghajanova; Rytkonen", 0), combos.get("Rytkonen; Stadtmauer", 0), combos.get("Aghajanova; Rytkonen; Stadtmauer", 0)
    three = sorted(best_c.index[best_c.n_mature == 3])
    n_orth = int(cor.strict_1to1_orthologs)
    sens_tab = S.loc["corrected_stadtmauer_table_genes_only"]
    assert int(sens_tab.consensus_genes) == len(best_c) - len(text_in)

    prev_methods = open(args.prev_methods).read()
    fmt = lambda x: f"{x:.1f}"
    words = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five'}
    prose = {'opposite direction': 'changes in the opposite direction', 'same direction': 'changes in the same direction',
             'tested, not called': 'is tested but not called', 'not testable': 'is not testable'}
    sci = lambda x: f"{x:.1e}".replace('e-0', ' × 10^-').replace('e-', ' × 10^-')
    st = manifest["stadtmauer"]
    no_call = sorted(st["table_genes_without_a_call"])
    reps = [
        ("The supplementary tables give values for 29 genes. Most other gene-level results are reported only in the Results text, which names the genes in each k-means cluster of the bulk data and states the cluster's direction under decidualizing treatment; we transcribed directions for 37 such genes and checked each against the text. For one tabulated gene, KLF4, the two treatment arms disagreed and the direction was set to up by judgment [state basis]. Genes shown only in single-cell figures, without a stated direction, were not used. In total, 67 genes.",
         f"Supplementary Tables S2 and S3 give mean expression (TPM) under each treatment and the FDR of each treatment versus base medium for {st['table_genes']} genes. For these genes, a treatment was counted as up- or downregulating a gene only if it met the authors' own criteria (FDR < 10^-6, fold change > 1.5 and expression > 3 TPM); {st['table_genes_with_a_call']} genes met them under at least one treatment, and the other {len(no_call)} ({', '.join(no_call)}) received no call. Most other gene-level results are reported only in the Results text, which names the genes in each k-means cluster of the bulk data and states the cluster's direction under decidualizing treatment; we transcribed directions for {st['text_genes']} such genes and checked each against the text. Genes shown only in single-cell figures, without a stated direction, were not used. In total, {T['Stadtmauer']['genes_listed']} genes."),
        ("genes with conflicting calls were excluded from that study's vote (none in Aghajanova et al., 45 in Rytkönen et al., 3 in Stadtmauer and Wagner).",
         f"genes with conflicting calls were excluded from that study's vote ({T['Rytkonen']['conflicting_excluded']} in Rytkönen et al., {T['Aghajanova']['conflicting_excluded'] or 'none'} in Aghajanova et al. and {T['Stadtmauer']['conflicting_excluded'] or 'none'} in Stadtmauer and Wagner)."),
        ("Of 147 genes called by two or more studies, 123 agreed (93 up, 30 down) and 24 had conflicting directions and were excluded. The 123 comprise 77 genes supported by Aghajanova and Rytkönen, 39 by Rytkönen and Stadtmauer, and 7 by all three (ABI3BP, FOXO1, HAND2, IGFBP1, OSR2, PRL, ZBTB16).",
         f"Of {T['vote']['genes_with_2plus_clean_calls']} genes called by two or more studies, {T['vote']['agreeing']} agreed ({int(cor.up)} up, {int(cor.down)} down) and {T['vote']['conflicting']} had conflicting directions and were excluded. The {T['vote']['agreeing']} comprise {n_ag_ry} genes supported by Aghajanova and Rytkönen, {n_ry_st} by Rytkönen and Stadtmauer, and {n_all} by all three ({', '.join(three)})."),
        ("114 of the 123 consensus genes had a one-to-one mouse ortholog.", f"{n_orth} of the {len(best_c)} consensus genes had a one-to-one mouse ortholog."),
        ("49 of the 114 orthologs were called in X-Mens: 19 concordant and 30 discordant; 54 were tested but not called and 11 were not testable",
         f"{int(cor.xmens_called)} of the {n_orth} orthologs were called in X-Mens: {int(cor.same_direction)} concordant and {int(cor.opposite_direction)} discordant ({fmt(cor.pct_same_direction)}% concordant; {fmt(cor.pct_expected_by_chance)}% expected by chance); {int(cor.tested_not_called)} were tested but not called and {int(cor.not_testable)} were not testable"),
        ("1. **Aghajanova timepoints.** The 6-h and 48-h lists are early progesterone responses rather than decidualized cells. Restricting Aghajanova et al. to 14 days removed one gene (NPR3, listed only at 48 h) and left OSR2 supported by two studies instead of three. Neither change affected the X-Mens comparison: NPR3 was not called in X-Mens, and OSR2 remains in the list.",
         f"1. **Aghajanova timepoints.** The 6-h and 48-h lists are early progesterone responses rather than decidualized cells. Restricting Aghajanova et al. to 14 days removed {words.get(len(a14_removed), len(a14_removed))} gene{'s' if len(a14_removed) != 1 else ''} ({', '.join(a14_removed)}, listed only at 48 h) and left {', '.join(a14_fewer)} supported by two studies instead of three. " + (f"The X-Mens comparison was unchanged ({int(cor.same_direction)} of {int(cor.xmens_called)} called genes concordant)." if (int(S.loc['corrected_aghajanova_14d_only'].same_direction), int(S.loc['corrected_aghajanova_14d_only'].xmens_called)) == (int(cor.same_direction), int(cor.xmens_called)) else f"X-Mens concordance was {int(S.loc['corrected_aghajanova_14d_only'].same_direction)} of {int(S.loc['corrected_aghajanova_14d_only'].xmens_called)} genes ({fmt(S.loc['corrected_aghajanova_14d_only'].pct_same_direction)}%), against {int(cor.same_direction)} of {int(cor.xmens_called)} for the full list.")),
        ("2. **Transcribed Stadtmauer directions.** Twenty-two consensus genes depend on Stadtmauer and Wagner directions taken from the Results text or set by judgment (KLF4); all are also supported by Rytkönen et al. Excluding them left 101 genes, for which X-Mens called 36 orthologs: 11 (31%) changed in the same direction (39% expected by chance), compared with 19 of 49 (39%; 45% expected) for the full list.",
         f"2. **Transcribed Stadtmauer directions.** Of the consensus genes, {len(text_in)} depend on Stadtmauer and Wagner directions taken from the Results text; each is otherwise supported only by Rytkönen et al. Rebuilding the list from the tabulated Stadtmauer genes alone gave {int(sens_tab.consensus_genes)} genes, for which X-Mens called {int(sens_tab.xmens_called)} orthologs: {int(sens_tab.same_direction)} ({fmt(sens_tab.pct_same_direction)}%) changed in the same direction ({fmt(sens_tab.pct_expected_by_chance)}% expected by chance), compared with {int(cor.same_direction)} of {int(cor.xmens_called)} ({fmt(cor.pct_same_direction)}%; {fmt(cor.pct_expected_by_chance)}% expected) for the full list.\n"
         f"3. **Reading of the Stadtmauer expression filter.** We read \"expressed at greater than 3 TPM\" as above 3 TPM in either condition compared. Reading it as the mean of the two conditions removed one or both treatment calls for {words.get(len(tpm_arm_genes), len(tpm_arm_genes))} genes ({', '.join(tpm_arm_genes)}); the list was {'unchanged' if tpm_same_list else 'changed'}, and {', '.join(tpm_fewer) if tpm_fewer else 'no gene'} {'was' if len(tpm_fewer) == 1 else 'were'} then supported by two studies instead of three.\n"
         f"4. **SOX4.** The Stadtmauer and Wagner text describes SOX4 as elevated under MPA + cAMP, but its tabulated FDR ({sci(sox_fdr)}) does not meet the authors' threshold, so it was not counted. Counting it would add SOX4 (also upregulated in Rytkönen et al.) to the list; in X-Mens it {prose.get(sox_status, sox_status)} (log2 fold change {sox_lfc:.2f})."),
        ("3. **Full differentially expressed gene sets.**", "5. **Full differentially expressed gene sets.**"),
        (' Lucas et al. 2020 [ref] was used for annotation only, because its contrast captures the onset of decidualization rather than the decidualized state measured by the other three.', ""),   # Lucas 2020 is carried only as an annotation and affects no result, so the paper does not mention it
        ("the consensus-building repository (methods.md, build_consensus.py, derive_signatures.py, scripts/aghajanova_de.R, scripts/rytkonen_de.py) and the three source papers.",
         "the consensus-building repository (methods.md, build_consensus.py, derive_signatures.py, scripts/aghajanova_de.R, scripts/rytkonen_de.py), the three source papers, and the corrected consensus (rebuild_consensus123_stadtmauer_thresholds.py, 2026-10-22)."),
    ]
    methods = methods_text(prev_methods, {"replacements": reps})
    (outdir / "METHODS_IVD_comparison_draft.md").write_text(methods)

    label = f"X-Mens vs human consensus in vitro ({len(best_c)} genes, {n_orth} with mouse ortholog)"
    summary_row = {"label": label, "chance": float(cor.pct_expected_by_chance),
                   "note_consensus": (f"Human consensus in vitro list: {len(best_c)} genes changed in the same direction in at least 2 of 3 human in vitro "
                                      f"decidualization studies (Aghajanova, Rytkönen, Stadtmauer); {n_orth} have a one-to-one mouse ortholog. "
                                      "Full list on the \"consensus 123 genes\" sheet."),
                   "note_xmens": (f"X-Mens vs human consensus: of the {n_orth} orthologs, {int(cor.xmens_called)} are called by X-Mens (adjusted p < 0.05, "
                                  "fold change > 2x); % same direction = share changed the same way as in vitro.")}
    note = ("Consensus list revised 2026-10-22: Stadtmauer genes with tabulated values are counted only where a treatment meets the authors' criteria "
            "(FDR < 1e-6, fold change > 1.5, > 3 TPM vs base medium). This removed KLF4 and TSHZ1 and added ADAMTS5 and CLU; X-Mens changed from "
            f"{int(rel.same_direction)} of {int(rel.xmens_called)} to {int(cor.same_direction)} of {int(cor.xmens_called)} genes in the same direction.")
    wb_check = update_workbook(args.prev_workbook, outdir / "IVD_comparisons_with_genes.xlsx", per_gene, job, summary_row, note)
    assert (wb_check["genes_compared"], wb_check["same_direction"]) == (int(cor.xmens_called), int(cor.same_direction)) and wb_check["other_summary_rows_unchanged"], wb_check
    summ = pd.read_csv(args.prev_summary)
    i = summ.index[summ.comparison.str.startswith("X-Mens vs human consensus")]
    assert len(i) == 1
    summ.loc[i, ["comparison", "genes compared", "same direction", "% same direction", "% expected by chance"]] = [
        label, wb_check["genes_compared"], wb_check["same_direction"], wb_check["pct_same"], float(cor.pct_expected_by_chance)]
    summ.to_csv(outdir / "SUMMARY_IVD_comparisons.csv", index=False)
    print(json.dumps({"workbook_row2": wb_check, "released": rel.to_dict(), "corrected": cor.to_dict(),
                      "aghajanova_14d": {"removed": a14_removed, "fewer_studies": a14_fewer}, "text_genes_in_list": len(text_in),
                      "tpm_mean_reading": {"arm_changes": tpm_arm_genes, "list_unchanged": tpm_same_list, "fewer_studies": tpm_fewer},
                      "SOX4_if_counted": {"xmens": sox_status, "lfc": round(sox_lfc, 2), "fdr_cAMP": sox_fdr}}, indent=1, default=str))


if __name__ == "__main__":
    main()
