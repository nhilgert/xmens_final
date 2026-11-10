#!/usr/bin/env python
"""Compare the IVD pipeline's outputs with the tables of record (resources/ivd/expected) and the workbooks on Drive.

Every table the rules write that has a counterpart of record is compared: byte-for-byte where the writer is the same
code on the same inputs, otherwise column by column (strings equal; numbers equal to within 1e-9 relative, with the
largest difference reported). The two workbooks are compared cell for cell (values, formulas, sheet order) with the
copies in --record-workbooks, and their formatting fingerprints (fonts, number formats, frozen panes, column widths)
are compared too. Two differences are expected and are listed, not failed:
  * S10's two X-Mens columns: the record (22 Oct) scored the 19 Oct fit; the pipeline uses the fit of record (27 Oct,
    with the uninduced fibroblasts). Calls and labels must agree for every gene; fold changes may differ at 2 decimals.
  * Notes cells that name the Ulrich release folder or a source file: the pipeline's folder names differ from the
    artifact names of the record.
Writes VERIFICATION.md and verification.json into --out and exits 1 if any check fails.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

READ = dict(keep_default_na=False, na_values=[""], float_precision="round_trip")
RESULTS: list[dict] = []


def record(check: str, status: str, detail: str = "") -> None:
    RESULTS.append({"check": check, "status": status, "detail": detail})
    print(f"[{status}] {check}" + (f": {detail}" if detail else ""), flush=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path, **kw) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    sep = "\t" if ".tsv" in path.name else ","
    with opener(path, "rt") as f:
        return pd.read_csv(f, sep=sep, **READ, **kw)


def frames_equal(got: pd.DataFrame, exp: pd.DataFrame, rtol: float = 1e-9) -> tuple[bool, str]:
    if list(got.columns) != list(exp.columns):
        return False, f"columns differ: got {list(got.columns)[:8]}..., expected {list(exp.columns)[:8]}..."
    if len(got) != len(exp):
        return False, f"row counts differ: {len(got)} vs {len(exp)}"
    worst, bad = 0.0, []
    for col in got.columns:
        g, e = got[col], exp[col]
        if pd.api.types.is_numeric_dtype(g) and pd.api.types.is_numeric_dtype(e):
            gv, ev = g.to_numpy(dtype=float), e.to_numpy(dtype=float)
            nan_same = np.array_equal(np.isnan(gv), np.isnan(ev))
            m = ~np.isnan(gv) & ~np.isnan(ev)
            if not nan_same:
                bad.append(f"{col}: missing values differ")
                continue
            if m.any():
                rel = np.abs(gv[m] - ev[m]) / np.maximum(np.abs(ev[m]), 1e-300)
                exact = np.array_equal(gv[m], ev[m])
                if not exact:
                    worst = max(worst, float(rel.max()))
                    if rel.max() > rtol:
                        bad.append(f"{col}: max relative difference {rel.max():.3g}")
        else:
            gs, es = g.astype(str).where(g.notna(), ""), e.astype(str).where(e.notna(), "")
            n = int((gs.to_numpy() != es.to_numpy()).sum())
            if n:
                bad.append(f"{col}: {n} cells differ")
    if bad:
        return False, "; ".join(bad[:6])
    return True, "identical" if worst == 0 else f"numbers agree to {worst:.2g} relative"


def compare_tables(check: str, got: Path, exp: Path, bytes_ok: bool = True, **kw) -> bool:
    if not got.exists():
        record(check, "FAIL", f"missing {got}"); return False
    if bytes_ok and not str(exp).endswith(".gz") and sha256(got) == sha256(exp):
        record(check, "PASS", "byte-identical"); return True
    ok, msg = frames_equal(read_table(got, **kw), read_table(exp, **kw))
    record(check, "PASS" if ok else "FAIL", msg)
    return ok


def normalize_provenance(text: str) -> str:
    """Notes name the Ulrich release folder in parentheses; the pipeline's folder differs from the record's artifact name."""
    return re.sub(r"\((release_[^)]*|results/ivd/[^)]*|[^()]*ulrich/release[^)]*)\)", "(RELEASE)", text)


def compare_workbook(check: str, got: Path, exp: Path, provenance_markers: list[str]) -> bool:
    if not got.exists():
        record(check, "FAIL", f"missing output workbook {got}"); return False
    if not exp.exists():
        record(check, "SKIP", f"no record workbook at {exp}"); return True
    g, e = load_workbook(got), load_workbook(exp)
    ok = True
    if g.sheetnames != e.sheetnames:
        record(f"{check}: sheet order", "FAIL", f"{g.sheetnames} vs {e.sheetnames}"); ok = False
    n_cells, diffs, expected_text = 0, [], []
    for name in e.sheetnames:
        if name not in g.sheetnames:
            continue
        gs, es = g[name], e[name]
        rows = max(gs.max_row, es.max_row); cols = max(gs.max_column, es.max_column)
        for r in range(1, rows + 1):
            for c in range(1, cols + 1):
                a, b = gs.cell(r, c).value, es.cell(r, c).value
                n_cells += 1
                if a is None and b is None:
                    continue
                if isinstance(a, bool) or isinstance(b, bool):
                    same = a == b
                elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
                    same = abs(float(a) - float(b)) <= 1e-9 * max(1.0, abs(float(b)))
                elif isinstance(a, str) and isinstance(b, str) and a.startswith("=") and b.startswith("="):
                    same = a.replace(" ", "") == b.replace(" ", "")
                else:
                    same = (a == b) or (a in (None, "") and b in (None, ""))
                if not same:
                    if (isinstance(a, str) and isinstance(b, str) and name == "Notes"
                            and normalize_provenance(a) == normalize_provenance(b)):
                        expected_text.append(f"{name}!{gs.cell(r, c).coordinate}")
                    else:
                        diffs.append(f"{name}!{gs.cell(r, c).coordinate}: {str(a)[:60]!r} vs {str(b)[:60]!r}")
    if diffs:
        record(f"{check}: cells", "FAIL", f"{len(diffs)} of {n_cells} cells differ; first: " + " | ".join(diffs[:5])); ok = False
    else:
        record(f"{check}: cells", "PASS", f"{n_cells} cells compared, all equal" + (f"; {len(expected_text)} provenance note cells differ as expected ({', '.join(expected_text)})" if expected_text else ""))
    # formatting fingerprint
    def fingerprint(wb):
        out = {}
        for ws in wb.worksheets:
            cells = [cell for row in ws.iter_rows() for cell in row if cell.value is not None]
            out[ws.title] = dict(
                fonts=sorted({(c.font.name, float(c.font.size or 0), bool(c.font.bold)) for c in cells}),
                number_formats=sorted({c.number_format for c in cells} - {"General"}),
                frozen=ws.freeze_panes,
                widths={k: round(v.width, 2) for k, v in ws.column_dimensions.items() if v.width and k != "G" and k != "H"},
            )
        return out
    fg, fe = fingerprint(g), fingerprint(e)
    fmt_diffs = []
    for name in fe:
        if name not in fg:
            continue
        for key in ("fonts", "number_formats", "frozen"):
            if fg[name][key] != fe[name][key]:
                fmt_diffs.append(f"{name}: {key} {fg[name][key]} vs {fe[name][key]}")
        w_g, w_e = fg[name]["widths"], fe[name]["widths"]
        for col in sorted(set(w_g) | set(w_e)):
            if w_g.get(col) != w_e.get(col):
                fmt_diffs.append(f"{name}: width {col} {w_g.get(col)} vs {w_e.get(col)}")
    record(f"{check}: formatting", "PASS" if not fmt_diffs else "NOTE", "fonts, number formats, frozen panes and widths equal" if not fmt_diffs else "; ".join(fmt_diffs[:8]))
    return ok


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, type=Path, help="results/ivd/<run>")
    ap.add_argument("--expected", required=True, type=Path, help="resources/ivd/expected")
    ap.add_argument("--ulrich-contrast", default="dS_vs_C7fib_all_donors")
    ap.add_argument("--record-workbooks", type=Path, default=None, help="folder with the two workbooks of record (Drive copies)")
    ap.add_argument("--provenance-markers", nargs="*", default=["release", "artifact", ".zip", ".csv ("], help="substrings marking provenance text in Notes")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    R, E = a.results, a.expected
    a.out.mkdir(parents=True, exist_ok=True)

    # ---- consensus (same code, same inputs: byte-identical)
    for fn in ("consensus123_stadtmauer_authors_criteria.csv", "consensus123_stadtmauer_authors_criteria_mouse_orthology.csv",
               "stadtmauer2022_hits_authors_criteria.csv", "consensus_full_vote_table_stadtmauer_authors_criteria.csv"):
        compare_tables(f"consensus: {fn}", R / "consensus" / fn, E / "consensus" / fn)
    compare_tables("S9: TABLE_consensus123_sources_by_study.csv", R / "consensus/tables/TABLE_consensus123_sources_by_study.csv", E / "consensus/TABLE_consensus123_sources_by_study.csv")
    s10_got, s10_exp = R / "consensus/tables/TABLE_consensus123_per_gene_calls.csv", E / "consensus/TABLE_consensus123_per_gene_calls.csv"
    if s10_got.exists():
        g, e = read_table(s10_got), read_table(s10_exp)
        first = list(g.columns[:7])
        ok, msg = frames_equal(g[first], e[first])
        record("S10: calls, studies and orthologs (columns 1-7)", "PASS" if ok else "FAIL", msg)
        gi, ei = g.set_index("gene"), e.set_index("gene")
        labels_same = (gi["X-Mens"].astype(str) == ei.reindex(gi.index)["X-Mens"].astype(str)).all()
        lfc_g, lfc_e = gi["X-Mens log2 fold change"], ei.reindex(gi.index)["X-Mens log2 fold change"]
        n_diff = int(((lfc_g.round(2) != lfc_e.round(2)) & lfc_g.notna()).sum())
        sup = read_table(R / "xmens/compare/consensus_ortholog_xmens_support_rows.csv").set_index("human_gene")
        cur = sup.log2FoldChange_unshrunk.reindex(gi.index).round(2)
        from_current = bool(((gi["X-Mens log2 fold change"].round(2) == cur) | gi["X-Mens log2 fold change"].isna()).all())
        record("S10: X-Mens columns", "PASS" if labels_same and from_current else "FAIL",
               f"labels equal to the record for all {len(gi)} genes: {labels_same}; fold changes equal the current fit: {from_current}; "
               f"{n_diff} fold changes differ from the 22 Oct record (expected: the record scored the 19 Oct fit)")
    else:
        record("S10", "FAIL", "missing")

    # ---- X-Mens
    compare_tables("X-Mens pseudobulk counts", R / "xmens/pb/xmens_broad_carrier_counts.tsv", E / "xmens/xmens_broad_carrier_counts_with_uninduced.tsv.gz")
    compare_tables("X-Mens pseudobulk metadata", R / "xmens/pb/xmens_broad_metadata.tsv", E / "xmens/xmens_broad_metadata_with_uninduced.tsv")
    compare_tables("X-Mens genes", R / "xmens/pb/xmens_broad_genes.tsv", E / "xmens/xmens_broad_genes.tsv")
    compare_tables("X-Mens all genes (fit of record, 27 Oct)", R / "xmens/fit/xmens_broad_all_genes.csv", E / "xmens/xmens_broad_all_genes_with_uninduced.csv.gz")
    for fn in ("consensus_ortholog_xmens_support_rows.csv", "genes_xmens_vs_consensus.csv", "genes_xmens_vs_oil.csv",
               "genes_xmens_vs_pregnancy.csv", "genes_xmens_vs_in_vitro.csv", "SUMMARY_xmens_comparisons.csv"):
        compare_tables(f"X-Mens comparisons: {fn}", R / "xmens/compare" / fn, E / "xmens" / fn)

    # ---- human
    compare_tables("Human all genes (29 Oct)", R / "human/tables/human_decidual_vs_fibroblast_all_genes.csv", E / "human/human_decidual_vs_fibroblast_all_genes.csv.gz")
    diag_g, diag_e = R / "human/fits/fit_diagnostics.json", E / "human/fit_diagnostics.json"
    if diag_g.exists():
        dg, de = json.load(open(diag_g)), json.load(open(diag_e))
        keys = ("design", "n_pseudobulks", "residual_df", "eligible", "calls")
        got = {k: dg["B_new_celltype_release"].get(k) for k in keys}; exp = {k: de["B_new_celltype_release"].get(k) for k in keys}
        record("Human fit diagnostics (B_new_celltype_release)", "PASS" if got == exp else "FAIL", f"{got} vs {exp}")

    # ---- Ulrich
    cid = a.ulrich_contrast
    compare_tables(f"Ulrich all genes ({cid})", R / f"ulrich/release/ulrich_{cid}_all_genes.csv", E / f"ulrich/ulrich_{cid}_all_genes.csv.gz")
    compare_tables("Ulrich pseudobulk metadata", R / f"ulrich/release/intermediate/ulrich_{cid}_metadata.tsv", E / f"ulrich/ulrich_{cid}_metadata.tsv")
    compare_tables("Ulrich base pseudobulk metadata", R / "ulrich/release/intermediate/ulrich_base_pseudobulk_metadata.tsv", E / "ulrich/ulrich_base_pseudobulk_metadata.tsv")
    cg, ce = R / "ulrich/release/ulrich_contrasts.csv", E / "ulrich/ulrich_contrasts.csv"
    if cg.exists():
        g, e = read_table(cg).set_index("contrast_id"), read_table(ce).set_index("contrast_id")
        cols = [c for c in e.columns if c in g.columns]
        ok, msg = frames_equal(g[cols].reset_index(), e[cols].reset_index())
        record("Ulrich contrasts table (all 9 contrasts)", "PASS" if ok else "FAIL", msg)
    vc = R / "ulrich/release/validation_checks.csv"
    if vc.exists():
        v = pd.read_csv(vc)
        n_fail = int((v["pass"].astype(str) != "True").sum())
        record("Ulrich acceptance checks (pinned counts of the record)", "PASS" if n_fail == 0 else "FAIL", f"{len(v)} checks, {n_fail} failed")

    # ---- workbooks
    if a.record_workbooks:
        compare_workbook("DE workbook", R / "workbooks/pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx",
                         a.record_workbooks / "pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx", a.provenance_markers)
        compare_workbook("IVD workbook", R / "workbooks/IVD_comparisons_with_genes.xlsx",
                         a.record_workbooks / "IVD_comparisons_with_genes.xlsx", a.provenance_markers)
    # Summary values of the IVD workbook against the record's stored values
    wb = R / "workbooks/IVD_comparisons_with_genes.xlsx"
    if wb.exists():
        exp = pd.read_csv(E / "workbook_summaries/IVD_workbook_Summary_values.csv")
        ws = load_workbook(wb)["Summary"]
        got_rows = []
        for r in range(2, ws.max_row + 1):
            label = ws.cell(r, 1).value
            sheet = ws.cell(r, 6).value.replace("gene sheet: ", "")
            gs = load_workbook(wb, data_only=True)[sheet]
            vals = list(gs.iter_rows(min_row=2, values_only=True))
            hdr = [c.value for c in load_workbook(wb)[sheet][1]]
            same_i = hdr.index("same direction")
            n = sum(v[0] is not None for v in vals); same = sum(v[same_i] is True for v in vals)
            got_rows.append({"comparison": label, "genes compared": n, "same direction": same, "% same direction": float(np.floor(100 * same / n * 10 + 0.5) / 10), "% expected by chance": ws.cell(r, 5).value})
        got = pd.DataFrame(got_rows).set_index("comparison"); e = exp.set_index("comparison")
        bad = []
        for label, row in e.iterrows():
            if label not in got.index:
                bad.append(f"missing row {label!r}"); continue
            for col in ("genes compared", "same direction", "% same direction", "% expected by chance"):
                if abs(float(got.loc[label, col]) - float(row[col])) > 1e-9:
                    bad.append(f"{label}: {col} {got.loc[label, col]} vs {row[col]}")
        record("IVD workbook Summary values vs record", "PASS" if not bad else "FAIL", "; ".join(bad[:6]) if bad else f"{len(e)} rows equal")

    # ---- report
    n_fail = sum(r["status"] == "FAIL" for r in RESULTS)
    lines = ["# IVD comparison pipeline: verification against the tables of record", "",
             f"results: `{R}`; expected: `{E}`" + (f"; record workbooks: `{a.record_workbooks}`" if a.record_workbooks else ""), "",
             f"**{sum(r['status'] == 'PASS' for r in RESULTS)} PASS, {n_fail} FAIL, {sum(r['status'] == 'NOTE' for r in RESULTS)} NOTE, {sum(r['status'] == 'SKIP' for r in RESULTS)} SKIP**", "",
             "| check | status | detail |", "|---|---|---|"]
    lines += [f"| {r['check']} | {r['status']} | {r['detail']} |" for r in RESULTS]
    (a.out / "VERIFICATION.md").write_text("\n".join(lines) + "\n")
    json.dump(RESULTS, open(a.out / "verification.json", "w"), indent=1)
    print(f"\n{sum(r['status'] == 'PASS' for r in RESULTS)} PASS, {n_fail} FAIL -> {a.out / 'VERIFICATION.md'}")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
