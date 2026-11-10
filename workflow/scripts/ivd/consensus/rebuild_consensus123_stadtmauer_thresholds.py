#!/usr/bin/env python3
"""Rebuild the 123-gene human in vitro decidualization consensus with the Stadtmauer & Wagner
(2022) table genes called by the authors' own differential-expression criteria, and re-score
X-Mens against the corrected list.

Why
---
derive_signatures.py (human_ivd_repro_final_defense_20260812, lines 98-99, 181 and 186) gave
every gene in Stadtmauer Supplementary Tables S2 and S3 a direction for each decidualizing
treatment from the sign of (treatment TPM - max(base TPM, 3)). The FDR and fold change printed
in the same tables were never used, so:
  * KLF4 and TSHZ1, which do not change significantly under either treatment, received a
    Stadtmauer vote. KLF4's two arms disagreed only because of the 3-TPM floor, and that
    disagreement was then set to "up" by hand (papers/stadtmauer/stadtmauer_manual_overlay.csv,
    line 12, "boundary FC~1 curation call").
  * ADAMTS5 and CLU, significantly up with MPA + cAMP and unchanged with MPA + PGE2, were
    counted as conflicting (up in one arm, "down" in the other) and excluded.

What changes
------------
A tabulated gene receives a direction for a treatment arm only if that arm meets the paper's
criteria versus base medium (Stadtmauer & Wagner 2022, Biol Reprod 106:155, p. 159):
FDR < 1e-6, average fold change > 1.5 (treatment/base for up, base/treatment for down) and
expression > 3 TPM (implemented as: in at least one of the two conditions compared; the
alternative reading, mean of the two > 3 TPM, is checked and reported). Genes reported only in
the Results text keep their transcribed directions (the other 37 rows of the overlay file); the
KLF4 judgment row is dropped because KLF4 has table values.

Nothing else changes. The Aghajanova and Rytkonen signatures are used as they are; the vote
(build_consensus.build), the orthology derivation (derive_consensus123_orthology.py) and the
X-Mens scoring (run_pipeline.build_consensus_ortholog_support) are imported or run unchanged
from the two repositories, which are only read.

Checks (the script stops if any fails)
--------------------------------------
  1. The table parser here, applied with the ORIGINAL rule, reproduces every table-derived row
     of signatures/stadtmauer2022_hits.csv (and KLF4's pre-override conflict).
  2. The unchanged vote on the original signatures reproduces resources/consensus123.csv.
  3. The unchanged orthology script on the released list reproduces
     resources/consensus123_mouse_orthology.csv byte for byte. (The script accepts only 123-gene
     lists; sensitivity lists of other lengths are run in padded 123-gene batches, see
     derive_orthology.)
  4. The unchanged scoring on the released list reproduces the verified release's
     consensus_ortholog_xmens_broad_support_rows.csv.

Usage
-----
  python rebuild_consensus123_stadtmauer_thresholds.py \\
      --ivd-repo /n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro_final_defense_20260812 \\
      --pipeline-repo /n/lab_storage/eddy_lab/users/nhilgert/cagri_ivd_end_to_end_20260815 \\
      --xmens xmens_broad_all_genes.csv \\
      --released-support consensus_ortholog_xmens_broad_support_rows.csv \\
      --out out
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import html
import importlib.util
import io
import json
import math
import platform
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd

# ---- Stadtmauer & Wagner 2022, p. 159: differential-expression criteria ----------------------
STAD_FDR_LT = 1e-6
STAD_FC_GT = 1.5
STAD_TPM_GT = 3.0
# Supplementary Tables S2/S3, 0-based columns: 0 gene, 1 cluster, 2 TPM base, 3 TPM MPA+PGE2,
# 4 TPM MPA+cAMP, 5 TPM MPA, 6 TPM PGE2, 7 FDR base vs MPA+PGE2, 8 FDR base vs MPA+cAMP,
# 9 FDR MPA+PGE2 vs MPA+cAMP. Only the two decidualizing arms vote, as in the original parser.
ARMS = {"MPA+PGE2": (3, 7), "MPA+cAMP": (4, 8)}
HEADER_FIRST_CELL = ("GENE NAME", "GENE")
ORIGINAL_FLOOR = 3.0          # derive_signatures.py line 98; used only to reproduce the old calls
TABLE_SOURCE_ORIGINAL = "Stadtmauer docx (baseline-floored FC)"
TABLE_SOURCE_NEW = "Stadtmauer supplementary table (authors' criteria)"
TEXT_SOURCE = "Stadtmauer Results text"
JUDGMENT_REASON = "boundary"  # the KLF4 overlay row's reason starts with this word

SIGNATURES = {
    "Aghajanova": "aghajanova2011_invitro_decidualization_nonendo.csv",
    "Rytkonen": "rytkonen2019_complete_signature.csv",
    "Stadtmauer": "stadtmauer2022_hits.csv",
    "Lucas": "lucas2020_S2vsS3_decidualization_onset.csv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def repo_paths(ivd_repo: Path, pipeline_repo: Path) -> dict[str, Path]:
    p = {
        "build_consensus": ivd_repo / "build_consensus.py",
        "derive_signatures": ivd_repo / "derive_signatures.py",
        "table_s2": ivd_repo / "papers/stadtmauer/tables2.docx",
        "table_s3": ivd_repo / "papers/stadtmauer/tables3.docx",
        "overlay": ivd_repo / "papers/stadtmauer/stadtmauer_manual_overlay.csv",
        "released_consensus": pipeline_repo / "resources/consensus123.csv",
        "released_orthology": pipeline_repo / "resources/consensus123_mouse_orthology.csv",
        "released_heca_mapping": pipeline_repo / "resources/consensus123_heca_mapping.csv",
        "orthology_script": pipeline_repo / "derive_consensus123_orthology.py",
        "run_pipeline": pipeline_repo / "run_pipeline.py",
        "pipeline_config": pipeline_repo / "config.yaml",
    }
    for label, fn in SIGNATURES.items():
        p[f"sig_{label}"] = ivd_repo / "signatures" / fn
    return p


# ---- Step 1: Stadtmauer tables --------------------------------------------------------------
def docx_rows(path: Path):
    """Yield the text of every table row in a .docx file (all tables, in document order)."""
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8")
    for table_i, table in enumerate(re.findall(r"<w:tbl>.*?</w:tbl>", xml, flags=re.S), start=1):
        for row_i, tr in enumerate(re.findall(r"<w:tr[ >].*?</w:tr>", table, flags=re.S), start=1):
            cells = ["".join(re.findall(r"<w:t[^>]*>([^<]*)</w:t>", tc))
                     for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, flags=re.S)]
            yield table_i, row_i, [html.unescape(c).strip() for c in cells]


def tpm_value(text: str) -> float:
    return float(re.split(r"±|\+/-", text)[0].replace(",", "").strip())


def fdr_value(text: str) -> float:
    return float(text.replace("<", "").replace(",", "").strip())


def parse_stadtmauer_tables(paths: dict[str, Path], norm) -> pd.DataFrame:
    rows = []
    for label in ("table_s2", "table_s3"):
        for table_i, row_i, c in docx_rows(paths[label]):
            if not c or c[0].upper() in HEADER_FIRST_CELL:
                continue
            if len(c) < 10:
                raise SystemExit(f"{paths[label].name} table {table_i} row {row_i}: expected 10 cells, got {c}")
            rec = {"gene": norm(c[0]), "gene_as_printed": c[0], "table": "S2" if label == "table_s2" else "S3",
                   "docx_table": table_i, "docx_row": row_i, "cluster": c[1], "tpm_base": tpm_value(c[2])}
            for arm, (ti, fi) in ARMS.items():
                rec[f"tpm_{arm}"] = tpm_value(c[ti])
                rec[f"fdr_{arm}"] = fdr_value(c[fi])
            rows.append(rec)
    tab = pd.DataFrame(rows)
    if tab.gene.duplicated().any():
        raise SystemExit(f"duplicated genes in the Stadtmauer tables: {tab.gene[tab.gene.duplicated()].tolist()}")
    return tab


def original_arm_call(base: float, treat: float) -> str:
    """derive_signatures.py lines 181 and 186: sign of treatment vs max(base, 3)."""
    return "up" if treat > max(base, ORIGINAL_FLOOR) else "down"


def fold_up(base: float, treat: float) -> float:
    if base > 0:
        return treat / base
    return math.inf if treat > 0 else 1.0


def paper_arm_call(base: float, treat: float, fdr: float, tpm_rule: str = "either") -> tuple[str | None, str]:
    """Authors' criteria for one arm vs base medium. Returns (direction or None, reason)."""
    if not fdr < STAD_FDR_LT:
        return None, f"FDR {fdr:.2g} not < {STAD_FDR_LT:g}"
    expressed = max(base, treat) > STAD_TPM_GT if tpm_rule == "either" else (base + treat) / 2 > STAD_TPM_GT
    if not expressed:
        return None, f"not > {STAD_TPM_GT:g} TPM ({tpm_rule})"
    fc = fold_up(base, treat)
    if fc > STAD_FC_GT:
        return "up", f"{fc:.2f}-fold up, FDR {fdr:.2g}"
    if fc > 0 and 1 / fc > STAD_FC_GT:
        return "down", f"{1 / fc:.2f}-fold down, FDR {fdr:.2g}"
    if fc == 0:
        return "down", f"absent after treatment, FDR {fdr:.2g}"
    return None, f"fold change {fc:.2f} within 1/{STAD_FC_GT:g}-{STAD_FC_GT:g}"


def call_tables(tab: pd.DataFrame) -> pd.DataFrame:
    tab = tab.copy()
    for arm in ARMS:
        tab[f"original_call_{arm}"] = [original_arm_call(b, t) for b, t in zip(tab.tpm_base, tab[f"tpm_{arm}"])]
        calls = [paper_arm_call(b, t, f) for b, t, f in zip(tab.tpm_base, tab[f"tpm_{arm}"], tab[f"fdr_{arm}"])]
        tab[f"fold_change_{arm}"] = [round(fold_up(b, t), 3) for b, t in zip(tab.tpm_base, tab[f"tpm_{arm}"])]
        tab[f"paper_call_{arm}"] = [d or "" for d, _ in calls]
        tab[f"paper_reason_{arm}"] = [r for _, r in calls]
        tab[f"paper_call_meanTPM_{arm}"] = [paper_arm_call(b, t, f, "mean")[0] or ""
                                            for b, t, f in zip(tab.tpm_base, tab[f"tpm_{arm}"], tab[f"fdr_{arm}"])]
    pc = tab[[f"paper_call_{a}" for a in ARMS]]
    tab["new_study_call"] = [("/".join(sorted(set(r) - {""})) or "no call") for r in pc.itertuples(index=False)]
    oc = tab[[f"original_call_{a}" for a in ARMS]]
    tab["original_study_call"] = ["/".join(sorted(set(r))) for r in oc.itertuples(index=False)]
    return tab


def hits_rows_from_tables(tab: pd.DataFrame, which: str) -> list[dict]:
    """which: 'original' (old rule), 'paper' (authors' criteria) or 'paper_meanTPM' (TPM read as a mean)."""
    col = {"original": "original_call_", "paper": "paper_call_", "paper_meanTPM": "paper_call_meanTPM_"}[which]
    rows = []
    for rec in tab.to_dict("records"):
        dirs = sorted({rec[f"{col}{arm}"] for arm in ARMS} - {""})
        src = TABLE_SOURCE_ORIGINAL if which == "original" else TABLE_SOURCE_NEW
        rows += [{"gene": rec["gene"], "direction": d, "source": src} for d in dirs]
    return rows


def text_rows_from_overlay(overlay: pd.DataFrame, norm) -> tuple[list[dict], pd.DataFrame]:
    keep = overlay[~overlay.reason.astype(str).str.startswith(JUDGMENT_REASON)]
    dropped = overlay[overlay.reason.astype(str).str.startswith(JUDGMENT_REASON)]
    rows = []
    for rec in keep.to_dict("records"):
        d = str(rec["direction"]).strip().lower()
        g = norm(rec["gene"])
        rows += ([{"gene": g, "direction": "up", "source": TEXT_SOURCE}, {"gene": g, "direction": "down", "source": TEXT_SOURCE}]
                 if d == "biphasic" else [{"gene": g, "direction": d, "source": TEXT_SOURCE}])
    return rows, dropped


# ---- Step 2: vote ----------------------------------------------------------------------------
def vote(bc, paths: dict[str, Path], workdir: Path, stad_rows: list[dict] | None,
         aghajanova_timepoints: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    workdir.mkdir(parents=True, exist_ok=True)
    for label, fn in SIGNATURES.items():
        shutil.copy(paths[f"sig_{label}"], workdir / fn)
    if stad_rows is not None:
        pd.DataFrame(stad_rows, columns=["gene", "direction", "source"]).to_csv(workdir / SIGNATURES["Stadtmauer"], index=False)
    if aghajanova_timepoints is not None:
        a = pd.read_csv(workdir / SIGNATURES["Aghajanova"])
        a[a.timepoint.isin(aghajanova_timepoints)].to_csv(workdir / SIGNATURES["Aghajanova"], index=False)
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        best, full = bc.build(str(workdir))
    # build_consensus.build orders genes by tier; the released resource is sorted by gene, which is
    # also the row order the orthology script and the scoring use. Sort to match before writing.
    best = best.sort_values("gene", kind="stable").reset_index(drop=True)
    best.to_csv(workdir / "best_hits.csv", index=False)
    full.to_csv(workdir / "consensus_full_table.csv", index=False)
    return pd.read_csv(workdir / "best_hits.csv"), pd.read_csv(workdir / "consensus_full_table.csv"), log.getvalue()


# ---- Steps 3-4: orthology and X-Mens scoring ----------------------------------------------------
ORTHOLOGY_BATCH = 123   # derive_consensus123_orthology.py stops on any list that is not exactly 123 genes (its line 69)


def run_orthology_script(paths: dict[str, Path], consensus: pd.DataFrame, workdir: Path) -> tuple[Path, str]:
    """Run derive_consensus123_orthology.py unchanged in a scratch tree (its REPO is its own folder)."""
    (workdir / "resources").mkdir(parents=True, exist_ok=True)
    shutil.copy(paths["orthology_script"], workdir / paths["orthology_script"].name)
    shutil.copy(paths["released_heca_mapping"], workdir / "resources/consensus123_heca_mapping.csv")
    consensus.to_csv(workdir / "resources/consensus123.csv", index=False)
    proc = subprocess.run([sys.executable, paths["orthology_script"].name], cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"orthology derivation failed in {workdir}:\n{proc.stdout}\n{proc.stderr}")
    out_csv = workdir / "resources/consensus123_mouse_orthology.csv"
    if not out_csv.exists():
        raise SystemExit(f"orthology script exited 0 but wrote no table in {workdir}:\n{proc.stdout}\n{proc.stderr}")
    return out_csv, proc.stdout + proc.stderr


def derive_orthology(paths: dict[str, Path], consensus: pd.DataFrame, released: pd.DataFrame, workdir: Path) -> tuple[Path, str]:
    """Orthology rows for a gene list, from the unchanged script.

    The script resolves each gene independently of the others but accepts only 123-gene lists, so a
    list of another length is run in 123-gene batches padded with genes from the released list, and
    only the requested genes' rows are kept. Padding rows are checked against the released table.
    """
    if len(consensus) == ORTHOLOGY_BATCH:
        return run_orthology_script(paths, consensus, workdir)
    released_orth = pd.read_csv(paths["released_orthology"], dtype=str, keep_default_na=False).set_index("human_gene")
    parts, logs = [], []
    for i in range(0, len(consensus), ORTHOLOGY_BATCH):
        assigned = consensus.iloc[i:i + ORTHOLOGY_BATCH]
        pad = released[~released.gene.isin(assigned.gene)].head(ORTHOLOGY_BATCH - len(assigned))
        batch = pd.concat([assigned, pad]).sort_values("gene", kind="stable").reset_index(drop=True)
        assert len(batch) == ORTHOLOGY_BATCH and batch.gene.is_unique
        out_csv, log = run_orthology_script(paths, batch, workdir / f"batch{i // ORTHOLOGY_BATCH + 1}")
        o = pd.read_csv(out_csv, dtype=str, keep_default_na=False)
        pad_rows = o[o.human_gene.isin(pad.gene)].set_index("human_gene")
        if not pad_rows.equals(released_orth.loc[pad_rows.index, pad_rows.columns]):
            raise SystemExit(f"padding genes resolved differently from the released orthology table in {workdir}")
        parts.append(o[o.human_gene.isin(assigned.gene)])
        logs.append(log)
    orth = pd.concat(parts).set_index("human_gene").loc[consensus.gene.tolist()].reset_index()
    out = workdir / "consensus_mouse_orthology.csv"
    orth.to_csv(out, index=False)
    return out, "\n".join(logs)


def score(RP, consensus_csv: Path, orthology_csv: Path, xmens: pd.DataFrame, calls: dict) -> pd.DataFrame:
    orth = pd.read_csv(orthology_csv)
    n_strict = int(orth.strict_one_to_one.map(RP.bool_value).sum())
    return RP.build_consensus_ortholog_support(consensus_csv, orthology_csv, xmens, calls, n_strict)


def summarize(name: str, best: pd.DataFrame, support: pd.DataFrame) -> dict:
    called = support[support.support_status.isin(["supported", "reversed"])]
    p_cons = (called.consensus_direction == "up").mean() if len(called) else float("nan")
    p_x = (called.log2FoldChange_unshrunk > 0).mean() if len(called) else float("nan")
    vc = support.support_status.value_counts()
    return {"variant": name, "consensus_genes": len(best), "up": int((best.consensus_direction == "up").sum()),
            "down": int((best.consensus_direction == "down").sum()), "strict_1to1_orthologs": len(support),
            "xmens_called": len(called), "same_direction": int(vc.get("supported", 0)),
            "opposite_direction": int(vc.get("reversed", 0)), "tested_not_called": int(vc.get("tested_not_called", 0)),
            "not_testable": int(vc.get("ineligible_or_unavailable", 0)),
            "pct_same_direction": round(100 * vc.get("supported", 0) / len(called), 1) if len(called) else None,
            "pct_expected_by_chance": round(100 * (p_cons * p_x + (1 - p_cons) * (1 - p_x)), 1) if len(called) else None}


def study_tallies(bc, hits_dir: Path, best: pd.DataFrame, full: pd.DataFrame) -> dict:
    in_list = set(best.gene)
    out = {}
    for label in ("Aghajanova", "Rytkonen", "Stadtmauer"):
        v = pd.Series(bc.collapse_paper(str(hits_dir / SIGNATURES[label])))
        out[label] = {"genes_listed": len(v), "consistent_up": int((v == "up").sum()), "consistent_down": int((v == "down").sum()),
                      "conflicting_excluded": int((v == "biphasic").sum()), "genes_in_list": len(in_list & set(v.index))}
    gcol = full.columns[0]
    two = full[full.n_mature_clean >= 2]
    agree = two[(two.n_up == two.n_mature_clean) | (two.n_down == two.n_mature_clean)]
    out["vote"] = {"genes_with_2plus_clean_calls": len(two), "agreeing": len(agree), "conflicting": len(two) - len(agree),
                   "by_study_combination": best.mature_papers.value_counts().to_dict()}
    assert set(agree[gcol]) == in_list, "vote tally does not match the best-hits list"
    return out


def frames_equal(a: pd.DataFrame, b: pd.DataFrame) -> tuple[bool, str]:
    try:
        pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True), check_dtype=False, check_exact=False, rtol=1e-12)
        return True, "identical"
    except AssertionError as exc:
        return False, str(exc)[:600]


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ivd-repo", type=Path, required=True)
    ap.add_argument("--pipeline-repo", type=Path, required=True)
    ap.add_argument("--xmens", type=Path, required=True, help="xmens_broad_all_genes.csv from the verified release")
    ap.add_argument("--released-support", type=Path, required=True, help="consensus_ortholog_xmens_broad_support_rows.csv from the same release")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)

    out, work = args.out, args.out / "work"
    out.mkdir(parents=True, exist_ok=True)
    paths = repo_paths(args.ivd_repo, args.pipeline_repo)
    checks = {}
    manifest = {"inputs": {}, "checks": checks}
    for label, p in list(paths.items()) + [("xmens", args.xmens), ("released_support", args.released_support)]:
        manifest["inputs"][label] = {"path": str(p), "sha256": sha256(p)}

    bc = load_module("build_consensus_upstream", paths["build_consensus"])
    RP = load_module("run_pipeline_upstream", paths["run_pipeline"])
    import yaml
    calls = yaml.safe_load(open(paths["pipeline_config"]))["support_calls"]
    assert float(calls["padj_lt"]) == 0.05 and float(calls["abs_unshrunk_lfc_gt"]) == 1.0, calls
    manifest["environment"] = {"python": platform.python_version(), "pandas": pd.__version__, "support_calls": calls}

    # 1. Stadtmauer tables, reproduced with the original rule, then called with the authors' criteria
    tab = call_tables(parse_stadtmauer_tables(paths, bc.norm))
    original_hits = pd.read_csv(paths["sig_Stadtmauer"])
    original_hits["gene"] = original_hits.gene.map(bc.norm)
    overlay = pd.read_csv(paths["overlay"])
    overlay_genes = set(overlay.gene.map(bc.norm))
    expected = {(r["gene"], r["direction"]) for r in hits_rows_from_tables(tab, "original") if r["gene"] not in overlay_genes}
    observed = {(g, d) for g, d, s in original_hits[["gene", "direction", "source"]].itertuples(index=False) if s == TABLE_SOURCE_ORIGINAL}
    klf4_pre_override = sorted({tab.set_index("gene").loc["KLF4", f"original_call_{a}"] for a in ARMS})
    checks["1_table_parser_reproduces_original_calls"] = {"pass": expected == observed and klf4_pre_override == ["down", "up"],
                                                          "table_genes": len(tab), "rows_compared": len(observed),
                                                          "missing": sorted(observed - expected), "extra": sorted(expected - observed),
                                                          "KLF4_original_rule": klf4_pre_override}
    if not checks["1_table_parser_reproduces_original_calls"]["pass"]:
        raise SystemExit(f"check 1 failed: {checks['1_table_parser_reproduces_original_calls']}")
    text_rows, dropped = text_rows_from_overlay(overlay, bc.norm)
    text_genes = {r["gene"] for r in text_rows}
    both = sorted(text_genes & set(tab.gene))
    assert not both, f"genes with both table values and text-derived calls: {both}"
    assert dropped.gene.map(bc.norm).tolist() == ["KLF4"], dropped
    new_rows = text_rows + hits_rows_from_tables(tab, "paper")
    mean_rule_changes = [g for g, r in tab.set_index("gene").iterrows()
                         if any((r[f"paper_call_{a}"] or "") != (r[f"paper_call_meanTPM_{a}"] or "") for a in ARMS)]
    manifest["stadtmauer"] = {"table_genes": len(tab), "table_genes_with_a_call": int((tab.new_study_call != "no call").sum()),
                              "text_genes": len(text_genes), "dropped_overlay_rows": dropped.to_dict("records"),
                              "genes_whose_arm_call_changes_under_mean_TPM_reading": mean_rule_changes}
    tab.to_csv(out / "stadtmauer_tables_S2_S3_parsed_and_called.csv", index=False)
    pd.DataFrame(new_rows).to_csv(out / "stadtmauer2022_hits_authors_criteria.csv", index=False)

    # 2. Vote: released inputs must reproduce the released list; then the corrected inputs
    best_rel, full_rel, _ = vote(bc, paths, work / "vote_released", None)
    released = pd.read_csv(paths["released_consensus"])
    ok, msg = frames_equal(best_rel, released)
    checks["2_vote_reproduces_released_consensus123"] = {"pass": ok, "detail": msg, "genes": len(best_rel)}
    if not ok:
        raise SystemExit(f"check 2 failed: {msg}")
    variants = {"corrected": dict(stad_rows=new_rows),
                "corrected_aghajanova_14d_only": dict(stad_rows=new_rows, aghajanova_timepoints=["14d"]),
                "corrected_stadtmauer_table_genes_only": dict(stad_rows=hits_rows_from_tables(tab, "paper")),
                "corrected_TPM_filter_read_as_mean": dict(stad_rows=text_rows + hits_rows_from_tables(tab, "paper_meanTPM")),
                "corrected_plus_SOX4_from_text": dict(stad_rows=new_rows + [{"gene": "SOX4", "direction": "up", "source": TEXT_SOURCE}])}
    lists = {"released": (best_rel, full_rel)}
    for name, kw in variants.items():
        b, f, log = vote(bc, paths, work / f"vote_{name}", **kw)
        lists[name] = (b, f)
    manifest["tallies"] = {name: study_tallies(bc, work / f"vote_{name}", *lists[name]) for name in ("released", "corrected")}

    # 3-4. Orthology and scoring, released first (regression), then every variant
    xmens = pd.read_csv(args.xmens)
    summaries, supports = [], {}
    for name, (b, f) in lists.items():
        cons_csv = work / f"vote_{name}" / "best_hits.csv"
        orth_csv, orth_log = derive_orthology(paths, pd.read_csv(cons_csv), released, work / f"orthology_{name}")
        if name == "released":
            same_bytes = orth_csv.read_bytes() == paths["released_orthology"].read_bytes()
            checks["3_orthology_script_reproduces_released_table"] = {"pass": same_bytes, "log": orth_log.strip()[-400:]}
            if not same_bytes:
                raise SystemExit("check 3 failed: orthology table differs from resources/consensus123_mouse_orthology.csv")
        sup = score(RP, cons_csv, orth_csv, xmens, calls)
        RP.write_csv(sup, work / f"support_{name}.csv")
        sup = pd.read_csv(work / f"support_{name}.csv")
        if name == "released":
            ok, msg = frames_equal(sup, pd.read_csv(args.released_support))
            checks["4_scoring_reproduces_released_support_rows"] = {"pass": ok, "detail": msg}
            if not ok:
                raise SystemExit(f"check 4 failed: {msg}")
        supports[name] = (sup, orth_csv)
        summaries.append(summarize(name, b, sup))
    pd.DataFrame(summaries).to_csv(out / "summary_released_corrected_and_sensitivity.csv", index=False)

    # Outputs for the corrected list (new filenames; nothing in either repository is written)
    best_c, full_c = lists["corrected"]
    best_c.to_csv(out / "consensus123_stadtmauer_authors_criteria.csv", index=False)
    full_c.to_csv(out / "consensus_full_vote_table_stadtmauer_authors_criteria.csv", index=False)
    shutil.copy(supports["corrected"][1], out / "consensus123_stadtmauer_authors_criteria_mouse_orthology.csv")
    supports["corrected"][0].to_csv(out / "consensus_ortholog_xmens_broad_support_rows_stadtmauer_authors_criteria.csv", index=False)
    for name in variants:
        if name != "corrected":
            (out / "sensitivity").mkdir(exist_ok=True)
            lists[name][0].to_csv(out / "sensitivity" / f"consensus_{name}.csv", index=False)
            supports[name][0].to_csv(out / "sensitivity" / f"xmens_support_rows_{name}.csv", index=False)
    rel, cor = set(best_rel.gene), set(best_c.gene)
    sup_r = supports["released"][0].set_index("human_gene").support_status
    sup_c = supports["corrected"][0].set_index("human_gene").support_status
    fr, fc = full_rel.set_index(full_rel.columns[0]), full_c.set_index(full_c.columns[0])
    changes = [{"gene": g, "change": "removed" if g in rel else "added",
                "studies_released": fr.mature_papers.get(g, ""), "studies_corrected": fc.mature_papers.get(g, ""),
                "stadtmauer_table_call_released": tab.set_index("gene").original_study_call.get(g, ""),
                "stadtmauer_table_call_corrected": tab.set_index("gene").new_study_call.get(g, ""),
                "xmens_released": sup_r.get(g, ""), "xmens_corrected": sup_c.get(g, "")} for g in sorted(rel ^ cor)]
    pd.DataFrame(changes).to_csv(out / "consensus123_changes_vs_released.csv", index=False)
    manifest["changes"] = changes
    manifest["summary"] = summaries
    json.dump(manifest, open(out / "manifest.json", "w"), indent=2, default=str)
    print(json.dumps({"checks": {k: v["pass"] for k, v in checks.items()}, "changes": changes}, indent=1, default=str))
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
