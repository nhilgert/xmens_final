#!/usr/bin/env python
"""Build the 123-gene human in vitro decidualization consensus and its mouse orthology from the curated inputs.

Pipeline front end (2026-11-02) for the 22 Oct 2026 correction (rebuild_consensus123_stadtmauer_thresholds.py in this
folder, imported for its functions and not edited). Differences from running that script directly:
  * the Stadtmauer and Wagner table values are read from a CSV of the parsed Supplementary Tables S2 and S3
    (resources/ivd/consensus/stadtmauer_tables_S2_S3.csv: gene, table, cluster, TPM and FDR per arm) instead of the
    authors' .docx files, which are third-party supplements and are not tracked; the authors' criteria are applied by the
    same function (paper_arm_call, via call_tables);
  * the orthology is derived by the pipeline copy of derive_consensus123_orthology.py (command-line inputs);
  * the X-Mens scoring and sensitivity lists of the correction job are not repeated here: X-Mens is scored downstream by
    compare_xmens_with_references.py on the current fit.
Regression checks (all must pass):
  1. the table rows called with the ORIGINAL rule (sign of treatment minus max(base, 3 TPM)) reproduce every
     table-derived row of the original signatures/stadtmauer2022_hits.csv, and KLF4's pre-override conflict;
  2. the unchanged vote on the original four signatures reproduces released_consensus123.csv;
  3. the orthology script on the released list reproduces released_consensus123_mouse_orthology.csv byte for byte.

Writes into --out (file names as in the 22 Oct job, so build_outputs_stadtmauer_fix.py's table functions read them):
  consensus123_stadtmauer_authors_criteria.csv            the consensus list (gene, consensus_direction, tier, ...)
  consensus_full_vote_table_stadtmauer_authors_criteria.csv   every gene with a call in any study, with the vote tally
  consensus123_stadtmauer_authors_criteria_mouse_orthology.csv  the orthology table (114 strict one-to-one)
  stadtmauer_tables_S2_S3_parsed_and_called.csv           the 30 tabulated genes with both rules' calls
  stadtmauer2022_hits_authors_criteria.csv                the Stadtmauer study list (table + Results-text rows)
  work/vote_corrected/<signature>.csv                     the four study tables the vote read
  manifest.json                                            inputs (sha256), checks, tallies (S9 numbers)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--signatures", required=True, type=Path, help="folder with the four per-study signature CSVs")
    ap.add_argument("--stadtmauer-table", required=True, type=Path, help="parsed Supplementary Tables S2 and S3 (CSV)")
    ap.add_argument("--overlay", required=True, type=Path, help="stadtmauer_manual_overlay.csv (Results-text directions)")
    ap.add_argument("--released-consensus", required=True, type=Path)
    ap.add_argument("--released-orthology", required=True, type=Path)
    ap.add_argument("--heca-mapping", required=True, type=Path)
    ap.add_argument("--components", required=True, type=Path)
    ap.add_argument("--oma", required=True, type=Path)
    ap.add_argument("--annotation", required=True, type=Path)
    ap.add_argument("--alias", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    RB = load_module("rebuild_stadtmauer", HERE / "rebuild_consensus123_stadtmauer_thresholds.py")
    bc = load_module("build_consensus_upstream", HERE / "build_consensus.py")
    orthology_script = HERE / "derive_consensus123_orthology.py"

    out, work = args.out, args.out / "work"
    out.mkdir(parents=True, exist_ok=True)
    inputs = {
        "signatures_dir": args.signatures, "stadtmauer_table": args.stadtmauer_table, "overlay": args.overlay,
        "released_consensus": args.released_consensus, "released_orthology": args.released_orthology,
        "heca_mapping": args.heca_mapping, "components": args.components, "oma": args.oma,
        "annotation": args.annotation, "alias": args.alias,
        "build_consensus": HERE / "build_consensus.py", "rebuild_script": HERE / "rebuild_consensus123_stadtmauer_thresholds.py",
        "orthology_script": orthology_script,
    }
    sig_paths = {label: args.signatures / fn for label, fn in RB.SIGNATURES.items()}
    for label, p in sig_paths.items():
        inputs[f"sig_{label}"] = p
    manifest = {"inputs": {k: {"path": str(p), "sha256": RB.sha256(p)} for k, p in inputs.items() if p.is_file()},
                "environment": {"python": platform.python_version(), "pandas": pd.__version__}, "checks": {}}
    checks = manifest["checks"]

    # ---- 1. Stadtmauer tables from the parsed CSV, called with both rules ----------------------------------------
    tab = pd.read_csv(args.stadtmauer_table, keep_default_na=False)
    needed = {"gene", "table", "cluster", "tpm_base"} | {f"tpm_{a}" for a in RB.ARMS} | {f"fdr_{a}" for a in RB.ARMS}
    missing = needed - set(tab.columns)
    if missing:
        raise SystemExit(f"{args.stadtmauer_table} lacks columns {sorted(missing)}")
    tab = tab[["gene", "table", "cluster", "tpm_base"] + [c for a in RB.ARMS for c in (f"tpm_{a}", f"fdr_{a}")]].copy()
    tab["gene"] = tab.gene.map(bc.norm)
    if tab.gene.duplicated().any():
        raise SystemExit("duplicated genes in the Stadtmauer table")
    tab = RB.call_tables(tab)
    original_hits = pd.read_csv(sig_paths["Stadtmauer"])
    original_hits["gene"] = original_hits.gene.map(bc.norm)
    overlay = pd.read_csv(args.overlay)
    overlay_genes = set(overlay.gene.map(bc.norm))
    expected = {(r["gene"], r["direction"]) for r in RB.hits_rows_from_tables(tab, "original") if r["gene"] not in overlay_genes}
    observed = {(g, d) for g, d, s in original_hits[["gene", "direction", "source"]].itertuples(index=False) if s == RB.TABLE_SOURCE_ORIGINAL}
    klf4_pre = sorted({tab.set_index("gene").loc["KLF4", f"original_call_{a}"] for a in RB.ARMS}) if "KLF4" in set(tab.gene) else []
    checks["1_table_parser_reproduces_original_calls"] = {
        "pass": expected == observed and (klf4_pre == ["down", "up"] or "KLF4" not in overlay_genes),
        "table_genes": int(len(tab)), "rows_compared": len(observed),
        "missing": sorted(observed - expected), "extra": sorted(expected - observed), "KLF4_original_rule": klf4_pre}
    if not checks["1_table_parser_reproduces_original_calls"]["pass"]:
        raise SystemExit(f"check 1 failed: {checks['1_table_parser_reproduces_original_calls']}")
    text_rows, dropped = RB.text_rows_from_overlay(overlay, bc.norm)
    text_genes = {r["gene"] for r in text_rows}
    both = sorted(text_genes & set(tab.gene))
    if both:
        raise SystemExit(f"genes with both table values and text-derived calls: {both}")
    new_rows = text_rows + RB.hits_rows_from_tables(tab, "paper")
    manifest["stadtmauer"] = {
        "table_genes": int(len(tab)), "table_genes_with_a_call": int((tab.new_study_call != "no call").sum()),
        "table_genes_without_a_call": sorted(tab.loc[tab.new_study_call == "no call", "gene"]),
        "text_genes": len(text_genes), "dropped_overlay_rows": dropped.to_dict("records"),
        "genes_whose_arm_call_changes_under_mean_TPM_reading": [
            g for g, r in tab.set_index("gene").iterrows()
            if any((r[f"paper_call_{a}"] or "") != (r[f"paper_call_meanTPM_{a}"] or "") for a in RB.ARMS)],
    }
    tab.to_csv(out / "stadtmauer_tables_S2_S3_parsed_and_called.csv", index=False)
    pd.DataFrame(new_rows).to_csv(out / "stadtmauer2022_hits_authors_criteria.csv", index=False)

    # ---- 2. Vote: the released inputs must reproduce the released list; then the corrected inputs ----------------
    paths = {f"sig_{label}": p for label, p in sig_paths.items()}
    best_rel, full_rel, _ = RB.vote(bc, paths, work / "vote_released", None)
    released = pd.read_csv(args.released_consensus)
    ok, msg = RB.frames_equal(best_rel, released)
    checks["2_vote_reproduces_released_consensus123"] = {"pass": ok, "detail": msg, "genes": int(len(best_rel))}
    if not ok:
        raise SystemExit(f"check 2 failed: {msg}")
    best_c, full_c, _ = RB.vote(bc, paths, work / "vote_corrected", new_rows)
    manifest["tallies"] = {"released": RB.study_tallies(bc, work / "vote_released", best_rel, full_rel),
                           "corrected": RB.study_tallies(bc, work / "vote_corrected", best_c, full_c)}
    best_c.to_csv(out / "consensus123_stadtmauer_authors_criteria.csv", index=False)
    full_c.to_csv(out / "consensus_full_vote_table_stadtmauer_authors_criteria.csv", index=False)

    # ---- 3. Orthology: the released list must reproduce the released table byte for byte; then the corrected list --
    def derive(consensus_csv: Path, out_csv: Path) -> str:
        cmd = [sys.executable, str(orthology_script), "--consensus", str(consensus_csv), "--heca-mapping", str(args.heca_mapping),
               "--components", str(args.components), "--oma", str(args.oma), "--annotation", str(args.annotation),
               "--alias", str(args.alias), "--out", str(out_csv), "--allow-any-length"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit(f"orthology derivation failed:\n{proc.stdout}\n{proc.stderr}")
        return proc.stdout
    (work / "orthology_released").mkdir(parents=True, exist_ok=True)
    log_rel = derive(args.released_consensus, work / "orthology_released" / "consensus123_mouse_orthology.csv")
    same = (work / "orthology_released" / "consensus123_mouse_orthology.csv").read_bytes() == args.released_orthology.read_bytes()
    checks["3_orthology_script_reproduces_released_table"] = {"pass": bool(same), "log": log_rel.strip()[-400:]}
    if not same:
        raise SystemExit("check 3 failed: orthology table differs from the released resource")
    log_cor = derive(out / "consensus123_stadtmauer_authors_criteria.csv", out / "consensus123_stadtmauer_authors_criteria_mouse_orthology.csv")
    orth = pd.read_csv(out / "consensus123_stadtmauer_authors_criteria_mouse_orthology.csv", dtype=str, keep_default_na=False)
    manifest["orthology"] = {"genes": int(len(orth)), "strict_one_to_one": int((orth.strict_one_to_one == "True").sum()),
                             "renamed": orth.loc[orth.human_gene != orth.current_human_symbol, ["human_gene", "current_human_symbol"]].to_dict("records"),
                             "not_one_to_one": orth.loc[orth.strict_one_to_one != "True", "human_gene"].tolist(), "log": log_cor.strip()[-400:]}

    # ---- changes vs the released list ------------------------------------------------------------------------------
    rel, cor = set(best_rel.gene), set(best_c.gene)
    fr, fc = full_rel.set_index(full_rel.columns[0]), full_c.set_index(full_c.columns[0])
    ti = tab.set_index("gene")
    manifest["changes_vs_released"] = [
        {"gene": g, "change": "removed" if g in rel else "added", "studies_released": fr.mature_papers.get(g, ""),
         "studies_corrected": fc.mature_papers.get(g, ""), "stadtmauer_table_call_released": ti.original_study_call.get(g, ""),
         "stadtmauer_table_call_corrected": ti.new_study_call.get(g, "")} for g in sorted(rel ^ cor)]
    manifest["summary"] = {"consensus_genes": int(len(best_c)), "up": int((best_c.consensus_direction == "up").sum()),
                           "down": int((best_c.consensus_direction == "down").sum()),
                           "strict_1to1_orthologs": manifest["orthology"]["strict_one_to_one"]}
    json.dump(manifest, open(out / "manifest.json", "w"), indent=2, default=str)
    print(json.dumps({"checks": {k: v["pass"] for k, v in checks.items()}, "summary": manifest["summary"],
                      "stadtmauer": {k: manifest["stadtmauer"][k] for k in ("table_genes", "table_genes_with_a_call", "text_genes")},
                      "vote": manifest["tallies"]["corrected"]["vote"], "changes_vs_released": manifest["changes_vs_released"]},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
