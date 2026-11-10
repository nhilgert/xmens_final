#!/usr/bin/env python
"""Supplementary Tables S9 and S10 of the in vitro consensus, with the X-Mens columns from the current fit.

Wraps the two table functions of build_outputs_stadtmauer_fix.py (in this folder, imported and not edited):
  per_study_table  -> TABLE_consensus123_sources_by_study.csv  (S9: what each study contributed)
  per_gene_table   -> TABLE_consensus123_per_gene_calls.csv    (S10: each gene's calls, ortholog and X-Mens result)
per_gene_table reads the X-Mens support rows from the consensus job folder under the 22 Oct file name; this script
copies the support rows of the current X-Mens fit (compare_xmens_with_references.py) there first, so S10's last two
columns describe the fit of record instead of the 19 Oct fit the 22 Oct job scored.

Usage
  python make_consensus_tables.py --consensus-dir CONSENSUS_OUT --xmens-support consensus_ortholog_xmens_support_rows.csv
         --build-consensus build_consensus.py --out OUT_DIR
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
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
    ap.add_argument("--consensus-dir", required=True, type=Path, help="output folder of build_consensus123.py")
    ap.add_argument("--xmens-support", required=True, type=Path, help="consensus_ortholog_xmens_support_rows.csv of the current X-Mens fit")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    job = args.consensus_dir
    BO = load_module("build_outputs_stadtmauer_fix", HERE / "build_outputs_stadtmauer_fix.py")
    norm = BO.load_norm(HERE / "build_consensus.py")
    manifest = json.load(open(job / "manifest.json"))
    if not all(c["pass"] for c in manifest["checks"].values()):
        raise SystemExit(f"consensus job checks did not all pass: {manifest['checks']}")

    # the X-Mens support rows of the current fit, under the file name per_gene_table expects
    support = pd.read_csv(args.xmens_support)
    consensus = pd.read_csv(job / "consensus123_stadtmauer_authors_criteria.csv")
    orth = pd.read_csv(job / "consensus123_stadtmauer_authors_criteria_mouse_orthology.csv", dtype=str, keep_default_na=False)
    strict = set(orth.loc[orth.strict_one_to_one == "True", "human_gene"])
    if set(support.human_gene) != strict:
        raise SystemExit("the X-Mens support rows do not cover exactly the strict one-to-one orthologs of this consensus list")
    if not set(support.human_gene) <= set(consensus.gene):
        raise SystemExit("the X-Mens support rows name genes outside the consensus list")
    shutil.copy(args.xmens_support, job / "consensus_ortholog_xmens_broad_support_rows_stadtmauer_authors_criteria.csv")

    tab = pd.read_csv(job / "stadtmauer_tables_S2_S3_parsed_and_called.csv", keep_default_na=False)
    manifest["stadtmauer"]["table_genes_without_a_call"] = tab.loc[tab.new_study_call == "no call", "gene"].tolist()
    args.out.mkdir(parents=True, exist_ok=True)
    per_gene = BO.per_gene_table(job, norm)
    per_gene.to_csv(args.out / "TABLE_consensus123_per_gene_calls.csv", index=False)
    per_study = BO.per_study_table(manifest)
    per_study.to_csv(args.out / "TABLE_consensus123_sources_by_study.csv", index=False)
    print(per_study[["study", "genes_with_a_call", "consistent_up", "consistent_down", "conflicting_excluded", "genes_in_consensus"]].to_string(index=False))
    print("S10 X-Mens column:", per_gene["X-Mens"].value_counts().to_dict(), "| source_id of the support rows:", support.source_id.unique().tolist())


if __name__ == "__main__":
    main()
