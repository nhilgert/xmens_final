#!/usr/bin/env python3
"""Rerun the checks listed in README.md ("Checks run on 29 October 2026 for this snapshot").

1. Workbooks. run_workbooks.sh is run on INPUT_DIR into a temporary folder, and both workbooks are compared with the
   reference copies (the files on Drive): every cell value (formulas as text) and the formatting fingerprint from
   5_workbooks/build_house_style_workbooks.py. The reference files' SHA-256 must also equal the Drive files' values below.
2. X-Mens sheets. 2_xmens/update_ivd_workbook.py is run on version 4 of the IVD workbook with the X-Mens comparison
   tables and the revision note that version 5 carries; the result is compared with version 5 the same way.
3. Hammoud files. The SHA-256 of 4_hammoud/ulrich_ivd_pipeline.py and of the config must equal the values the release
   recorded (run_metadata.json, execution_sha256); those of the files in 4_hammoud/resources/ must equal the values
   pinned in the config; and every row of the release's validation_checks.csv must pass.

Usage
  python verify_bundle.py --inputs INPUT_DIR --reference-dir DIR --ivd-v4 V4.xlsx --ivd-v5 V5.xlsx --xmens-comparisons DIR

  --inputs             the step-5 input folder (see run_workbooks.sh)
  --reference-dir      folder holding the two workbooks to compare against (the files on Drive)
  --ivd-v4, --ivd-v5   IVD_comparisons_with_genes.xlsx, versions 4 and 5 (Claude Science artifact
                       c0e13e9f-9222-48cd-9a8a-e96f51d5ab5b)
  --xmens-comparisons  the X-Mens comparison tables (xmens_comparisons_with_uninduced.zip, unpacked)

Needs pandas, openpyxl and PyYAML. Prints one line per check and exits with status 1 if any check fails.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
# SHA-256 of the two workbooks in xmens_figures/DESeq2_DE_and_IVD_comparisons/ on Drive (29 October 2026)
DRIVE_SHA256 = {
    "pseudobulk_DE_decidual_vs_fibroblast_XMens_and_human.xlsx": "b6584047b5b4615a56fb2dc8c505324e41be94ae142d1ec7fe75a3003772dae5",
    "IVD_comparisons_with_genes.xlsx": "24b3a50898037a809f9166d9ed5383419dc3d22fcf905a9981cb906797ce9ce9",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cell_values(path):
    workbook = load_workbook(path)
    return {(sheet.title, cell.coordinate): cell.value
            for sheet in workbook for row in sheet.iter_rows() for cell in row if cell.value is not None}


def load_fingerprint():
    namespace = {}
    exec((HERE / "5_workbooks" / "build_house_style_workbooks.py").read_text(), namespace)
    return namespace["format_fingerprint"]


def compare(new, reference, fingerprint):
    a, b = cell_values(new), cell_values(reference)
    differing = sorted(key for key in set(a) | set(b) if a.get(key) != b.get(key))
    same_format = fingerprint(str(new)) == fingerprint(str(reference))
    detail = f"{len(a):,} vs {len(b):,} cells, {len(differing)} differ"
    if differing:
        detail += f" (first: {differing[:3]})"
    detail += f"; formatting {'identical' if same_format else 'differs'}"
    return not differing and same_format, detail


def check_workbooks(args, fingerprint, work):
    out = work / "workbooks"
    run = subprocess.run(["bash", str(HERE / "run_workbooks.sh"), str(args.inputs), str(out)],
                         capture_output=True, text=True, env=dict(os.environ, PYTHON=sys.executable))
    if run.returncode != 0:
        return [("workbooks: run_workbooks.sh", False, run.stderr.strip()[-500:])]
    results = []
    for name, drive_sha in DRIVE_SHA256.items():
        reference = Path(args.reference_dir) / name
        ok, detail = compare(out / name, reference, fingerprint)
        is_drive = sha256(reference) == drive_sha
        results.append((f"workbooks: {name}", ok and is_drive,
                        detail + ("" if is_drive else "; the reference is not the 29 October Drive file")))
    return results


def check_xmens_sheets(args, fingerprint, work):
    def notes(path):
        return [row[0] for row in load_workbook(path)["Notes"].iter_rows(values_only=True)]
    notes_v4, notes_v5 = notes(args.ivd_v4), notes(args.ivd_v5)
    if notes_v5[:-1] != notes_v4:
        return [("X-Mens sheets: version 4 -> version 5", False, "version 5's Notes are not version 4's plus one line")]
    out = work / "IVD_v5_rebuilt.xlsx"
    run = subprocess.run([sys.executable, str(HERE / "2_xmens" / "update_ivd_workbook.py"), "--previous", str(args.ivd_v4),
                          "--comparisons", str(args.xmens_comparisons), "--note", notes_v5[-1], "--out", str(out)],
                         capture_output=True, text=True)
    if run.returncode != 0:
        return [("X-Mens sheets: update_ivd_workbook.py", False, run.stderr.strip()[-500:])]
    ok, detail = compare(out, args.ivd_v5, fingerprint)
    return [("X-Mens sheets: version 4 -> version 5", ok, detail)]


def check_hammoud(args):
    results = []
    config_path = HERE / "4_hammoud" / "config_ulrich_dS_vs_C7fib_luteal_follicular.yaml"
    config = yaml.safe_load(config_path.read_text())
    release = Path(args.inputs) / "hammoud_release"
    recorded = json.loads((release / "run_metadata.json").read_text())["execution_sha256"]
    for label, path, expected in (("pipeline script", HERE / "4_hammoud" / "ulrich_ivd_pipeline.py", recorded["pipeline"]),
                                  ("config", config_path, recorded["config"])):
        results.append((f"Hammoud: {label} = SHA-256 recorded by the release", sha256(path) == expected, expected[:12]))
    for entry in config["resources"].values():
        if not str(entry.get("path", "")).startswith("resources/"):
            continue
        path = HERE / "4_hammoud" / entry["path"]
        results.append((f"Hammoud: {entry['path']} = SHA-256 pinned in the config",
                        path.exists() and sha256(path) == entry["sha256"], entry["sha256"][:12]))
    checks = pd.read_csv(release / "validation_checks.csv")
    passed = checks["pass"].astype(str).str.lower().eq("true")
    results.append(("Hammoud: the release's validation checks", bool(passed.all()), f"{int(passed.sum())} of {len(checks)} pass"))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--reference-dir", required=True, type=Path)
    parser.add_argument("--ivd-v4", required=True, type=Path)
    parser.add_argument("--ivd-v5", required=True, type=Path)
    parser.add_argument("--xmens-comparisons", required=True, type=Path)
    args = parser.parse_args()
    fingerprint = load_fingerprint()
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        results = check_workbooks(args, fingerprint, work) + check_xmens_sheets(args, fingerprint, work) + check_hammoud(args)
    width = max(len(name) for name, _, _ in results)
    for name, ok, detail in results:
        print(f"{'PASS' if ok else 'FAIL'}  {name:<{width}}  {detail}")
    failed = sum(not ok for _, ok, _ in results)
    print(f"{len(results) - failed} of {len(results)} checks pass")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
