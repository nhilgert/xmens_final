#!/usr/bin/env python3
"""Run tools/check_annotations.py on every staged object for one species."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


sm = globals()["snakemake"]
species = sm.params.species
repo = Path(sm.params.repo)
checker = repo / "tools/check_annotations.py"
objects = []
objects.extend((Path(path), "reannotated") for path in sm.input.reannotated)
objects.extend([
    (Path(sm.input.refined), "refined"),
    (Path(sm.input.final), "final"),
    (Path(sm.input.fibdec), "fibdec"),
    (Path(sm.input.deposit_all), "all_deposit"),
    (Path(sm.input.deposit_fibdec), "fibdec_deposit"),
])
Path(sm.output.report).parent.mkdir(parents=True, exist_ok=True)

records = []
for path, role in objects:
    fd, tmp_name = tempfile.mkstemp(prefix="xmens_annotation_check_", suffix=".json",
                                    dir=Path(sm.output.report).parent)
    os.close(fd)
    tmp = Path(tmp_name)
    command = [
        sys.executable, str(checker), "--repo", str(repo),
        "--h5ad", str(path), "--species", species, "--role", role,
        "--report", str(tmp), "--nonfatal-counts",
    ]
    result = subprocess.run(command, cwd=repo, text=True, capture_output=True)
    try:
        payload = json.loads(tmp.read_text())
    except Exception as exc:
        payload = {"path": str(path), "role": role, "status": "FAIL",
                   "schema_failures": [f"checker produced no readable report: {exc}"],
                   "stdout": result.stdout, "stderr": result.stderr}
    payload["checker_exit_code"] = result.returncode
    records.append(payload)
    tmp.unlink(missing_ok=True)

states = [record.get("status", "FAIL") for record in records]
overall = "FAIL" if "FAIL" in states else ("SUSPECT" if "SUSPECT" in states else "PASS")
output = {"species": species, "status": overall, "objects": records}
report = Path(sm.output.report)
report.parent.mkdir(parents=True, exist_ok=True)
tmp_report = report.with_name(f".{report.name}.{os.getpid()}.tmp")
tmp_report.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
os.replace(tmp_report, report)
print(json.dumps(output, indent=2, sort_keys=True))
