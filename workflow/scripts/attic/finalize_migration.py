#!/usr/bin/env python3
"""Aggregate staging validation without hiding SUSPECT/FAIL evidence."""
import json
import hashlib
import os
from pathlib import Path


sm = globals()["snakemake"]
reports = [json.loads(Path(path).read_text()) for path in sm.input.validation]
annotation_checks = [json.loads(Path(path).read_text()) for path in sm.input.annotation_checks]
exports = [json.loads(Path(path).read_text()) for path in sm.input.exports]
protected = json.loads(Path(sm.input.protected).read_text())
drafts_path = Path(sm.input.drafts)
drafts = json.loads(drafts_path.read_text())
draft_ok = (drafts.get("schema_version") == "xmens.annotation.v1" and
            drafts.get("status") == "DRAFT" and
            len(drafts.get("notebooks", [])) == 11)
states = ([item.get("status", "FAIL") for item in reports + annotation_checks + exports] +
          [protected.get("status", "FAIL"), "PASS" if draft_ok else "FAIL"])
overall = "FAIL" if "FAIL" in states else ("SUSPECT" if "SUSPECT" in states else "PASS")
payload = {
    "run_id": sm.params.run_id,
    "schema_version": "xmens.annotation.v1",
    "status": overall,
    "validation": reports,
    "annotation_checks": annotation_checks,
    "exports": exports,
    "protected_live_immutability": protected,
    "draft_notebooks": {
        "path": str(drafts_path),
        "status": drafts.get("status"),
        "count": len(drafts.get("notebooks", [])),
        "gate": "PASS" if draft_ok else "FAIL",
        "sha256": hashlib.sha256(drafts_path.read_bytes()).hexdigest(),
    },
    "samap": "SKIPPED (optional; live SAMap outputs untouched)",
}
status_path = Path(sm.output.status)
ready_path = Path(sm.output.ready)
status_path.parent.mkdir(parents=True, exist_ok=True)
tmp_status = status_path.with_name(f".{status_path.name}.{os.getpid()}.tmp")
tmp_ready = ready_path.with_name(f".{ready_path.name}.{os.getpid()}.tmp")
tmp_status.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
tmp_ready.write_text(f"{overall}\n")
os.replace(tmp_status, status_path)
os.replace(tmp_ready, ready_path)
print(json.dumps(payload, indent=2, sort_keys=True))
