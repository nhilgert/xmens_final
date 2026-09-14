#!/usr/bin/env python3
"""Prepare and verify backups for the v2 annotation promotion; never promote."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

RUN_ID = "xmens_annotation_v2_20260811"
OLD_RUN_ID = "xmens_annotation_v1_20260810"
PLAN_REL = Path("tools/promotion_plan_xmens_annotation_v1_20260810.json")
PLAN_FILE_SHA256 = "bcfc92d45cd46e944a608ce9c20bedecf020cf39818d507506fb1d5cf92edc27"
MCK = Path("/n/lab_storage/eddy_lab/Lab/mckinley_backup_2025-12-18/data")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fsync_file(path):
    with Path(path).open("rb+") as handle:
        os.fsync(handle.fileno())


def fsync_dir(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def metadata(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"not a regular non-symlink file: {path}")
    stat = path.stat()
    if stat.st_nlink != 1:
        raise ValueError(f"hard-linked file is not allowed: {path}")
    return {
        "sha256": sha256(path), "size": stat.st_size, "mode": f"{stat.st_mode & 0o7777:04o}",
        "uid": stat.st_uid, "gid": stat.st_gid, "mtime_ns": stat.st_mtime_ns,
        "inode": stat.st_ino, "device": stat.st_dev,
    }


def atomic_json(path, payload):
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    fsync_file(tmp)
    os.replace(tmp, path)
    fsync_dir(path.parent)


def require_v2_pass(repo):
    root = repo / f"data/interim/migrated/{RUN_ID}"
    ready = root / "MIGRATION_REVIEW_READY"
    status_path = root / "manifests/MIGRATION_STATUS.json"
    if ready.read_text() != "PASS\n":
        raise ValueError(f"{ready} is not exact PASS")
    status = json.loads(status_path.read_text())
    if status.get("run_id") != RUN_ID or status.get("status") != "PASS":
        raise ValueError(f"v2 migration status is not PASS: {status_path}")
    for section in ("validation", "annotation_checks", "exports"):
        if not status.get(section) or any(item.get("status") != "PASS" for item in status[section]):
            raise ValueError(f"v2 status section is not wholly PASS: {section}")
    protected = status.get("protected_live_immutability", {})
    if protected.get("status") != "PASS" or protected.get("changes") != []:
        raise ValueError("protected-live validation is not exact PASS/changes=[]")
    pins = {}
    for report in status["annotation_checks"]:
        for record in report.get("objects", []):
            path = Path(record["path"])
            if not path.is_absolute():
                path = repo / path
            path = path.resolve()
            if (record.get("status") != "PASS" or record.get("checker_exit_code") != 0
                    or not record.get("sha256") or record.get("size") is None):
                raise ValueError(f"annotation source is not cryptographically PASS: {path}")
            pins[str(path)] = {"sha256": record["sha256"], "size": record["size"]}
    if len(pins) != 21:
        raise ValueError(f"expected 21 unique certified staged sources; found {len(pins)}")
    return root, ready, status_path, pins


def build_rows(repo, plan):
    rows = []
    for mapping in sorted(plan["plan"]["mappings"], key=lambda row: row["order"]):
        source_rel = mapping["source"]["path"].replace(OLD_RUN_ID, RUN_ID)
        source = repo / source_rel
        destination = repo / mapping["destination"]["path"]
        rows.append({
            "physical_order": len(rows) + 1,
            "logical_order": mapping["order"],
            "stage": mapping["stage"],
            "kind": "reviewed_destination",
            "source": str(source),
            "destination": str(destination),
            "expected_old": mapping["destination"]["old"],
        })
        if mapping["stage"] == "deposit":
            rows.append({
                "physical_order": len(rows) + 1,
                "logical_order": mapping["order"],
                "stage": "deposit_mirror",
                "kind": "mckinley_mirror",
                "source": str(source),
                "copy_source_after_promotion": str(destination),
                "destination": str(MCK / destination.name),
                "expected_old": {
                    "sha256": mapping["destination"]["old"]["sha256"],
                    "size": mapping["destination"]["old"]["size"],
                },
            })
    if len(rows) != 25 or [row["physical_order"] for row in rows] != list(range(1, 26)):
        raise ValueError("physical promotion mapping must contain exactly 25 ordered rows")
    if len({row["destination"] for row in rows}) != 25:
        raise ValueError("promotion destinations are not unique")
    return rows


parser = argparse.ArgumentParser()
parser.add_argument("--repo", required=True)
parser.add_argument("--backup-dir", required=True)
parser.add_argument("--execute", action="store_true")
args = parser.parse_args()
repo = Path(args.repo).resolve()
expected_repo = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
if repo != expected_repo:
    raise SystemExit(f"repo must be exactly {expected_repo}")
backup_dir = Path(args.backup_dir)
if not backup_dir.is_absolute():
    backup_dir = repo / backup_dir
backup_dir = backup_dir.resolve(strict=False)
expected_parent = (repo / "data/interim").resolve()
if backup_dir.parent != expected_parent or not backup_dir.name.startswith("_promotion_backup_"):
    raise SystemExit(f"backup directory must be a new direct child of {expected_parent}")
if backup_dir.exists():
    raise FileExistsError(f"refusing to reuse backup directory: {backup_dir}")

root, ready_path, status_path, certified_sources = require_v2_pass(repo)
plan_path = repo / PLAN_REL
if sha256(plan_path) != PLAN_FILE_SHA256:
    raise ValueError("reviewed promotion plan changed")
plan = json.loads(plan_path.read_text())
rows = build_rows(repo, plan)

for row in rows:
    source, destination = Path(row["source"]), Path(row["destination"])
    source_meta = metadata(source)
    certified = certified_sources.get(str(source.resolve()))
    if certified is None:
        raise ValueError(f"staged source is absent from v2 certification: {source}")
    for key in ("sha256", "size"):
        if source_meta[key] != certified[key]:
            raise ValueError(f"staged source differs from v2 certification ({key}): {source}")
    old_meta = metadata(destination)
    expected = row["expected_old"]
    for key in ("sha256", "size"):
        if old_meta[key] != expected[key]:
            raise ValueError(f"live destination drift for {destination}: {key}")
    if row["kind"] == "reviewed_destination":
        for key in ("mode", "uid", "gid", "mtime_ns", "inode"):
            if old_meta[key] != expected[key]:
                raise ValueError(f"reviewed destination metadata drift for {destination}: {key}")
    if not os.access(destination.parent, os.W_OK):
        raise PermissionError(f"destination parent is not writable: {destination.parent}")
    if source_meta["device"] != old_meta["device"]:
        raise ValueError(f"source/destination are not on one filesystem: {source} -> {destination}")
    row["source_observed"] = source_meta
    row["source_certified"] = certified
    row["old_observed"] = old_meta
    row["dry_run"] = "would promote"

required_backup_bytes = sum(row["old_observed"]["size"] for row in rows)
free_bytes = shutil.disk_usage(backup_dir.parent).free
required_with_margin = required_backup_bytes + 5 * 1024 ** 3
if free_bytes < required_with_margin:
    raise OSError(
        f"insufficient free space: free={free_bytes}, required={required_with_margin}")

dry_payload = {
    "status": "DRY_RUN_PASS",
    "run_id": RUN_ID,
    "repo": str(repo),
    "migration_status": str(status_path),
    "migration_status_sha256": sha256(status_path),
    "migration_ready": str(ready_path),
    "migration_ready_sha256": sha256(ready_path),
    "reviewed_plan": str(plan_path),
    "reviewed_plan_sha256": sha256(plan_path),
    "backup_dir": str(backup_dir),
    "backup_bytes_required": required_backup_bytes,
    "free_bytes_observed": free_bytes,
    "required_with_5gib_margin": required_with_margin,
    "rows": rows,
}
print(json.dumps(dry_payload, indent=2, sort_keys=True))
if not args.execute:
    raise SystemExit(0)

backup_dir.mkdir(mode=0o755)
fsync_dir(backup_dir.parent)
try:
    for row in rows:
        destination = Path(row["destination"])
        group = "mckinley" if row["kind"] == "mckinley_mirror" else "repo"
        backup = backup_dir / group / f"{row['physical_order']:02d}_{destination.name}.old"
        backup.parent.mkdir(parents=True, exist_ok=True)
        if backup.exists():
            raise FileExistsError(backup)
        shutil.copyfile(destination, backup)
        shutil.copystat(destination, backup)
        fsync_file(backup)
        fsync_dir(backup.parent)
        observed = metadata(backup)
        if observed["sha256"] != row["old_observed"]["sha256"] or observed["size"] != row["old_observed"]["size"]:
            raise ValueError(f"backup verification failed: {backup}")
        row["backup"] = str(backup)
        row["backup_observed"] = observed
        print(f"BACKUP {row['physical_order']:02d}/25 {destination} -> {backup}", flush=True)

    # Re-hash all three sides after copying so the receipt binds a stable state.
    for row in rows:
        source_now = metadata(row["source"])
        destination_now = metadata(row["destination"])
        backup_now = metadata(row["backup"])
        if source_now["sha256"] != row["source_observed"]["sha256"]:
            raise ValueError(f"staged source changed during backup: {row['source']}")
        if destination_now["sha256"] != row["old_observed"]["sha256"]:
            raise ValueError(f"live destination changed during backup: {row['destination']}")
        if backup_now["sha256"] != row["old_observed"]["sha256"]:
            raise ValueError(f"backup changed after copy: {row['backup']}")
    if sha256(status_path) != dry_payload["migration_status_sha256"]:
        raise ValueError("v2 migration status changed during backup preparation")
    if sha256(ready_path) != dry_payload["migration_ready_sha256"]:
        raise ValueError("v2 readiness marker changed during backup preparation")

    payload = dict(dry_payload)
    payload["status"] = "PREPARED_NOT_PROMOTED"
    payload["rows"] = rows
    manifest_path = backup_dir / "BACKUP_MANIFEST.json"
    atomic_json(manifest_path, payload)
    tsv_path = backup_dir / "BACKUP_MANIFEST.tsv"
    with tsv_path.open("w", newline="") as handle:
        columns = ["physical_order", "logical_order", "stage", "kind", "source",
                   "destination", "backup", "old_sha256", "new_sha256", "dry_run"]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                **{key: row[key] for key in columns if key in row},
                "old_sha256": row["old_observed"]["sha256"],
                "new_sha256": row["source_observed"]["sha256"],
            })
    fsync_file(tsv_path)
    receipt = backup_dir / "PREPARED_NOT_PROMOTED"
    receipt.write_text("PREPARED_NOT_PROMOTED\n")
    fsync_file(receipt)
    fsync_dir(backup_dir)
    print(f"PASS: verified 25 independent backups; no live file was changed.\n{manifest_path}")
except Exception:
    incomplete = backup_dir / "PREPARATION_INCOMPLETE"
    incomplete.write_text("PREPARATION_INCOMPLETE\n")
    fsync_file(incomplete)
    fsync_dir(backup_dir)
    raise
