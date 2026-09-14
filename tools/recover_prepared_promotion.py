#!/usr/bin/env python3
"""Idempotently restore all 25 old destinations after an interrupted v2 promotion."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import sys
from datetime import datetime, timezone

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


def atomic_json(path, payload):
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    if tmp.exists():
        raise FileExistsError(tmp)
    with tmp.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    fsync_dir(path.parent)


def atomic_bytes(path, data):
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    if tmp.exists():
        raise FileExistsError(tmp)
    with tmp.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
    fsync_dir(path.parent)


def append_journal(path, payload):
    path = Path(path)
    first = not path.exists()
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    if first:
        fsync_dir(path.parent)


def check_regular(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_nlink != 1:
        raise ValueError(f"not a regular independent file: {path}")


def old_metadata_matches(path, old):
    stat = Path(path).stat()
    return (
        stat.st_size == old["size"]
        and f"{stat.st_mode & 0o7777:04o}" == old["mode"]
        and stat.st_uid == old["uid"]
        and stat.st_gid == old["gid"]
        and stat.st_mtime_ns == old["mtime_ns"]
    )


def classify(row):
    destination = Path(row["destination"])
    actual = sha256(destination)
    if actual == row["old_observed"]["sha256"]:
        return "old" if old_metadata_matches(destination, row["old_observed"]) else "old_metadata_drift"
    if actual == row["source_observed"]["sha256"]:
        return "promoted"
    return "conflict"


def copy_restore(row, path):
    path = Path(path)
    if path.exists() or path.is_symlink():
        raise FileExistsError(path)
    shutil.copy2(row["backup"], path)
    os.chmod(path, int(row["old_observed"]["mode"], 8))
    os.chown(path, row["old_observed"]["uid"], row["old_observed"]["gid"])
    os.utime(path, ns=(row["old_observed"]["mtime_ns"], row["old_observed"]["mtime_ns"]))
    fsync_file(path)
    fsync_dir(path.parent)
    if (sha256(path) != row["old_observed"]["sha256"]
            or not old_metadata_matches(path, row["old_observed"])):
        raise ValueError(f"restore temp failed verification: {path}")


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def signal_as_exception(signum, _frame):
    raise KeyboardInterrupt(f"received signal {signum}")


parser = argparse.ArgumentParser()
parser.add_argument("--repo", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--expected-manifest-sha256", required=True)
args = parser.parse_args()
repo = Path(args.repo).resolve()
expected_repo = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
if repo != expected_repo or Path.cwd().resolve() != expected_repo:
    raise SystemExit(f"run from exactly {expected_repo}")
manifest_path = Path(args.manifest).resolve()
backup_dir = manifest_path.parent
if (manifest_path.name != "BACKUP_MANIFEST.json"
        or backup_dir.parent != (repo / "data/interim").resolve()
        or not backup_dir.name.startswith("_promotion_backup_")):
    raise ValueError("manifest path is not one exact promotion backup directory")
if sha256(manifest_path) != args.expected_manifest_sha256:
    raise ValueError("prepared manifest hash changed")
manifest = json.loads(manifest_path.read_text())
if manifest.get("status") != "PREPARED_NOT_PROMOTED" or manifest.get("run_id") != RUN_ID:
    raise ValueError("prepared manifest status/run is wrong")
if len(manifest.get("rows", [])) != 25:
    raise ValueError("prepared manifest must contain exactly 25 rows")
recovery = backup_dir / "RECOVERY_REQUIRED"
commit_marker = backup_dir / "PROMOTION_COMMITTED"
if commit_marker.exists():
    raise ValueError("promotion has an authoritative commit marker; automatic rollback is forbidden")
if not recovery.is_file():
    raise FileNotFoundError(f"recovery marker is absent: {recovery}")

plan_path = repo / PLAN_REL
if sha256(plan_path) != PLAN_FILE_SHA256:
    raise ValueError("reviewed promotion plan changed")
plan = json.loads(plan_path.read_text())
expected_rows = []
for mapping in sorted(plan["plan"]["mappings"], key=lambda row: row["order"]):
    source = repo / mapping["source"]["path"].replace(OLD_RUN_ID, RUN_ID)
    destination = repo / mapping["destination"]["path"]
    expected_rows.append({
        "physical_order": len(expected_rows) + 1,
        "source": str(source), "destination": str(destination),
        "expected_old": mapping["destination"]["old"],
    })
    if mapping["stage"] == "deposit":
        expected_rows.append({
            "physical_order": len(expected_rows) + 1,
            "source": str(source), "destination": str(MCK / destination.name),
            "expected_old": {
                "sha256": mapping["destination"]["old"]["sha256"],
                "size": mapping["destination"]["old"]["size"],
            },
        })
if len(expected_rows) != 25:
    raise ValueError("reviewed mapping did not expand to exactly 25 physical rows")
for row, expected in zip(manifest["rows"], expected_rows):
    for key in ("physical_order", "source", "destination"):
        if row.get(key) != expected[key]:
            raise ValueError(f"prepared mapping differs from reviewed mapping: {key}")
    for key, value in expected["expected_old"].items():
        if row.get("old_observed", {}).get(key) != value:
            raise ValueError(f"prepared old-state claim differs from reviewed plan: {key}")
    backup = Path(row["backup"]).resolve()
    if backup_dir not in backup.parents:
        raise ValueError(f"backup escapes prepared directory: {backup}")
    check_regular(backup)
    if sha256(backup) != row["old_observed"]["sha256"]:
        raise ValueError(f"backup hash changed: {backup}")
    check_regular(row["destination"])

states = [classify(row) for row in manifest["rows"]]
if "conflict" in states:
    raise ValueError("at least one destination matches neither reviewed old nor certified promoted SHA")
if not sys.stdin.isatty():
    raise SystemExit("interactive terminal required")
token = f"ROLLBACK:{RUN_ID}:{args.expected_manifest_sha256}"
entered = input(f"Type exactly {token}\n> ")
if entered != token:
    raise SystemExit("ABORT: confirmation did not match")

lock_path = repo / "data/interim/.annotation_promotion.lock"
with lock_path.open("a+") as lock:
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    # The pre-prompt marker check is advisory. Recheck under the shared
    # promotion lock so a just-committed promotion can never be rolled back.
    if commit_marker.exists():
        raise ValueError("promotion committed while recovery waited; rollback is forbidden")
    if not recovery.is_file():
        raise FileNotFoundError("recovery marker disappeared before the lock was acquired")
    if sha256(manifest_path) != args.expected_manifest_sha256:
        raise ValueError("prepared manifest changed before recovery acquired the lock")
    states = [classify(row) for row in manifest["rows"]]
    if "conflict" in states:
        raise ValueError("destination conflict appeared before recovery")
    restore_temps = []
    old_handlers = {}
    try:
        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            old_handlers[sig] = signal.signal(sig, signal_as_exception)
        missing_bytes = 0
        planned_temps = []
        initial_states = {row["destination"]: classify(row) for row in manifest["rows"]}
        for row in manifest["rows"]:
            destination = Path(row["destination"])
            temp = destination.with_name(
                f".{destination.name}.manual_restore.{args.expected_manifest_sha256[:16]}.tmp")
            if initial_states[row["destination"]] == "old":
                if temp.exists() or temp.is_symlink():
                    temp.unlink()
                    fsync_dir(temp.parent)
                continue
            if temp.exists() or temp.is_symlink():
                valid = False
                try:
                    check_regular(temp)
                    valid = (sha256(temp) == row["old_observed"]["sha256"]
                             and old_metadata_matches(temp, row["old_observed"]))
                except Exception:
                    valid = False
                if not valid:
                    temp.unlink()
                    fsync_dir(temp.parent)
            if not temp.exists():
                missing_bytes += row["old_observed"]["size"]
            planned_temps.append((row, temp))
        if shutil.disk_usage(repo / "data/interim").free < missing_bytes + 5 * 1024 ** 3:
            raise OSError("insufficient free space to complete the deterministic rollback set")
        for row, temp in planned_temps:
            destination = Path(row["destination"])
            restore_temps.append((row, temp))
            if not temp.exists():
                copy_restore(row, temp)
            print(f"PREPARED RESTORE {row['physical_order']:02d}/25 {destination}", flush=True)
        states = [classify(row) for row in manifest["rows"]]
        if "conflict" in states:
            raise ValueError("destination conflict appeared during recovery preparation")
        planned_destinations = {row["destination"] for row, _temp in restore_temps}
        for row, state in zip(manifest["rows"], states):
            if state != "old" and row["destination"] not in planned_destinations:
                raise ValueError(f"previously-old destination changed during recovery prep: {row['destination']}")
        for row, temp in restore_temps:
            check_regular(temp)
            if (sha256(temp) != row["old_observed"]["sha256"]
                    or not old_metadata_matches(temp, row["old_observed"])):
                raise ValueError(f"restore temp changed during recovery prep: {temp}")
        journal = backup_dir / "RECOVERY_JOURNAL.jsonl"
        state_path = backup_dir / "RECOVERY_STATE.json"
        rollback_marker = backup_dir / "RECOVERY_ROLLED_BACK"
        rollback_payload = (
            "PASS\n"
            f"run_id={RUN_ID}\n"
            f"prepared_manifest_sha256={args.expected_manifest_sha256}\n"
        ).encode()
        atomic_json(state_path, {"status": "ROLLING_BACK", "utc": utc_now()})
        try:
            for row, temp in reversed(restore_temps):
                destination = Path(row["destination"])
                state = classify(row)
                if state == "conflict":
                    raise ValueError(f"recovery conflict: {destination}")
                if state != "old":
                    check_regular(temp)
                    if (sha256(temp) != row["old_observed"]["sha256"]
                            or not old_metadata_matches(temp, row["old_observed"])):
                        raise ValueError(f"pending restore temp drift: {temp}")
                    append_journal(journal, {
                        "event": "RESTORE_INTENT", "destination": str(destination),
                        "utc": utc_now(),
                    })
                    os.replace(temp, destination)
                    fsync_dir(destination.parent)
                    append_journal(journal, {
                        "event": "RESTORE_DONE", "destination": str(destination),
                        "utc": utc_now(),
                    })
            for row in manifest["rows"]:
                check_regular(row["destination"])
            if any(classify(row) != "old" for row in manifest["rows"]):
                raise ValueError("final recovery verification failed")
            atomic_json(state_path, {"status": "ROLLED_BACK", "utc": utc_now()})
            atomic_json(backup_dir / "ROLLBACK_RESULT.json", {
                "status": "PASS", "run_id": RUN_ID,
                "prepared_manifest_sha256": args.expected_manifest_sha256,
                "n_restored": 25, "completed_utc": utc_now(),
            })
            atomic_bytes(recovery, rollback_payload)
            os.replace(recovery, rollback_marker)
            fsync_dir(backup_dir)
            print("PASS: all 25 destinations match the reviewed old state.")
        except BaseException as exc:
            # As in promotion, a handled signal can arrive after the marker
            # rename but before the following Python instruction. The exact,
            # durable rollback marker is authoritative in that narrow window.
            if rollback_marker.exists() and rollback_marker.read_bytes() == rollback_payload:
                fsync_dir(backup_dir)
                print("Rollback marker is authoritative; all old hashes were already verified.",
                      file=sys.stderr)
            else:
                if not recovery.exists():
                    atomic_bytes(
                        recovery,
                        ("RECOVERY REQUIRED: rollback marker was not authoritative.\n"
                         f"run_id={RUN_ID}\n").encode(),
                    )
                atomic_json(state_path, {
                    "status": "RECOVERY_REQUIRED", "error": repr(exc), "utc": utc_now()})
                print("Recovery was interrupted; rerun this same reviewed recovery command.",
                      file=sys.stderr)
                raise
    finally:
        for _row, temp in restore_temps:
            try:
                temp.unlink()
            except FileNotFoundError:
                pass
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
