#!/usr/bin/env python3
"""Human-only application of one reviewed PREPARED_NOT_PROMOTED receipt."""
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


def check_regular(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_nlink != 1:
        raise ValueError(f"not a regular independent file: {path}")


def utc_now():
    return datetime.now(timezone.utc).isoformat()


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


def atomic_json(path, payload):
    atomic_bytes(path, (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode())


def append_journal(path, payload):
    path = Path(path)
    first = not path.exists()
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    if first:
        fsync_dir(path.parent)


def copy_temp(source, destination, mode, uid, gid, mtime_ns=None):
    destination = Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if destination.parent.is_symlink() or not destination.parent.is_dir():
        raise ValueError(f"unsafe temporary-file parent: {destination.parent}")
    shutil.copy2(source, destination)
    os.chmod(destination, int(mode, 8))
    os.chown(destination, uid, gid)
    if mtime_ns is not None:
        os.utime(destination, ns=(mtime_ns, mtime_ns))
    fsync_file(destination)
    fsync_dir(destination.parent)


def cleanup(paths):
    for path in paths:
        try:
            Path(path).unlink()
            fsync_dir(Path(path).parent)
        except FileNotFoundError:
            pass


def old_metadata_matches(path, old):
    stat = Path(path).stat()
    return (
        stat.st_size == old["size"]
        and f"{stat.st_mode & 0o7777:04o}" == old["mode"]
        and stat.st_uid == old["uid"]
        and stat.st_gid == old["gid"]
        and stat.st_mtime_ns == old["mtime_ns"]
    )


def classify_destination(row):
    destination = Path(row["destination"])
    actual = sha256(destination)
    if actual == row["old_observed"]["sha256"]:
        return "old" if old_metadata_matches(destination, row["old_observed"]) else "old_metadata_drift"
    if actual == row["source_observed"]["sha256"]:
        return "promoted"
    return "conflict"


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
expected_backup_parent = (repo / "data/interim").resolve()
if (manifest_path.name != "BACKUP_MANIFEST.json"
        or manifest_path.parent.parent != expected_backup_parent
        or not manifest_path.parent.name.startswith("_promotion_backup_")):
    raise ValueError("manifest must be BACKUP_MANIFEST.json in one promotion backup dir")
if sha256(manifest_path) != args.expected_manifest_sha256:
    raise ValueError("prepared manifest hash changed")
manifest = json.loads(manifest_path.read_text())
if manifest.get("status") != "PREPARED_NOT_PROMOTED" or manifest.get("run_id") != RUN_ID:
    raise ValueError("prepared manifest status/run is wrong")
if len(manifest.get("rows", [])) != 25:
    raise ValueError("prepared manifest does not contain 25 rows")
backup_dir = manifest_path.parent
if (backup_dir / "PREPARATION_INCOMPLETE").exists():
    raise ValueError("backup preparation is marked incomplete")
if (backup_dir / "PREPARED_NOT_PROMOTED").read_text() != "PREPARED_NOT_PROMOTED\n":
    raise ValueError("prepared receipt is missing")
if (backup_dir / "PROMOTION_RESULT.json").exists():
    raise FileExistsError("this prepared receipt was already promoted")

plan_path = repo / PLAN_REL
if sha256(plan_path) != PLAN_FILE_SHA256:
    raise ValueError("reviewed promotion plan changed")
plan = json.loads(plan_path.read_text())
expected_rows = []
for mapping in sorted(plan["plan"]["mappings"], key=lambda row: row["order"]):
    source = repo / mapping["source"]["path"].replace(OLD_RUN_ID, RUN_ID)
    destination = repo / mapping["destination"]["path"]
    expected_rows.append({
        "physical_order": len(expected_rows) + 1, "logical_order": mapping["order"],
        "stage": mapping["stage"], "kind": "reviewed_destination",
        "source": str(source), "destination": str(destination),
        "expected_old": mapping["destination"]["old"],
    })
    if mapping["stage"] == "deposit":
        expected_rows.append({
            "physical_order": len(expected_rows) + 1, "logical_order": mapping["order"],
            "stage": "deposit_mirror", "kind": "mckinley_mirror",
            "source": str(source), "destination": str(MCK / destination.name),
            "copy_source_after_promotion": str(destination),
            "expected_old": {
                "sha256": mapping["destination"]["old"]["sha256"],
                "size": mapping["destination"]["old"]["size"],
            },
        })
for observed, expected in zip(manifest["rows"], expected_rows):
    for key, value in expected.items():
        if key == "expected_old":
            continue
        if observed.get(key) != value:
            raise ValueError(f"prepared mapping differs from reviewed mapping ({key})")
    for key, value in expected["expected_old"].items():
        if observed.get("old_observed", {}).get(key) != value:
            raise ValueError(f"prepared old-state claim differs from reviewed plan ({key})")
    backup = Path(observed["backup"]).resolve()
    if backup_dir not in backup.parents:
        raise ValueError(f"backup path escapes prepared directory: {backup}")
if [row["physical_order"] for row in manifest["rows"]] != list(range(1, 26)):
    raise ValueError("prepared physical order is not exactly 1..25")
if len({row["destination"] for row in manifest["rows"]}) != 25:
    raise ValueError("prepared destinations are not unique")
if not sys.stdin.isatty():
    raise SystemExit("interactive terminal required")
token = f"PROMOTE:{RUN_ID}:{args.expected_manifest_sha256}"
print("This will atomically replace 21 live files and four McKinley mirrors.")
print(f"Verified backups remain at: {backup_dir}")
entered = input(f"Type exactly {token}\n> ")
if entered != token:
    raise SystemExit("ABORT: confirmation did not match")

lock_path = repo / "data/interim/.annotation_promotion.lock"
with lock_path.open("a+") as lock:
    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    journal = backup_dir / "PROMOTION_JOURNAL.jsonl"
    state_path = backup_dir / "PROMOTION_STATE.json"
    recovery = backup_dir / "RECOVERY_REQUIRED"
    commit_marker = backup_dir / "PROMOTION_COMMITTED"
    result_path = backup_dir / "PROMOTION_RESULT.json"
    live_manifest_path = backup_dir / "PROMOTED_LIVE_MANIFEST.json"
    # A hard stop while only preparing same-directory temps cannot have
    # changed live data (RECOVERY_REQUIRED is created before the first replace).
    # Reclaim only exact manifest-scoped temps from that safe pre-mutation
    # state, making the human command idempotently retryable.
    if state_path.exists() and not any(
            path.exists() for path in (recovery, commit_marker, result_path, live_manifest_path)):
        previous = json.loads(state_path.read_text())
        if previous.get("status") != "PREPARING_TEMPS":
            raise ValueError(f"unresolved promotion state: {previous}")
        suffix = args.expected_manifest_sha256[:16]
        for row in manifest["rows"]:
            destination = Path(row["destination"])
            for label in ("promote", "restore"):
                temp = destination.with_name(f".{destination.name}.{label}.{suffix}.tmp")
                if temp.exists() or temp.is_symlink():
                    temp.unlink()
                    fsync_dir(temp.parent)
        state_path.unlink()
        journal.unlink(missing_ok=True)
        fsync_dir(backup_dir)
    for path in (journal, state_path, recovery, commit_marker, result_path, live_manifest_path):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"refusing to reuse promotion transaction artifact: {path}")

    # Full preflight before the first mutation.
    if sha256(manifest["migration_status"]) != manifest["migration_status_sha256"]:
        raise ValueError("certified v2 migration status changed")
    if sha256(manifest["migration_ready"]) != manifest["migration_ready_sha256"]:
        raise ValueError("certified v2 readiness marker changed")
    status = json.loads(Path(manifest["migration_status"]).read_text())
    if status.get("run_id") != RUN_ID or status.get("status") != "PASS":
        raise ValueError("certified v2 status no longer has the exact PASS contract")
    for section in ("validation", "annotation_checks", "exports"):
        if (not status.get(section)
                or any(item.get("status") != "PASS" for item in status[section])):
            raise ValueError(f"certified v2 status section is not wholly PASS: {section}")
    protected = status.get("protected_live_immutability", {})
    if protected.get("status") != "PASS" or protected.get("changes") != []:
        raise ValueError("certified v2 protected-live result is not PASS/changes=[]")
    certified_sources = {}
    for report in status["annotation_checks"]:
        for record in report.get("objects", []):
            path = Path(record["path"])
            if not path.is_absolute():
                path = repo / path
            path = str(path.resolve())
            if (record.get("status") != "PASS" or record.get("checker_exit_code") != 0
                    or not record.get("sha256") or record.get("size") is None):
                raise ValueError(f"certified staged object is not cryptographically PASS: {path}")
            certified_sources[path] = {
                "sha256": record["sha256"], "size": record["size"]}
    if len(certified_sources) != 21:
        raise ValueError(f"certified v2 status must bind exactly 21 staged objects; got {len(certified_sources)}")
    for row in manifest["rows"]:
        source = Path(row["source"])
        destination = Path(row["destination"])
        backup = Path(row["backup"])
        for path in (source, destination, backup):
            check_regular(path)
        if sha256(source) != row["source_observed"]["sha256"]:
            raise ValueError(f"staged source changed: {source}")
        if row["source_observed"]["sha256"] != row["source_certified"]["sha256"]:
            raise ValueError(f"staged source is not the certified v2 artifact: {source}")
        certified = certified_sources.get(str(source.resolve()))
        if certified != {
                "sha256": row["source_observed"]["sha256"],
                "size": row["source_observed"]["size"]}:
            raise ValueError(f"staged source differs from certified v2 status: {source}")
        if sha256(destination) != row["old_observed"]["sha256"]:
            raise ValueError(f"live destination changed: {destination}")
        if sha256(backup) != row["old_observed"]["sha256"]:
            raise ValueError(f"backup changed: {backup}")
        source_stat = source.stat()
        destination_stat = destination.stat()
        if source_stat.st_size != row["source_observed"]["size"]:
            raise ValueError(f"staged source size changed: {source}")
        if source_stat.st_dev != destination_stat.st_dev:
            raise ValueError(f"source and destination are on different filesystems: {source}")
        for key, actual in (
            ("size", destination_stat.st_size),
            ("mode", f"{destination_stat.st_mode & 0o7777:04o}"),
            ("uid", destination_stat.st_uid),
            ("gid", destination_stat.st_gid),
            ("mtime_ns", destination_stat.st_mtime_ns),
            ("inode", destination_stat.st_ino),
        ):
            if actual != row["old_observed"][key]:
                raise ValueError(f"live destination metadata changed ({key}): {destination}")
        if (f"{source_stat.st_mode & 0o7777:04o}" != row["old_observed"]["mode"]
                or source_stat.st_uid != row["old_observed"]["uid"]
                or source_stat.st_gid != row["old_observed"]["gid"]):
            raise ValueError(f"staged source ownership/mode is not promotion-safe: {source}")

    old_total = sum(row["old_observed"]["size"] for row in manifest["rows"])
    new_total = sum(row["source_observed"]["size"] for row in manifest["rows"])
    free = shutil.disk_usage(repo / "data/interim").free
    # Reserve the primary restore set + promotion set and a second complete
    # restore set for the reviewed hard-crash recovery helper.
    required = 2 * old_total + new_total + 5 * 1024 ** 3
    if free < required:
        raise OSError(f"insufficient space for install/rollback temps: need={required}, free={free}")

    prepared = []
    created = []
    mutation_started = False
    committed = False
    commit_payload = None
    old_handlers = {}
    try:
        atomic_json(state_path, {
            "status": "PREPARING_TEMPS", "run_id": RUN_ID,
            "prepared_manifest_sha256": args.expected_manifest_sha256,
            "utc": utc_now(),
        })
        for row in manifest["rows"]:
            destination = Path(row["destination"])
            suffix = args.expected_manifest_sha256[:16]
            promote_tmp = destination.with_name(
                f".{destination.name}.promote.{suffix}.tmp")
            restore_tmp = destination.with_name(
                f".{destination.name}.restore.{suffix}.tmp")
            created.extend((promote_tmp, restore_tmp))
            copy_temp(
                row["source"], promote_tmp, row["old_observed"]["mode"],
                row["old_observed"]["uid"], row["old_observed"]["gid"])
            copy_temp(
                row["backup"], restore_tmp, row["old_observed"]["mode"],
                row["old_observed"]["uid"], row["old_observed"]["gid"],
                row["old_observed"]["mtime_ns"])
            if sha256(promote_tmp) != row["source_observed"]["sha256"]:
                raise ValueError(f"promotion temp failed verification: {promote_tmp}")
            if (sha256(restore_tmp) != row["old_observed"]["sha256"]
                    or not old_metadata_matches(restore_tmp, row["old_observed"])):
                raise ValueError(f"rollback temp failed verification: {restore_tmp}")
            prepared.append({"row": row, "promote_tmp": promote_tmp,
                             "restore_tmp": restore_tmp})
            print(f"PREPARED {row['physical_order']:02d}/25 {destination}", flush=True)

        # Revalidate every trust input after the lengthy copy phase and before
        # creating the durable recovery marker or mutating a live pathname.
        if sha256(manifest_path) != args.expected_manifest_sha256:
            raise ValueError("prepared manifest changed during install preparation")
        if sha256(manifest["migration_status"]) != manifest["migration_status_sha256"]:
            raise ValueError("certified v2 migration status changed during preparation")
        if sha256(manifest["migration_ready"]) != manifest["migration_ready_sha256"]:
            raise ValueError("certified v2 readiness marker changed during preparation")
        for item in prepared:
            row = item["row"]
            if sha256(row["source"]) != row["source_observed"]["sha256"]:
                raise ValueError(f"staged source changed during preparation: {row['source']}")
            if classify_destination(row) != "old":
                raise ValueError(f"live destination changed during preparation: {row['destination']}")
            if sha256(row["backup"]) != row["old_observed"]["sha256"]:
                raise ValueError(f"backup changed during preparation: {row['backup']}")
            if sha256(item["promote_tmp"]) != row["source_observed"]["sha256"]:
                raise ValueError(f"promotion temp changed during preparation: {item['promote_tmp']}")
            if (sha256(item["restore_tmp"]) != row["old_observed"]["sha256"]
                    or not old_metadata_matches(item["restore_tmp"], row["old_observed"])):
                raise ValueError(f"rollback temp changed during preparation: {item['restore_tmp']}")

        recovery_text = (
            "RECOVERY REQUIRED until PROMOTION_COMMITTED exists.\n"
            f"run_id={RUN_ID}\n"
            f"prepared_manifest_sha256={args.expected_manifest_sha256}\n"
        ).encode()
        atomic_bytes(recovery, recovery_text)
        mutation_started = True
        atomic_json(state_path, {"status": "APPLYING", "utc": utc_now()})
        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            old_handlers[sig] = signal.signal(sig, signal_as_exception)

        for item in prepared:
            row = item["row"]
            destination = Path(row["destination"])
            # A non-cooperating writer can change a late destination while
            # earlier rows are being installed; verify it immediately before
            # its atomic replacement.
            if classify_destination(row) != "old":
                raise ValueError(f"pending live destination drift: {destination}")
            check_regular(item["promote_tmp"])
            if sha256(item["promote_tmp"]) != row["source_observed"]["sha256"]:
                raise ValueError(f"pending promotion temp drift: {item['promote_tmp']}")
            promote_stat = item["promote_tmp"].stat()
            if (promote_stat.st_size != row["source_observed"]["size"]
                    or f"{promote_stat.st_mode & 0o7777:04o}" != row["old_observed"]["mode"]
                    or promote_stat.st_uid != row["old_observed"]["uid"]
                    or promote_stat.st_gid != row["old_observed"]["gid"]):
                raise ValueError(f"pending promotion temp metadata drift: {item['promote_tmp']}")
            append_journal(journal, {
                "event": "INSTALL_INTENT", "physical_order": row["physical_order"],
                "destination": str(destination), "utc": utc_now(),
            })
            os.replace(item["promote_tmp"], destination)
            fsync_dir(destination.parent)
            if sha256(destination) != row["source_observed"]["sha256"]:
                raise ValueError(f"installed hash failed: {destination}")
            append_journal(journal, {
                "event": "INSTALL_DONE", "physical_order": row["physical_order"],
                "destination": str(destination), "utc": utc_now(),
            })
            print(f"INSTALLED {row['physical_order']:02d}/25 {destination}", flush=True)

        atomic_json(state_path, {"status": "VERIFYING", "utc": utc_now()})
        live_records = []
        for row in manifest["rows"]:
            destination = Path(row["destination"])
            check_regular(destination)
            actual = sha256(destination)
            stat = destination.stat()
            if (actual != row["source_observed"]["sha256"]
                    or stat.st_size != row["source_observed"]["size"]
                    or f"{stat.st_mode & 0o7777:04o}" != row["old_observed"]["mode"]
                    or stat.st_uid != row["old_observed"]["uid"]
                    or stat.st_gid != row["old_observed"]["gid"]):
                raise ValueError(f"final promoted state failed: {destination}")
            live_records.append({
                "physical_order": row["physical_order"], "destination": str(destination),
                "sha256": actual, "size": stat.st_size,
                "mode": f"{stat.st_mode & 0o7777:04o}",
                "uid": stat.st_uid, "gid": stat.st_gid,
            })
        live_manifest = {
            "status": "VERIFIED_NOT_COMMITTED", "run_id": RUN_ID,
            "prepared_manifest_sha256": args.expected_manifest_sha256,
            "destinations": live_records, "verified_utc": utc_now(),
        }
        atomic_json(live_manifest_path, live_manifest)
        live_manifest_sha = sha256(live_manifest_path)
        atomic_json(state_path, {
            "status": "COMMIT_READY", "live_manifest_sha256": live_manifest_sha,
            "utc": utc_now(),
        })

        commit_payload = (
            "PASS\n"
            f"run_id={RUN_ID}\n"
            f"prepared_manifest_sha256={args.expected_manifest_sha256}\n"
            f"live_manifest_sha256={live_manifest_sha}\n"
        ).encode()
        atomic_bytes(recovery, commit_payload)
        os.replace(recovery, commit_marker)
        fsync_dir(backup_dir)
        committed = True
        postcommit_warnings = []
        try:
            atomic_json(result_path, {
                "status": "PASS", "run_id": RUN_ID,
                "prepared_manifest_sha256": args.expected_manifest_sha256,
                "live_manifest_sha256": live_manifest_sha,
                "n_installed": 25, "downstream_jobs": "NONE",
                "completed_utc": utc_now(),
            })
            atomic_json(state_path, {"status": "COMMITTED", "utc": utc_now()})
        except BaseException as exc:
            postcommit_warnings.append(repr(exc))
        try:
            cleanup([item["restore_tmp"] for item in prepared])
        except BaseException as exc:
            postcommit_warnings.append(f"restore-temp cleanup: {exc!r}")
        print("PASS: promotion complete; no downstream job was run.")
        for warning in postcommit_warnings:
            print(f"POST-COMMIT WARNING (do not roll back): {warning}", file=sys.stderr)
    except BaseException as promotion_error:
        # Resolve the signal/os.replace ambiguity. Once the exact durable
        # commit marker exists, it is authoritative and rollback is forbidden.
        if mutation_started and not committed and commit_marker.exists():
            if commit_payload is None or commit_marker.read_bytes() != commit_payload:
                atomic_json(state_path, {
                    "status": "RECOVERY_REQUIRED",
                    "error": "invalid commit marker after commit-window exception",
                    "utc": utc_now(),
                })
                raise RuntimeError("ambiguous promotion commit marker; manual recovery required")
            fsync_dir(backup_dir)
            committed = True
            print("Commit marker is authoritative; refusing automatic rollback.", file=sys.stderr)
        if mutation_started and not committed:
            append_journal(journal, {"event": "AUTO_ROLLBACK_BEGIN", "utc": utc_now()})
            rollback_error = None
            try:
                for item in reversed(prepared):
                    row = item["row"]
                    destination = Path(row["destination"])
                    state = classify_destination(row)
                    if state == "conflict":
                        raise RuntimeError(f"rollback conflict: {destination}")
                    if state in {"promoted", "old_metadata_drift"}:
                        restore_tmp = item["restore_tmp"]
                        if (not restore_tmp.exists()
                                or sha256(restore_tmp) != row["old_observed"]["sha256"]):
                            restore_tmp = destination.with_name(
                                f".{destination.name}.emergency_restore.{os.getpid()}.tmp")
                            copy_temp(
                                row["backup"], restore_tmp, row["old_observed"]["mode"],
                                row["old_observed"]["uid"], row["old_observed"]["gid"],
                                row["old_observed"]["mtime_ns"])
                        if (sha256(restore_tmp) != row["old_observed"]["sha256"]
                                or not old_metadata_matches(restore_tmp, row["old_observed"])):
                            raise RuntimeError(f"emergency rollback temp failed: {restore_tmp}")
                        check_regular(restore_tmp)
                        append_journal(journal, {
                            "event": "RESTORE_INTENT", "destination": str(destination),
                            "utc": utc_now(),
                        })
                        os.replace(restore_tmp, destination)
                        fsync_dir(destination.parent)
                        append_journal(journal, {
                            "event": "RESTORE_DONE", "destination": str(destination),
                            "utc": utc_now(),
                        })
                for row in manifest["rows"]:
                    check_regular(row["destination"])
                    if classify_destination(row) != "old":
                        raise RuntimeError(f"rollback verification failed: {row['destination']}")
                atomic_json(state_path, {
                    "status": "ROLLED_BACK", "error": repr(promotion_error),
                    "utc": utc_now(),
                })
                atomic_json(result_path, {
                    "status": "ROLLED_BACK", "run_id": RUN_ID,
                    "reason": repr(promotion_error), "utc": utc_now(),
                })
                atomic_bytes(backup_dir / "AUTO_ROLLBACK_COMPLETE", b"PASS\n")
                recovery.unlink(missing_ok=True)
                fsync_dir(backup_dir)
                try:
                    cleanup(created)
                except BaseException as cleanup_error:
                    print(f"ROLLBACK CLEANUP WARNING: {cleanup_error!r}", file=sys.stderr)
                print("Promotion failed; all 25 old files were restored.", file=sys.stderr)
            except BaseException as exc:
                rollback_error = exc
                atomic_json(state_path, {
                    "status": "RECOVERY_REQUIRED", "promotion_error": repr(promotion_error),
                    "rollback_error": repr(exc), "utc": utc_now(),
                })
                print(f"CRITICAL: automatic rollback incomplete: {exc!r}", file=sys.stderr)
            if rollback_error is not None:
                raise rollback_error from promotion_error
        elif not mutation_started:
            cleanup(created)
        raise
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
