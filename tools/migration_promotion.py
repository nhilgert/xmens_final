#!/usr/bin/env python3
"""Manifest-bound promotion and rollback for xmens.annotation.v1.

This module is invoked only through promote_migration.sh and
rollback_migration.sh.  Its default modes are read-only.  Mutation requires an
explicit action, a unique promotion id, and a typed confirmation string.
"""

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import signal
import stat
import subprocess
import sys
from datetime import datetime, timezone


RUN_ID = "xmens_annotation_v1_20260810"
REPO = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
PLAN_PATH = REPO / "tools/promotion_plan_xmens_annotation_v1_20260810.json"
EXPECTED_PLAN_SHA256 = "5c0dda6041034be2a6d941efe2430b028d143eaa0b1e40b055164972d0503ad3"
REVIEWED_BASE_GIT_HEAD = "2cdc6a50ef2c84c9c3340a86fafff94e95929fa6"
BACKUP_BASE = REPO / "data/interim/migration_backups" / RUN_ID
LOCK_PATH = REPO / "data/interim/.annotation_promotion.lock"
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$")
PHASES = ["reannotated", "refined", "final", "fibdec", "deposit"]
CRITICAL_PATHS = {
    "tools/migration_promotion.py",
    "tools/promote_migration.sh",
    "tools/rollback_migration.sh",
    "tools/promotion_plan_xmens_annotation_v1_20260810.json",
}
ALLOWED_POST_REVIEW_COMMITS = CRITICAL_PATHS | {"MIGRATION_VERIFY_REPORT.md"}


class SafetyError(RuntimeError):
    pass


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path, chunk=16 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb", buffering=0) as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def canonical_bytes(value):
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def fsync_dir(path):
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def fsync_file(path):
    with open(path, "rb+", buffering=0) as handle:
        os.fsync(handle.fileno())


def atomic_bytes(path, payload, mode=0o644):
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    if tmp.exists():
        raise SafetyError(f"temporary path already exists: {tmp}")
    with open(tmp, "xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.chmod(tmp, mode)
    fsync_file(tmp)
    os.replace(tmp, path)
    fsync_dir(path.parent)


def atomic_json(path, value):
    atomic_bytes(path, json.dumps(value, indent=2, sort_keys=True).encode() + b"\n")


def append_journal(path, event):
    existed = Path(path).exists()
    line = json.dumps(event, sort_keys=True, separators=(",", ":")) + "\n"
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())
    if not existed:
        fsync_dir(Path(path).parent)


def require_relative(path_text):
    pure = PurePosixPath(path_text)
    if pure.is_absolute() or ".." in pure.parts or str(pure) != path_text:
        raise SafetyError(f"unsafe relative path: {path_text!r}")
    return REPO / pure


def reject_symlink_components(path, allow_missing_leaf=False):
    path = Path(path)
    try:
        rel = path.relative_to(REPO)
    except ValueError as exc:
        raise SafetyError(f"path escapes repository: {path}") from exc
    cursor = REPO
    if os.path.islink(REPO):
        raise SafetyError(f"repository root is a symlink: {REPO}")
    for index, part in enumerate(rel.parts):
        cursor = cursor / part
        if allow_missing_leaf and index == len(rel.parts) - 1 and not cursor.exists():
            break
        try:
            info = os.lstat(cursor)
        except FileNotFoundError as exc:
            raise SafetyError(f"missing path component: {cursor}") from exc
        if stat.S_ISLNK(info.st_mode):
            raise SafetyError(f"symlink is forbidden: {cursor}")


def regular_single_link(path):
    reject_symlink_components(path)
    info = os.lstat(path)
    if not stat.S_ISREG(info.st_mode):
        raise SafetyError(f"not a regular file: {path}")
    if info.st_nlink != 1:
        raise SafetyError(f"hard-linked file is forbidden (nlink={info.st_nlink}): {path}")
    return info


def assert_path_contract(plan):
    mappings = sorted(plan["mappings"], key=lambda item: item["order"])
    if len(mappings) != 21 or [m["order"] for m in mappings] != list(range(1, 22)):
        raise SafetyError("plan must contain exactly ordinals 1..21")
    if [m["stage"] for m in mappings] != sorted(
        [m["stage"] for m in mappings], key=PHASES.index
    ):
        raise SafetyError("plan is not in dependency order")
    sources = [m["source"]["path"] for m in mappings]
    destinations = [m["destination"]["path"] for m in mappings]
    if len(set(sources)) != 21 or len(set(destinations)) != 21:
        raise SafetyError("plan contains duplicate sources or destinations")

    source_prefixes = (
        f"data/interim/migrated/{RUN_ID}/",
        f"dataverse/objects_migrated/{RUN_ID}/",
    )
    exact_interim = {
        "data/interim/refined_human.h5ad",
        "data/interim/refined_mouse.h5ad",
        "data/interim/final_human.h5ad",
        "data/interim/final_mouse.h5ad",
        "data/interim/final_human_fib_dec.h5ad",
        "data/interim/final_mouse_fib_dec.h5ad",
    }
    for mapping in mappings:
        source = mapping["source"]["path"]
        destination = mapping["destination"]["path"]
        require_relative(source)
        require_relative(destination)
        if not source.startswith(source_prefixes):
            raise SafetyError(f"source is outside the staged roots: {source}")
        destination_ok = (
            destination.startswith("data/interim/reannotated/")
            and destination.count("/") == 3
        ) or destination in exact_interim or (
            destination.startswith("dataverse/objects/")
            and destination.count("/") == 2
        )
        if not destination_ok:
            raise SafetyError(f"destination is outside the literal live roots: {destination}")
    return mappings


def load_plan(path):
    path = Path(path)
    regular_single_link(path)
    audit = json.loads(path.read_text())
    plan = audit.get("plan")
    if not isinstance(plan, dict):
        raise SafetyError("plan wrapper does not contain an object named 'plan'")
    digest = hashlib.sha256(canonical_bytes(plan)).hexdigest()
    if digest != EXPECTED_PLAN_SHA256:
        raise SafetyError(f"plan digest mismatch: {digest}")
    if audit.get("canonical_plan_sha256") != EXPECTED_PLAN_SHA256:
        raise SafetyError("plan wrapper advertises the wrong canonical digest")
    if audit.get("reviewed_repository_head") != REVIEWED_BASE_GIT_HEAD:
        raise SafetyError("plan wrapper advertises the wrong reviewed repository head")
    if plan.get("plan_schema_version") != "xmens.promotion-plan.v1":
        raise SafetyError("unsupported plan schema")
    if plan.get("run_id") != RUN_ID or plan.get("repo_root") != str(REPO):
        raise SafetyError("plan run id or repository root mismatch")
    return plan, assert_path_contract(plan)


def load_saved_plan(path):
    """Load the immutable canonical plan copied into a recovery directory."""
    path = Path(path)
    regular_single_link(path)
    plan = json.loads(path.read_text())
    digest = hashlib.sha256(canonical_bytes(plan)).hexdigest()
    if digest != EXPECTED_PLAN_SHA256:
        raise SafetyError(f"saved recovery plan digest mismatch: {digest}")
    if plan.get("plan_schema_version") != "xmens.promotion-plan.v1":
        raise SafetyError("unsupported saved recovery plan schema")
    if plan.get("run_id") != RUN_ID or plan.get("repo_root") != str(REPO):
        raise SafetyError("saved recovery plan run id or repository root mismatch")
    return plan, assert_path_contract(plan)


def verify_status_pins(plan):
    for pin in plan["status_pins"]:
        path = require_relative(pin["path"])
        info = regular_single_link(path)
        if info.st_size != pin["size"] or sha256_file(path) != pin["sha256"]:
            raise SafetyError(f"pinned status/input changed: {pin['path']}")

    ready = REPO / f"data/interim/migrated/{RUN_ID}/MIGRATION_REVIEW_READY"
    if ready.read_bytes() != b"PASS\n":
        raise SafetyError("MIGRATION_REVIEW_READY is not exactly PASS\\n")

    manifests = REPO / f"data/interim/migrated/{RUN_ID}/manifests"
    json_names = [
        "MIGRATION_STATUS.json",
        "validation_human.json",
        "validation_mouse.json",
        "check_annotations_human.json",
        "check_annotations_mouse.json",
        "export_human.json",
        "export_mouse.json",
        "protected_live_after.json",
    ]
    for name in json_names:
        payload = json.loads((manifests / name).read_text())
        if payload.get("status") != "PASS":
            raise SafetyError(f"{name} does not report PASS")
        if payload.get("failures") not in (None, []):
            raise SafetyError(f"{name} contains failures")
        if payload.get("suspects") not in (None, []):
            raise SafetyError(f"{name} contains suspects")
        if name == "protected_live_after.json" and payload.get("changes") != []:
            raise SafetyError("protected_live_after.json contains changes")


def verify_inventory(mappings):
    expected = {require_relative(m["source"]["path"]) for m in mappings}
    roots = [
        REPO / f"data/interim/migrated/{RUN_ID}",
        REPO / f"dataverse/objects_migrated/{RUN_ID}",
    ]
    observed = set()
    for root in roots:
        reject_symlink_components(root)
        for path in root.rglob("*"):
            if ".tmp" in path.name:
                raise SafetyError(f"temporary remnant in staging: {path}")
            if path.suffix == ".h5ad":
                observed.add(path)
    if observed != expected:
        missing = sorted(str(p.relative_to(REPO)) for p in expected - observed)
        extra = sorted(str(p.relative_to(REPO)) for p in observed - expected)
        raise SafetyError(f"staged H5AD inventory mismatch; missing={missing}, extra={extra}")


def load_baseline(plan):
    path = require_relative(plan["protected_baseline"]["path"])
    payload = json.loads(path.read_text())
    protected = {item["path"]: item for item in payload["protected"]}
    if len(protected) != 26 or len(payload["protected"]) != 26:
        raise SafetyError("protected baseline must contain 26 unique objects")
    return protected


def verify_sources_and_old_live(mappings, protected, strict_metadata=True):
    for mapping in mappings:
        source = require_relative(mapping["source"]["path"])
        destination = require_relative(mapping["destination"]["path"])
        src_info = regular_single_link(source)
        dst_info = regular_single_link(destination)
        expected_source = mapping["source"]
        expected_old = mapping["destination"]["old"]
        if src_info.st_size != expected_source["size"]:
            raise SafetyError(f"source size drift: {source}")
        if sha256_file(source) != expected_source["sha256"]:
            raise SafetyError(f"source SHA-256 drift: {source}")
        if dst_info.st_size != expected_old["size"]:
            raise SafetyError(f"live size drift: {destination}")
        if sha256_file(destination) != expected_old["sha256"]:
            raise SafetyError(f"live SHA-256 drift: {destination}")
        baseline = protected.get(mapping["destination"]["path"])
        if baseline is None:
            raise SafetyError(f"destination absent from protected baseline: {destination}")
        for key in ("sha256", "size", "inode", "mtime_ns"):
            if baseline[key] != expected_old[key]:
                raise SafetyError(f"plan/baseline {key} mismatch: {destination}")
        if strict_metadata:
            actual = {
                "inode": dst_info.st_ino,
                "mtime_ns": dst_info.st_mtime_ns,
                "mode": f"{stat.S_IMODE(dst_info.st_mode):04o}",
                "uid": dst_info.st_uid,
                "gid": dst_info.st_gid,
            }
            for key, value in actual.items():
                if value != expected_old[key]:
                    raise SafetyError(f"live metadata drift ({key}): {destination}")


def non_targets(mappings, protected):
    destinations = {m["destination"]["path"] for m in mappings}
    remainder = {path: item for path, item in protected.items() if path not in destinations}
    if len(remainder) != 5:
        raise SafetyError(f"expected five protected non-targets, found {len(remainder)}")
    return remainder


def verify_non_targets(items):
    for relative, expected in items.items():
        path = require_relative(relative)
        info = regular_single_link(path)
        if info.st_size != expected["size"] or sha256_file(path) != expected["sha256"]:
            raise SafetyError(f"protected non-target changed: {relative}")


def current_git_head():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
    ).strip()


def verify_reviewed_code(expected_git_head):
    head = current_git_head()
    if head != expected_git_head:
        raise SafetyError(f"Git HEAD changed: expected {expected_git_head}, found {head}")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", REVIEWED_BASE_GIT_HEAD, head],
        cwd=REPO,
    )
    if ancestor.returncode != 0:
        raise SafetyError(f"reviewed base {REVIEWED_BASE_GIT_HEAD} is not an ancestor of {head}")
    changed = set(
        subprocess.check_output(
            ["git", "diff", "--name-only", f"{REVIEWED_BASE_GIT_HEAD}..{head}"],
            cwd=REPO,
            text=True,
        ).splitlines()
    )
    unexpected = sorted(changed - ALLOWED_POST_REVIEW_COMMITS)
    if unexpected:
        raise SafetyError(f"unexpected committed paths since reviewed base: {unexpected}")
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--", *sorted(CRITICAL_PATHS)],
        cwd=REPO,
        text=True,
    ).strip()
    if dirty:
        raise SafetyError(f"critical promotion files are not committed/clean:\n{dirty}")


def preflight(plan, mappings, expected_git_head=None, strict_live=True):
    if Path.cwd().resolve() != REPO:
        raise SafetyError(f"run from exactly {REPO}")
    if expected_git_head is not None:
        verify_reviewed_code(expected_git_head)
    verify_status_pins(plan)
    verify_inventory(mappings)
    protected = load_baseline(plan)
    if strict_live:
        verify_sources_and_old_live(mappings, protected, strict_metadata=True)
    targets = non_targets(mappings, protected)
    verify_non_targets(targets)
    blockers = incomplete_transactions()
    if blockers:
        raise SafetyError(f"unresolved promotion/recovery transaction(s): {blockers}")
    return protected, targets


def ensure_devices_and_space(mappings):
    repo_dev = os.stat(REPO).st_dev
    old_total = sum(m["destination"]["old"]["size"] for m in mappings)
    new_total = sum(m["source"]["size"] for m in mappings)
    for mapping in mappings:
        source = require_relative(mapping["source"]["path"])
        destination = require_relative(mapping["destination"]["path"])
        if os.stat(source).st_dev != repo_dev or os.stat(destination.parent).st_dev != repo_dev:
            raise SafetyError(f"source/destination is not on repository filesystem: {destination}")
    # Reserve backup + primary restore temps + promote temps, plus one more
    # complete old set so a post-crash manual rollback can pre-stage safely
    # without relying on any surviving primary restore temp.
    required = old_total * 3 + new_total + 5 * 1024**3
    free = shutil.disk_usage(REPO).free
    if free < required:
        raise SafetyError(f"insufficient free space: need {required}, have {free}")
    return old_total, new_total, required, free


def copy_independent(source, destination, mode=None, uid=None, gid=None):
    destination = Path(destination)
    if destination.exists() or destination.is_symlink():
        raise SafetyError(f"refusing to overwrite prepared file: {destination}")
    reject_symlink_components(destination.parent)
    subprocess.run(
        [
            "cp",
            "--reflink=never",
            "--sparse=always",
            "--preserve=mode,timestamps",
            "--",
            str(source),
            str(destination),
        ],
        check=True,
    )
    if mode is not None:
        os.chmod(destination, mode)
    if uid is not None and gid is not None:
        os.chown(destination, uid, gid)
    fsync_file(destination)
    fsync_dir(destination.parent)


def acquire_lock():
    reject_symlink_components(LOCK_PATH.parent)
    handle = open(LOCK_PATH, "a+")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        handle.close()
        raise SafetyError("another promotion/rollback process holds the lock") from exc
    return handle


def incomplete_transactions():
    if not BACKUP_BASE.exists():
        return []
    reject_symlink_components(BACKUP_BASE)
    blockers = set(BACKUP_BASE.glob("*/RECOVERY_REQUIRED"))
    incomplete_statuses = {
        "PREPARING",
        "PREPARED",
        "APPLYING",
        "VERIFYING",
        "COMMIT_READY",
        "RECOVERY_REQUIRED",
        "ROLLING_BACK",
    }
    for state_path in BACKUP_BASE.glob("*/STATE.json"):
        try:
            status_value = json.loads(state_path.read_text()).get("status")
        except Exception:
            blockers.add(state_path)
            continue
        if status_value in incomplete_statuses:
            blockers.add(state_path)
    return sorted(blockers)


def create_backup_base():
    parent = REPO / "data/interim"
    reject_symlink_components(parent)
    for directory in (parent / "migration_backups", BACKUP_BASE):
        if directory.exists():
            reject_symlink_components(directory)
        else:
            directory.mkdir()
            fsync_dir(directory.parent)


def restorable_metadata_matches(destination, old):
    info = os.lstat(destination)
    return (
        f"{stat.S_IMODE(info.st_mode):04o}" == old["mode"]
        and info.st_uid == old["uid"]
        and info.st_gid == old["gid"]
        and info.st_mtime_ns == old["mtime_ns"]
    )


def classify(destination, old, new_sha):
    actual = sha256_file(destination)
    if actual == old["sha256"]:
        return "old" if restorable_metadata_matches(destination, old) else "old_metadata_drift"
    if actual == new_sha:
        return "promoted"
    return "conflict"


def verify_pending_destination_old(mapping):
    destination = require_relative(mapping["destination"]["path"])
    old = mapping["destination"]["old"]
    info = regular_single_link(destination)
    actual = {
        "size": info.st_size,
        "inode": info.st_ino,
        "mtime_ns": info.st_mtime_ns,
        "mode": f"{stat.S_IMODE(info.st_mode):04o}",
        "uid": info.st_uid,
        "gid": info.st_gid,
    }
    for key, value in actual.items():
        if value != old[key]:
            raise SafetyError(f"pending destination drift ({key}): {destination}")
    if sha256_file(destination) != old["sha256"]:
        raise SafetyError(f"pending destination SHA-256 drift: {destination}")


def observed_record(relative):
    path = require_relative(relative)
    info = regular_single_link(path)
    return {
        "path": relative,
        "sha256": sha256_file(path),
        "size": info.st_size,
        "inode": info.st_ino,
        "mtime_ns": info.st_mtime_ns,
        "mode": f"{stat.S_IMODE(info.st_mode):04o}",
        "uid": info.st_uid,
        "gid": info.st_gid,
    }


def cleanup_paths(paths):
    for path in paths:
        path = Path(path)
        try:
            path.unlink()
            fsync_dir(path.parent)
        except FileNotFoundError:
            pass


def prepare_backup_and_temps(mappings, promotion_dir, promotion_id):
    files_dir = promotion_dir / "files"
    files_dir.mkdir(parents=True, exist_ok=False)
    fsync_dir(files_dir.parent)
    prepared = []
    created_temps = []
    try:
        for mapping in mappings:
            print(
                f"Preparing {mapping['order']:02d}/21 {mapping['destination']['path']}",
                flush=True,
            )
            source = require_relative(mapping["source"]["path"])
            destination = require_relative(mapping["destination"]["path"])
            old = mapping["destination"]["old"]
            backup = files_dir / mapping["destination"]["path"]
            backup.parent.mkdir(parents=True, exist_ok=True)
            copy_independent(destination, backup)
            if sha256_file(backup) != old["sha256"]:
                raise SafetyError(f"backup verification failed: {backup}")

            promote_tmp = destination.parent / (
                f".{destination.name}.promote.{promotion_id}.{os.getpid()}.tmp"
            )
            restore_tmp = destination.parent / (
                f".{destination.name}.restore.{promotion_id}.{os.getpid()}.tmp"
            )
            created_temps.extend((promote_tmp, restore_tmp))
            mode = int(old["mode"], 8)
            copy_independent(source, promote_tmp, mode, old["uid"], old["gid"])
            copy_independent(backup, restore_tmp, mode, old["uid"], old["gid"])
            if sha256_file(promote_tmp) != mapping["source"]["sha256"]:
                raise SafetyError(f"promotion temp verification failed: {promote_tmp}")
            if sha256_file(restore_tmp) != old["sha256"]:
                raise SafetyError(f"rollback temp verification failed: {restore_tmp}")
            prepared.append(
                {
                    "mapping": mapping,
                    "backup": backup,
                    "promote_tmp": promote_tmp,
                    "restore_tmp": restore_tmp,
                }
            )
    except BaseException:
        cleanup_paths(created_temps)
        raise
    return prepared


def auto_rollback(prepared, journal, state_path, recovery, protected_non_targets):
    append_journal(journal, {"event": "AUTO_ROLLBACK_BEGIN", "utc": utc_now()})
    conflicts = []
    for item in reversed(prepared):
        mapping = item["mapping"]
        destination = require_relative(mapping["destination"]["path"])
        old = mapping["destination"]["old"]
        old_sha = old["sha256"]
        new_sha = mapping["source"]["sha256"]
        state = classify(destination, old, new_sha)
        if state == "conflict":
            conflicts.append(str(destination.relative_to(REPO)))
            append_journal(
                journal,
                {"event": "RESTORE_CONFLICT", "path": str(destination.relative_to(REPO)), "utc": utc_now()},
            )
            continue
        if state in {"promoted", "old_metadata_drift"}:
            restore_tmp = item["restore_tmp"]
            if not restore_tmp.exists() or sha256_file(restore_tmp) != old_sha:
                restore_tmp = destination.parent / (
                    f".{destination.name}.auto-restore.{os.getpid()}.tmp"
                )
                copy_independent(
                    item["backup"],
                    restore_tmp,
                    int(old["mode"], 8),
                    old["uid"],
                    old["gid"],
                )
            append_journal(
                journal,
                {"event": "RESTORE_INTENT", "path": str(destination.relative_to(REPO)), "utc": utc_now()},
            )
            os.replace(restore_tmp, destination)
            fsync_dir(destination.parent)
            if sha256_file(destination) != old_sha:
                raise SafetyError(f"automatic rollback verification failed: {destination}")
            append_journal(
                journal,
                {"event": "RESTORE_DONE", "path": str(destination.relative_to(REPO)), "utc": utc_now()},
            )
    verify_non_targets(protected_non_targets)
    if conflicts:
        raise SafetyError(f"automatic rollback preserved conflicting live paths: {conflicts}")
    for item in prepared:
        mapping = item["mapping"]
        destination = require_relative(mapping["destination"]["path"])
        if classify(
            destination,
            mapping["destination"]["old"],
            mapping["source"]["sha256"],
        ) != "old":
            raise SafetyError(f"automatic rollback final verification failed: {destination}")
    atomic_json(state_path, {"status": "ROLLED_BACK", "utc": utc_now()})
    atomic_json(
        recovery.parent / "PROMOTION_RESULT.json",
        {"status": "ROLLED_BACK", "reason": "automatic rollback", "utc": utc_now()},
    )
    atomic_bytes(recovery.with_name("AUTO_ROLLBACK_COMPLETE"), b"PASS\n")
    recovery.unlink(missing_ok=True)
    fsync_dir(recovery.parent)


def do_promote(plan, mappings, args):
    if not ID_RE.fullmatch(args.promotion_id):
        raise SafetyError("invalid promotion id")
    if args.confirm != f"PROMOTE:{RUN_ID}":
        raise SafetyError(f"confirmation must be PROMOTE:{RUN_ID}")
    if not re.fullmatch(r"[0-9a-f]{40}", args.expected_git_head):
        raise SafetyError("--expected-git-head must be a full lowercase SHA-1")
    lock = acquire_lock()
    try:
        protected, targets = preflight(
            plan, mappings, expected_git_head=args.expected_git_head, strict_live=True
        )
        old_total, new_total, required, free = ensure_devices_and_space(mappings)
        promotion_dir = BACKUP_BASE / args.promotion_id
        create_backup_base()
        promotion_dir.mkdir(exist_ok=False)
        fsync_dir(promotion_dir.parent)
        state_path = promotion_dir / "STATE.json"
        journal = promotion_dir / "JOURNAL.jsonl"
        recovery = promotion_dir / "RECOVERY_REQUIRED"
        commit_marker = promotion_dir / "PROMOTION_COMMITTED"
        atomic_json(
            state_path,
            {
                "status": "PREPARING",
                "promotion_id": args.promotion_id,
                "plan_sha256": EXPECTED_PLAN_SHA256,
                "git_head": current_git_head(),
                "utc": utc_now(),
            },
        )
        recovery_tools = promotion_dir / "recovery_tools"
        critical_records = []
        try:
            atomic_bytes(promotion_dir / "PROMOTION_PLAN.json", canonical_bytes(plan))
            recovery_tools.mkdir()
            for relative in sorted(CRITICAL_PATHS):
                source = require_relative(relative)
                regular_single_link(source)
                destination = recovery_tools / Path(relative).name
                copy_independent(source, destination)
                critical_records.append(
                    {
                        "path": relative,
                        "recovery_copy": str(destination.relative_to(REPO)),
                        "sha256": sha256_file(destination),
                        "size": destination.stat().st_size,
                    }
                )
        except BaseException:
            atomic_json(
                state_path,
                {"status": "ABORTED_BEFORE_MUTATION", "utc": utc_now()},
            )
            if recovery_tools.exists():
                shutil.rmtree(recovery_tools)
                fsync_dir(promotion_dir)
            raise

        prepared = []
        mutation_started = False
        committed = False
        try:
            print("Preparing and verifying 21 independent backups and 42 install/restore temps...", flush=True)
            prepared = prepare_backup_and_temps(mappings, promotion_dir, args.promotion_id)
            backup_manifest = {
                "schema_version": "xmens.backup-manifest.v1",
                "promotion_id": args.promotion_id,
                "plan_sha256": EXPECTED_PLAN_SHA256,
                "created_utc": utc_now(),
                "critical_recovery_tools": critical_records,
                "protected_non_targets": [
                    {
                        "path": relative,
                        "sha256": targets[relative]["sha256"],
                        "size": targets[relative]["size"],
                    }
                    for relative in sorted(targets)
                ],
                "entries": [
                    {
                        "order": item["mapping"]["order"],
                        "source": item["mapping"]["source"]["path"],
                        "destination": item["mapping"]["destination"]["path"],
                        "backup": str(item["backup"].relative_to(REPO)),
                        "old_sha256": item["mapping"]["destination"]["old"]["sha256"],
                        "promoted_sha256": item["mapping"]["source"]["sha256"],
                        "old_mode": item["mapping"]["destination"]["old"]["mode"],
                        "old_uid": item["mapping"]["destination"]["old"]["uid"],
                        "old_gid": item["mapping"]["destination"]["old"]["gid"],
                    }
                    for item in prepared
                ],
            }
            backup_manifest_path = promotion_dir / "BACKUP_MANIFEST.json"
            atomic_json(backup_manifest_path, backup_manifest)
            backup_manifest_sha = sha256_file(backup_manifest_path)

            # Catch a concurrent writer after all lengthy copy preparation.
            verify_reviewed_code(args.expected_git_head)
            verify_status_pins(plan)
            verify_inventory(mappings)
            for record in critical_records:
                if sha256_file(require_relative(record["path"])) != record["sha256"]:
                    raise SafetyError(f"critical promotion file changed: {record['path']}")
            verify_sources_and_old_live(mappings, protected, strict_metadata=True)
            verify_non_targets(targets)
            atomic_json(
                state_path,
                {
                    "status": "PREPARED",
                    "promotion_id": args.promotion_id,
                    "plan_sha256": EXPECTED_PLAN_SHA256,
                    "backup_manifest_sha256": backup_manifest_sha,
                    "space": {
                        "old_bytes": old_total,
                        "new_bytes": new_total,
                        "required_bytes": required,
                        "free_bytes": free,
                    },
                    "utc": utc_now(),
                },
            )
            atomic_bytes(
                recovery,
                (
                    f"RECOVERY REQUIRED until PROMOTION_COMMITTED exists.\n"
                    f"promotion_id={args.promotion_id}\n"
                    f"backup_manifest_sha256={backup_manifest_sha}\n"
                ).encode(),
            )
            mutation_started = True
            atomic_json(state_path, {"status": "APPLYING", "utc": utc_now()})
            for item in prepared:
                mapping = item["mapping"]
                destination = require_relative(mapping["destination"]["path"])
                print(
                    f"Installing {mapping['order']:02d}/21 {mapping['destination']['path']}",
                    flush=True,
                )
                # Recheck immediately before replacement; the earlier complete
                # pass may have taken minutes on these large objects.
                verify_pending_destination_old(mapping)
                append_journal(
                    journal,
                    {"event": "INSTALL_INTENT", "order": mapping["order"], "path": mapping["destination"]["path"], "utc": utc_now()},
                )
                os.replace(item["promote_tmp"], destination)
                fsync_dir(destination.parent)
                if sha256_file(destination) != mapping["source"]["sha256"]:
                    raise SafetyError(f"installed SHA-256 mismatch: {destination}")
                append_journal(
                    journal,
                    {"event": "INSTALL_DONE", "order": mapping["order"], "path": mapping["destination"]["path"], "utc": utc_now()},
                )

            atomic_json(state_path, {"status": "VERIFYING", "utc": utc_now()})
            for mapping in mappings:
                source = require_relative(mapping["source"]["path"])
                destination = require_relative(mapping["destination"]["path"])
                if sha256_file(source) != mapping["source"]["sha256"]:
                    raise SafetyError(f"staged source changed during promotion: {source}")
                if sha256_file(destination) != mapping["source"]["sha256"]:
                    raise SafetyError(f"post-promotion verification failed: {destination}")
                info = regular_single_link(destination)
                old = mapping["destination"]["old"]
                if (
                    f"{stat.S_IMODE(info.st_mode):04o}" != old["mode"]
                    or info.st_uid != old["uid"]
                    or info.st_gid != old["gid"]
                ):
                    raise SafetyError(f"post-promotion mode/owner mismatch: {destination}")
            verify_non_targets(targets)
            target_records = []
            for mapping in mappings:
                record = observed_record(mapping["destination"]["path"])
                old = mapping["destination"]["old"]
                if (
                    record["sha256"] != mapping["source"]["sha256"]
                    or record["size"] != mapping["source"]["size"]
                    or record["mode"] != old["mode"]
                    or record["uid"] != old["uid"]
                    or record["gid"] != old["gid"]
                ):
                    raise SafetyError(
                        f"final captured target does not match plan: {record['path']}"
                    )
                target_records.append(record)
            non_target_records = []
            for relative in sorted(targets):
                record = observed_record(relative)
                expected = targets[relative]
                if (
                    record["sha256"] != expected["sha256"]
                    or record["size"] != expected["size"]
                ):
                    raise SafetyError(
                        f"final captured protected non-target changed: {relative}"
                    )
                non_target_records.append(record)
            live_manifest = {
                "status": "VERIFIED_NOT_COMMITTED",
                "promotion_id": args.promotion_id,
                "plan_sha256": EXPECTED_PLAN_SHA256,
                "backup_manifest_sha256": backup_manifest_sha,
                "targets": target_records,
                "protected_non_targets": non_target_records,
                "verified_utc": utc_now(),
            }
            live_manifest_path = promotion_dir / "PROMOTED_LIVE_MANIFEST.json"
            atomic_json(live_manifest_path, live_manifest)
            live_manifest_sha = sha256_file(live_manifest_path)
            atomic_json(
                state_path,
                {
                    "status": "COMMIT_READY",
                    "live_manifest_sha256": live_manifest_sha,
                    "utc": utc_now(),
                },
            )

            # The atomic marker rename is the irreversible commit boundary.
            # Before it, every failure rolls back. After it, no failure does.
            atomic_bytes(
                recovery,
                (
                    "PASS\n"
                    f"plan_sha256={EXPECTED_PLAN_SHA256}\n"
                    f"backup_manifest_sha256={backup_manifest_sha}\n"
                    f"live_manifest_sha256={live_manifest_sha}\n"
                ).encode(),
            )
            os.replace(recovery, commit_marker)
            committed = True
            postcommit_warnings = []
            try:
                fsync_dir(promotion_dir)
            except BaseException as commit_fsync_error:
                postcommit_warnings.append(f"commit-marker directory fsync: {commit_fsync_error!r}")
            try:
                result = {
                    "status": "PASS",
                    "promotion_id": args.promotion_id,
                    "plan_sha256": EXPECTED_PLAN_SHA256,
                    "backup_manifest_sha256": backup_manifest_sha,
                    "live_manifest_sha256": live_manifest_sha,
                    "n_promoted": 21,
                    "n_protected_non_targets_unchanged": 5,
                    "completed_utc": utc_now(),
                }
                atomic_json(promotion_dir / "PROMOTION_RESULT.json", result)
                atomic_json(state_path, {"status": "COMMITTED", "utc": utc_now()})
            except BaseException as finalize_error:
                postcommit_warnings.append(f"metadata finalization: {finalize_error!r}")
            try:
                cleanup_paths([item["restore_tmp"] for item in prepared])
            except BaseException as cleanup_error:
                postcommit_warnings.append(f"restore-temp cleanup: {cleanup_error!r}")
            print("PROMOTION COMPLETE: 21/21 destinations match the reviewed staged SHA-256.")
            print(f"Recovery directory: {promotion_dir.relative_to(REPO)}")
            print(f"Backup manifest SHA-256: {backup_manifest_sha}")
            for warning in postcommit_warnings:
                print(f"POST-COMMIT WARNING (do not auto-rollback): {warning}", file=sys.stderr)
            print("No downstream job was run. Do not use unrestricted samap_all.")
        except BaseException:
            # Resolve the narrow signal/syscall ambiguity: a handled signal can
            # arrive after os.replace() has committed but before the Python flag
            # assignment. A valid commit marker is authoritative and forbids
            # rollback even if the in-memory flag is still false.
            if mutation_started and not committed and commit_marker.exists():
                expected_marker = (
                    "PASS\n"
                    f"plan_sha256={EXPECTED_PLAN_SHA256}\n"
                    f"backup_manifest_sha256={backup_manifest_sha}\n"
                    f"live_manifest_sha256={live_manifest_sha}\n"
                ).encode()
                if commit_marker.read_bytes() != expected_marker:
                    atomic_bytes(
                        recovery,
                        b"RECOVERY REQUIRED: ambiguous promotion commit marker.\n",
                    )
                    atomic_json(
                        state_path,
                        {
                            "status": "RECOVERY_REQUIRED",
                            "error": "invalid commit marker after commit-window exception",
                            "utc": utc_now(),
                        },
                    )
                    raise SafetyError(
                        f"commit marker exists but failed validation: {commit_marker}"
                    )
                fsync_dir(promotion_dir)
                committed = True
                print(
                    "Promotion commit marker is authoritative; refusing automatic rollback.",
                    file=sys.stderr,
                )
            if mutation_started and not committed:
                try:
                    auto_rollback(prepared, journal, state_path, recovery, targets)
                    cleanup_paths(
                        [item["promote_tmp"] for item in prepared]
                        + [item["restore_tmp"] for item in prepared]
                    )
                    print("Promotion failed; automatic rollback restored all 21 old hashes.", file=sys.stderr)
                except BaseException as rollback_error:
                    atomic_json(
                        state_path,
                        {"status": "RECOVERY_REQUIRED", "error": repr(rollback_error), "utc": utc_now()},
                    )
                    print(
                        f"CRITICAL: automatic rollback did not complete; preserve {recovery}",
                        file=sys.stderr,
                    )
            elif not mutation_started:
                cleanup_paths(
                    [item["promote_tmp"] for item in prepared]
                    + [item["restore_tmp"] for item in prepared]
                )
                atomic_json(
                    state_path,
                    {"status": "ABORTED_BEFORE_MUTATION", "utc": utc_now()},
                )
                files_dir = promotion_dir / "files"
                if files_dir.exists():
                    try:
                        shutil.rmtree(files_dir)
                        fsync_dir(promotion_dir)
                    except BaseException as cleanup_error:
                        print(
                            f"PRE-MUTATION CLEANUP WARNING: {cleanup_error!r}",
                            file=sys.stderr,
                        )
            raise
    finally:
        lock.close()


def load_backup_manifest(promotion_dir, expected_sha):
    path = promotion_dir / "BACKUP_MANIFEST.json"
    regular_single_link(path)
    actual = sha256_file(path)
    if actual != expected_sha:
        raise SafetyError(f"backup manifest SHA-256 mismatch: {actual}")
    payload = json.loads(path.read_text())
    if payload.get("plan_sha256") != EXPECTED_PLAN_SHA256:
        raise SafetyError("backup manifest is bound to another plan")
    if len(payload.get("entries", [])) != 21:
        raise SafetyError("backup manifest must contain 21 entries")
    return payload


def verify_backup_rows(mappings, manifest):
    by_order = {item["order"]: item for item in manifest["entries"]}
    if len(by_order) != 21:
        raise SafetyError("backup manifest contains duplicate ordinals")
    rows = []
    for mapping in mappings:
        item = by_order.get(mapping["order"])
        if item is None:
            raise SafetyError("backup manifest is missing an ordinal")
        expected = {
            "source": mapping["source"]["path"],
            "destination": mapping["destination"]["path"],
            "old_sha256": mapping["destination"]["old"]["sha256"],
            "promoted_sha256": mapping["source"]["sha256"],
        }
        for key, value in expected.items():
            if item.get(key) != value:
                raise SafetyError(f"backup manifest mapping mismatch at {mapping['order']}: {key}")
        backup = require_relative(item["backup"])
        info = regular_single_link(backup)
        if info.st_size != mapping["destination"]["old"]["size"]:
            raise SafetyError(f"backup size mismatch: {backup}")
        if sha256_file(backup) != mapping["destination"]["old"]["sha256"]:
            raise SafetyError(f"backup SHA-256 mismatch: {backup}")
        rows.append((mapping, item, backup))
    return rows


def rollback_preflight(plan, mappings, args):
    if not ID_RE.fullmatch(args.promotion_id):
        raise SafetyError("invalid promotion id")
    if args.confirm != f"ROLLBACK:{args.promotion_id}":
        raise SafetyError(f"confirmation must be ROLLBACK:{args.promotion_id}")
    if not re.fullmatch(r"[0-9a-f]{64}", args.expected_backup_manifest_sha256):
        raise SafetyError("expected backup manifest SHA-256 must be 64 lowercase hex characters")
    promotion_dir = BACKUP_BASE / args.promotion_id
    reject_symlink_components(promotion_dir)
    manifest = load_backup_manifest(promotion_dir, args.expected_backup_manifest_sha256)
    critical = {
        item["path"]: item for item in manifest.get("critical_recovery_tools", [])
    }
    helper_record = critical.get("tools/migration_promotion.py")
    if helper_record is None:
        raise SafetyError("backup manifest does not bind the recovery helper")
    executing_helper = Path(__file__).resolve()
    if sha256_file(executing_helper) != helper_record["sha256"]:
        raise SafetyError(
            "executing recovery helper differs from the version preserved by promotion; "
            "invoke the recovery_tools copy"
        )
    rows = verify_backup_rows(mappings, manifest)
    target_rows = manifest.get("protected_non_targets", [])
    targets = {item["path"]: item for item in target_rows}
    if len(target_rows) != 5 or len(targets) != 5:
        raise SafetyError("backup manifest must bind exactly five protected non-targets")
    verify_non_targets(targets)
    states = []
    for mapping, item, backup in rows:
        destination = require_relative(mapping["destination"]["path"])
        regular_single_link(destination)
        state = classify(
            destination,
            mapping["destination"]["old"],
            mapping["source"]["sha256"],
        )
        if state == "conflict":
            raise SafetyError(f"rollback conflict; destination is neither old nor promoted: {destination}")
        states.append((mapping, item, backup, destination, state))
    return promotion_dir, states, targets


def do_rollback(plan, mappings, args):
    lock = acquire_lock()
    try:
        promotion_dir, states, targets = rollback_preflight(plan, mappings, args)
        prepared = []
        try:
            # Prepare every restore before the first mutation. This also handles
            # a destination changing from old to promoted during preparation.
            for mapping, item, backup, destination, state in states:
                old = mapping["destination"]["old"]
                restore_tmp = destination.parent / (
                    f".{destination.name}.manual-restore.{args.promotion_id}.{os.getpid()}.tmp"
                )
                copy_independent(
                    backup, restore_tmp, int(old["mode"], 8), old["uid"], old["gid"]
                )
                if sha256_file(restore_tmp) != old["sha256"]:
                    raise SafetyError(f"manual rollback temp verification failed: {restore_tmp}")
                prepared.append((mapping, destination, restore_tmp))
        except BaseException:
            cleanup_paths([row[2] for row in prepared])
            raise

        # Catch changes made while the restore files were being prepared.
        refreshed = {}
        for mapping, item, backup, destination, state in states:
            refreshed[mapping["order"]] = classify(
                destination,
                mapping["destination"]["old"],
                mapping["source"]["sha256"],
            )
            if refreshed[mapping["order"]] == "conflict":
                cleanup_paths([row[2] for row in prepared])
                raise SafetyError(f"rollback conflict after preparation: {destination}")
        verify_non_targets(targets)

        recovery = promotion_dir / "RECOVERY_REQUIRED"
        journal = promotion_dir / "ROLLBACK_JOURNAL.jsonl"
        state_path = promotion_dir / "ROLLBACK_STATE.json"
        atomic_bytes(
            recovery,
            f"Manual rollback in progress: {args.promotion_id}\n".encode(),
        )
        atomic_json(state_path, {"status": "ROLLING_BACK", "utc": utc_now()})
        by_order = {mapping["order"]: (mapping, destination, tmp) for mapping, destination, tmp in prepared}
        for mapping, item, backup, destination, state in reversed(states):
            current = classify(
                destination,
                mapping["destination"]["old"],
                mapping["source"]["sha256"],
            )
            if current == "old":
                append_journal(
                    journal,
                    {"event": "RESTORE_ALREADY_OLD", "order": mapping["order"], "path": mapping["destination"]["path"], "utc": utc_now()},
                )
                continue
            if current not in {"promoted", "old_metadata_drift"}:
                raise SafetyError(f"rollback conflict during apply: {destination}")
            restore_tmp = by_order[mapping["order"]][2]
            append_journal(
                journal,
                {"event": "RESTORE_INTENT", "order": mapping["order"], "path": mapping["destination"]["path"], "utc": utc_now()},
            )
            os.replace(restore_tmp, destination)
            fsync_dir(destination.parent)
            if sha256_file(destination) != mapping["destination"]["old"]["sha256"]:
                raise SafetyError(f"manual rollback verification failed: {destination}")
            append_journal(
                journal,
                {"event": "RESTORE_DONE", "order": mapping["order"], "path": mapping["destination"]["path"], "utc": utc_now()},
            )
        for mapping, item, backup, destination, state in states:
            if classify(
                destination,
                mapping["destination"]["old"],
                mapping["source"]["sha256"],
            ) != "old":
                raise SafetyError(f"final rollback verification failed: {destination}")
        verify_non_targets(targets)
        rollback_target_records = []
        for mapping, item, backup, destination, state in states:
            record = observed_record(mapping["destination"]["path"])
            old = mapping["destination"]["old"]
            if (
                record["sha256"] != old["sha256"]
                or record["size"] != old["size"]
                or record["mode"] != old["mode"]
                or record["uid"] != old["uid"]
                or record["gid"] != old["gid"]
                or record["mtime_ns"] != old["mtime_ns"]
            ):
                raise SafetyError(f"final captured rollback target mismatch: {record['path']}")
            rollback_target_records.append(record)
        rollback_non_target_records = []
        for relative in sorted(targets):
            record = observed_record(relative)
            expected = targets[relative]
            if (
                record["sha256"] != expected["sha256"]
                or record["size"] != expected["size"]
            ):
                raise SafetyError(f"final captured rollback non-target changed: {relative}")
            rollback_non_target_records.append(record)
        rollback_live_manifest = {
            "status": "VERIFIED_OLD",
            "promotion_id": args.promotion_id,
            "plan_sha256": EXPECTED_PLAN_SHA256,
            "targets": rollback_target_records,
            "protected_non_targets": rollback_non_target_records,
            "verified_utc": utc_now(),
        }
        rollback_manifest_path = promotion_dir / "ROLLBACK_LIVE_MANIFEST.json"
        atomic_json(rollback_manifest_path, rollback_live_manifest)
        rollback_manifest_sha = sha256_file(rollback_manifest_path)
        result = {
            "status": "PASS",
            "promotion_id": args.promotion_id,
            "plan_sha256": EXPECTED_PLAN_SHA256,
            "backup_manifest_sha256": args.expected_backup_manifest_sha256,
            "rollback_live_manifest_sha256": rollback_manifest_sha,
            "n_restored_or_already_old": 21,
            "completed_utc": utc_now(),
        }
        atomic_json(promotion_dir / "ROLLBACK_RESULT.json", result)
        atomic_json(state_path, {"status": "ROLLED_BACK", "utc": utc_now()})
        # Preserve historical commit evidence but remove its claim to be the
        # current state. STATE.json and PROMOTION_RESULT.json become authoritative.
        committed_marker = promotion_dir / "PROMOTION_COMMITTED"
        if committed_marker.exists():
            os.replace(committed_marker, promotion_dir / "PROMOTION_COMMITTED_HISTORICAL")
        promotion_result = promotion_dir / "PROMOTION_RESULT.json"
        historical_result = promotion_dir / "PROMOTION_RESULT_BEFORE_ROLLBACK.json"
        if promotion_result.exists() and not historical_result.exists():
            os.replace(promotion_result, historical_result)
        atomic_json(
            promotion_result,
            {
                "status": "ROLLED_BACK",
                "promotion_id": args.promotion_id,
                "rollback_live_manifest_sha256": rollback_manifest_sha,
                "completed_utc": utc_now(),
            },
        )
        atomic_json(promotion_dir / "STATE.json", {"status": "ROLLED_BACK", "utc": utc_now()})
        atomic_bytes(promotion_dir / "ROLLBACK_COMPLETE", b"PASS\n")
        recovery.unlink(missing_ok=True)
        fsync_dir(promotion_dir)
        try:
            cleanup_paths([row[2] for row in prepared])
        except BaseException as cleanup_error:
            print(f"POST-ROLLBACK TEMP CLEANUP WARNING: {cleanup_error!r}", file=sys.stderr)
        print("ROLLBACK COMPLETE: 21/21 destinations match the reviewed old SHA-256.")
        print(f"Recovery copies retained in: {promotion_dir.relative_to(REPO)}")
    finally:
        lock.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("promote", "rollback"))
    parser.add_argument("--plan", type=Path, default=PLAN_PATH)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--execute", action="store_true")
    parser.add_argument("--promotion-id")
    parser.add_argument("--expected-git-head")
    parser.add_argument("--expected-backup-manifest-sha256")
    parser.add_argument("--confirm")
    return parser.parse_args()


def main():
    args = parse_args()
    os.chdir(REPO)
    if args.action == "rollback" and args.promotion_id:
        if not ID_RE.fullmatch(args.promotion_id):
            raise SafetyError("invalid promotion id")
        saved_plan = BACKUP_BASE / args.promotion_id / "PROMOTION_PLAN.json"
        plan, mappings = load_saved_plan(saved_plan)
    else:
        plan, mappings = load_plan(args.plan)
    if args.action == "promote" and args.expected_backup_manifest_sha256 is not None:
        raise SafetyError("--expected-backup-manifest-sha256 applies only to rollback")
    if args.action == "rollback" and args.expected_git_head is not None:
        raise SafetyError("--expected-git-head applies only to promotion")
    if not args.execute:
        if args.action == "promote":
            protected, targets = preflight(
                plan,
                mappings,
                expected_git_head=args.expected_git_head,
                strict_live=True,
            )
            old_total, new_total, required, free = ensure_devices_and_space(mappings)
            print("READ-ONLY PROMOTION PREFLIGHT PASS")
            print(f"plan_sha256={EXPECTED_PLAN_SHA256}")
            print(f"git_head={current_git_head()}")
            print(f"objects=21 old_bytes={old_total} staged_bytes={new_total}")
            print(f"protected_non_targets={len(targets)} free_bytes={free} required_bytes={required}")
            print("No file was written. Re-run with reviewed execute arguments to mutate.")
        else:
            if not args.promotion_id or not args.expected_backup_manifest_sha256:
                raise SafetyError(
                    "rollback dry-run requires --promotion-id and --expected-backup-manifest-sha256"
                )
            if args.confirm is None:
                args.confirm = f"ROLLBACK:{args.promotion_id}"
            promotion_dir, states, targets = rollback_preflight(plan, mappings, args)
            counts = {"old": 0, "old_metadata_drift": 0, "promoted": 0}
            for row in states:
                counts[row[4]] += 1
            print("READ-ONLY ROLLBACK PREFLIGHT PASS")
            print(f"promotion_dir={promotion_dir.relative_to(REPO)}")
            print(
                f"old={counts['old']} old_metadata_drift={counts['old_metadata_drift']} "
                f"promoted={counts['promoted']} conflicts=0"
            )
            print(f"protected_non_targets={len(targets)}")
            print("No file was written. Re-run with --execute and typed confirmation to mutate.")
        return

    if args.action == "promote":
        required = (args.promotion_id, args.expected_git_head, args.confirm)
        if any(value is None for value in required):
            raise SafetyError(
                "promotion execute requires --promotion-id, --expected-git-head, and --confirm"
            )
        do_promote(plan, mappings, args)
    else:
        required = (
            args.promotion_id,
            args.expected_backup_manifest_sha256,
            args.confirm,
        )
        if any(value is None for value in required):
            raise SafetyError(
                "rollback execute requires --promotion-id, "
                "--expected-backup-manifest-sha256, and --confirm"
            )
        do_rollback(plan, mappings, args)


if __name__ == "__main__":
    try:
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            signal.signal(signum, lambda sig, frame: (_ for _ in ()).throw(KeyboardInterrupt()))
        main()
    except (SafetyError, subprocess.CalledProcessError, KeyboardInterrupt) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
