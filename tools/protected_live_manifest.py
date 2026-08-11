#!/usr/bin/env python3
"""Capture or verify immutable live-object filesystem metadata."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def protected_paths(repo: Path):
    fixed = [
        repo / f"data/interim/refined_{species}.h5ad" for species in ("human", "mouse")
    ] + [
        repo / f"data/interim/final_{species}.h5ad" for species in ("human", "mouse")
    ] + [
        repo / f"data/interim/final_{species}_fib_dec.h5ad" for species in ("human", "mouse")
    ]
    globbed = (
        list((repo / "data/interim/reannotated").glob("*.h5ad")) +
        list((repo / "results/samap_transfer").glob("**/*")) +
        list((repo / "dataverse/objects").glob("*.h5ad"))
    )
    return sorted({path for path in fixed + globbed if path.is_file()})


def state(repo: Path, hash_content: bool = False):
    rows = []
    for path in protected_paths(repo):
        stat = path.stat()
        row = {
            "path": str(path.relative_to(repo)), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns, "inode": stat.st_ino,
        }
        if hash_content:
            row["sha256"] = sha256(path)
        rows.append(row)
    return rows


ap = argparse.ArgumentParser()
ap.add_argument("--repo", default=".")
ap.add_argument("--capture", default=None)
ap.add_argument("--verify", default=None)
ap.add_argument("--report", default=None)
ap.add_argument("--nonfatal", action="store_true")
ap.add_argument("--hash-content", action="store_true")
args = ap.parse_args()
repo = Path(args.repo).resolve()
if bool(args.capture) == bool(args.verify):
    raise SystemExit("choose exactly one of --capture or --verify")

if args.capture:
    now = state(repo, hash_content=args.hash_content)
    small = [repo / "resources/annotation_vocabulary.yaml",
             repo / "resources/lineage_spec.yaml"]
    small += sorted((repo / "resources/annotations").glob("*.csv"))
    payload = {
        "git_head": subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "protected": now,
        "small_input_sha256": {
            str(path.relative_to(repo)): sha256(path) for path in small if path.is_file()
        },
    }
    out = Path(args.capture)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise FileExistsError(f"refusing to overwrite {out}")
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
else:
    before = json.loads(Path(args.verify).read_text())
    hash_content = args.hash_content or any("sha256" in row for row in before["protected"])
    now = state(repo, hash_content=hash_content)
    old = {row["path"]: row for row in before["protected"]}
    current = {row["path"]: row for row in now}
    changes = []
    for path in sorted(set(old) | set(current)):
        if old.get(path) != current.get(path):
            changes.append({"path": path, "before": old.get(path), "after": current.get(path)})
    payload = {"status": "PASS" if not changes else "FAIL", "changes": changes,
               "n_protected": len(current)}
    if args.report:
        out = Path(args.report)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    raise SystemExit(0 if (not changes or args.nonfatal) else 1)
