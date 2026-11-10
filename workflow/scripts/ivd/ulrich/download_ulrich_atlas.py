#!/usr/bin/env python
"""Fetch the Ulrich et al. 2024 stromal atlas (CZ CELLxGENE dataset c93622aa-a356-4e29-b45a-c0a532c65d10) and check its SHA-256.

The file is 1,283,567,526 bytes; its SHA-256 is pinned in config (ivd.ulrich.sha256) and in the Ulrich DE config.
If --mirror names an existing copy with the right checksum (the lab's netscratch copy), it is copied instead of
downloaded; otherwise the file is downloaded with curl, resuming a partial download. Either way the checksum is
verified before the output is put in place.

Usage: python download_ulrich_atlas.py --url URL --sha256 HEX --out PATH [--mirror PATH]
"""
import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path


def sha256(path: Path, chunk: int = 16 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--url", required=True)
    ap.add_argument("--sha256", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--mirror", type=Path, default=None, help="optional local copy to use instead of downloading")
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    partial = args.out.with_suffix(args.out.suffix + ".partial")

    if args.mirror and args.mirror.is_file():
        got = sha256(args.mirror)
        if got == args.sha256:
            print(f"mirror {args.mirror} has the pinned checksum; copying", flush=True)
            shutil.copyfile(args.mirror, partial)
        else:
            print(f"mirror {args.mirror} has sha256 {got}, not the pinned {args.sha256}; downloading instead", flush=True)
    if not partial.exists():
        cmd = ["curl", "-fL", "--retry", "5", "--retry-delay", "30", "-C", "-", "-o", str(partial), args.url]
        print(" ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
    got = sha256(partial)
    if got != args.sha256:
        partial.unlink(missing_ok=True)
        sys.exit(f"checksum mismatch: got {got}, expected {args.sha256}")
    partial.rename(args.out)
    print(f"wrote {args.out} ({args.out.stat().st_size:,} bytes), sha256 {got}")


if __name__ == "__main__":
    main()
