#!/usr/bin/env python
"""Write the Ulrich DE config for this checkout from the tracked template.

The template (config_ulrich_dS_vs_C7fib.template.yaml) is the config of the 28 Oct 2026 release with every path
replaced by a placeholder; this script fills them in so the pipeline's own resource and object checks (SHA-256 pins)
run against the repository copies:
  {ROOT}        repository root (absolute)
  {ULRICH_H5AD} the downloaded atlas object
  {UPSTREAM}    folder holding the pinned copy of run_pipeline.py
  {PYTHON}      interpreter the fit runs in (recorded only; the version checks are in the config's runtime block)
  {OUTPUT}      output folder (recorded only; the driver takes --out)

Usage: python render_ulrich_config.py --template T.yaml --root ROOT --ulrich-h5ad PATH --upstream DIR --python PY --output DIR --out CONFIG.yaml
"""
import argparse
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("template", "root", "ulrich-h5ad", "upstream", "python", "output", "out"):
        ap.add_argument(f"--{name}", required=True)
    a = ap.parse_args()
    text = Path(a.template).read_text()
    subs = {"{ROOT}": str(Path(a.root).resolve()), "{ULRICH_H5AD}": str(Path(a.ulrich_h5ad).resolve()),
            "{UPSTREAM}": str(Path(a.upstream).resolve()), "{PYTHON}": a.python, "{OUTPUT}": str(Path(a.output).resolve())}
    for key, value in subs.items():
        text = text.replace(key, value)
    left = [k for k in subs if k in text]
    if left:
        raise SystemExit(f"unreplaced placeholders in the rendered config: {left}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(text)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
