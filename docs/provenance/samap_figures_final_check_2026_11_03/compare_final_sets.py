"""Compare every output of a SAMap figure set with the same file of an earlier render of that set.

Used to check the set drawn from the deposited final objects by the committed rules (tag final, 3 Nov 2026) against
the scratch render of 2 Nov 2026 behind the Drive figures. Each file under the old results folder is matched by its
path under the new one (the graph block moved from results/samap_figures/graph_block to <tag>/graph_block) and
compared by type: PNG pixels; SVG text after removing dates and matplotlib's hashed ids, plus the embedded rasters
decoded to pixels; CSV and parquet values; the sparse graph block entry by entry; other text byte by byte.
PDF and provenance JSON carry run dates and paths and are listed without comparison.
Writes one JSON report (--out) and prints its summary.
"""
import argparse
import base64
import difflib
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from PIL import Image


Image.MAX_IMAGE_PIXELS = None

# Old path prefix -> new path prefix, for files that changed folder between the two renders.
MOVED = {"samap_figures/graph_block/": "samap_figures/final/graph_block/"}


def png_report(old: Path, new: Path) -> dict:
    a = np.asarray(Image.open(old).convert("RGBA"))
    b = np.asarray(Image.open(new).convert("RGBA"))
    if a.shape != b.shape:
        return {"identical": False, "shape_old": list(a.shape), "shape_new": list(b.shape)}
    differ = np.any(a != b, axis=-1)
    report = {"identical": not differ.any(), "n_pixels_differ": int(differ.sum())}
    if differ.any():
        ys, xs = np.nonzero(differ)
        # Bounding box of the changed pixels as fractions of width and height, origin top left.
        report["changed_box"] = [round(xs.min() / a.shape[1], 3), round(ys.min() / a.shape[0], 3),
                                 round(xs.max() / a.shape[1], 3), round(ys.max() / a.shape[0], 3)]
    return report


def svg_parts(path: Path) -> tuple[list[str], list[np.ndarray]]:
    """SVG lines with dates and hashed ids normalized, and the embedded PNG rasters as pixel arrays."""
    text = path.read_text()
    rasters = []

    def take_raster(match):
        rasters.append(np.asarray(Image.open(io.BytesIO(base64.b64decode(match.group(1)))).convert("RGBA")))
        return f"data:image/png;base64,RASTER{len(rasters) - 1}"

    text = re.sub(r"data:image/png;base64,([A-Za-z0-9+/=\s]+)", take_raster, text)
    text = re.sub(r"<dc:date>.*?</dc:date>", "<dc:date/>", text, flags=re.S)

    # matplotlib hashes element ids per run; number them in order of first appearance instead.
    ids = list(dict.fromkeys(re.findall(r'id="([^"]+)"', text)))
    for ident in sorted(ids, key=len, reverse=True):
        token = f"ID{ids.index(ident)}"
        text = text.replace(f'id="{ident}"', f'id="{token}"').replace(f"#{ident}", f"#{token}")
    return text.splitlines(), rasters


def svg_report(old: Path, new: Path) -> dict:
    old_lines, old_rasters = svg_parts(old)
    new_lines, new_rasters = svg_parts(new)
    changed = [line[:200] for line in difflib.unified_diff(old_lines, new_lines, lineterm="", n=0)
               if line[:1] in "+-" and not line.startswith(("+++", "---"))]
    rasters_equal = len(old_rasters) == len(new_rasters) and all(
        a.shape == b.shape and np.array_equal(a, b) for a, b in zip(old_rasters, new_rasters))
    return {"identical": not changed and rasters_equal, "rasters_identical": rasters_equal,
            "n_rasters": len(new_rasters), "changed_lines": changed[:10], "n_changed_lines": len(changed)}


def frame_report(old: pd.DataFrame, new: pd.DataFrame) -> dict:
    report = {"rows_old": len(old), "rows_new": len(new),
              "columns_only_old": sorted(set(old.columns) - set(new.columns)),
              "columns_only_new": sorted(set(new.columns) - set(old.columns)), "columns_differ": {}}
    if len(old) == len(new) and old.index.equals(new.index):
        for column in [c for c in old.columns if c in new.columns]:
            a, b = old[column], new[column]
            if a.dtype.kind in "fiub" and b.dtype.kind in "fiub":
                same = np.array_equal(a.to_numpy(float), b.to_numpy(float), equal_nan=True)
            else:
                same = bool((a.astype(str).to_numpy() == b.astype(str).to_numpy()).all())
            if not same:
                report["columns_differ"][column] = int((a.astype(str).to_numpy() != b.astype(str).to_numpy()).sum())
    else:
        report["index_differs"] = True
    report["identical"] = (not report["columns_only_old"] and not report["columns_only_new"]
                           and not report["columns_differ"] and "index_differs" not in report)
    return report


def edge_share_variants_report(old: Path, new: Path) -> dict:
    """The variants table: rows of the label sets both renders computed must agree; missing sets are listed."""
    a, b = pd.read_csv(old), pd.read_csv(new)
    keys = ["graph", "labels", "direction", "source_compartment"]
    shared = a.merge(b, on=keys, suffixes=("_old", "_new"))
    values = [c for c in a.columns if c not in keys]
    gaps = {c: float(np.nanmax(np.abs(shared[f"{c}_old"] - shared[f"{c}_new"]))) for c in values}
    return {"identical": False if len(a) != len(b) else max(gaps.values()) == 0.0,
            "label_sets_old": sorted(a["labels"].unique()), "label_sets_new": sorted(b["labels"].unique()),
            "rows_shared": len(shared), "max_abs_diff_shared_rows": max(gaps.values())}


def compare(old: Path, new: Path) -> dict:
    suffix = old.suffix
    if suffix == ".png":
        return png_report(old, new)
    if suffix == ".svg":
        return svg_report(old, new)
    if old.name == "samap_compartment_edge_share_all.csv":
        return edge_share_variants_report(old, new)
    if suffix == ".csv":
        return frame_report(pd.read_csv(old, index_col=0), pd.read_csv(new, index_col=0))
    if suffix == ".parquet":
        return frame_report(pd.read_parquet(old), pd.read_parquet(new))
    if suffix == ".npz":
        a, b = sp.load_npz(old).tocsr(), sp.load_npz(new).tocsr()
        same = a.shape == b.shape and (a != b).nnz == 0
        return {"identical": bool(same), "shape": list(a.shape), "nnz": int(a.nnz)}
    if suffix in (".pdf", ".json"):
        return {"identical": None, "note": "not compared (run dates and paths)"}
    return {"identical": old.read_bytes() == new.read_bytes()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--old", type=Path, required=True, help="results/ folder of the earlier render")
    parser.add_argument("--new", type=Path, required=True, help="results/ folder of the new render")
    parser.add_argument("--subdirs", nargs="+", required=True, help="folders under results/ to compare")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    report = {}
    for subdir in args.subdirs:
        for old_file in sorted((args.old / subdir).rglob("*")):
            if not old_file.is_file():
                continue
            relative = str(old_file.relative_to(args.old))
            moved = next((relative.replace(a, b, 1) for a, b in MOVED.items() if relative.startswith(a)), relative)
            new_file = args.new / moved
            if not new_file.exists():
                report[relative] = {"identical": False, "missing_in_new": moved}
                continue
            report[relative] = compare(old_file, new_file)

    args.out.write_text(json.dumps(report, indent=1, default=str) + "\n")
    summary = {
        "files": len(report),
        "identical": sorted(k for k, v in report.items() if v.get("identical") is True),
        "not_compared": sorted(k for k, v in report.items() if v.get("identical") is None),
        "different": {k: v for k, v in report.items() if v.get("identical") is False},
    }
    print(json.dumps({"files": summary["files"], "n_identical": len(summary["identical"]),
                      "n_not_compared": len(summary["not_compared"]), "different": summary["different"]},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
