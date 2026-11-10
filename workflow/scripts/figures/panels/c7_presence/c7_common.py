"""Inputs, checks, style and provenance shared by the C7-fibroblast figure scripts (rules in
workflow/rules/41_fig_c7_presence.smk).

Imported by the scripts in this directory; Python puts a script's own directory on the import path, so no path setup
is needed. Nothing is read or written on import.
"""
import datetime
import hashlib
import json
import platform
import sys
from importlib.metadata import version
from pathlib import Path

import h5py
import matplotlib as mpl
import numpy as np
import pandas as pd
from anndata.io import read_elem
from matplotlib import font_manager


# Inputs ----

def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_sha256(path: Path, expected: str, role: str) -> str:
    """The file's sha256, or exit when it is not the one config/c7_presence.yaml pins."""
    got = sha256_of(Path(path))
    if got != expected:
        sys.exit(f"{role}: {path} has sha256 {got}; config/c7_presence.yaml pins {expected}")
    return got


def read_obs(path: Path, columns=None) -> pd.DataFrame:
    """obs of an .h5ad (h5py, read-only; the matrix is never opened), index as str."""
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
    obs.index = obs.index.astype(str)
    if columns is not None:
        missing = sorted(set(columns) - set(obs.columns))
        if missing:
            raise KeyError(f"{path} obs lacks columns {missing}")
        obs = obs[list(columns)]
    return obs


def read_obsm(path: Path, key: str) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        return np.asarray(read_elem(handle["obsm"][key]))


def gene_columns(path: Path, genes, element: str = "X") -> dict:
    """Dense columns of a CSR matrix (X or layers/<name>) for the named genes, in one pass over the indices."""
    with h5py.File(path, "r") as handle:
        var_index = read_elem(handle["var"]).index.astype(str)
        group = handle[element]
        if not isinstance(group, h5py.Group) or group.attrs.get("encoding-type") != "csr_matrix":
            raise TypeError(f"{path}:{element} is {group.attrs.get('encoding-type')}, expected csr_matrix")
        n = int(group.attrs["shape"][0])
        indptr = group["indptr"][:]
        indices = group["indices"][:]
        data = group["data"][:]
    rows = np.repeat(np.arange(n), np.diff(indptr))
    out = {}
    for gene in genes:
        hits = np.flatnonzero(var_index == gene)
        if len(hits) != 1:
            raise KeyError(f"{path}: gene {gene} found {len(hits)} times in var")
        on = indices == int(hits[0])
        column = np.zeros(n, dtype=float)
        column[rows[on]] = data[on]
        out[gene] = column
    return out


def transfer(C, values: np.ndarray, source: np.ndarray) -> np.ndarray:
    """Connectivity-weighted mean of the source cells' values over each cell's neighbours in graph C.

    The operator of the pseudotime transfer (transfer_label_staining.project): numerator C @ value, denominator
    C @ indicator, over the source cells that have a value; NaN where no neighbour has one.
    """
    has = source & np.isfinite(values)
    num = C.dot(np.where(has, values, 0.0))
    den = C.dot(has.astype(float))
    out = np.full(C.shape[0], np.nan)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


# Style ----

# Helvetica Neue on macOS (the version of record); on the cluster Nimbus Sans, URW's Helvetica clone with the same
# metrics (installed with the ghostscript base-35 fonts), so the layout matches; DejaVu Sans last. SVG text stays
# live and names Helvetica Neue first.
FONT_FAMILY = ["Helvetica Neue", "Helvetica", "Nimbus Sans", "DejaVu Sans"]
MACOS_HELVETICA_NEUE = Path("/System/Library/Fonts/HelveticaNeue.ttc")
# Fixed SVG element ids and no creation dates, so a rerun on the same inputs writes the same bytes.
SAVE_RC = {"svg.fonttype": "none", "svg.hashsalt": "xmens", "pdf.fonttype": 42}
RASTER_DPI = 600   # scatter layers rasterized inside the PDF and SVG
PNG_DPI = 300


def set_font() -> str:
    """Register Helvetica Neue where it exists and set the family list; returns the face that will be used."""
    if MACOS_HELVETICA_NEUE.exists():
        font_manager.fontManager.addfont(str(MACOS_HELVETICA_NEUE))
    mpl.rcParams.update({"font.family": FONT_FAMILY, **SAVE_RC})
    return font_manager.FontProperties(family=FONT_FAMILY).get_name()


def save_figure(fig, paths: dict) -> None:
    """PDF and SVG (live text, rasterized dots) and a PNG preview, tight bounding box, white page."""
    with mpl.rc_context(SAVE_RC):
        fig.savefig(paths["pdf"], dpi=RASTER_DPI, bbox_inches="tight", metadata={"CreationDate": None})
        fig.savefig(paths["svg"], dpi=RASTER_DPI, bbox_inches="tight", metadata={"Date": None})
        fig.savefig(paths["png"], dpi=PNG_DPI, bbox_inches="tight")


# Provenance ----

def file_record(path: Path, role: str, sha256: str = None) -> dict:
    """Path, size, mtime and sha256 of one input file (sha256 passed in when already computed)."""
    path = Path(path)
    stat = path.stat()
    return {
        "role": role,
        "path": str(path.resolve()),
        "bytes": stat.st_size,
        "mtime": datetime.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
        "sha256": sha256 or sha256_of(path),
    }


def base_record(script: Path, inputs: dict, outputs: dict, known_sha256: dict = None) -> dict:
    """Outputs, script and its md5, the command, inputs with checksums, host and versions."""
    known_sha256 = known_sha256 or {}
    versions = {"python": platform.python_version()}
    for package in ("numpy", "pandas", "scipy", "anndata", "scanpy", "matplotlib", "igraph", "leidenalg", "h5py"):
        try:
            versions[package] = version(package)
        except Exception:
            pass
    return {
        "outputs": {kind: str(Path(path).resolve()) for kind, path in outputs.items()},
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "host": platform.node(),
        "script": str(Path(script).resolve()),
        "script_md5": hashlib.md5(Path(script).read_bytes()).hexdigest(),
        "command": " ".join(sys.argv),
        "inputs": [file_record(Path(path), role, known_sha256.get(role)) for role, path in inputs.items()],
        "versions": versions,
    }


def write_json(path: Path, record: dict) -> None:
    Path(path).write_text(json.dumps(record, indent=1, default=str))


class Checks:
    """Named comparisons of a recomputed value with the one it must equal; fail() exits when any failed."""

    def __init__(self, log=print):
        self.rows, self.log = [], log

    def equal(self, claim, got, want):
        ok = got == want
        self.rows.append({"claim": claim, "got": str(got), "want": str(want), "ok": bool(ok)})
        self.log(f"[{'PASS' if ok else 'FAIL'}] {claim}: {got} vs {want}")
        return ok

    def within(self, claim, got, tolerance):
        ok = bool(np.isfinite(got) and abs(float(got)) <= tolerance)
        self.rows.append({"claim": claim, "got": f"{float(got):.3g}", "want": f"<= {tolerance:g}", "ok": ok})
        self.log(f"[{'PASS' if ok else 'FAIL'}] {claim}: {float(got):.3g} (tolerance {tolerance:g})")
        return ok

    def failed(self):
        return [row for row in self.rows if not row["ok"]]

    def fail_if_any(self, where):
        bad = self.failed()
        if bad:
            sys.exit(f"{where}: {len(bad)} check(s) failed: " + "; ".join(row["claim"] for row in bad))
