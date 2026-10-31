"""Style, inputs and provenance shared by the cross-species panels of Fig. 6 (6A, 6B, 6D).

Imported by joint_umap_panel.py and spatial_efib_dec_panel.py in this directory; Python puts a script's own
directory on the import path, so no path setup is needed. Nothing is read or written on import.
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


# Style ----

# Font list of the Fig. 6C notebook, so 6A-6C share one face: Helvetica Neue on macOS (the version of record),
# otherwise DejaVu Sans (the cluster has no Helvetica). The provenance record names the face used.
FONT_FAMILY = ["Helvetica Neue", "Helvetica", "DejaVu Sans"]
MACOS_HELVETICA_NEUE = Path("/System/Library/Fonts/HelveticaNeue.ttc")
TEXT = "white"
BACKGROUND = "black"

# Colors are the manuscript's (legends of Fig. 5B, 6A, 6B and the dots of 6D, read off the figure images).
# The release objects carry no usable palette: the human uns['label_colors'] holds matplotlib defaults
# (#ff7f0e, #2ca02c) and the mouse objects carry none.
CELL_TYPE_COLORS = {"endometrial_fibroblast": "#CC79A7", "decidual_cell": "#00B4C5"}
CELL_TYPE_LABELS = {"endometrial_fibroblast": "eFib", "decidual_cell": "dec"}
SPECIES_COLORS = {"mm": "#E6782A", "hs": "#2A6FDB"}
SPECIES_LABELS = {"mm": "mouse", "hs": "human"}

# Scatter layers are rasterized at RASTER_DPI inside the PDF and SVG (text and axes stay vector, as in 6C);
# the PNG preview is written at PNG_DPI.
RASTER_DPI = 600
PNG_DPI = 300

# Fixed SVG element ids and no creation dates, so a rerun on the same inputs writes the same bytes.
SAVE_RC = {"svg.fonttype": "none", "svg.hashsalt": "xmens", "pdf.fonttype": 42}


def set_font() -> str:
    """Register Helvetica Neue where it exists and set the family list; returns the face that will be used."""
    if MACOS_HELVETICA_NEUE.exists():
        font_manager.fontManager.addfont(str(MACOS_HELVETICA_NEUE))
    mpl.rcParams.update({"font.family": FONT_FAMILY, **SAVE_RC})
    return font_manager.FontProperties(family=FONT_FAMILY).get_name()


def save_panel(fig, paths: dict) -> None:
    """Write the panel as PDF and SVG (live text, rasterized dots) and as a PNG preview, on black."""
    with mpl.rc_context(SAVE_RC):
        fig.savefig(paths["pdf"], dpi=RASTER_DPI, facecolor=BACKGROUND, metadata={"CreationDate": None})
        fig.savefig(paths["svg"], dpi=RASTER_DPI, facecolor=BACKGROUND, metadata={"Date": None})
        fig.savefig(paths["png"], dpi=PNG_DPI, facecolor=BACKGROUND)


# Inputs ----

def load_release_obs(path: Path, columns: list) -> pd.DataFrame:
    """Selected obs columns of a release .h5ad, read with h5py in read-only mode; the matrix is never opened."""
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
    missing = sorted(set(columns) - set(obs.columns))
    if missing:
        raise KeyError(f"{path} obs lacks columns {missing}")
    return obs[columns]


# Provenance ----

def file_record(path: Path, role: str) -> dict:
    """Path, size, mtime and sha256 of one input file."""
    stat = path.stat()
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return {
        "role": role,
        "path": str(path.resolve()),
        "bytes": stat.st_size,
        "mtime": datetime.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
        "sha256": digest.hexdigest(),
    }


def base_record(script: Path, inputs: dict, outputs: dict, font: str) -> dict:
    """What every provenance record holds: outputs, script and its md5, inputs with checksums, font and versions."""
    return {
        "outputs": {kind: str(Path(path).resolve()) for kind, path in outputs.items()},
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "host": platform.node(),
        "script": str(script.resolve()),
        "script_md5": hashlib.md5(script.read_bytes()).hexdigest(),
        "command": " ".join(sys.argv),
        "inputs": [file_record(Path(path), role) for role, path in inputs.items()],
        "font": font,
        "versions": {
            "python": platform.python_version(),
            "matplotlib": mpl.__version__,
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "anndata": version("anndata"),
        },
    }


def write_json(path: Path, record: dict) -> None:
    Path(path).write_text(json.dumps(record, indent=1, default=str))
