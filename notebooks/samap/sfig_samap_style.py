"""House style, 12-column grid and output files of the SAMap supplementary figure (panels a-g).

Imported by the panel notebooks sfig_samap_<letter>_*.ipynb and by sfig_samap_compose.ipynb in this directory, so
every panel is drawn at its exact grid size and the composite places it where it was designed to sit. Nothing is
read or written on import; save_panel() and write_json() write the paths they are given.
"""
import datetime
import hashlib
import json
import platform
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

# Helvetica Neue where installed (the version of record was drawn with it on macOS), else the first installed
# fallback. The provenance record names the font actually used.
FONT_PREFERENCE = ["Helvetica Neue", "Helvetica", "Nimbus Sans", "Arial", "Liberation Sans", "DejaVu Sans"]

# Compartment colors (panels d and e). Blue and orange stay apart under the common color-vision deficiencies
# and do not collide with viridis, which carries pseudotime in the other panels.
FIBROBLAST = "#3a78b8"
DECIDUAL = "#e0836b"
COMPARTMENT_COLORS = {"endometrial_fibroblast": FIBROBLAST, "decidual_cell": DECIDUAL}

# Pseudotime on one fixed 0-1 viridis scale with values as stored, never re-ranked (as in Fig. 6 and the methods).
PSEUDOTIME_CMAP = mpl.colormaps["viridis"]
PSEUDOTIME_NORM = mpl.colors.Normalize(0, 1)

# Grid arithmetic is done in 300-dpi pixels, the resolution of the panel PNGs; the 600-dpi files re-render the
# same figures.
DPI = 300
PX_PER_MM = DPI / 25.4
GUTTER_MM = 4
TITLE_INSET_MM = 5.2  # from a panel's left edge to its title, leaving room for the letter the composite adds

# 180 mm wide and 12 columns; every panel spans whole columns of one row. Row heights in mm.
LAYOUT = {
    "claim": (
        "Mouse pseudotime carried across the SAMap alignment recovers human pseudotime ordering, "
        "strongly within decidual cells, and mouse pseudotime shuffled before the transfer does not."
    ),
    "width_mm": 180,
    "ncol": 12,
    "row_heights_mm": [60, 55, 84],
    "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 4, "chart_family": "embedding",
         "message": "mouse source values"},
        {"letter": "b", "row": 0, "col": 4, "colspan": 4, "chart_family": "embedding",
         "message": "transferred values"},
        {"letter": "c", "row": 0, "col": 8, "colspan": 4, "chart_family": "embedding",
         "message": "human values"},
        {"letter": "d", "row": 1, "col": 0, "colspan": 4, "chart_family": "100% stacked bars",
         "message": "edges link matching compartments"},
        {"letter": "e", "row": 1, "col": 4, "colspan": 8, "chart_family": "scatter + binned median",
         "message": "order preserved by compartment"},
        {"letter": "f", "row": 2, "col": 0, "colspan": 7, "chart_family": "spatial maps",
         "message": "shuffles keep at most the compartment split"},
        {"letter": "g", "row": 2, "col": 7, "colspan": 5, "chart_family": "null strips",
         "message": "observed vs 1,000 shuffles"},
    ],
}


# Style ----

def choose_font() -> str:
    """First installed family of FONT_PREFERENCE."""
    installed = {entry.name for entry in fm.fontManager.ttflist}
    return next(name for name in FONT_PREFERENCE if name in installed)


def font_families(font: str) -> list:
    """rcParams font.family: Helvetica Neue first, then `font`, then DejaVu Sans.

    SVG text names this list, so it opens in Helvetica Neue wherever that is installed. Where it is not, matplotlib
    lays out and rasterizes the text with `font`; DejaVu Sans supplies glyphs both lack.
    """
    return list(dict.fromkeys(["Helvetica Neue", font, "DejaVu Sans"]))


def apply_style(font: str) -> None:
    """Role-mapped size ladder: titles and axis labels 9 pt, legends 8 pt, ticks 7 pt."""
    mpl.rcParams.update({
        "font.family": font_families(font),
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.grid": False,
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


# Grid ----

def grid_geometry(dpi: int = DPI, gutter_mm: float = GUTTER_MM) -> dict:
    """Figure width, column width, gutter, row heights and row tops, in whole pixels at `dpi`."""
    px_per_mm = dpi / 25.4
    width = int(LAYOUT["width_mm"] * px_per_mm)
    ncol = LAYOUT["ncol"]
    gutter = int(gutter_mm * px_per_mm)
    col_width = (width - gutter * (ncol - 1)) // ncol
    row_heights = [int(h * px_per_mm) for h in LAYOUT["row_heights_mm"]]
    row_tops = [sum(row_heights[:i]) + gutter * i for i in range(len(row_heights))]
    return {"width": width, "col_width": col_width, "gutter": gutter, "row_heights": row_heights,
            "row_tops": row_tops, "height": row_tops[-1] + row_heights[-1]}


def panel_spec(letter: str) -> dict:
    """The LAYOUT entry of one panel."""
    for spec in LAYOUT["panels"]:
        if spec["letter"] == letter:
            return spec
    raise KeyError(f"no panel {letter!r} in LAYOUT")


def panel_size_px(letter: str, dpi: int = DPI, gutter_mm: float = GUTTER_MM) -> tuple:
    """(width, height) of a panel in pixels at `dpi`: whole columns plus the gutters between them."""
    grid = grid_geometry(dpi, gutter_mm)
    spec = panel_spec(letter)
    cols, rows, row = spec["colspan"], spec.get("rowspan", 1), spec["row"]
    width = grid["col_width"] * cols + grid["gutter"] * (cols - 1)
    height = sum(grid["row_heights"][row:row + rows]) + grid["gutter"] * (rows - 1)
    return width, height


def panel_origin_px(letter: str, dpi: int = DPI, gutter_mm: float = GUTTER_MM) -> tuple:
    """(x, y) of a panel's top-left corner in the composite, in pixels at `dpi`."""
    grid = grid_geometry(dpi, gutter_mm)
    spec = panel_spec(letter)
    return spec["col"] * (grid["col_width"] + grid["gutter"]), grid["row_tops"][spec["row"]]


def new_panel(letter: str):
    """Empty figure of the panel's exact grid size; returns (fig, width_mm, height_mm)."""
    width, height = panel_size_px(letter)
    fig = plt.figure(figsize=(width / DPI, height / DPI), dpi=DPI)
    return fig, width / PX_PER_MM, height / PX_PER_MM


def panel_title(fig, width_mm: float, height_mm: float, text: str) -> None:
    """Title at the panel's top left, 1 mm down, inset for the letter."""
    fig.text(TITLE_INSET_MM / width_mm, 1 - 1.0 / height_mm, text, ha="left", va="top", fontsize=9)


# Files ----

def find_repo(repo) -> Path:
    """The repository root: `repo` if given, else the git top level above the working directory."""
    if repo is not None:
        return Path(repo)
    git_root = ["git", "rev-parse", "--show-toplevel"]
    return Path(subprocess.run(git_root, capture_output=True, text=True, check=True).stdout.strip())


def panel_paths(out_dir: Path, letter: str) -> dict:
    """The four files a panel notebook writes."""
    stem = out_dir / f"sfig_samap_panel_{letter}"
    return {
        "png": Path(f"{stem}.png"),
        "png_600dpi": Path(f"{stem}_600dpi.png"),
        "svg": Path(f"{stem}.svg"),
        "provenance": Path(f"{stem}.provenance.json"),
    }


def refuse_overwrite(paths, force: bool) -> None:
    """Stop before any existing output is replaced, unless force."""
    existing = [str(path) for path in paths if Path(path).exists()]
    if existing and not force:
        raise FileExistsError(f"Refusing to overwrite {existing}; set FORCE = True")


def save_panel(fig, letter: str, paths: dict) -> None:
    """Write the panel as a 300-dpi PNG, an SVG with live text and a 600-dpi PNG, at its exact grid size."""
    # 'standard' bounding box: 'tight' would trim the panel and break its fit in the grid.
    with mpl.rc_context({"savefig.bbox": "standard"}):
        fig.savefig(paths["png"], dpi=DPI, transparent=True)
        with mpl.rc_context({"svg.fonttype": "none"}):
            fig.savefig(paths["svg"], dpi=600, transparent=True)
        fig.savefig(paths["png_600dpi"], dpi=600, transparent=True)

    size = Image.open(paths["png"]).size
    if size != panel_size_px(letter):
        raise ValueError(f"panel {letter} was saved at {size} px; its grid size is {panel_size_px(letter)} px")


def check_same_transfer(example_transfers: pd.DataFrame, human_cells: pd.DataFrame) -> None:
    """Fail unless the permutation run transferred the same pseudotime as the transfer of record."""
    observed = example_transfers["observed"].reindex(human_cells.index)
    if not np.allclose(observed, human_cells["pt_transferred"], atol=1e-9, equal_nan=True):
        raise ValueError("example_transfers.csv 'observed' differs from human_cells.csv 'pt_transferred': "
                         "the permutation outputs and the tables come from different runs")


# Provenance ----

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path, role: str) -> dict:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "role": role,
        "bytes": stat.st_size,
        "mtime": datetime.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
        "sha256": sha256(path),
    }


def git_revision(path: Path):
    """Commit of the checkout holding path, and whether path matches it; None outside git."""
    git = ["git", "-C", str(path.parent)]
    try:
        commit = subprocess.run(git + ["rev-parse", "HEAD"], capture_output=True, text=True, check=True)
        status = subprocess.run(git + ["status", "--porcelain", "--", path.name],
                                capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    return {"commit": commit.stdout.strip(), "notebook_matches_commit": status.stdout.strip() == ""}


def provenance(notebook: Path, parameters: dict, inputs: dict, outputs, font: str) -> dict:
    """What every record shares: outputs, notebook and commit, parameters, input files, font and versions."""
    return {
        "outputs": [str(Path(path).resolve()) for path in outputs],
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "host": platform.node(),
        "notebook": str(notebook.resolve()),
        "notebook_md5": hashlib.md5(notebook.read_bytes()).hexdigest(),
        "git": git_revision(notebook),
        "parameters": parameters,
        "inputs": [file_record(Path(path), role) for role, path in inputs.items()],
        "font": font,
        "versions": {
            "python": platform.python_version(),
            "matplotlib": mpl.__version__,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
    }


def write_json(path: Path, record: dict) -> None:
    path.write_text(json.dumps(record, indent=1, default=str))
