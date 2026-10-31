"""Vocabulary, colours, display transform and output helpers shared by the human whole-tissue panels.

Imported by the panel_*.py scripts of the rules in workflow/rules/91_panels_human_tissue.smk
(Fig. 5B, 5C, S25, S26). Reads only the plot table written by extract_human_tissue.py.
"""

import datetime
import hashlib
import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.text

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager


# Vocabulary ----

# Cell classes (obs cell_type) in the order of the Fig. 5B legend.
CLASS_ORDER = [
    "endometrial_fibroblast",
    "decidual_cell",
    "epithelial_cell",
    "immune_cell",
    "endothelial_cell",
    "perivascular_cell",
    "smooth_muscle_cell",
]

# Abbreviations of the Fig. 5 legend. S25 writes perivascular cells as "pvc" (its own legend).
ABBREVIATIONS = {
    "endometrial_fibroblast": "eFib",
    "decidual_cell": "dec",
    "epithelial_cell": "epi",
    "immune_cell": "imm",
    "endothelial_cell": "ec",
    "perivascular_cell": "pv",
    "smooth_muscle_cell": "smc",
}

# Row labels of S26, worded as in the published supplement.
FULL_NAMES = {
    "epithelial_cell": "epithelial cells",
    "endometrial_fibroblast": "endometrial\nfibroblasts",
    "decidual_cell": "decidual cells",
    "immune_cell": "immune cells",
    "endothelial_cell": "endothelial cells",
    "perivascular_cell": "perivascular cells",
    "smooth_muscle_cell": "smooth muscle\ncells",
}

# obs timepoint values in display order: left/right in 5B and S26, top/bottom in the S25 table.
SAMPLES = ["secretory", "menstrual"]


# Colours ----

# Class colours sampled from the legend dots of the published Fig. 5B. The release objects carry no
# class palette (uns['label_colors'] holds two scanpy default colours), so the manuscript is the reference.
CLASS_COLORS = {
    "endometrial_fibroblast": "#CD79A7",
    "decidual_cell": "#00B4C5",
    "epithelial_cell": "#D4B900",
    "immune_cell": "#D55E00",
    "endothelial_cell": "#5B7FE6",
    "perivascular_cell": "#8F63CA",
    "smooth_muscle_cell": "#4DAF4A",
}

# Sample colours of the published S25 violins; uns['timepoint_colors'] holds matplotlib defaults.
SAMPLE_COLORS = {"secretory": "#EAB073", "menstrual": "#1E6A37"}

COLOR_SOURCE = (
    "manuscript: class colours sampled from the Fig. 5B legend dots, sample colours from the S25 violins; "
    "the release uns palettes are scanpy/matplotlib defaults and are not used"
)


# Spatial display ----

# Display orientation of each section, registered against the published Fig. 5B: the colour-segmented
# panel was compared with the nuclei under all eight right-angle rotations and mirrors, and both sections
# match best mirrored, the menstrual one as a plain swap of the x and y axes.
# Values: (mirror x before rotating, counterclockwise rotation in degrees).
DISPLAY_TRANSFORM = {
    "secretory": (True, 90),
    "menstrual": (True, 270),
}

SCALE_BAR_UM = 1000

# Spatial panels draw the nuclei that passed spatial QC; embedding outliers are hidden as the release
# rule prescribes (uns['embedding_outlier_rule']: retained in analyses, omitted from embedding and
# spatial plots).
SPATIAL_KEPT = "kept"

SEED = 0  # shuffled draw order, so no class or library systematically covers another


logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("panels_human_tissue")


def set_style() -> None:
    """House fonts (Helvetica, Nimbus Sans on the cluster) and editable text in PDF and SVG."""
    helvetica = Path("/System/Library/Fonts/HelveticaNeue.ttc")
    if helvetica.exists():
        font_manager.fontManager.addfont(str(helvetica))
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

    # Type sizes follow a three-step ladder: labels 7 pt, secondary 6.5 pt, ticks 6 pt.
    plt.rcParams.update({
        "font.family": ["Helvetica Neue", "Helvetica", "Nimbus Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "legend.fontsize": 6.5,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        # A fixed salt makes the SVG element ids, and so the files, identical between runs.
        "svg.hashsalt": "panels_human_tissue",
    })


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def load_cells(path: Path) -> pd.DataFrame:
    """Plot table of one release; fails on a cell class the panels have no colour or label for."""
    cells = pd.read_parquet(path)
    cells["cell_type"] = cells["cell_type"].astype(str)
    cells["timepoint"] = cells["timepoint"].astype(str)

    unknown_classes = set(cells["cell_type"]) - set(CLASS_ORDER)
    unknown_samples = set(cells["timepoint"]) - set(SAMPLES)
    if unknown_classes or unknown_samples:
        raise ValueError(f"unexpected cell_type {sorted(unknown_classes)} or timepoint {sorted(unknown_samples)}")
    return cells


def spatial_cells(cells: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Nuclei drawn in the spatial panels, with display coordinates x_plot/y_plot in um.

    Returns the drawn nuclei in a seeded random order and per-sample counts of every exclusion.
    """
    is_kept = cells["spatial_qc_reason"].astype(str) == SPATIAL_KEPT
    has_xy = cells[["x_um_qc", "y_um_qc"]].notna().all(axis=1)
    if (is_kept != has_xy).any():
        raise ValueError("spatial_qc_reason 'kept' and non-missing x_um_qc/y_um_qc disagree")

    is_outlier = cells["embedding_outlier"].astype(bool)
    drawn = cells[is_kept & ~is_outlier].copy()

    counts = {}
    for sample in SAMPLES:
        in_sample = cells["timepoint"] == sample
        counts[sample] = {
            "nuclei_in_release": int(in_sample.sum()),
            "spatial_qc_reason": cells.loc[in_sample, "spatial_qc_reason"].astype(str).value_counts().to_dict(),
            "kept_hidden_embedding_outlier": int((in_sample & is_kept & is_outlier).sum()),
            "drawn": int((in_sample & is_kept & ~is_outlier).sum()),
        }

    # Each section is centred on its bounding box, then mirrored and rotated for display.
    drawn["x_plot"] = np.nan
    drawn["y_plot"] = np.nan
    for sample, (mirror, degrees_ccw) in DISPLAY_TRANSFORM.items():
        in_sample = drawn["timepoint"] == sample
        x_um = drawn.loc[in_sample, "x_um_qc"].to_numpy()
        y_um = drawn.loc[in_sample, "y_um_qc"].to_numpy()
        x_c = x_um - (x_um.min() + x_um.max()) / 2
        y_c = y_um - (y_um.min() + y_um.max()) / 2
        if mirror:
            x_c = -x_c
        theta = np.deg2rad(degrees_ccw)
        drawn.loc[in_sample, "x_plot"] = x_c * np.cos(theta) - y_c * np.sin(theta)
        drawn.loc[in_sample, "y_plot"] = x_c * np.sin(theta) + y_c * np.cos(theta)

    order = np.random.default_rng(SEED).permutation(len(drawn))
    return drawn.iloc[order], counts


def shared_half_width_um(drawn: pd.DataFrame, pad: float = 1.03) -> float:
    """Half the side of a square that holds every section, so all maps share one spatial scale."""
    spans = drawn.groupby("timepoint")[["x_plot", "y_plot"]].agg(np.ptp).to_numpy()
    return float(spans.max() / 2 * pad)


def draw_scale_bar(ax, half_width_um: float, color: str, label: bool = True) -> None:
    """1 mm bar in the lower-left corner of a map axis."""
    x0 = -half_width_um * 0.94
    y0 = -half_width_um * 0.95
    ax.plot([x0, x0 + SCALE_BAR_UM], [y0, y0], color=color, lw=1.0, solid_capstyle="butt")
    if label:
        ax.text(x0 + SCALE_BAR_UM / 2, y0 + half_width_um * 0.03, "1 mm", color=color, fontsize=6, ha="center", va="bottom")


def check_text_layout(fig) -> None:
    """Fail when a label leaves the figure or two labels overlap; font metrics differ between the Mac
    (Helvetica Neue) and the cluster (Nimbus Sans), so the layout is checked where it is rendered."""
    # Legends place their text at draw time, so the layout is computed before measuring.
    fig.draw_without_rendering()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in texts]

    tolerance_px = 1.0
    outside = [
        text for text, box in boxes
        if box.x0 < -tolerance_px or box.y0 < -tolerance_px
        or box.x1 > fig.bbox.x1 + tolerance_px or box.y1 > fig.bbox.y1 + tolerance_px
    ]
    overlapping = [
        (text_a, text_b)
        for i, (text_a, box_a) in enumerate(boxes)
        for text_b, box_b in boxes[i + 1:]
        if box_a.overlaps(box_b)
    ]
    if outside or overlapping:
        raise ValueError(f"text layout: outside the figure {outside}; overlapping {overlapping}")


def save_panel(fig, stem: Path, facecolor: str) -> list[str]:
    """PDF and SVG keep text editable with points rasterized at 600 dpi; the PNG is 300 dpi."""
    check_text_layout(fig)
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    # Creation dates are left out, so a re-run with the same inputs writes byte-identical files.
    formats = (
        ("pdf", 600, {"CreationDate": None}),
        ("svg", 600, {"Date": None}),
        ("png", 300, None),
    )
    for ext, dpi, metadata in formats:
        path = stem.parent / f"{stem.name}.{ext}"
        fig.savefig(path, dpi=dpi, facecolor=facecolor, metadata=metadata)
        paths.append(str(path))
    plt.close(fig)
    return paths


def resolved_font() -> dict:
    """Font file matplotlib uses for the configured family list (Helvetica Neue on a Mac, Nimbus Sans on the cluster)."""
    path = font_manager.findfont(font_manager.FontProperties(family=plt.rcParams["font.family"]))
    return {"requested": list(plt.rcParams["font.family"]), "resolved": font_manager.get_font(path).family_name, "file": path}


def base_provenance(panel: str, release: str, cells_path: Path, extract_provenance_path: Path, scripts: list[Path]) -> dict:
    extract = json.loads(Path(extract_provenance_path).read_text())
    return {
        "panel": panel,
        "release": release,
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "release_object": extract["input"],
        "release_manifest": extract["release_manifest"],
        "plot_table": {"path": str(Path(cells_path).resolve()), "sha256": sha256(cells_path)},
        "scripts": [{"path": str(Path(s).resolve()), "sha256": sha256(s)} for s in scripts],
        "font": resolved_font(),
        "colors": {"source": COLOR_SOURCE},
    }


def write_json(path: Path, payload: dict) -> None:
    Path(path).write_text(json.dumps(payload, indent=1, default=str))
