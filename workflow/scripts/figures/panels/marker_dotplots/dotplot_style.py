"""Fonts, saving and provenance shared by the marker dot plot scripts (91_panels_marker_dotplots.smk).

Every plot is saved as SVG with live text (Helvetica Neue named first, for Illustrator) and a 300 dpi PNG
preview; the review plots of a species are also drawn into one reading PDF. SVG ids and file metadata are
pinned so that a re-render of the same inputs is byte-identical.
"""

import hashlib
import json
from pathlib import Path

import matplotlib
from matplotlib import font_manager
from matplotlib.text import Text


# Helvetica Neue as in the manuscript where installed. The cluster has none and lays text out with Nimbus Sans
# (Helvetica metrics); the SVG still names Helvetica Neue first, so Illustrator shows it.
FONT_FAMILY = ["Helvetica Neue", "Helvetica", "Arial", "Nimbus Sans", "DejaVu Sans"]

# Base size of the earlier dot plots: scanpy draws every label 'small' (8.33 pt) on a 10 pt base.
BASE_FONT_SIZE = 10

# Fixed salt for the ids matplotlib writes into SVG files; with the date left out, renders are reproducible.
SVG_HASH_SALT = "panels_marker_dotplots"
PNG_DPI = 300

# PDF metadata without a creation date; matplotlib's PDF writer adds no document id, so pages are reproducible.
PDF_METADATA = {"CreationDate": None}

# Font family and file settings of every plot; sizes stay with each plot type.
FONT_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": FONT_FAMILY,
    # Live text in SVG and embedded TrueType in PDF, so every label stays editable.
    "svg.fonttype": "none",
    "svg.hashsalt": SVG_HASH_SALT,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def resolved_font() -> str:
    return Path(font_manager.findfont(font_manager.FontProperties(family=FONT_FAMILY))).name


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 24), b""):
            digest.update(chunk)
    return digest.hexdigest()


# Rotated text ----

# Text rotated by 90 degrees and aligned by its rotated box is placed with the metrics of the font that drew
# it; Illustrator then shifts every label by the width difference between Nimbus Sans and Helvetica Neue.
# Anchored at the end that touches the axis, the SVG carries a text-anchor and the label sits right in any
# font. Maps the box alignment of a 90-degree text to the equivalent anchor alignment.
ANCHOR_HA = {"bottom": "left", "baseline": "left", "top": "right", "center": "center", "center_baseline": "center"}
ANCHOR_VA = {"left": "top", "right": "bottom", "center": "center"}


def anchor_rotated_text(fig) -> int:
    """Re-anchor every visible 90-degree text of the figure; returns how many were changed."""
    changed = 0
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text() or text.get_rotation_mode() == "anchor":
            continue
        if round(text.get_rotation()) % 360 != 90:
            continue
        horizontal, vertical = text.get_horizontalalignment(), text.get_verticalalignment()
        text.set_rotation_mode("anchor")
        text.set_horizontalalignment(ANCHOR_HA[vertical])
        text.set_verticalalignment(ANCHOR_VA[horizontal])
        changed += 1
    return changed


# Saving ----

def save_svg_png(fig, out_prefix: Path, bbox_inches="tight") -> list:
    """Write <prefix>.svg and <prefix>.png (300 dpi); returns the paths."""
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    svg, png = Path(f"{out_prefix}.svg"), Path(f"{out_prefix}.png")
    fig.savefig(svg, bbox_inches=bbox_inches, metadata={"Date": None})
    fig.savefig(png, dpi=PNG_DPI, bbox_inches=bbox_inches, metadata={"Software": None})
    return [svg, png]


def write_provenance(path: Path, record: dict) -> None:
    record = {
        **record,
        "matplotlib": matplotlib.__version__,
        "font_used_for_layout": resolved_font(),
        "svg_font_family": FONT_FAMILY,
    }
    Path(path).write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n")
