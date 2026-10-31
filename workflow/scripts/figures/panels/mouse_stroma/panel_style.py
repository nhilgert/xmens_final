"""Shared style, vocabulary and helpers for the mouse stroma supplementary panels (Figs. S15, S16, S17, S22).

Imported by the plot_*.py scripts next to it (rules: workflow/rules/91_panels_mouse_stroma.smk). Holds the
type sizes, the stage, class and subtype vocabularies with their manuscript labels, the subtype colours, the
physical-section map of the mouse Slide-tags pucks, the nucleus filters for spatial panels (turned by the shared
section rotation table, figures/section_rotation.py), and the save/provenance helpers. Reads and writes nothing
itself.
"""
import datetime
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# The rotation reader is shared by every mouse spatial panel and lives in workflow/scripts/figures/.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import section_rotation  # noqa: E402

# Type ladder (pt), as figure-style: titles and labels, legend text, tick labels.
SIZE_TITLE = 8
SIZE_LEGEND = 7
SIZE_TICK = 6

# Live SVG text names Helvetica Neue first (the manuscript font); the cluster has none of the first two,
# so PNG and PDF fall back to Nimbus Sans, a metric-compatible Helvetica clone.
FONT_FAMILY = ["Helvetica Neue", "Helvetica", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
PNG_DPI = 300

# Release timepoints (vocabulary v5+, +0.5 dpi convention) in time order, with the manuscript labels.
STAGE_LABELS = {
    "uninduced": "uninduced",
    "1.5_dpi_early": "1.5 dpi early",
    "1.5_dpi_late": "1.5 dpi late",
    "2.5_dpi": "2.5 dpi",
    "3.5_dpi": "3.5 dpi",
    "4.5_dpi_bleeding": "4.5 dpi menstrual",
}
STAGES = list(STAGE_LABELS)

# Major classes: release obs['coarse_refined'] -> manuscript abbreviation (Fig. 3 and S12 legends).
# 'vascular cells' are blood and (from v9) lymphatic endothelium; v9 pericytes are 'perivascular cells'.
CLASS_ABBREV = {
    "epithelial cells": "epi",
    "endometrial fibroblasts": "eFib",
    "decidual cells": "dec",
    "immune cells": "imm",
    "vascular cells": "vsc",
    "smooth muscle cells": "smc",
    "perivascular cells": "pv",
    "perimetrium cells": "pm",
}

# Fibroblast/decidual subtypes as stored in the release (vocabulary v4+ codes).
SUBTYPES = ["eFib1", "eFib2", "eFib2i", "eFib3", "eFibP", "dec1", "dec2", "dec3", "decAp"]

# The mouse release objects carry no uns['label_colors'], so the subtype colours are read off the manuscript:
# sampled from the solid bars of Fig. S17B (the same hues as Fig. 4A).
SUBTYPE_COLORS = {
    "eFib1": "#A77207",
    "eFib2": "#00401D",
    "eFib2i": "#127D35",
    "eFib3": "#4694C5",
    "eFibP": "#A0E1B0",
    "dec1": "#762A83",
    "dec2": "#1E1EE2",
    "dec3": "#B11A2D",
    "decAp": "#FCCC3F",
}

# Physical sections: obs['tissue_section'] -> section id. Lanes loaded from one section share its coordinate
# frame and are drawn together (as in plot_pt_spatial_sections.py). The uninduced puck holds two tissue
# pieces; piece B is the section shown in the manuscript (Figs. 4A, S16, S22), so it represents uninduced.
SECTION_OF = {
    "uninduced_mck8a__A": "uninduced_A",
    "uninduced_mck8b__A": "uninduced_A",
    "uninduced_mck8a__B": "uninduced_B",
    "uninduced_mck8b__B": "uninduced_B",
    "1_dpi_early_mck1__A": "1.5_dpi_early",
    "1_dpi_late_mck12a__A": "1.5_dpi_late",
    "1_dpi_late_mck12b__A": "1.5_dpi_late",
    "2_dpi_mck15__A": "2.5_dpi",
    "3_dpi_mck10a__A": "3.5_dpi",
    "3_dpi_mck10b__A": "3.5_dpi",
    "4_dpi_bleeding_mck11a__A": "4.5_dpi_bleeding",
}
SECTION_OF_STAGE = {
    "uninduced": "uninduced_B",
    "1.5_dpi_early": "1.5_dpi_early",
    "1.5_dpi_late": "1.5_dpi_late",
    "2.5_dpi": "2.5_dpi",
    "3.5_dpi": "3.5_dpi",
    "4.5_dpi_bleeding": "4.5_dpi_bleeding",
}

# Display rotation of each stage: resources/figures/mouse_section_rotation.yaml, the table every mouse spatial
# panel reads through figures/section_rotation.py (passed in by the spatial scripts; see placed_nuclei).

# Spatial QC value of nuclei whose x_um_qc/y_um_qc are trusted (release obs['spatial_qc_reason']).
SPATIAL_QC_KEPT = "kept"

# Release labels still awaiting sign-off; panels that include such nuclei say so in their provenance.
PENDING_LABELS = {"neutrophil": "mouse neutrophils, new in v9: label pending PI sign-off"}


def apply_style(dark: bool) -> str:
    """Set the shared rcParams; return the foreground colour (white on black for spatial panels)."""
    foreground = "white" if dark else "black"
    background = "black" if dark else "white"
    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": SIZE_TITLE,
        "axes.titlesize": SIZE_TITLE,
        "axes.labelsize": SIZE_TITLE,
        "legend.fontsize": SIZE_LEGEND,
        "xtick.labelsize": SIZE_TICK,
        "ytick.labelsize": SIZE_TICK,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.facecolor": background,
        "axes.facecolor": background,
        "savefig.facecolor": background,
        "text.color": foreground,
        "axes.labelcolor": foreground,
        "axes.edgecolor": foreground,
        "xtick.color": foreground,
        "ytick.color": foreground,
    })
    return foreground


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_table(table_path: Path, table_prov_path: Path) -> tuple[pd.DataFrame, dict]:
    """Per-nucleus plot table and the input record (table sha256 plus the release object it came from)."""
    table = pd.read_csv(table_path, index_col="cell", low_memory=False)
    table["embedding_outlier"] = table["embedding_outlier"].astype(str).str.lower().map({"true": True, "false": False})
    if table["embedding_outlier"].isna().any():
        raise ValueError(f"{table_path}: embedding_outlier holds values other than True/False")
    table_prov = json.loads(Path(table_prov_path).read_text())
    record = {
        "table": {"path": str(Path(table_path).resolve()), "sha256": sha256_of(table_path)},
        "table_provenance": str(Path(table_prov_path).resolve()),
        "release_object": table_prov["input"],
        "expression": table_prov.get("expression_source"),
        "x_normalization": table_prov.get("x_normalization"),
    }
    return table, record


def placed_nuclei(table: pd.DataFrame, stages: list[str], angles: dict) -> tuple[pd.DataFrame, dict]:
    """Nuclei drawn in spatial panels of the given stages, with display coordinates x_plot/y_plot.

    Spatial panels show nuclei with spatial_qc_reason 'kept' (x_um_qc/y_um_qc set) in the stage's section.
    Embedding outliers are omitted, as the release rule says (uns['embedding_outlier_rule']: retained in all
    analyses, omitted from embedding and spatial plots). angles: display rotation per stage
    (section_rotation.load_rotation_table). Each stage's whole array (all lanes, both uninduced pieces) turns
    as one rigid body before the drawn nuclei are picked from it. Returns the nuclei and per-stage counts of
    each filter.
    """
    counts = {}
    frames = []
    for stage in stages:
        in_stage = table["timepoint"] == stage
        is_kept = in_stage & (table["spatial_qc_reason"] == SPATIAL_QC_KEPT)
        section = table["tissue_section"].map(SECTION_OF)
        if section[is_kept].isna().any():
            unknown = sorted(table.loc[is_kept & section.isna(), "tissue_section"].astype(str).unique())
            raise ValueError(f"kept nuclei of {stage} with an unmapped tissue_section: {unknown}")
        in_section = is_kept & (section == SECTION_OF_STAGE[stage])
        drawn = in_section & ~table["embedding_outlier"]

        # Turn the stage's whole array as one body, then pick the drawn nuclei from it (table order kept).
        array = table.loc[is_kept].copy()
        array["x_plot"], array["y_plot"] = section_rotation.rotate_by_timepoint(
            array["x_um_qc"], array["y_um_qc"], array["timepoint"], angles
        )
        stage_cells = array[drawn.loc[array.index].to_numpy()].copy()
        stage_cells["section"] = SECTION_OF_STAGE[stage]
        frames.append(stage_cells)

        counts[stage] = {
            "section": SECTION_OF_STAGE[stage],
            "nuclei_in_stage": int(in_stage.sum()),
            "spatial_qc_kept": int(is_kept.sum()),
            "kept_in_other_piece": int((is_kept & ~in_section).sum()),
            "embedding_outliers_hidden": int((in_section & table["embedding_outlier"]).sum()),
            "nuclei_drawn": int(drawn.sum()),
            "rotation_deg_ccw": angles[stage],
        }
    return pd.concat(frames), counts


def pending_label_notes(cell_subtypes: pd.Series) -> dict:
    """Counts of nuclei carrying a label that awaits sign-off, with the reason."""
    counts = cell_subtypes.value_counts()
    return {label: {"n": int(counts[label]), "note": note} for label, note in PENDING_LABELS.items() if label in counts}


def marker_area(n_nuclei: int, panel_in: float, coverage: float = 0.45) -> float:
    """Scatter marker area (pt^2) so markers cover a similar share of a panel whatever the section's size."""
    return float(np.clip(coverage * (panel_in * 72) ** 2 / n_nuclei, 0.15, 6.0))


def dot_area(fraction: np.ndarray, dot_max: float, largest_area: float) -> np.ndarray:
    """Dot-plot marker area (pt^2): scanpy's dotplot rule, area = (fraction / dot_max)^1.5 * largest area.

    The manuscript dot plots (Figs. S12B, S15A, S17C) were drawn with this rule; dot_max is the largest
    fraction rounded up to the next 0.1, as scanpy does.
    """
    return (np.clip(fraction / dot_max, 0, 1) ** 1.5) * largest_area


def scanpy_dot_max(max_fraction: float) -> float:
    return float(np.ceil(max_fraction * 10) / 10)


def mixed_style_text(fig, x: float, y: float, parts: list[tuple[str, str]], rotation: float = 0, **kwargs) -> None:
    """Draw consecutive text pieces with their own font style (e.g. an italic gene name) as live text.

    parts: (text, fontstyle) pairs. x, y in figure coordinates are the start of the first piece; pieces
    advance along the reading direction (right for rotation 0, up for rotation 90).
    """
    renderer = fig.canvas.get_renderer()
    x_px, y_px = fig.transFigure.transform((x, y))
    label_id = f"mixed-label-{len(fig.texts)}"  # shared by the pieces so the overlap check skips them
    for text, style in parts:
        artist = fig.text(
            *fig.transFigure.inverted().transform((x_px, y_px)),
            text,
            fontstyle=style,
            rotation=rotation,
            rotation_mode="anchor",
            ha="left",
            va="baseline",
            gid=label_id,
            **kwargs,
        )
        extent = artist.get_window_extent(renderer)
        if rotation == 90:
            y_px += extent.height
        else:
            x_px += extent.width


def text_overlaps(fig) -> int:
    """Count overlapping drawn text boxes (figure-style 9.1); written to provenance, checked by eye.

    Skipped by construction: tick labels of axes that are not drawn, pieces of one mixed-style label (they
    abut), and pairs of tick labels set at an angle other than 0 or 90 degrees, whose axis-aligned boxes
    overlap while the glyphs do not.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    hidden_ticks = set()
    angled_ticks = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tick in axis.get_major_ticks():
                for label in (tick.label1, tick.label2):
                    if not ax.axison or not label.get_visible():
                        hidden_ticks.add(id(label))
                    elif label.get_rotation() % 90:
                        angled_ticks.add(id(label))
    texts = [
        t for t in fig.findobj(matplotlib.text.Text)
        if t.get_visible() and t.get_text().strip() and id(t) not in hidden_ticks
    ]
    boxes = [t.get_window_extent(renderer) for t in texts]
    n_overlaps = 0
    for i, (text_a, box_a) in enumerate(zip(texts, boxes)):
        for text_b, box_b in zip(texts[i + 1:], boxes[i + 1:]):
            same_label = text_a.get_gid() is not None and text_a.get_gid() == text_b.get_gid()
            both_angled = id(text_a) in angled_ticks and id(text_b) in angled_ticks
            if box_a.overlaps(box_b) and not (same_label or both_angled):
                n_overlaps += 1
    return n_overlaps


def save_panel(fig, stem: Path, provenance: dict) -> None:
    """Write <stem>.pdf/.svg/.png (300 dpi) and <stem>.provenance.json with each file's sha256."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    provenance["text_overlaps"] = text_overlaps(fig)
    files = {}
    for ext in ("pdf", "svg", "png"):
        path = stem.parent / f"{stem.name}.{ext}"
        # Scatter layers are rasterized at the PNG resolution inside PDF/SVG; text and axes stay vector.
        fig.savefig(path, dpi=PNG_DPI)
        files[path.name] = sha256_of(path)
    plt.close(fig)
    provenance["written"] = datetime.datetime.now().isoformat(timespec="seconds")
    provenance["matplotlib"] = matplotlib.__version__
    provenance["font_family"] = FONT_FAMILY
    provenance["files"] = files
    (stem.parent / f"{stem.name}.provenance.json").write_text(json.dumps(provenance, indent=1, default=str))
    print(f"wrote {stem}.{{pdf,svg,png,provenance.json}}")
