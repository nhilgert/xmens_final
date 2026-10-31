"""Palettes, labels and helpers shared by the mouse whole-tissue panels (Fig. 3A-D, Figs. S12-S13).

Imported by the plot_*.py scripts in this directory, which rule panels_mouse_tissue_panel
(workflow/rules/91_panels_mouse_tissue.smk) runs once per panel and release version. The XMens
release objects carry no cell-type colours (uns holds timepoint_colors only, matplotlib's defaults),
so every colour below is the manuscript's, read off the published figure images; colours of the
classes added in vocabulary v9 are marked as new.
"""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import font_manager
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


# Stages ----

STAGE_ORDER = ["uninduced", "1.5_dpi_early", "1.5_dpi_late", "2.5_dpi", "3.5_dpi", "4.5_dpi_bleeding"]

# Each panel family names the stages as its manuscript panel does.
BAR_STAGE_LABELS = {
    "uninduced": "uninduced",
    "1.5_dpi_early": "1.5dpi early",
    "1.5_dpi_late": "1.5dpi late",
    "2.5_dpi": "2.5dpi",
    "3.5_dpi": "3.5dpi",
    "4.5_dpi_bleeding": "4.5dpi\nmenstrual",
}
TABLE_STAGE_LABELS = {
    "uninduced": "uninduced",
    "1.5_dpi_early": "1.5 dpi early",
    "1.5_dpi_late": "1.5 dpi late",
    "2.5_dpi": "2.5 dpi",
    "3.5_dpi": "3.5 dpi",
    "4.5_dpi_bleeding": "4.5 dpi (mens)",
}
VIOLIN_STAGE_LABELS = {
    "uninduced": "uninduced",
    "1.5_dpi_early": "1.5 dpi\nearly",
    "1.5_dpi_late": "1.5 dpi\nlate",
    "2.5_dpi": "2.5 dpi",
    "3.5_dpi": "3.5 dpi",
    "4.5_dpi_bleeding": "4.5 dpi\nmenstrual",
}
GRID_STAGE_LABELS = {
    "uninduced": "uninduced",
    "1.5_dpi_early": "1.5 dpi early",
    "1.5_dpi_late": "1.5 dpi late",
    "2.5_dpi": "2.5 dpi",
    "3.5_dpi": "3.5 dpi",
    "4.5_dpi_bleeding": "4.5 dpi menstrual",
}

# Stage colours of the Fig. S12C-D violins. The release's uns['timepoint_colors'] are matplotlib's
# tab10 defaults and differ, so the manuscript's fills are used.
STAGE_COLORS = {
    "uninduced": "#871844",
    "1.5_dpi_early": "#dc7a5c",
    "1.5_dpi_late": "#dfeaa5",
    "2.5_dpi": "#71b6a2",
    "3.5_dpi": "#f1d99b",
    "4.5_dpi_bleeding": "#685a99",
}


# Major cell classes ----

# Stacking order of Fig. 3B. Colours and abbreviations are Fig. 3A-B's; perivascular_cell (v9) takes
# Fig. 5's colour and "pv" for human perivascular cells, so the class reads the same in both species.
CLASS_ORDER = [
    "endometrial_fibroblast",
    "decidual_cell",
    "epithelial_cell",
    "immune_cell",
    "endothelial_cell",
    "perivascular_cell",
    "smooth_muscle_cell",
    "perimetrial_cell",
]
CLASS_LABELS = {
    "endometrial_fibroblast": "eFib",
    "decidual_cell": "dec",
    "epithelial_cell": "epi",
    "immune_cell": "imm",
    "endothelial_cell": "vsc",
    "perivascular_cell": "pv",
    "smooth_muscle_cell": "smc",
    "perimetrial_cell": "pm",
}
CLASS_COLORS = {
    "endometrial_fibroblast": "#cc79a7",
    "decidual_cell": "#00b4c5",
    "epithelial_cell": "#d4b900",
    "immune_cell": "#d55e00",
    "endothelial_cell": "#5b7fe6",
    "perivascular_cell": "#8e63c9",
    "smooth_muscle_cell": "#4daf4a",
    "perimetrial_cell": "#8c8c8c",
}

# Row order of Fig. S12B and column order of Fig. S13 (the manuscript's order, with the v9
# perivascular class placed after the vasculature it lines).
SUPP_CLASS_ORDER = [
    "epithelial_cell",
    "endometrial_fibroblast",
    "decidual_cell",
    "endothelial_cell",
    "perivascular_cell",
    "immune_cell",
    "smooth_muscle_cell",
    "perimetrial_cell",
]
DOTPLOT_CLASS_LABELS = {**CLASS_LABELS, "smooth_muscle_cell": "SMC"}
GRID_CLASS_TITLES = {
    "epithelial_cell": "epithelial cells",
    "endometrial_fibroblast": "endometrial\nfibroblasts",
    "decidual_cell": "decidual cells",
    "endothelial_cell": "vasculature cells",
    "perivascular_cell": "perivascular\ncells",
    "immune_cell": "immune cells",
    "smooth_muscle_cell": "smooth muscle\ncells",
    "perimetrial_cell": "perimetrial cells",
}

# Fig. S12B marker genes in the manuscript's column order. Rgs5 marks the v9 perivascular class
# (pericytes): detected in 79% of them and in at most 7% of any other class in release 2026_10_30.
DOTPLOT_GENES = [
    "Cdh1", "Pax8", "Pdgfra", "Wnt4", "Prl8a2", "Bmp2", "Cdh3", "Flt1", "Pecam1",
    "Ptprc", "Acta2", "Myh11", "Sox6", "Pgr", "Esr1",
]
PERIVASCULAR_GENES = ["Rgs5"]
PERIVASCULAR_GENES_AFTER = "Pecam1"


# Immune subtypes ----

# Stacking order and colours of Fig. 3D, matched to the manuscript bars by their v8 fractions.
# neutrophil and plasma_cell are v9 labels with new colours: neutrophils a dark red beside the
# myeloid oranges, plasma cells a light green beside the B cells they derive from.
IMMUNE_ORDER = [
    "macrophage",
    "cycling_macrophage",
    "monocyte",
    "dendritic_cell",
    "neutrophil",
    "t_cell",
    "gamma_delta_t_cell",
    "b_cell",
    "plasma_cell",
    "cytotoxic_unk_cell",
    "unk_cell",
    "cycling_unk_cell",
]
IMMUNE_LABELS = {
    "macrophage": "Mφ",
    "cycling_macrophage": "cycling Mφ",
    "monocyte": "Mono",
    "dendritic_cell": "DC",
    "neutrophil": "Neu",
    "t_cell": "T",
    "gamma_delta_t_cell": "γδ T",
    "b_cell": "B",
    "plasma_cell": "plasma",
    "cytotoxic_unk_cell": "cytotoxic uNK",
    "unk_cell": "uNK",
    "cycling_unk_cell": "cycling uNK",
}
IMMUNE_COLORS = {
    "macrophage": "#d55e00",
    "cycling_macrophage": "#f28e5c",
    "monocyte": "#c9a227",
    "dendritic_cell": "#8c5a2b",
    "neutrophil": "#b2182b",
    "t_cell": "#0072b2",
    "gamma_delta_t_cell": "#56b4e9",
    "b_cell": "#009e73",
    "plasma_cell": "#a6d96a",
    "cytotoxic_unk_cell": "#cc79a7",
    "unk_cell": "#7b5ea7",
    "cycling_unk_cell": "#b39ddb",
}
NEW_IN_V9 = ["perivascular_cell", "neutrophil", "plasma_cell"]
PENDING_SIGN_OFF = {"neutrophil": "label approved 30 Oct 2026 for vocabulary v9; PI sign-off pending"}


# Spatial grid (Fig. S13) ----

GRID_HIGHLIGHT = "#ff00ff"
GRID_OTHERS = "#b0b0b0"


# Matplotlib settings ----

# Helvetica as in the manuscript where installed; the cluster has none and falls back to Nimbus
# Sans or DejaVu Sans. Font sizes: 8 pt for labels and titles, 7 pt for ticks and legends.
FONT_FAMILY = ["Helvetica Neue", "Helvetica", "Arial", "Nimbus Sans", "DejaVu Sans"]
FONT_SIZE = 8
FONT_SIZE_SMALL = 7


def apply_panel_style() -> None:
    """The figure-style rcParams (open frame, outward ticks, embedded fonts), text kept as text."""
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": FONT_FAMILY,
        "font.size": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "legend.fontsize": FONT_SIZE_SMALL,
        "xtick.labelsize": FONT_SIZE_SMALL,
        "ytick.labelsize": FONT_SIZE_SMALL,
        "axes.linewidth": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        # Type-42 fonts in PDF and live text in SVG, so Illustrator can edit every label.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def resolved_font() -> str:
    return Path(font_manager.findfont(font_manager.FontProperties(family=FONT_FAMILY))).name


# Command line, inputs and provenance ----

def panel_arguments(description: str) -> argparse.ArgumentParser:
    """Arguments every panel script takes; scripts add their own on top."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--tables-dir", required=True, type=Path, help="panels_mouse_tissue_extract output")
    parser.add_argument("--version", required=True, choices=["v8", "v9"])
    parser.add_argument("--seed", required=True, type=int, help="config seed; orders overplotted points")
    parser.add_argument("--out-prefix", required=True, type=Path)
    return parser


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 24), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_cells(tables_dir: Path) -> pd.DataFrame:
    return pd.read_csv(tables_dir / "cells.csv.gz", index_col="obs_name", low_memory=False)


def base_provenance(panel: str, args: argparse.Namespace, inputs: list, script: str) -> dict:
    """Provenance fields every panel records: inputs with sha256 and the release they came from."""
    extract = json.loads((args.tables_dir / "extract_provenance.json").read_text())
    return {
        "panel": panel,
        "version": args.version,
        "release": extract["release"],
        "vocabulary_schema": extract["vocabulary_schema"],
        "release_inputs": extract["inputs"],
        "plot_inputs": {str(p): sha256_of(p) for p in inputs},
        "script": script,
        "script_sha256": sha256_of(Path(script)),
        "style_sha256": sha256_of(Path(__file__)),
        "seed": args.seed,
    }


def save_panel(fig, out_prefix: Path, provenance: dict, facecolor: str = "white") -> None:
    """Write <prefix>.pdf, .svg, .png (300 dpi) and <prefix>.provenance.json."""
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # Scatter layers are rasterised inside the PDF and SVG at 600 dpi; text and axes stay vector.
    fig.savefig(f"{out_prefix}.pdf", dpi=600, facecolor=facecolor, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{out_prefix}.svg", dpi=600, facecolor=facecolor, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{out_prefix}.png", dpi=300, facecolor=facecolor, bbox_inches="tight", pad_inches=0.02)

    record = {
        **provenance,
        "outputs": [f"{out_prefix}.{ext}" for ext in ("pdf", "svg", "png")],
        "matplotlib": matplotlib.__version__,
        "font": resolved_font(),
        "rendered_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    Path(f"{out_prefix}.provenance.json").write_text(json.dumps(record, indent=2, ensure_ascii=False))


# Direct labels on embeddings ----

def largest_island(points: np.ndarray, link_px: float) -> np.ndarray:
    """Points of the largest island: groups of points joined by gaps shorter than link_px."""
    pairs = cKDTree(points).query_pairs(r=link_px, output_type="ndarray")
    n = len(points)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, island_of = connected_components(graph, directed=False)
    return points[island_of == np.bincount(island_of).argmax()]


# Label rings beyond the island's edge, in display pixels. Labels on the inner rings sit flush
# with their island, as in the manuscript; labels pushed to the outer rings get a leader line.
ADJACENT_RINGS_PX = (0.0, 3.0, 6.0)
LEADER_RINGS_PX = (12.0, 18.0, 26.0, 36.0)
LEADER_COLOR = "#606060"
LEADER_WIDTH_PT = 0.4

# Labels keep this gap from each other (box growth factors), so neighbours never read as one word.
LABEL_GAP_SCALE = (1.25, 1.4)

# An adjacent label must sit at least this much closer to its own group than to any other (px).
AMBIGUITY_MARGIN_PX = 2.0


def _island_core(points_px: np.ndarray, link_px: float):
    """Median and dense core (90% nearest the median) of a group's largest island."""
    island = largest_island(points_px, link_px)
    centre = np.median(island, axis=0)
    distance = np.hypot(*(island - centre).T)
    return centre, island[distance <= np.quantile(distance, 0.9)]


def _distance_to_box(points_px: np.ndarray, box) -> float:
    """Smallest distance from any of the points to the rectangle (0 when a point is inside)."""
    if len(points_px) == 0:
        return np.inf
    dx = np.maximum(np.maximum(box.x0 - points_px[:, 0], 0), points_px[:, 0] - box.x1)
    dy = np.maximum(np.maximum(box.y0 - points_px[:, 1], 0), points_px[:, 1] - box.y1)
    return float(np.hypot(dx, dy).min())


def _segment_hits_box(start: np.ndarray, end: np.ndarray, box, n_samples: int = 32) -> bool:
    t = np.linspace(0, 1, n_samples)[:, None]
    points = start + t * (end - start)
    return bool((
        (points[:, 0] >= box.x0) & (points[:, 0] <= box.x1)
        & (points[:, 1] >= box.y0) & (points[:, 1] <= box.y1)
    ).any())


def _points_near_segment(points_px: np.ndarray, start: np.ndarray, end: np.ndarray, radius_px: float) -> int:
    along = end - start
    t = np.clip(((points_px - start) @ along) / max(along @ along, 1e-9), 0, 1)
    nearest = start + t[:, None] * along
    return int((np.hypot(*(points_px - nearest).T) <= radius_px).sum())


def _segments_cross(a1: np.ndarray, a2: np.ndarray, b1: np.ndarray, b2: np.ndarray) -> bool:
    """Whether segments a1-a2 and b1-b2 properly cross (touching ends do not count)."""
    def side(p, q, r):
        return np.sign((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
    return side(a1, a2, b1) * side(a1, a2, b2) < 0 and side(b1, b2, a1) * side(b1, b2, a2) < 0


def _plan_labels(ax, points_px, groups, placing_order, labels, cores, pad_px, angles_deg, bounds, keep_out):
    """Greedy placement in the given order; returns (total cost, placements) or None if stuck."""
    renderer = ax.figure.canvas.get_renderer()
    to_data = ax.transData.inverted().transform
    placed_boxes = []
    placed_leaders = []
    placements = []
    total = 0.0

    for group in placing_order:
        centre, core = cores[group]
        others = points_px[groups != group]
        best = None
        for ring_px in ADJACENT_RINGS_PX + LEADER_RINGS_PX:
            has_leader = ring_px in LEADER_RINGS_PX
            for angle in np.linspace(0, 2 * np.pi, 16, endpoint=False):
                direction = np.array([np.cos(angle), np.sin(angle)])
                preferred = group in angles_deg and np.isclose(
                    np.cos(angle - np.deg2rad(angles_deg[group])), 1.0
                )
                off_preference = group in angles_deg and not preferred
                edge_px = centre + ((core - centre) @ direction).max() * direction
                anchor_px = edge_px + (pad_px + ring_px) * direction
                ha = "left" if direction[0] > 0.38 else "right" if direction[0] < -0.38 else "center"
                va = "bottom" if direction[1] > 0.38 else "top" if direction[1] < -0.38 else "center"

                text = ax.text(*to_data(anchor_px), labels[group], ha=ha, va=va)
                box = text.get_window_extent(renderer)
                text.remove()

                inside_bounds = (
                    box.x0 >= bounds.x0 and box.x1 <= bounds.x1
                    and box.y0 >= bounds.y0 and box.y1 <= bounds.y1
                )
                spaced_box = box.expanded(*LABEL_GAP_SCALE)
                if not inside_bounds or any(spaced_box.overlaps(other) for other in placed_boxes + keep_out):
                    continue
                # Leaders and labels never cross each other.
                if any(_segment_hits_box(start, end, spaced_box) for start, end in placed_leaders):
                    continue
                if has_leader and any(_segment_hits_box(anchor_px, edge_px, other) for other in placed_boxes + keep_out):
                    continue
                if has_leader and any(_segments_cross(anchor_px, edge_px, start, end) for start, end in placed_leaders):
                    continue
                covered = int((
                    (points_px[:, 0] >= box.x0) & (points_px[:, 0] <= box.x1)
                    & (points_px[:, 1] >= box.y0) & (points_px[:, 1] <= box.y1)
                ).sum())
                # The label reads as the name of whichever group lies closest to its box; a
                # direction set by hand was judged readable on the rendered panel.
                misleading = (not has_leader) and (not preferred) and (
                    _distance_to_box(others, box)
                    < _distance_to_box(points_px[groups == group], box) + AMBIGUITY_MARGIN_PX
                )
                crossed = _points_near_segment(others, anchor_px, edge_px, 2.0) if has_leader else 0
                cost = (
                    2 * covered
                    + 2 * crossed
                    + 100 * misleading
                    + 30 * has_leader
                    + 50 * off_preference
                    + 0.01 * ring_px
                    + 1e-3 * np.hypot(*(anchor_px - centre))
                )
                if best is None or cost < best["cost"]:
                    best = {
                        "group": group,
                        "cost": cost,
                        "anchor_px": anchor_px,
                        "edge_px": edge_px,
                        "ha": ha,
                        "va": va,
                        "box": spaced_box,
                        "covered": covered,
                        "crossed": crossed,
                        "misleading": bool(misleading),
                        "leader": has_leader,
                    }
        if best is None:
            return None
        placed_boxes.append(best["box"])
        if best["leader"]:
            placed_leaders.append((best["anchor_px"], best["edge_px"]))
        placements.append(best)
        total += best["cost"]
    return total, placements


def place_direct_labels(ax, xy: np.ndarray, groups: np.ndarray, order: list, labels: dict,
                        angles_deg: dict, n_orders: int, seed: int,
                        pad_px: float = 3.0, link_px: float = 6.0,
                        bounds=None, keep_out=()) -> dict:
    """Write one label per group beside its largest island, where it reads as that island's name.

    Candidates sit just outside the island's dense core (the 90% of its points nearest its
    median) in 16 directions, on the adjacent rings and then on the leader rings. Boxes that leave
    the figure, crowd a placed label or touch a placed leader, and leaders that cross a placed
    label, are ruled out. A candidate costs 2 per plotted point under
    its box and 2 per other group's point its leader passes; an adjacent label whose own group is
    not at least AMBIGUITY_MARGIN_PX closer to its box than every other group could name the wrong
    island and costs 100 more, while a leader label costs 30 more (the leader names the island). Ties go to the
    inner ring and the shorter distance. A group listed in angles_deg (degrees counter-clockwise
    from east, one of the 16 directions) prefers that direction: it carries no wrong-island cost
    there and 50 more anywhere else. Labels are planned
    greedily in n_orders placing orders (largest group first, smallest first, then seeded random
    orders) and the cheapest plan is drawn. Labels stay inside bounds (a display-space Bbox;
    default the whole figure) and off every keep_out Bbox (display space; e.g. an inset drawn over
    the axes), and no leader crosses a keep_out Bbox or another leader. Returns the chosen anchor (data coordinates),
    alignment, leader and cost terms per group, for the provenance record.
    """
    to_data = ax.transData.inverted().transform
    points_px = ax.transData.transform(xy)
    present = [g for g in order if (groups == g).any()]
    cores = {g: _island_core(points_px[groups == g], link_px) for g in present}
    by_size = sorted(present, key=lambda g: -(groups == g).sum())
    rng = np.random.default_rng(seed)
    placing_orders = [by_size, by_size[::-1]]
    placing_orders += [list(rng.permutation(by_size)) for _ in range(n_orders - 2)]

    plans = [
        _plan_labels(ax, points_px, groups, placing_order, labels, cores, pad_px, angles_deg,
                     bounds=ax.figure.bbox if bounds is None else bounds, keep_out=list(keep_out))
        for placing_order in placing_orders
    ]
    plans = [plan for plan in plans if plan is not None]
    if not plans:
        raise RuntimeError("no label plan fits; enlarge the figure")
    _, placements = min(plans, key=lambda plan: plan[0])

    chosen = {}
    for best in placements:
        group = best["group"]
        anchor = to_data(best["anchor_px"])
        ax.text(*anchor, labels[group], ha=best["ha"], va=best["va"])
        if best["leader"]:
            # Thin grey leader from the label to the island's edge, starting just off the text.
            gap = best["edge_px"] - best["anchor_px"]
            start_px = best["anchor_px"] + 2.0 * gap / np.hypot(*gap)
            line = to_data(np.vstack([start_px, best["edge_px"]]))
            ax.plot(line[:, 0], line[:, 1], color=LEADER_COLOR, linewidth=LEADER_WIDTH_PT, solid_capstyle="butt")
        chosen[group] = {
            "label": labels[group],
            "anchor": [round(float(v), 3) for v in anchor],
            "ha": best["ha"],
            "va": best["va"],
            "leader_line": best["leader"],
            "points_under_label": best["covered"],
            "leader_points_crossed": best["crossed"],
            "nearer_another_group": best["misleading"],
            "direction_set_by_hand": group in angles_deg,
        }
    return chosen
