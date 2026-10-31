"""Fig. 5E: correlation of SCENIC regulon activity with pseudotime, human (x) against X-Mens (y).

Reads the regulon table staged in the repository (one row per transcription factor with a regulon in both
species; columns human_tf, mouse_tf, rho_human, rho_mouse; written by scenic_quadrant_rho.py) and checks its
sha256 against the value pinned in config. Writes <out>.svg (live text, Helvetica Neue named first),
<out>.png (300 dpi preview) and <out>.provenance.json. Drawn as the manuscript panel: axes crossing at the
origin, dotted lines at |rho| = threshold, quadrant labels in the corners, the manuscript's regulons named.
Rule: panels_programs_5e (workflow/rules/91_panels_programs.smk).
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.markers import MarkerStyle
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, TextArea, VPacker
from matplotlib.patches import Rectangle
from matplotlib.path import Path as MplPath

# Type sizes (pt), measured on the manuscript panel at the 180 mm page width: title, axis names and
# quadrant labels in bold; regulon names and tick labels one step down.
SIZE_LABEL = 7
SIZE_NAME = 6
CAP_HEIGHT_EM = 0.72  # Helvetica cap height; names are centred on their point by their capitals

# Live SVG text names Helvetica Neue first (the manuscript font); the cluster has neither Helvetica, so the
# PNG falls back to Nimbus Sans, a metric-compatible clone.
FONT_FAMILY = ["Helvetica Neue", "Helvetica", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"]
PNG_DPI = 300
SVG_HASH_SALT = "panel_5E"  # fixed SVG element ids: the same table gives the same bytes

# Geometry (inches) of the manuscript panel: a 2.94 in square plot, the title 0.24 in to its left and the
# x-axis name running 11 pt past the arrow tip.
FIG_IN = (3.5, 3.0)
AXES_LEFT_IN = 0.36
AXES_BOTTOM_IN = 0.03
AXES_SIZE_IN = 2.94
TITLE_X_IN = 0.07  # centre of the rotated title
AXIS_NAME_OVERHANG_PT = 11
AXIS_END = 0.95  # rho at the arrow tips
TICKS = [-0.8, -0.2, 0.2, 0.8]
TICK_LENGTH_PT = 3.2  # tick marks cross the axis line, half on each side
TICK_PAD_PT = 2.5  # tick label to axis line
LINE_WIDTH_PT = 0.7
ARROWHEAD_PT = 6.5  # marker size of the arrowheads at both ends of each axis
THRESHOLD_LINE_END = AXIS_END - 0.02  # dotted lines stop just short of the arrow tips

# Quadrant-label chevrons (pt): three per word, each 3 pt wide and 1.1 pt tall, stacked over one cap height.
CHEVRON_WIDTH_PT = 3.0
CHEVRON_RISE_PT = 1.1
CHEVRON_LINE_PT = 0.75

# Classes: coloured when |rho| >= threshold in both species, by the signs of the two rho; grey otherwise.
# Hues sampled from the manuscript panel; down-in-human/up-in-mouse has no regulon there or here, so its
# colour only marks the quadrant label if a regulon ever lands in it.
CLASS_COLORS = {
    "up_both": "#B2182B",
    "down_both": "#2166AC",
    "up_human_down_mouse": "#E08214",
    "down_human_up_mouse": "#5E4FA2",
}
GREY_POINT = "#CCCCCC"  # regulons below the threshold in at least one species
GREY_NAME = "#767676"  # names of named regulons below the threshold: grey, yet 4.5:1 on white
GREY_THRESHOLD = "#BDBDBD"
GREY_LEADER = "#808080"

# Quadrant of each class (sign of human rho, sign of mouse rho) and the corner its label sits in.
QUADRANTS = {
    "up_both": (1, 1, "upper right"),
    "down_human_up_mouse": (-1, 1, "upper left"),
    "down_both": (-1, -1, "lower left"),
    "up_human_down_mouse": (1, -1, "lower right"),
}

# Point areas (pt^2), from the manuscript dot diameters: 2.7 pt coloured, 2.0 pt grey.
AREA_COLOURED = 7.3
AREA_GREY = 4.0

# Name placement: offset of the name's anchor from its point (pt; +x right, +y up) and whether a leader
# line joins them. The side follows the sign of dx (right of the point by default, as in the manuscript).
# Names inside crowded groups sit outside the group with a leader, as JUN and WT1 do in the manuscript.
# Offsets were chosen by eye on the rendered panel with at least 2.5 pt between a name and any other point.
DEFAULT_OFFSET_PT = (3.2, 0.0, False)
NAME_OFFSETS_PT = {
    "WT1": (9.0, 5.0, True),  # right side blocked by FOSL2, upper left by ELF1
    "JUN": (-3.2, 0.0, False),  # right side blocked by XBP1
    "JUNB": (5.0, 9.0, True),  # right side blocked by ETV6, left by ELF2
    "KLF6": (3.2, -2.0, False),  # clears ETV6 above
    "BPTF": (-3.2, 0.0, False),  # left, as in the manuscript
    "IRF1": (3.2, -7.0, True),  # right side blocked by ELK4, left by NFIC
    "ETS1": (3.2, 2.5, False),  # ETS1 and FOXO1 are 4 pt apart: names spread up and down
    "FOXO1": (3.2, -2.5, False),
}


def sha256_of(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def classify(rho_human: np.ndarray, rho_mouse: np.ndarray, threshold: float) -> np.ndarray:
    """Class of each regulon: the quadrant if |rho| >= threshold in both species, else 'below'."""
    strong = (np.abs(rho_human) >= threshold) & (np.abs(rho_mouse) >= threshold)
    conditions = [
        strong & (rho_human > 0) & (rho_mouse > 0),
        strong & (rho_human < 0) & (rho_mouse < 0),
        strong & (rho_human > 0) & (rho_mouse < 0),
        strong & (rho_human < 0) & (rho_mouse > 0),
    ]
    return np.select(conditions, ["up_both", "down_both", "up_human_down_mouse", "down_human_up_mouse"], "below")


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.size": SIZE_NAME,
        "svg.fonttype": "none",
        "svg.hashsalt": SVG_HASH_SALT,
        "lines.solid_capstyle": "butt",
    })


def stealth_head() -> MplPath:
    """Swept-back arrowhead pointing +x with its tip at the origin, as on the manuscript axes."""
    vertices = [(0, 0), (-1.0, 0.5), (-0.72, 0), (-1.0, -0.5), (0, 0)]
    codes = [MplPath.MOVETO, MplPath.LINETO, MplPath.LINETO, MplPath.LINETO, MplPath.CLOSEPOLY]
    return MplPath(vertices, codes)


def draw_axes_cross(ax) -> list:
    """Both axes as double-headed arrows through the origin, with ticks crossing the line.

    Returns the tick labels, so mask_under_tick_labels can clear the threshold lines beneath them.
    """
    head = stealth_head()
    arrows = [
        ((-AXIS_END, 0), (AXIS_END, 0), 0),  # x axis: heads point to 0 and 180 degrees
        ((0, -AXIS_END), (0, AXIS_END), 90),  # y axis: heads point to 90 and 270 degrees
    ]
    for start, end, angle in arrows:
        ax.plot(*zip(start, end), color="black", lw=LINE_WIDTH_PT, zorder=3, clip_on=False)
        for tip, tip_angle in ((end, angle), (start, angle + 180)):
            ax.plot(*tip, marker=MarkerStyle(head).rotated(deg=tip_angle), markersize=ARROWHEAD_PT,
                    color="black", markeredgewidth=0, zorder=3, clip_on=False)

    # Tick marks and labels: x labels below the axis, y labels right of it, as in the manuscript. Labels
    # are drawn above the points: regulons cluster near rho = 0 in human, so some grey points fall
    # under the y labels and stay visible around the digits.
    tick_labels = []
    for tick in TICKS:
        label = f"{tick:.1f}".replace("-", "\u2212")
        ax.plot(tick, 0, marker="|", markersize=TICK_LENGTH_PT, markeredgewidth=LINE_WIDTH_PT,
                color="black", zorder=3, clip_on=False)
        ax.plot(0, tick, marker="_", markersize=TICK_LENGTH_PT, markeredgewidth=LINE_WIDTH_PT,
                color="black", zorder=3, clip_on=False)
        tick_labels.append(ax.annotate(label, (tick, 0), xytext=(0, -TICK_PAD_PT), textcoords="offset points",
                                       ha="center", va="top", fontsize=SIZE_NAME, zorder=4.5))
        tick_labels.append(ax.annotate(label, (0, tick), xytext=(TICK_PAD_PT, 0), textcoords="offset points",
                                       ha="left", va="center_baseline", fontsize=SIZE_NAME, zorder=4.5))
    return tick_labels


def mask_under_tick_labels(ax, tick_labels: list) -> None:
    """White boxes under the tick labels, below the points: the dotted lines at +-threshold would
    otherwise run through the digits (0.2 reads as 0:2), while every point stays drawn."""
    renderer = ax.figure.canvas.get_renderer()
    to_data = ax.transData.inverted()
    for label in tick_labels:
        box = label.get_window_extent(renderer)
        (x0, y0), (x1, y1) = to_data.transform([(box.x0, box.y0), (box.x1, box.y1)])
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="white", edgecolor="none", zorder=2))


def chevrons(direction: int, color: str) -> DrawingArea:
    """Three stacked chevrons (up for direction 1, down for -1), one cap height tall, before a word."""
    width, height, n = CHEVRON_WIDTH_PT, SIZE_LABEL * CAP_HEIGHT_EM, 3
    pitch = (height - CHEVRON_RISE_PT) / (n - 1)
    area = DrawingArea(width, height)
    for k in range(n):
        base = k * pitch
        tip, ends = (base + CHEVRON_RISE_PT, base) if direction > 0 else (base, base + CHEVRON_RISE_PT)
        area.add_artist(Line2D([0.2, width / 2, width - 0.2], [ends, tip, ends], color=color,
                               lw=CHEVRON_LINE_PT, solid_joinstyle="miter", solid_capstyle="butt"))
    return area


def quadrant_label(ax, signs: tuple[int, int], corner: str, color: str) -> None:
    """Chevrons and species names in one corner: up or down in human (top row) and mouse (bottom row)."""
    rows = []
    for species, sign in zip(("human", "mouse"), signs):
        word = TextArea(species, textprops=dict(color=color, fontsize=SIZE_LABEL, fontweight="bold"))
        rows.append(HPacker(children=[chevrons(sign, color), word], align="baseline", pad=0, sep=1.2))
    block = VPacker(children=rows, align="left", pad=0, sep=0.6)
    anchor = {"upper right": (1, 1), "upper left": (0, 1), "lower left": (0, 0), "lower right": (1, 0)}[corner]
    box = AnchoredOffsetbox(loc=corner, child=block, bbox_to_anchor=anchor, bbox_transform=ax.transAxes,
                            frameon=False, pad=0, borderpad=0)
    ax.add_artist(box)


def place_names(ax, table: pd.DataFrame, named: list[str], classes: np.ndarray) -> list[dict]:
    """Name each listed regulon beside its point; returns the placements for the provenance."""
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    points_px = ax.transData.transform(table[["rho_human", "rho_mouse"]].to_numpy())
    pt = fig.dpi / 72
    placements = []
    for name in named:
        row = int(np.flatnonzero(table["human_tf"].to_numpy() == name)[0])
        dx, dy, has_leader = NAME_OFFSETS_PT.get(name, DEFAULT_OFFSET_PT)
        is_below = classes[row] == "below"
        color = GREY_NAME if is_below else CLASS_COLORS[classes[row]]

        # Anchor on the baseline, shifted by half the cap height so the capitals centre on (dx, dy).
        ha = "left" if dx > 0 else "right" if dx < 0 else "center"
        text = ax.annotate(name, table.iloc[row][["rho_human", "rho_mouse"]].to_numpy(float),
                           xytext=(dx, dy - SIZE_NAME * CAP_HEIGHT_EM / 2), textcoords="offset points",
                           ha=ha, va="baseline", fontsize=SIZE_NAME, color=color, zorder=5)
        box = text.get_window_extent(renderer)

        # Leader: from just outside the point to the nearest edge of the name's box.
        leader = None
        if has_leader:
            start = points_px[row]
            end = np.clip(start, (box.x0 - 0.8 * pt, box.y0 - 0.8 * pt), (box.x1 + 0.8 * pt, box.y1 + 0.8 * pt))
            unit = (end - start) / np.hypot(*(end - start))
            start = start + unit * 2.2 * pt
            line = ax.transData.inverted().transform([start, end])
            ax.plot(line[:, 0], line[:, 1], color=GREY_LEADER, lw=0.5, zorder=4)
            leader = line.round(4).tolist()
        placements.append({
            "name": name,
            "class": str(classes[row]),
            "below_threshold": bool(is_below),
            "rho_human": float(table.iloc[row]["rho_human"]),
            "rho_mouse": float(table.iloc[row]["rho_mouse"]),
            "offset_pt": [dx, dy],
            "leader_data_xy": leader,
        })
    return placements


def layout_checks(fig, ax, table: pd.DataFrame, named: list[str]) -> dict:
    """Overlaps between drawn texts, and points lying under a name or other text (figure-style 9.1).

    Names must clear every point; tick labels and axis names are drawn above the points, so points under
    them are listed for the record.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    boxes = [t.get_window_extent(renderer) for t in texts]
    text_overlaps = [
        (texts[i].get_text(), texts[j].get_text())
        for i in range(len(texts)) for j in range(i + 1, len(texts)) if boxes[i].overlaps(boxes[j])
    ]
    points_px = ax.transData.transform(table[["rho_human", "rho_mouse"]].to_numpy())
    radius_px = np.sqrt(AREA_COLOURED) / 2 * fig.dpi / 72  # the larger marker, for every point
    covered = []
    for text, box in zip(texts, boxes):
        near_x = np.clip(points_px[:, 0], box.x0, box.x1)
        near_y = np.clip(points_px[:, 1], box.y0, box.y1)
        hit = np.hypot(points_px[:, 0] - near_x, points_px[:, 1] - near_y) < radius_px
        covered += [(text.get_text(), table["human_tf"].iloc[i]) for i in np.flatnonzero(hit)]
    figure_box = fig.bbox
    outside = [t.get_text() for t, b in zip(texts, boxes)
               if b.x0 < figure_box.x0 or b.x1 > figure_box.x1 or b.y0 < figure_box.y0 or b.y1 > figure_box.y1]
    return {
        "text_overlaps": text_overlaps,
        "points_under_names": [pair for pair in covered if pair[0] in named],
        "points_under_axis_text": [pair for pair in covered if pair[0] not in named],
        "text_outside_figure": outside,
    }


def save(fig, stem: Path, provenance: dict) -> None:
    """Write <stem>.svg and .png byte-reproducibly (no dates, fixed SVG ids), then the provenance."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    files = {}
    for ext, metadata in (("svg", {"Date": None}), ("png", None)):
        path = stem.parent / f"{stem.name}.{ext}"
        fig.savefig(path, dpi=PNG_DPI, metadata=metadata)
        files[path.name] = sha256_of(path)
    provenance["files"] = files
    provenance["written"] = datetime.datetime.now().isoformat(timespec="seconds")
    provenance["matplotlib"] = matplotlib.__version__
    provenance["font_family"] = FONT_FAMILY
    (stem.parent / f"{stem.name}.provenance.json").write_text(json.dumps(provenance, indent=1))
    print(f"wrote {stem}.{{svg,png,provenance.json}}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--table-sha256", required=True, help="sha256 pinned in config")
    parser.add_argument("--table-provenance", required=True, type=Path, help="upstream provenance, staged")
    parser.add_argument("--source-of-record", required=True, help="where the staged table was copied from")
    parser.add_argument("--threshold", required=True, type=float)
    parser.add_argument("--named", required=True, help="comma-separated human TF symbols to name")
    parser.add_argument("--out", required=True, type=Path, help="output stem")
    args = parser.parse_args()

    # Load ----

    observed_sha = sha256_of(args.table)
    if observed_sha != args.table_sha256:
        raise ValueError(f"{args.table}: sha256 {observed_sha} differs from the pinned {args.table_sha256}")
    table = pd.read_csv(args.table)
    named = args.named.split(",")
    missing = [name for name in named if name not in set(table["human_tf"])]
    if missing:
        raise KeyError(f"named regulons not in {args.table.name}: {missing}")
    classes = classify(table["rho_human"].to_numpy(), table["rho_mouse"].to_numpy(), args.threshold)
    counts = {c: int((classes == c).sum()) for c in [*CLASS_COLORS, "below"]}
    print(f"{len(table)} regulons; " + ", ".join(f"{c} {n}" for c, n in counts.items()))

    # Figure ----

    apply_style()
    fig = plt.figure(figsize=FIG_IN)
    width_in, height_in = FIG_IN
    ax = fig.add_axes([
        AXES_LEFT_IN / width_in,
        AXES_BOTTOM_IN / height_in,
        AXES_SIZE_IN / width_in,
        AXES_SIZE_IN / height_in,
    ])
    ax.set_xlim(-AXIS_END, AXIS_END)
    ax.set_ylim(-AXIS_END, AXIS_END)
    ax.set_axis_off()

    # Dotted lines at |rho| = threshold, the colouring rule, across the whole plot as in the manuscript.
    dotted = dict(color=GREY_THRESHOLD, lw=0.7, linestyle=(0, (1, 2)), zorder=1)
    for value in (-args.threshold, args.threshold):
        ax.plot([-THRESHOLD_LINE_END, THRESHOLD_LINE_END], [value, value], **dotted)
        ax.plot([value, value], [-THRESHOLD_LINE_END, THRESHOLD_LINE_END], **dotted)
    tick_labels = draw_axes_cross(ax)

    # One point per regulon at (rho human, rho mouse): grey below the threshold, else the quadrant colour.
    is_below = classes == "below"
    ax.scatter(table.loc[is_below, "rho_human"], table.loc[is_below, "rho_mouse"], s=AREA_GREY,
               color=GREY_POINT, linewidths=0, zorder=4)
    for class_name, color in CLASS_COLORS.items():
        chosen = classes == class_name
        ax.scatter(table.loc[chosen, "rho_human"], table.loc[chosen, "rho_mouse"], s=AREA_COLOURED,
                   color=color, linewidths=0, zorder=4)

    # Quadrant labels in the corners; an empty quadrant's label is grey, as in the manuscript.
    for class_name, (sign_human, sign_mouse, corner) in QUADRANTS.items():
        color = CLASS_COLORS[class_name] if counts[class_name] else GREY_POINT
        quadrant_label(ax, (sign_human, sign_mouse), corner, color)

    # Axis names at the positive arrow tips, and the panel title on the left, as in the manuscript.
    ax.annotate("human (\u03c1)", (AXIS_END, 0), xytext=(AXIS_NAME_OVERHANG_PT, 3), textcoords="offset points",
                ha="right", va="bottom", fontsize=SIZE_LABEL, fontweight="bold", annotation_clip=False, zorder=5)
    ax.annotate("mouse (\u03c1)", (0, AXIS_END), xytext=(-3, 0), textcoords="offset points", ha="right",
                va="bottom", rotation=90, rotation_mode="anchor", fontsize=SIZE_LABEL, fontweight="bold",
                annotation_clip=False, zorder=5)
    fig.text(TITLE_X_IN / width_in, (AXES_BOTTOM_IN + AXES_SIZE_IN / 2) / height_in,
             "SCENIC regulons of X-Mens and human", rotation=90, ha="center", va="center",
             fontsize=SIZE_LABEL, fontweight="bold")

    fig.canvas.draw()
    mask_under_tick_labels(ax, tick_labels)
    placements = place_names(ax, table, named, classes)
    checks = layout_checks(fig, ax, table, named)
    print("layout checks:", checks)

    # Save ----

    sign_agree = counts["up_both"] + counts["down_both"]
    n_coloured = sum(counts[c] for c in CLASS_COLORS)
    provenance = {
        "panel": "5E",
        "script": Path(__file__).name,
        "script_sha256": sha256_of(Path(__file__)),
        "inputs": {
            "table": {"path": str(args.table.resolve()), "sha256": observed_sha, "matches_pinned": True},
            "table_source_of_record": args.source_of_record,
            "upstream_provenance": json.loads(args.table_provenance.read_text()),
        },
        "unit": "one point per transcription factor with a floor-6 SCENIC regulon in both species (Ensembl one-to-one pair)",
        "x": "Spearman rho of regulon activity with rank-normalized pseudotime (pt_harmonized_rank), human",
        "y": "the same, X-Mens (mouse)",
        "threshold": args.threshold,
        "n_regulons": int(len(table)),
        "class_counts": counts,
        "sign_agreement": f"{sign_agree} of {n_coloured} coloured regulons agree in sign",
        "named": placements,
        "named_below_threshold": [p["name"] for p in placements if p["below_threshold"]],
        "label_column_of_table": "not used; the named regulons come from config (manuscript panel)",
        "colours": CLASS_COLORS,
        "grey": {"points": GREY_POINT, "names_below_threshold": GREY_NAME},
        "layout_checks": checks,
    }
    save(fig, args.out, provenance)
    plt.close(fig)


if __name__ == "__main__":
    main()
