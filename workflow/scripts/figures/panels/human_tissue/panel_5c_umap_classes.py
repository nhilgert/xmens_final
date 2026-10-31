"""Fig. 5C: UMAP of all human Slide-tags nuclei coloured by cell class, with direct labels, on white.

Called by rule `panels_human_tissue_5c`, once per release.
Reads:  the plot table and extraction provenance of one release (extract_human_tissue.py).
Writes: <out-stem>.pdf/.svg/.png and <out-stem>.provenance.json.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

import panel_common as pc


# Geometry (in): 1.95 in tall as panel C of the published Fig. 5; the width follows the embedding's
# aspect so UMAP distances are not stretched along one axis.
HEIGHT_IN = 1.95
PAD_FRACTION = 0.06  # whitespace around the embedding, room for labels at the edges

# Marker area (pt^2): about 68k nuclei, small enough that the cluster outlines stay sharp.
MARKER_AREA = 0.35

LABEL_SIZE = 7

# Label placement: the largest connected cluster of each class is found on a grid of GRID_UNITS
# UMAP units; the label goes to the nearest spot along 24 directions whose text box holds at most
# MAX_POINTS_UNDER_LABEL nuclei of any class. A leader line is drawn only when that spot is
# further than LEADER_MIN_GAP_PT from the cluster.
GRID_UNITS = 0.25
N_DIRECTIONS = 24
MAX_POINTS_UNDER_LABEL = 3
LEADER_MIN_GAP_PT = 6
LABEL_PAD_PT = 1.5  # clearance kept around each label, from nuclei and from other labels


def main_cluster(xy: np.ndarray) -> np.ndarray:
    """Boolean mask of the points in the largest 8-connected grid component of one class."""
    cell_index = np.floor((xy - xy.min(axis=0)) / GRID_UNITS).astype(int)
    occupancy = np.zeros(cell_index.max(axis=0) + 1, dtype=int)
    np.add.at(occupancy, (cell_index[:, 0], cell_index[:, 1]), 1)

    components, _ = ndimage.label(occupancy > 0, structure=np.ones((3, 3)))
    point_component = components[cell_index[:, 0], cell_index[:, 1]]
    largest = np.bincount(point_component[point_component > 0]).argmax()
    return point_component == largest


def label_half_size(ax, text: str) -> np.ndarray:
    """Half width and half height of a label in data units; the axes scale is linear, so any position gives the same size."""
    artist = ax.text(0, 0, text, fontsize=LABEL_SIZE, ha="center", va="center")
    box = artist.get_window_extent(ax.figure.canvas.get_renderer()).transformed(ax.transData.inverted())
    artist.remove()
    return np.array([box.width / 2, box.height / 2])


def points_per_data_unit(ax) -> float:
    """Typographic points per UMAP unit (equal aspect, so the same along x and y)."""
    pixels = ax.transData.transform([1, 0])[0] - ax.transData.transform([0, 0])[0]
    return pixels * 72 / ax.figure.dpi


def make_point_counter(xy: np.ndarray):
    """Function counting the points inside an [x0, y0, x1, y1] box, from points sorted by x."""
    order = np.argsort(xy[:, 0])
    xs = xy[order, 0]
    ys = xy[order, 1]

    def count(box: np.ndarray) -> int:
        lo, hi = np.searchsorted(xs, [box[0], box[2]])
        y = ys[lo:hi]
        return int(((y > box[1]) & (y < box[3])).sum())

    return count


def boxes_overlap(a: np.ndarray, b: np.ndarray) -> bool:
    return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])


def place_labels(ax, cells) -> dict:
    """Place one label per class next to its largest cluster, in whitespace; returns the placements."""
    count_points = make_point_counter(cells[["umap_1", "umap_2"]].to_numpy())
    x_lim, y_lim = ax.get_xlim(), ax.get_ylim()
    max_distance = max(np.ptp(x_lim), np.ptp(y_lim))
    step = 0.4 * GRID_UNITS
    data_units_per_pt = 1 / points_per_data_unit(ax)
    placed_boxes = []
    placements = {}

    for cell_class in pc.CLASS_ORDER:
        xy = cells.loc[cells["cell_type"] == cell_class, ["umap_1", "umap_2"]].to_numpy()
        cluster = xy[main_cluster(xy)]
        center = np.median(cluster, axis=0)
        text = pc.ABBREVIATIONS[cell_class]
        half = label_half_size(ax, text) + LABEL_PAD_PT * data_units_per_pt
        best = None

        # Walk outward from the cluster centre along each direction until the label box is clear.
        for angle in np.linspace(0, 2 * np.pi, N_DIRECTIONS, endpoint=False):
            direction = np.array([np.cos(angle), np.sin(angle)])
            for distance in np.arange(0, max_distance, step):
                position = center + distance * direction
                box = np.concatenate([position - half, position + half])
                inside_axes = box[0] > x_lim[0] and box[2] < x_lim[1] and box[1] > y_lim[0] and box[3] < y_lim[1]
                if not inside_axes:
                    break
                is_clear = count_points(box) <= MAX_POINTS_UNDER_LABEL
                if is_clear and not any(boxes_overlap(box, b) for b in placed_boxes):
                    if best is None or distance < best["distance"]:
                        best = {"distance": distance, "position": position, "box": box}
                    break

        if best is None:
            raise ValueError(f"no clear label position for {cell_class}")
        placed_boxes.append(best["box"])

        # Nearest point of the cluster, for the leader-line decision and its end point.
        gaps = np.linalg.norm(cluster - best["position"], axis=1)
        placements[cell_class] = {
            "text": text,
            "position": best["position"],
            "nearest_cluster_point": cluster[gaps.argmin()],
            "cluster_size": int(len(cluster)),
            "class_size": int(len(xy)),
        }
    return placements


def draw_labels(ax, placements: dict) -> None:
    # Labels in black as in the published panel; a thin leader joins a label to its cluster only when
    # the nearest clear spot is detached from it.
    points_per_unit = points_per_data_unit(ax)
    for placement in placements.values():
        ax.text(*placement["position"], placement["text"], fontsize=LABEL_SIZE, ha="center", va="center")
        gap_pt = np.linalg.norm(placement["position"] - placement["nearest_cluster_point"]) * points_per_unit
        placement["leader"] = bool(gap_pt > LEADER_MIN_GAP_PT + LABEL_SIZE)
        if placement["leader"]:
            ax.annotate(
                "",
                xy=placement["nearest_cluster_point"],
                xytext=placement["position"],
                arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "black", "shrinkA": 5, "shrinkB": 1},
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--extract-provenance", required=True, type=Path)
    parser.add_argument("--release", required=True)
    parser.add_argument("--out-stem", required=True, type=Path)
    args = parser.parse_args()

    pc.set_style()

    # Load: one point per nucleus; embedding outliers hidden (release rule) ----
    cells = pc.load_cells(args.cells)
    is_outlier = cells["embedding_outlier"].astype(bool)
    drawn = cells[~is_outlier]
    pc.log.info(f"{args.release}: {len(cells):,} nuclei, {int(is_outlier.sum())} embedding outliers hidden, {len(drawn):,} drawn")

    order = np.random.default_rng(pc.SEED).permutation(len(drawn))
    drawn = drawn.iloc[order]

    # Figure: equal aspect, width from the embedding's extent ----
    x_span = np.ptp(drawn["umap_1"])
    y_span = np.ptp(drawn["umap_2"])
    pad = PAD_FRACTION * max(x_span, y_span)
    width_in = HEIGHT_IN * (x_span + 2 * pad) / (y_span + 2 * pad)
    fig = plt.figure(figsize=(width_in, HEIGHT_IN), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])

    # Points: colour = cell class, shuffled draw order so no class hides another ----
    ax.scatter(
        drawn["umap_1"],
        drawn["umap_2"],
        c=drawn["cell_type"].map(pc.CLASS_COLORS),
        s=MARKER_AREA,
        linewidths=0,
        rasterized=True,
    )
    ax.set_xlim(drawn["umap_1"].min() - pad, drawn["umap_1"].max() + pad)
    ax.set_ylim(drawn["umap_2"].min() - pad, drawn["umap_2"].max() + pad)
    ax.set_aspect("equal")
    ax.axis("off")

    # Direct labels instead of a legend, as in the published panel ----
    placements = place_labels(ax, drawn)
    draw_labels(ax, placements)

    outputs = pc.save_panel(fig, args.out_stem, facecolor="white")

    # Provenance ----
    scripts = [Path(__file__), Path(pc.__file__)]
    provenance = pc.base_provenance("5C", args.release, args.cells, args.extract_provenance, scripts)
    provenance.update({
        "embedding": "obsm X_umap_scVI of the release object",
        "filters": ["obs embedding_outlier hidden (release rule: omitted from embedding and spatial plots)"],
        "n_nuclei_release": int(len(cells)),
        "n_embedding_outliers_hidden": int(is_outlier.sum()),
        "n_drawn": int(len(drawn)),
        "drawn_by_class": drawn["cell_type"].value_counts().reindex(pc.CLASS_ORDER).to_dict(),
        "labels": {
            p["text"]: {
                "x": round(float(p["position"][0]), 3),
                "y": round(float(p["position"][1]), 3),
                "leader": p["leader"],
                "nuclei_in_labelled_cluster": p["cluster_size"],
                "nuclei_in_class": p["class_size"],
            }
            for p in placements.values()
        },
        "figure_size_in": [round(width_in, 3), HEIGHT_IN],
        "marker_area_pt2": MARKER_AREA,
        "draw_order_seed": pc.SEED,
        "outputs": outputs,
    })
    provenance["colors"]["classes"] = {pc.ABBREVIATIONS[c]: pc.CLASS_COLORS[c] for c in pc.CLASS_ORDER}
    pc.write_json(Path(f"{args.out_stem}.provenance.json"), provenance)


if __name__ == "__main__":
    main()
