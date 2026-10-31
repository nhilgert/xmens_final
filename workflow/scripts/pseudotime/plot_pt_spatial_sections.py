"""Rank pseudotime on the tissue: one black-background viridis panel per physical section, plus an overview.

Reads the fib/dec release object (.h5ad, obs only) or its plot table (.parquet) for one species and writes to --outdir:
  pt_spatial_<species>_<section>.{png,pdf,svg}       one panel per physical section, own scale bar
  pt_spatial_<species>_all_sections.{png,pdf,svg}    every section on one page, one shared spatial scale
  pt_spatial_<species>_cells.csv                     every nucleus considered, with each display decision
  pt_spatial_<species>_provenance.json               input sha256, parameters, per-section counts
With --no-filter the spatial discontinuity filter is skipped and every file name gets the suffix _unfiltered.
Mouse sections are turned for display by --section-rotation (resources/figures/mouse_section_rotation.yaml, the table
every mouse spatial panel reads through workflow/scripts/figures/section_rotation.py): each timepoint's array turns
as one rigid body, both uninduced pieces together. Human sections are drawn unrotated. No rule runs this script.

usage: plot_pt_spatial_sections.py --input OBJ --species {mouse,human} --outdir DIR [--pt-col pt_harmonized_rank]
                                   [--section-rotation resources/figures/mouse_section_rotation.yaml]   (mouse: required)
"""
import argparse
import datetime
import hashlib
import json
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.spatial import cKDTree

# The mouse rotation reader is shared by every mouse spatial panel and lives in workflow/scripts/figures/.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
import section_rotation  # noqa: E402


# Physical sections ----
# One Slide-tags tissue section can be spread over several libraries (10x lanes loaded from the same nuclei
# suspension). Libraries of one section share a coordinate frame: the nearest nucleus from another library of the
# section is as close as the nearest nucleus of the same library (mouse mck10a/b 38.6 vs 40.2 um; human mck6a/b/c
# and mck5a/b/c 19.6-26.1 vs 19.4-22.8 um), so they are drawn together in native coordinates. The mouse uninduced
# puck holds two separate tissue pieces (A, B; about 2 mm apart), drawn as two sections.
# Key: obs['tissue_section'] value -> (section id used in file names, panel title). Keys carry the sample IDs
# (original timepoint names); ids and titles use the released timepoint names (vocabulary v5).
SECTIONS = {
    "mouse": {
        "uninduced_mck8a__A": ("uninduced_A", "uninduced, piece A"),
        "uninduced_mck8b__A": ("uninduced_A", "uninduced, piece A"),
        "uninduced_mck8a__B": ("uninduced_B", "uninduced, piece B"),
        "uninduced_mck8b__B": ("uninduced_B", "uninduced, piece B"),
        "1_dpi_early_mck1__A": ("1.5_dpi_early", "1.5 dpi early"),
        "1_dpi_late_mck12a__A": ("1.5_dpi_late", "1.5 dpi late"),
        "1_dpi_late_mck12b__A": ("1.5_dpi_late", "1.5 dpi late"),
        "2_dpi_mck15__A": ("2.5_dpi", "2.5 dpi"),
        "3_dpi_mck10a__A": ("3.5_dpi", "3.5 dpi"),
        "3_dpi_mck10b__A": ("3.5_dpi", "3.5 dpi"),
        "4_dpi_bleeding_mck11a__A": ("4.5_dpi_bleeding", "4.5 dpi bleeding"),
    },
    "human": {
        "secretory_mck6a__A": ("secretory", "secretory (decidualized)"),
        "secretory_mck6b__A": ("secretory", "secretory (decidualized)"),
        "secretory_mck6c__A": ("secretory", "secretory (decidualized)"),
        "menstrual_mck5a__A": ("menstrual", "menstrual"),
        "menstrual_mck5b__A": ("menstrual", "menstrual"),
        "menstrual_mck5c__A": ("menstrual", "menstrual"),
    },
}
OVERVIEW_COLUMNS = {"mouse": 4, "human": 2}

# Spatial discontinuity filter: a nucleus is hidden when |its pseudotime - median pseudotime of its K nearest
# spatial neighbors| is above the Q-quantile of that difference within its section, i.e. the most discordant 10%.
FILTER_K = 5
FILTER_Q = 0.90

SCALE_BAR_UM = 1000
SEED = 0  # draw order is shuffled so no library or region systematically covers another

# Type sizes (pt) and panel geometry (in); white on black throughout.
SIZE_TITLE, SIZE_LABEL, SIZE_TICK = 8, 7, 6
PANEL_IN = 2.4
HEADER_IN, FOOTER_IN = 0.55, 0.6  # single panel: title band above the map, colorbar band below
OVERVIEW_PANEL_IN = 1.9
FOREGROUND = "white"
CMAP = "viridis"  # perceptually uniform; rank pseudotime 0 (root) to 1 (most mature)

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("pt_spatial")
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def set_fonts() -> None:
    helvetica = Path("/System/Library/Fonts/HelveticaNeue.ttc")
    if helvetica.exists():
        font_manager.fontManager.addfont(str(helvetica))
    plt.rcParams.update({
        "font.family": ["Helvetica Neue", "Helvetica", "Nimbus Sans", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })


def load_obs(path: Path, columns: list[str]) -> pd.DataFrame:
    """obs columns from a release .h5ad (read backed, matrix untouched) or from its plot-table parquet."""
    if path.suffix == ".h5ad":
        import anndata

        adata = anndata.read_h5ad(path, backed="r")
        obs = adata.obs[columns].copy()
        adata.file.close()
        return obs
    return pd.read_parquet(path, columns=columns)


def discontinuity_residual(xy_um: np.ndarray, pt: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """|pt - median pt of the k nearest spatial neighbors| and that neighbor median, per nucleus.

    A few nuclei from different libraries of one section share an identical position (2 pairs in mouse, 1 in
    human), so each nucleus is removed from its own neighbor list by index rather than by rank.
    """
    _, neighbors = cKDTree(xy_um).query(xy_um, k=k + 1)
    is_other = neighbors != np.arange(len(xy_um))[:, None]
    if not (is_other.sum(axis=1) >= k).all():
        raise ValueError(f"more than {k} nuclei share one position; neighbor sets are undefined")
    others = np.array([row[keep][:k] for row, keep in zip(neighbors, is_other)])
    neighbor_median = np.median(pt[others], axis=1)
    return np.abs(pt - neighbor_median), neighbor_median


def section_rotation_of(section_id: str, species: str, angles: dict, obs: pd.DataFrame, cells: pd.DataFrame) -> float:
    """Display rotation (degrees CCW) of one physical section: its timepoint's angle (mouse) or 0 (human)."""
    if species != "mouse":
        return 0.0
    timepoints = set(obs.loc[cells["section"] == section_id, "timepoint"].astype(str))
    if len(timepoints) != 1:
        raise ValueError(f"section {section_id} spans timepoints {sorted(timepoints)}")
    return angles[timepoints.pop()]


def marker_area(n_nuclei: int, panel_in: float) -> float:
    """Scatter marker area (pt^2) so markers cover a similar share of the panel whatever the section's size."""
    return float(np.clip(0.45 * (panel_in * 72) ** 2 / n_nuclei, 0.4, 6.0))


def draw_section(ax, cells: pd.DataFrame, half_width_um: float, panel_in: float, norm: Normalize) -> None:
    # Position = tissue coordinates (rotated for display), color = rank pseudotime on one fixed 0-1 scale so
    # colors compare across sections; equal aspect keeps the tissue undistorted.
    ax.set_facecolor("black")
    ax.scatter(
        cells["x_plot"],
        cells["y_plot"],
        c=cells["pt"],
        cmap=CMAP,
        norm=norm,
        s=marker_area(len(cells), panel_in),
        linewidths=0,
        rasterized=True,
    )
    cx = (cells["x_plot"].min() + cells["x_plot"].max()) / 2
    cy = (cells["y_plot"].min() + cells["y_plot"].max()) / 2
    ax.set_xlim(cx - half_width_um, cx + half_width_um)
    ax.set_ylim(cy - half_width_um, cy + half_width_um)
    ax.set_aspect("equal")
    ax.axis("off")


def add_scale_bar(ax, half_width_um: float, label: bool = True) -> None:
    x0, x1 = ax.get_xlim()
    y0, _ = ax.get_ylim()
    bar_x = x0 + 0.05 * (x1 - x0)
    bar_y = y0 + 0.04 * (x1 - x0)
    ax.plot([bar_x, bar_x + SCALE_BAR_UM], [bar_y, bar_y], color=FOREGROUND, lw=1.2, solid_capstyle="butt")
    if label:
        ax.text(bar_x + SCALE_BAR_UM / 2, bar_y + 0.02 * (x1 - x0), "1 mm", color=FOREGROUND,
                fontsize=SIZE_TICK, ha="center", va="bottom")


def add_colorbar(fig, rect: list[float], norm: Normalize) -> None:
    cax = fig.add_axes(rect)
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(colors=FOREGROUND, labelsize=SIZE_TICK, length=2, width=0.5, pad=1.5)
    colorbar.ax.set_xticklabels(["0", "0.5", "1"])
    cax.set_title("pseudotime (rank)", color=FOREGROUND, fontsize=SIZE_LABEL, pad=2)


def save(fig, stem: Path) -> list[str]:
    paths = []
    for ext in ("png", "pdf", "svg"):
        # Not stem.with_suffix(): section ids such as 1.5_dpi_late contain a dot, which with_suffix takes for an
        # extension (1.5_dpi_early and 1.5_dpi_late both became pt_spatial_mouse_1.*, one overwriting the other).
        path = stem.parent / f"{stem.name}.{ext}"
        fig.savefig(path, dpi=600, facecolor="black")
        paths.append(str(path))
    plt.close(fig)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--species", required=True, choices=sorted(SECTIONS))
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--pt-col", default="pt_harmonized_rank")
    parser.add_argument("--no-filter", action="store_true")
    parser.add_argument("--section-rotation", type=Path,
                        help="mouse only (required): resources/figures/mouse_section_rotation.yaml")
    args = parser.parse_args()
    if args.species == "mouse" and args.section_rotation is None:
        parser.error("--section-rotation is required for mouse sections")

    args.outdir.mkdir(parents=True, exist_ok=True)
    suffix = "_unfiltered" if args.no_filter else ""
    set_fonts()
    sections = SECTIONS[args.species]

    # Load: one row per nucleus of the fib/dec object ----
    columns = [args.pt_col, "x_um_qc", "y_um_qc", "tissue_section", "embedding_outlier"]
    obs = load_obs(args.input, columns + (["timepoint"] if args.species == "mouse" else []))
    obs["tissue_section"] = obs["tissue_section"].astype(object).where(obs["tissue_section"].notna(), None)
    obs.loc[obs["tissue_section"].astype(str) == "nan", "tissue_section"] = None
    unknown = set(obs["tissue_section"].dropna()) - set(sections)
    if unknown:
        raise ValueError(f"tissue_section values without a section assignment: {sorted(unknown)}")
    log.info(f"{args.species}: {len(obs):,} nuclei in {args.input.name}")

    # Display rules, in order; each nucleus records the first rule it fails ----
    # Pseudotime NaN = mouse decAp (excluded from the trajectory); no QC'd coordinates = not placed
    # on the tissue; embedding_outlier = release display filter (kept in statistics, never drawn).
    cells = pd.DataFrame({
        "tissue_section": obs["tissue_section"],
        "section": obs["tissue_section"].map(lambda s: sections[s][0] if s else None),
        "x_um": obs["x_um_qc"].astype(float),
        "y_um": obs["y_um_qc"].astype(float),
        "pt": obs[args.pt_col].astype(float),
        "is_embedding_outlier": obs["embedding_outlier"].astype(bool),
    }, index=obs.index)
    cells["excluded_because"] = None
    rules = [
        ("no pseudotime", cells["pt"].isna()),
        ("no QC'd coordinates", cells[["x_um", "y_um"]].isna().any(axis=1) | cells["section"].isna()),
        ("embedding outlier", cells["is_embedding_outlier"]),
    ]
    for reason, fails in rules:
        cells.loc[fails & cells["excluded_because"].isna(), "excluded_because"] = reason

    # Discontinuity filter (native coordinates) and display coordinates, per physical section ----
    cells["neighbor_median_pt"] = np.nan
    cells["abs_residual"] = np.nan
    cells["x_plot"] = np.nan
    cells["y_plot"] = np.nan

    # Display position of every placed nucleus. Mouse: each timepoint's whole array (all libraries, both uninduced
    # pieces, outliers included) turns as one rigid body by the shared table's angle. Human: as placed. Stored in
    # x_plot/y_plot below for the nuclei that enter the display filter, as before.
    placed = cells["section"].notna() & cells[["x_um", "y_um"]].notna().all(axis=1)
    display_xy = cells.loc[placed, ["x_um", "y_um"]].copy()
    angles = {}
    if args.species == "mouse":
        angles = section_rotation.load_rotation_table(args.section_rotation)
        display_xy["x_um"], display_xy["y_um"] = section_rotation.rotate_by_timepoint(
            display_xy["x_um"], display_xy["y_um"], obs.loc[placed, "timepoint"].astype(str), angles
        )
    per_section = {}
    for section_id in dict.fromkeys(s for s, _ in sections.values()):
        in_section = (cells["section"] == section_id) & cells["excluded_because"].isna()
        xy_um = cells.loc[in_section, ["x_um", "y_um"]].to_numpy()
        pt = cells.loc[in_section, "pt"].to_numpy()

        residual, neighbor_median = discontinuity_residual(xy_um, pt, FILTER_K)
        cut = float(np.quantile(residual, FILTER_Q))
        cells.loc[in_section, "abs_residual"] = residual
        cells.loc[in_section, "neighbor_median_pt"] = neighbor_median
        if not args.no_filter:
            hidden = in_section.copy()
            hidden.loc[in_section] = residual > cut
            cells.loc[hidden, "excluded_because"] = "spatial discontinuity"

        section_index = cells.index[in_section]
        cells.loc[in_section, "x_plot"] = display_xy.loc[section_index, "x_um"].to_numpy()
        cells.loc[in_section, "y_plot"] = display_xy.loc[section_index, "y_um"].to_numpy()

        all_in_section = cells["section"] == section_id
        per_section[section_id] = dict(
            title=next(t for s, t in sections.values() if s == section_id),
            libraries=sorted({ts.split("__")[0].rsplit("_", 1)[-1] for ts, (s, _) in sections.items() if s == section_id}),
            nuclei_in_section=int(all_in_section.sum()),
            excluded_no_pseudotime=int((all_in_section & (cells["excluded_because"] == "no pseudotime")).sum()),
            excluded_embedding_outlier=int((all_in_section & (cells["excluded_because"] == "embedding outlier")).sum()),
            nuclei_filter_input=int(in_section.sum()),
            discontinuity_cut=round(cut, 4),
            excluded_discontinuity=int((all_in_section & (cells["excluded_because"] == "spatial discontinuity")).sum()),
            rotation_deg_ccw=section_rotation_of(section_id, args.species, angles, obs, cells),
        )
        drawn = all_in_section & cells["excluded_because"].isna()
        per_section[section_id]["nuclei_drawn"] = int(drawn.sum())
        per_section[section_id]["median_pt_drawn"] = round(float(cells.loc[drawn, "pt"].median()), 3)
        log.info(f"  {section_id}: {per_section[section_id]}")

    cells["is_drawn"] = cells["excluded_because"].isna()
    order = np.random.default_rng(SEED).permutation(len(cells))
    drawn_cells = cells.iloc[order][cells["is_drawn"].iloc[order].to_numpy()]
    norm = Normalize(0, 1)
    outputs = []

    # One panel per section: own spatial scale, explicit section name, libraries and counts ----
    for section_id, info in per_section.items():
        section_cells = drawn_cells[drawn_cells["section"] == section_id]
        half_width = max(np.ptp(section_cells["x_plot"]), np.ptp(section_cells["y_plot"])) / 2 * 1.05
        fig_h = HEADER_IN + PANEL_IN + FOOTER_IN
        fig = plt.figure(figsize=(PANEL_IN, fig_h), facecolor="black")
        ax = fig.add_axes([0, FOOTER_IN / fig_h, 1, PANEL_IN / fig_h])
        draw_section(ax, section_cells, half_width, PANEL_IN, norm)
        add_scale_bar(ax, half_width)
        fig.text(0.5, 1 - 0.1 / fig_h, f"{args.species.capitalize()} {info['title']}", color=FOREGROUND,
                 fontsize=SIZE_TITLE, ha="center", va="top")
        fig.text(0.5, 1 - 0.3 / fig_h, f"libraries {' + '.join(info['libraries'])}  ·  n = {info['nuclei_drawn']:,} nuclei",
                 color=FOREGROUND, fontsize=SIZE_TICK, ha="center", va="top")
        add_colorbar(fig, [0.3, 0.24 / fig_h, 0.4, 0.06 / fig_h], norm)
        outputs += save(fig, args.outdir / f"pt_spatial_{args.species}_{section_id}{suffix}")

    # Overview: every section at one shared spatial scale, so section sizes compare directly ----
    n_columns = OVERVIEW_COLUMNS[args.species]
    n_rows = int(np.ceil(len(per_section) / n_columns))
    half_width = max(
        max(np.ptp(g["x_plot"]), np.ptp(g["y_plot"])) / 2 for _, g in drawn_cells.groupby("section")
    ) * 1.05
    title_in, colorbar_in = 0.35, 0.55
    fig_w = n_columns * OVERVIEW_PANEL_IN
    fig_h = n_rows * (OVERVIEW_PANEL_IN + title_in) + colorbar_in
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="black")
    for i, (section_id, info) in enumerate(per_section.items()):
        row, column = divmod(i, n_columns)
        top_in = row * (OVERVIEW_PANEL_IN + title_in) + title_in
        ax = fig.add_axes([column * OVERVIEW_PANEL_IN / fig_w, 1 - (top_in + OVERVIEW_PANEL_IN) / fig_h,
                           OVERVIEW_PANEL_IN / fig_w, OVERVIEW_PANEL_IN / fig_h])
        draw_section(ax, drawn_cells[drawn_cells["section"] == section_id], half_width, OVERVIEW_PANEL_IN, norm)
        ax.set_title(f"{info['title']}\n{' + '.join(info['libraries'])}  ·  n = {info['nuclei_drawn']:,}",
                     color=FOREGROUND, fontsize=SIZE_LABEL, pad=2)
        if i == 0:
            add_scale_bar(ax, half_width)
    add_colorbar(fig, [0.5 - 0.6 / fig_w, 0.32 * colorbar_in / fig_h, 1.2 / fig_w, 0.08 / fig_h], norm)
    outputs += save(fig, args.outdir / f"pt_spatial_{args.species}_all_sections{suffix}")

    # Per-nucleus table and provenance ----
    table_path = args.outdir / f"pt_spatial_{args.species}_cells{suffix}.csv"
    cells.rename_axis("cell").to_csv(table_path)
    provenance = dict(
        written=datetime.datetime.now().isoformat(timespec="seconds"),
        input=dict(path=str(args.input.resolve()), sha256=sha256(args.input)),
        section_rotation=(section_rotation.table_record(args.section_rotation) if args.species == "mouse"
                          else "human sections drawn unrotated"),
        pt_column=args.pt_col,
        filter=None if args.no_filter else dict(k=FILTER_K, quantile=FILTER_Q, scope="per physical section, native coordinates"),
        display_rules=[r for r, _ in rules],
        sections=per_section,
        seed=SEED,
        outputs=outputs + [str(table_path)],
    )
    (args.outdir / f"pt_spatial_{args.species}_provenance{suffix}.json").write_text(json.dumps(provenance, indent=1))
    log.info(f"wrote {len(outputs)} figure files and {table_path.name}")


if __name__ == "__main__":
    main()
