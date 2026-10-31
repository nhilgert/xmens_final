"""Display rotation of the mouse Slide-tags sections: the one reader of resources/figures/mouse_section_rotation.yaml.

Every mouse spatial panel turns its sections through this module: Fig. S13 (panels/mouse_tissue/plot_s13_spatial.py),
Figs. S15B, S16 and S22 (panels/mouse_stroma/panel_style.py), the mouse row of Fig. 6D
(panels/cross_species/spatial_efib_dec_panel.py), the S24 spatial maps (zone_modules/fig_zone_modules.py) and the
pseudotime section panels (pseudotime/plot_pt_spatial_sections.py).

Convention: counter-clockwise degrees in data coordinates, drawn with y increasing upwards. All nuclei of one
timepoint turn as one rigid body about one centre, the midpoint of their bounding box, so the libraries of a section
and the two uninduced pieces keep their placement relative to each other. Callers pass every placed nucleus of the
timepoint (the whole array) and select what they draw afterwards; a crop defined in raw coordinates is applied to
the raw coordinates, so it selects the same nuclei whatever the angle. Display only: callers store the result in
new columns, and the coordinates of record are never changed.
"""
import hashlib
from pathlib import Path

import numpy as np
import yaml

# Release timepoints (obs['timepoint']); the table must give an angle for each, so no section is drawn unrotated
# by omission.
TIMEPOINTS = ("uninduced", "1.5_dpi_early", "1.5_dpi_late", "2.5_dpi", "3.5_dpi", "4.5_dpi_bleeding")

CONVENTION = ("counter-clockwise degrees in data coordinates, y up; each timepoint's array rotated as one rigid body "
              "about the midpoint of its bounding box; no mirror; display only")


def load_rotation_table(path) -> dict:
    """Counter-clockwise display rotation (degrees) per release timepoint."""
    spec = yaml.safe_load(Path(path).read_text())
    angles = {str(timepoint): float(degrees) for timepoint, degrees in spec["rotation_deg_ccw"].items()}
    missing = [t for t in TIMEPOINTS if t not in angles]
    unknown = sorted(set(angles) - set(TIMEPOINTS))
    if missing or unknown:
        raise ValueError(f"{path}: the table must list each of {TIMEPOINTS}; missing {missing}, unknown {unknown}")
    return angles


def table_record(path) -> dict:
    """Provenance entry for a panel: where the angles came from, their sha256 and the angles themselves."""
    return {
        "table": str(path),
        "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "rotation_deg_ccw": load_rotation_table(path),
        "convention": CONVENTION,
    }


def rotate_about(xy: np.ndarray, degrees_ccw: float, centre: np.ndarray) -> np.ndarray:
    """Points (n x 2) turned counter-clockwise, y up, about centre."""
    theta = np.deg2rad(degrees_ccw)
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    return (xy - centre) @ rotation.T + centre


def array_centre(xy: np.ndarray) -> np.ndarray:
    """Pivot of one timepoint's array: the midpoint of the bounding box of its nuclei."""
    return (xy.min(axis=0) + xy.max(axis=0)) / 2


def rotate_by_timepoint(x, y, timepoints, angles: dict) -> tuple:
    """Display coordinates (x, y arrays in input order): each timepoint's nuclei turned by its angle about one centre.

    Pass placed nuclei only; a nucleus without coordinates has no place to turn from.
    """
    xy = np.column_stack([np.asarray(x, dtype=float), np.asarray(y, dtype=float)])
    if np.isnan(xy).any():
        raise ValueError("nuclei without coordinates cannot be rotated; pass placed nuclei only")
    timepoints = np.asarray(timepoints).astype(str)
    rotated = np.empty_like(xy)
    for timepoint in np.unique(timepoints):
        if timepoint not in angles:
            raise KeyError(f"no display rotation for timepoint {timepoint!r}")
        rows = timepoints == timepoint
        rotated[rows] = rotate_about(xy[rows], angles[timepoint], array_centre(xy[rows]))
    return rotated[:, 0], rotated[:, 1]
