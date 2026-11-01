# workflow/rules/91_panels_versions.smk
# Release versions the figure panels draw, read from the release interface of config/config.yaml:
#   releases               version tag -> release folder under data/deliverables/ (v9 -> 2026_10_30, ...)
#   deliverables.version   tag of the release of record; deliverables.release is its folder
# and from config/panels.yaml (read here):
#   panels.versions_drawn  tags whose panels the *_all targets ask for
#   panels.python          interpreter of every panel rule; panels.cache_dir for matplotlib and numba caches
#
# Tags after the release of record are skipped: their release is not cut yet. The release of record is the live
# tag: its panels follow the code, as any Snakemake output does. Earlier tags are frozen: panel_targets() asks only
# for their missing files, so an edit made for a new release never redraws, and never rewrites, an earlier panel.
# Families whose files carried no version tag before v10 (cross-species, SAMap, zone modules, human epithelium,
# marker dot plots, mouse stroma S16-S22, review notebooks) draw the live tag only; their earlier files keep the
# paths they were written to.
#
# Included by every panel rule file; Snakemake reads it once.

import os

configfile: "config/panels.yaml"

PANEL_RELEASES = config["releases"]
PANEL_LIVE_TAG = config["deliverables"]["version"]
PANEL_LIVE_RELEASE = config["deliverables"]["release"]
PANEL_PY = config["panels"]["python"]
PANEL_CACHE = config["panels"]["cache_dir"]

if PANEL_RELEASES.get(PANEL_LIVE_TAG) != PANEL_LIVE_RELEASE:
    raise ValueError(
        f"config deliverables: version {PANEL_LIVE_TAG!r} is release {PANEL_RELEASES.get(PANEL_LIVE_TAG)!r} in "
        f"config releases, but deliverables.release is {PANEL_LIVE_RELEASE!r}"
    )

# Tags in release order (the order of config releases), up to and including the release of record.
_PANEL_RANK = {tag: rank for rank, tag in enumerate(PANEL_RELEASES)}
PANEL_TAGS = [
    tag for tag in config["panels"]["versions_drawn"] if _PANEL_RANK[tag] <= _PANEL_RANK[PANEL_LIVE_TAG]
]
PANEL_FROZEN_TAGS = [tag for tag in PANEL_TAGS if tag != PANEL_LIVE_TAG]

# Matplotlib and numba caches default to HOME, which is full; HDF5 locking fails on lab_storage reads.
PANEL_ENV = (
    "export XDG_CACHE_HOME={c} MPLCONFIGDIR={c}/mpl NUMBA_CACHE_DIR={c}/numba HDF5_USE_FILE_LOCKING=FALSE; "
    "mkdir -p {c}/mpl {c}/numba"
).format(c=PANEL_CACHE)


def panel_release_dir(tag):
    """data/deliverables/<release folder> of a version tag."""
    return f"data/deliverables/{PANEL_RELEASES[tag]}"


def panel_targets(pattern, **wildcards):
    """pattern expanded over the drawn tags ({tag}) and the other wildcards; for a frozen tag only missing files."""
    paths = []
    for tag in PANEL_TAGS:
        tag_paths = expand(pattern, tag=tag, **wildcards)
        if tag in PANEL_FROZEN_TAGS:
            tag_paths = [path for path in tag_paths if not os.path.exists(path)]
        paths += tag_paths
    return paths


def panel_tags_up_to_live():
    """Every tag of config releases up to the release of record, oldest first."""
    return [tag for tag in PANEL_RELEASES if _PANEL_RANK[tag] <= _PANEL_RANK[PANEL_LIVE_TAG]]
