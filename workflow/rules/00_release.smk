# 00_release -- the release of record and the releases before it, for every rule file that names a release.
#
# config `releases` maps each version label to its release folder (data/deliverables/<folder>), and
# deliverables.version / deliverables.release name the release of record. Cutting the next release means setting
# those two keys; every reader takes release folders and version labels from the globals below, never from a
# typed date.
#
# A release is cut once its MANIFEST.json exists. A cut release is immutable (sha256 manifest, protected outputs),
# so the rules that build the release of record (04c2 round-2 annotation, 04g2, finalize_deliverables) are defined
# only while it is not cut yet. Once cut, its files are plain inputs of the DAG: a later code or config edit cannot
# schedule a rebuild into its folders, and a change to a cut release is made by cutting a new one.

import re
from pathlib import Path


def version_rank(version: str) -> int:
    """Sort key of a version label: v8 < v9 < v10."""
    if not re.fullmatch(r"v\d+", str(version)):
        raise ValueError(f"config releases: version labels are v<number>, not {version!r}")
    return int(str(version)[1:])


def release_is_cut(release: str) -> bool:
    """True once finalize_deliverables has written the release's manifest."""
    return Path(f"data/deliverables/{release}/MANIFEST.json").exists()


RELEASES = {str(version): str(folder) for version, folder in config["releases"].items()}
RELEASE_VERSION = str(config["deliverables"]["version"])
RELEASE = str(config["deliverables"]["release"])

for _version in RELEASES:
    version_rank(_version)
if RELEASE_VERSION not in RELEASES:
    raise ValueError(f"config deliverables.version {RELEASE_VERSION} is not a key of config releases {RELEASES}")
if RELEASES[RELEASE_VERSION] != RELEASE:
    raise ValueError(
        f"config deliverables.release is {RELEASE}, but config releases maps {RELEASE_VERSION} to "
        f"{RELEASES[RELEASE_VERSION]}: set both to the release being cut"
    )
if len(set(RELEASES.values())) != len(RELEASES):
    raise ValueError(f"config releases: two versions share a release folder: {RELEASES}")

RELEASE_DIR = f"data/deliverables/{RELEASE}"
RELEASE_MANIFEST = f"{RELEASE_DIR}/MANIFEST.json"
RELEASE_IS_CUT = release_is_cut(RELEASE)

# The release the release of record follows: the newest cut release of an earlier version (None for the first).
_earlier_cut = sorted(
    (v for v in RELEASES if version_rank(v) < version_rank(RELEASE_VERSION) and release_is_cut(RELEASES[v])),
    key=version_rank,
)
PREVIOUS_VERSION = _earlier_cut[-1] if _earlier_cut else None
PREVIOUS_RELEASE = RELEASES[PREVIOUS_VERSION] if PREVIOUS_VERSION else None

# Folder of the release's 04c2/04g2 intermediates: release-stamped, on lab_storage, unless the version recorded
# another folder when it was cut (release_round2.outdir_of_version).
_RR2_SETTINGS = config.get("release_round2") or {}
RELEASE_BUILD_DIR = (_RR2_SETTINGS.get("outdir_of_version") or {}).get(
    RELEASE_VERSION, str(_RR2_SETTINGS.get("outdir", "data/interim/release_{release}")).format(release=RELEASE)
)

# Every rule file that reads the release of record registers what it reads here ({reader: release folder}); the
# guard at the end of workflow/Snakefile fails the parse when a reader names another release.
RELEASE_READERS = {}
