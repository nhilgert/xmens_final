"""Paths of the released objects that the methods-number registry reads.

The release is the one config/config.yaml names under deliverables.release. Its directory, data/deliverables/<release>/,
holds MANIFEST.json and the objects; object paths are taken from the manifest's file names. From 2026_10_29b each
species also has an extended all-cells object (the primary nuclei plus the supplementary populations); the registry's
numbers come from the primary objects, and the extended ones are listed under their own tags. Used by collect_facts.py
and collect_facts2.py; writes nothing.
"""
import json

import yaml

# Registry tag of each released object, keyed by the object's name in MANIFEST.json.
OBJECT_TAGS = {
    "Cevrim_XMens_All_Cells": "mouse_all",
    "Cevrim_XMens_Only_Fib_Dec_Cells": "mouse_fibdec",
    "Cevrim_Human_All_Cells": "human_all",
    "Cevrim_Human_Only_Fib_Dec_Cells": "human_fibdec",
    "Cevrim_XMens_All_Cells_extended": "mouse_all_extended",       # 2026_10_29b on; not read by the registry
    "Cevrim_Human_All_Cells_extended": "human_all_extended",
}
REQUIRED_TAGS = ("mouse_all", "mouse_fibdec", "human_all", "human_fibdec")


def read_release(repo):
    """Release name, directory, manifest and object paths (keyed by registry tag) of the release config.yaml names."""
    config = yaml.safe_load(open(f"{repo}/config/config.yaml"))
    release = str(config["deliverables"]["release"])
    release_dir = f"{repo}/data/deliverables/{release}"
    manifest = json.load(open(f"{release_dir}/MANIFEST.json"))

    if manifest["release"] != release:
        raise ValueError(f"{release_dir}/MANIFEST.json describes release {manifest['release']}, config names {release}")

    paths = {OBJECT_TAGS[name]: f"{release_dir}/{entry['file']}" for name, entry in manifest["objects"].items()}
    missing = sorted(set(REQUIRED_TAGS) - set(paths))
    if missing:
        raise KeyError(f"release {release} lacks the objects {missing}")

    return {"release": release, "dir": release_dir, "manifest": manifest, "paths": paths}
