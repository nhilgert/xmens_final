"""The nuclei that the round-2 atlas table keeps, as a coordinates-only object for spatial_qc.py.

04c2 atlas_kept_coords_round2. Spatial QC (04f) reads the atlas-annotated object, i.e. every nucleus whose
whole-tissue cluster is not EXCLUDE in the atlas table, and uses only the coordinates, sample and timepoint. This
step selects the same nuclei under the round-2 atlas table (annotation_round2.atlas_decisions) from the July
whole-tissue fit (pinned by sha256) and writes those obs columns, plus the atlas cluster and lineage for the record.
With the July table instead of the round-2 table it selects exactly the nuclei of data/interim/annotated_{species}.h5ad.
"""
import hashlib
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "atlas"))
from apply_annotation import load_choices  # noqa: E402

snk = snakemake  # noqa: F821 (injected by Snakemake)
LOG = open(snk.log[0], "w")


def say(msg: str):
    print(f"[atlas_kept_coords:{snk.params.species}] {msg}", file=LOG, flush=True)


def sha256_of(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    species = str(snk.params.species)
    sqc = dict(snk.params.sqc)
    columns = [sqc["x_key"], sqc["y_key"], sqc["sample_key"], sqc["timepoint_key"]]
    found = sha256_of(snk.input.atlas)
    if found != str(snk.params.atlas_sha256):
        raise ValueError(f"{snk.input.atlas}: sha256 {found} != pinned {snk.params.atlas_sha256}")
    spec = yaml.safe_load(open(snk.input.spec))
    cluster_column = spec["atlas"][species]["res"]

    atlas = ad.read_h5ad(snk.input.atlas, backed="r")
    obs = atlas.obs[columns + [cluster_column]].copy()
    atlas.file.close()
    decisions = load_choices(snk.input.decisions)
    cluster = obs[cluster_column].astype(str)
    missing = sorted(set(cluster) - set(decisions))
    if missing:
        raise ValueError(f"{snk.input.decisions}: no row for atlas clusters {missing}")
    lineage = cluster.map(decisions)
    kept = (lineage != "EXCLUDE").to_numpy()
    say(f"{len(obs):,} nuclei in {Path(snk.input.atlas).name}; {int(kept.sum()):,} kept by {snk.input.decisions} "
        f"({int((~kept).sum()):,} in EXCLUDE clusters)")

    out = obs.loc[kept, columns].copy()
    out["atlas_cluster"] = cluster[kept].values
    out["lineage"] = pd.Categorical(lineage[kept].values)
    coords = ad.AnnData(obs=out)
    coords.obsm["spatial"] = out[[sqc["x_key"], sqc["y_key"]]].to_numpy(np.float64)
    path = Path(snk.output.h5ad)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp.h5ad")
    coords.write_h5ad(tmp)
    os.replace(tmp, path)
    say(f"wrote {coords.n_obs:,} nuclei -> {path}; lineages {out['lineage'].value_counts().to_dict()}")


main()
LOG.close()
