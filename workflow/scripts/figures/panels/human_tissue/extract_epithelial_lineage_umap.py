"""Extract the epithelial lineage fit's UMAP for the released human epithelial nuclei (Fig. S26C).

Called by rule `panels_human_epithelium_lineage_umap`, once per release.
Reads:  the human epithelial lineage fit (obs and obsm X_umap_scVI only, through h5py), and the plot
        table and uns export of the release written by extract_human_tissue.py.
Writes: lineage_umap.parquet  one row per released epithelial nucleus: lineage-fit UMAP and cluster
        provenance.json       fit sha256 (checked against the release's own record of it) and counts
"""

import argparse
import datetime
import hashlib
import json
import logging
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem


# The clustering on which the human epithelium was split into inner and outer layers
# (annotation vocabulary v7, content task human_epithelial_split), and its embedding.
CLUSTER_KEY = "leiden_scVI_3.0"
UMAP_KEY = "X_umap_scVI"

# The release records which lineage fit its epithelial labels come from, with the fit's sha256.
RELEASE_RECORD_KEY = "epithelial_reinstatement"

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("extract_epithelial_lineage_umap")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fit(path: Path) -> pd.DataFrame:
    """UMAP coordinates and cluster of every nucleus of the lineage fit, indexed by obs name."""
    with h5py.File(path, "r") as h5:
        obs = read_elem(h5["obs"])
        umap = h5["obsm"][UMAP_KEY][:]
    return pd.DataFrame(
        {
            "lineage_umap_1": umap[:, 0],
            "lineage_umap_2": umap[:, 1],
            "lineage_cluster": obs[CLUSTER_KEY].astype(str).to_numpy(),
        },
        index=obs.index.astype(str),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit", required=True, type=Path)
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--uns", required=True, type=Path)
    parser.add_argument("--release", required=True)
    parser.add_argument("--out-table", required=True, type=Path)
    parser.add_argument("--out-provenance", required=True, type=Path)
    args = parser.parse_args()

    # Fit identity: the bytes must be the fit the release says its epithelial clusters come from ----
    record = json.loads(args.uns.read_text())[RELEASE_RECORD_KEY]
    if record["cluster_column"] != CLUSTER_KEY:
        raise ValueError(f"release records cluster column {record['cluster_column']!r}, expected {CLUSTER_KEY!r}")
    fit_sha = sha256(args.fit)
    if fit_sha != record["clustering_sha256"]:
        raise ValueError(f"{args.fit} sha256 {fit_sha} differs from the release record {record['clustering_sha256']}")
    log.info(f"{args.fit.name}: sha256 matches uns['{RELEASE_RECORD_KEY}'] of release {args.release}")

    fit = read_fit(args.fit)

    # Every released epithelial nucleus must have a lineage-fit position; fit nuclei that the release
    # does not carry as epithelium (clusters labelled EXCLUDE) are counted and left out.
    cells = pd.read_parquet(args.cells, columns=["cell_type", "cell_subtype"])
    epithelial = cells[cells["cell_type"].astype(str) == "epithelial_cell"]
    missing = epithelial.index.difference(fit.index)
    if len(missing):
        raise ValueError(f"{len(missing)} released epithelial nuclei are not in {args.fit}: {list(missing[:5])}")
    table = fit.loc[epithelial.index]
    not_released = fit.index.difference(epithelial.index)
    log.info(
        f"{len(fit):,} fit nuclei; {len(table):,} released epithelial nuclei kept; "
        f"{len(not_released):,} fit nuclei not released as epithelium, left out"
    )

    # Cluster-to-subtype check: the layers were assigned per cluster, so each cluster should carry one subtype.
    subtype = epithelial["cell_subtype"].astype(str)
    cluster_subtypes = pd.crosstab(table["lineage_cluster"], subtype)
    subtypes_per_cluster = cluster_subtypes.gt(0).sum(axis=1)

    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    table.index.name = "cell"
    table.to_parquet(args.out_table)

    stat = args.fit.stat()
    provenance = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "release": args.release,
        "fit": {
            "path": str(args.fit.resolve()),
            "bytes": stat.st_size,
            "mtime": datetime.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
            "sha256": fit_sha,
            "sha256_matches_release_record": True,
            "release_record": record,
            "embedding": f"obsm {UMAP_KEY}, as stored",
            "cluster_column": CLUSTER_KEY,
        },
        "plot_table": {"path": str(args.cells.resolve()), "sha256": sha256(args.cells)},
        "script": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "n_fit_nuclei": int(len(fit)),
        "n_released_epithelial": int(len(table)),
        "released_epithelial_by_subtype": subtype.value_counts().to_dict(),
        "n_fit_nuclei_not_released": int(len(not_released)),
        "fit_nuclei_not_released_by_cluster": fit.loc[not_released, "lineage_cluster"].value_counts().to_dict(),
        "clusters_by_subtype": cluster_subtypes.idxmax(axis=1).value_counts().to_dict(),
        "clusters_with_more_than_one_subtype": cluster_subtypes.index[subtypes_per_cluster > 1].tolist(),
    }
    args.out_provenance.write_text(json.dumps(provenance, indent=1))
    log.info(f"wrote {args.out_table} ({len(table):,} rows) and {args.out_provenance}")


if __name__ == "__main__":
    main()
