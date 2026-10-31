"""Write the plot tables behind the mouse whole-tissue panels (Fig. 3A-D, Figs. S12-S13).

Called by rule panels_mouse_tissue_extract (workflow/rules/91_panels_mouse_tissue.smk), once per
release version. Reads the release's Cevrim_XMens_All_Cells.h5ad (obs, obsm X_umap_scVI, and X and
layers/counts in row blocks, never the whole matrix), the release MANIFEST.json, the immune lineage
fit drawn in Fig. 3C (obsm X_umap_scVI) and config/samples.tsv.
Writes to --out-dir: cells.csv.gz (one row per released nucleus), marker_stats.csv (cell class x
gene), lanes.csv (10X lanes per stage) and extract_provenance.json.
"""

import argparse
import hashlib
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from mouse_tissue_style import DOTPLOT_GENES, PERIVASCULAR_GENES


# Rows of the CSR matrices read per block: about 1 GB of layers/counts per block at most.
BLOCK_ROWS = 8192

OBS_COLUMNS = [
    "sample",
    "timepoint",
    "cell_type",
    "cell_subtype",
    "tissue_section",
    "spatial_qc_reason",
    "x_um_dbscan",
    "y_um_dbscan",
    "x_um_qc",
    "y_um_qc",
    "embedding_outlier",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", required=True, type=Path)
    parser.add_argument("--object-name", required=True)
    parser.add_argument("--immune-fit", required=True, type=Path)
    parser.add_argument("--immune-fit-sha256", default="", help="expected sha256; empty = record only")
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 24), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_against_manifest(h5ad: Path, manifest: dict, object_name: str) -> str:
    """Return the object's sha256 after checking it against the release MANIFEST."""
    expected = manifest["objects"][object_name]["sha256"]
    observed = sha256_of(h5ad)
    if observed != expected:
        raise ValueError(f"{h5ad}: sha256 {observed} differs from MANIFEST {expected}")
    return observed


def block_rows(indptr: np.ndarray):
    """Yield (start, stop, row index within block per stored value) for CSR row blocks."""
    n_rows = len(indptr) - 1
    for start in range(0, n_rows, BLOCK_ROWS):
        stop = min(start + BLOCK_ROWS, n_rows)
        row_of = np.repeat(np.arange(stop - start), np.diff(indptr[start:stop + 1]))
        yield start, stop, row_of


def count_metrics(counts: h5py.Group) -> pd.DataFrame:
    """Total counts and detected genes per nucleus from the raw counts layer."""
    if counts.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError(f"layers/counts is {counts.attrs.get('encoding-type')!r}, expected csr_matrix")
    indptr = counts["indptr"][:]
    n_rows = len(indptr) - 1
    total = np.zeros(n_rows)
    n_genes = np.zeros(n_rows, dtype=np.int64)
    for start, stop, row_of in block_rows(indptr):
        data = counts["data"][indptr[start]:indptr[stop]]
        total[start:stop] = np.bincount(row_of, weights=data, minlength=stop - start)
        n_genes[start:stop] = np.bincount(row_of[data > 0], minlength=stop - start)
    return pd.DataFrame({"total_counts": total, "n_genes": n_genes})


def gene_columns(matrix: h5py.Group, gene_idx: np.ndarray) -> np.ndarray:
    """Dense (n_obs x n_genes) values of the selected columns of a CSR matrix."""
    if matrix.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError(f"X is {matrix.attrs.get('encoding-type')!r}, expected csr_matrix")
    indptr = matrix["indptr"][:]
    column_slot = -np.ones(matrix.attrs["shape"][1], dtype=np.int64)
    column_slot[gene_idx] = np.arange(len(gene_idx))
    values = np.zeros((len(indptr) - 1, len(gene_idx)), dtype=np.float32)
    for start, stop, row_of in block_rows(indptr):
        lo, hi = indptr[start], indptr[stop]
        slots = column_slot[matrix["indices"][lo:hi]]
        hit = slots >= 0
        values[start + row_of[hit], slots[hit]] = matrix["data"][lo:hi][hit]
    return values


def marker_stats(expr_lognorm: np.ndarray, genes: list, cell_type: pd.Series) -> pd.DataFrame:
    """Fraction of nuclei detecting each gene and its mean log-normalised expression, per class.

    X of the release is log1p of normalised counts, so a value > 0 is the same as a raw count > 0.
    """
    expr = pd.DataFrame(expr_lognorm, columns=genes, index=cell_type.index)
    means = expr.groupby(cell_type.to_numpy()).mean()
    detected = (expr > 0).groupby(cell_type.to_numpy()).mean()
    stats = (
        means.stack()
        .rename("mean_lognorm")
        .to_frame()
        .join(detected.stack().rename("frac_detected"))
        .reset_index(names=["cell_type", "gene"])
    )
    stats["n_nuclei"] = stats["cell_type"].map(cell_type.value_counts())
    return stats


def immune_umap(fit_path: Path, immune_names: pd.Index) -> pd.DataFrame:
    """Fig. 3C coordinates: the immune lineage fit's UMAP for the released immune nuclei."""
    with h5py.File(fit_path, "r") as fit:
        fit_names = pd.Index(read_elem(fit["obs"]).index)
        coords = fit["obsm"]["X_umap_scVI"][:]
    missing = immune_names.difference(fit_names)
    if len(missing):
        raise ValueError(f"{len(missing)} released immune nuclei are not in {fit_path}: {list(missing[:5])}")
    table = pd.DataFrame(coords, index=fit_names, columns=["immune_umap_1", "immune_umap_2"])
    logging.info(
        "Immune fit %s: %d nuclei, %d released immune nuclei drawn, %d fit nuclei not released",
        fit_path,
        len(fit_names),
        len(immune_names),
        len(fit_names) - len(immune_names),
    )
    return table.loc[immune_names]


def lanes_per_stage(samples_path: Path, obs: pd.DataFrame) -> pd.DataFrame:
    """10X lanes per stage: one sample-sheet row per lane, mapped to the release timepoint."""
    sheet = pd.read_csv(samples_path, sep="\t")
    sheet = sheet[(sheet["species"] == "mouse") & (sheet["modality"] == "slidetags")]
    stages_of_sample = (
        obs.astype({"sample": str, "timepoint": str})
        .groupby("sample")["timepoint"]
        .agg(lambda s: sorted(s.unique()))
    )
    ambiguous = stages_of_sample[stages_of_sample.map(len) != 1]
    if len(ambiguous):
        raise ValueError(f"samples spanning several timepoints: {ambiguous.to_dict()}")
    unreleased = set(sheet["sample_id"]) - set(stages_of_sample.index)
    if unreleased:
        raise ValueError(f"sample-sheet lanes without released nuclei: {sorted(unreleased)}")
    sheet["timepoint"] = sheet["sample_id"].map(stages_of_sample.str[0])
    return (
        sheet.groupby("timepoint")["sample_id"]
        .agg(n_lanes="count", samples=lambda s: ";".join(sorted(s)))
        .reset_index()
    )


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Fig. S12B genes: the manuscript's markers and the v9 perivascular marker (drawn only when
    # the release has perivascular cells).
    genes = DOTPLOT_GENES + PERIVASCULAR_GENES

    # Release object, checked against its MANIFEST so the panels name the exact bytes they drew.
    h5ad = args.release_dir / f"{args.object_name}.h5ad"
    manifest_path = args.release_dir / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    object_sha256 = verify_against_manifest(h5ad, manifest, args.object_name)

    fit_sha256 = sha256_of(args.immune_fit)
    if args.immune_fit_sha256 and fit_sha256 != args.immune_fit_sha256:
        raise ValueError(f"{args.immune_fit}: sha256 {fit_sha256} differs from the pinned {args.immune_fit_sha256}")

    with h5py.File(h5ad, "r") as release:
        obs = read_elem(release["obs"])[OBS_COLUMNS]
        var_names = pd.Index(read_elem(release["var"]).index)
        umap = release["obsm"]["X_umap_scVI"][:]
        qc = count_metrics(release["layers"]["counts"])
        missing_genes = [g for g in genes if g not in var_names]
        if missing_genes:
            raise ValueError(f"marker genes not in var: {missing_genes}")
        expr = gene_columns(release["X"], var_names.get_indexer(genes))
    logging.info("Read %s: %d nuclei, %d genes", h5ad, len(obs), len(var_names))

    # One row per released nucleus. Spatially mapped = has Slide-tags coordinates at all, before
    # spatial QC: the definition behind the manuscript's S12A column.
    cells = obs.astype({c: str for c in ["sample", "timepoint", "cell_type", "cell_subtype", "spatial_qc_reason"]})
    cells["has_coordinates"] = cells["x_um_dbscan"].notna() & cells["y_um_dbscan"].notna()
    disagree = int((cells["has_coordinates"] != (cells["spatial_qc_reason"] != "unmapped")).sum())
    if disagree:
        raise ValueError(f"{disagree} nuclei: coordinates disagree with spatial_qc_reason 'unmapped'")
    cells["umap_1"] = umap[:, 0]
    cells["umap_2"] = umap[:, 1]
    cells["total_counts"] = qc["total_counts"].to_numpy()
    cells["n_genes"] = qc["n_genes"].to_numpy()

    immune_names = cells.index[cells["cell_type"] == "immune_cell"]
    cells = cells.join(immune_umap(args.immune_fit, immune_names))
    cells = cells.drop(columns=["x_um_dbscan", "y_um_dbscan"])
    cells.index.name = "obs_name"

    stats = marker_stats(expr, genes, cells["cell_type"])
    lanes = lanes_per_stage(args.samples, obs)

    cells.to_csv(args.out_dir / "cells.csv.gz")
    stats.to_csv(args.out_dir / "marker_stats.csv", index=False)
    lanes.to_csv(args.out_dir / "lanes.csv", index=False)

    provenance = {
        "rule": "panels_mouse_tissue_extract",
        "version": args.version,
        "release": manifest["release"],
        "vocabulary_schema": manifest.get("vocabulary_schema"),
        "inputs": {
            "release_object": {"path": str(h5ad), "sha256": object_sha256, "matches_manifest": True},
            "manifest": {"path": str(manifest_path), "sha256": sha256_of(manifest_path)},
            "immune_fit": {
                "path": str(args.immune_fit),
                "sha256": fit_sha256,
                "pinned_sha256_checked": bool(args.immune_fit_sha256),
                "embedding": "obsm X_umap_scVI",
            },
            "samples": {"path": str(args.samples), "sha256": sha256_of(args.samples)},
        },
        "embedding": "obsm X_umap_scVI of the release object, as stored",
        "expression": "dot-plot means from X (log1p of normalised counts); detection = X > 0",
        "qc_metrics": "total_counts and n_genes summed from layers/counts",
        "counts": {
            "n_nuclei": int(len(cells)),
            "n_embedding_outlier": int(cells["embedding_outlier"].sum()),
            "n_with_coordinates": int(cells["has_coordinates"].sum()),
            "spatial_qc_reason": cells["spatial_qc_reason"].value_counts().to_dict(),
            "per_timepoint": cells["timepoint"].value_counts().to_dict(),
            "per_cell_type": cells["cell_type"].value_counts().to_dict(),
            "immune_per_subtype": cells.loc[immune_names, "cell_subtype"].value_counts().to_dict(),
            "immune_nuclei_with_fit_umap": int(cells["immune_umap_1"].notna().sum()),
        },
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    (args.out_dir / "extract_provenance.json").write_text(json.dumps(provenance, indent=2, default=int))
    logging.info("Wrote cells.csv.gz (%d rows), marker_stats.csv, lanes.csv to %s", len(cells), args.out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    main()
