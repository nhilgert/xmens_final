#!/usr/bin/env python3
"""Step 1: pseudobulks of the X-Mens mouse object, for the decidual vs fibroblast DE.

Reads, without modifying:
  final_mouse.h5ad            the mouse object: obs and layers/counts
  xmens_feature_manifest.tsv  how the object's features map to Ensembl 116 genes
Both must be the files the release of 19 October 2026 used; their SHA-256 checksums are checked first.

Nuclei: every nucleus labelled "decidual cells" or "endometrial fibroblasts" in coarse_refined, so all decidual
subtypes (decidual_apoptotic included) and all fibroblast subtypes.

Pseudobulks: UMI counts summed over a timepoint's nuclei of one arm, across all of the timepoint's samples. Each
pseudobulk needs at least 20 nuclei (25 until 29 October 2026; the smallest has 1,170). Two designs (--design):
  with_uninduced  (default, adopted 27 October 2026) one decidual and one fibroblast pseudobulk for each of the five
                  induced timepoints, plus one fibroblast pseudobulk from the uninduced sample (9,270 nuclei): 11
                  pseudobulks. The uninduced sample's 64 decidual-labelled nuclei are left out; an uninduced uterus has
                  no decidual cells to compare.
  induced_only    the release of 19 October 2026: the five induced timepoints only, 10 pseudobulks. Kept so the
                  release can be rebuilt exactly.

Genes: the object's features are mapped to Ensembl 116 genes with the feature manifest. Features that the manifest
excludes, or that have no Ensembl gene, are dropped; features that map to the same gene are summed.

Writes into --out:
  xmens_broad_carrier_counts.tsv  summed UMIs, pseudobulks x genes (first column pseudobulk_id)
  xmens_broad_metadata.tsv        one row per pseudobulk: timepoint, group (arm), number of nuclei, samples
  xmens_broad_genes.tsv           one row per gene, in the column order of the counts: Ensembl gene ID, symbol,
                                  biotype, and the object's features that were summed into it
With --design induced_only, the first two have the same names and contents as the release's intermediate files.

Usage
  python mouse_pseudobulk_job.py --out OUT_DIR [--design with_uninduced|induced_only]
                                 [--object final_mouse.h5ad] [--manifest xmens_feature_manifest.tsv]
"""
import argparse
import hashlib
import os

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

DEFAULT_OBJECT = "/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/data/interim/final_mouse.h5ad"
DEFAULT_MANIFEST = ("/n/lab_storage/eddy_lab/users/nhilgert/cagri_ivd_end_to_end_20260815/resources/"
                    "xmens_feature_manifest.tsv")

# Checksums of the files the release used. Any other file stops the job.
EXPECTED_SHA256 = {
    "object": "7cfa0ebda34812fbfbc9ea4fbe23888e32997ad0bc8263634abc7f9f30ad7116",
    "manifest": "551649d3d30bee8a5a40515bd4b3c28e0426ff0fe1c6c1bbd97816733fbb79cf",
}

LABEL_COLUMN = "coarse_refined"
ARM = {"decidual cells": "decidual", "endometrial fibroblasts": "fibroblast"}
INDUCED_TIMEPOINTS = ["1_dpi_early", "1_dpi_late", "2_dpi", "3_dpi", "4_dpi_bleeding"]   # in order
UNINDUCED = "uninduced"
MIN_NUCLEI = 20   # the Hammoud (Ulrich) release's cell floor, also used for the human fit; 25 until 29 October 2026

# The pseudobulks of each design: (timepoint, arms), in output order
DESIGNS = {
    "with_uninduced": [(timepoint, ("decidual", "fibroblast")) for timepoint in INDUCED_TIMEPOINTS]
                      + [(UNINDUCED, ("fibroblast",))],
    "induced_only": [(timepoint, ("decidual", "fibroblast")) for timepoint in INDUCED_TIMEPOINTS],
}


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_true(value):
    """The manifest's flags may be booleans or text; a missing flag counts as False."""
    if value is None or pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def read_object(path):
    """obs, the feature names, and layers/counts as an integer sparse matrix (nuclei x features)."""
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        feature_names = read_elem(f["var"]).index.astype(str)
        counts = sp.csr_matrix(read_elem(f["layers/counts"]))
    assert obs.index.is_unique and feature_names.is_unique
    assert counts.shape == (len(obs), len(feature_names))
    if not np.issubdtype(counts.dtype, np.integer):
        assert np.array_equal(counts.data, np.round(counts.data)), "layers/counts holds non-integer values"
    assert counts.data.min() >= 0, "layers/counts holds negative values"
    return obs, feature_names, counts.astype(np.int64)


def map_features_to_genes(manifest, feature_names):
    """Which gene each of the object's features belongs to, and the table of genes.

    Returns gene_of_feature, the position of each feature's gene in the gene table (-1 for a dropped feature), and
    genes: one row per gene, in order of first appearance in the manifest, with the gene's symbol and biotype (left
    empty when the manifest gives none, or when the gene's features disagree) and the features summed into it."""
    manifest = manifest.sort_values("source_index", kind="stable").reset_index(drop=True)
    assert manifest.source_index.tolist() == list(range(len(manifest)))
    assert manifest.source_feature_id.astype(str).tolist() == list(feature_names), \
        "the manifest does not list the object's features in the object's order"

    kept = (~manifest.excluded.map(is_true)
            & manifest.canonical_id.notna()
            & manifest.canonical_id.astype(str).ne(""))
    kept_features = manifest[kept]
    gene_ids = list(dict.fromkeys(kept_features.canonical_id.astype(str)))      # in order of first appearance
    position_of_gene = {gene_id: position for position, gene_id in enumerate(gene_ids)}

    gene_of_feature = np.full(len(manifest), -1, dtype=np.int64)
    for feature_position, gene_id in zip(kept_features.source_index, kept_features.canonical_id.astype(str)):
        gene_of_feature[int(feature_position)] = position_of_gene[gene_id]

    rows = []
    for gene_id, features in kept_features.groupby("canonical_id", sort=False):
        symbols = [symbol for symbol in features.gene_symbol.dropna().astype(str).unique() if symbol]
        biotypes = [biotype for biotype in features.biotype.dropna().astype(str).unique() if biotype]
        rows.append({
            "ensembl_gene_id": str(gene_id),
            "gene_symbol": symbols[0] if len(symbols) == 1 else pd.NA,
            "biotype": biotypes[0] if len(biotypes) == 1 else pd.NA,
            "source_feature_ids": ";".join(features.source_feature_id.astype(str)),
            "source_feature_n": len(features),
        })
    genes = pd.DataFrame(rows)
    assert genes.ensembl_gene_id.tolist() == gene_ids
    return gene_of_feature, genes


def describe_pseudobulks(obs, design):
    """One row per pseudobulk (timepoint x arm) of the design, and the row numbers of each pseudobulk's nuclei."""
    arm = obs[LABEL_COLUMN].astype(str).map(ARM)          # missing for every other cell type
    timepoint = obs["timepoint"].astype(str)
    print(pd.crosstab(timepoint[arm.notna()], arm[arm.notna()]).to_string(), flush=True)

    descriptions, members = [], []
    for timepoint_name, arm_names in DESIGNS[design]:
        for arm_name in arm_names:
            rows = np.flatnonzero((timepoint == timepoint_name).to_numpy() & (arm == arm_name).to_numpy())
            samples = obs["sample"].iloc[rows]
            descriptions.append({
                "pseudobulk_id": f"timepoint={timepoint_name}__group={arm_name}",
                "timepoint": timepoint_name,
                "group": arm_name,
                "n_cells": len(rows),
                "n_samples": int(samples.nunique()),
                "samples": ";".join(sorted(samples.astype(str).unique())),
            })
            members.append(rows)
    pseudobulks = pd.DataFrame(descriptions)
    assert (pseudobulks.n_cells >= MIN_NUCLEI).all(), pseudobulks[["pseudobulk_id", "n_cells"]]
    return pseudobulks, members


def sum_pseudobulks(counts, members, gene_of_feature, n_genes):
    """UMIs summed over each pseudobulk's nuclei and over each gene's features (pseudobulks x genes)."""
    # membership[p, n] = 1 if nucleus n is in pseudobulk p
    pseudobulk_of_member = np.repeat(np.arange(len(members)), [len(rows) for rows in members])
    membership = sp.csr_matrix((np.ones(len(pseudobulk_of_member), dtype=np.int64),
                                (pseudobulk_of_member, np.concatenate(members))),
                               shape=(len(members), counts.shape[0]))
    # feature_to_gene[f, g] = 1 if feature f is summed into gene g; dropped features have no entry
    mapped = np.flatnonzero(gene_of_feature >= 0)
    feature_to_gene = sp.csr_matrix((np.ones(len(mapped), dtype=np.int64), (mapped, gene_of_feature[mapped])),
                                    shape=(counts.shape[1], n_genes))
    summed = (membership @ counts @ feature_to_gene).toarray()

    # Every UMI of a used nucleus in a mapped feature is counted exactly once.
    used = np.concatenate(members)
    assert summed.sum() == counts[used][:, mapped].sum()
    return summed


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True, help="output folder")
    parser.add_argument("--design", choices=sorted(DESIGNS), default="with_uninduced",
                        help="with_uninduced (default) or induced_only (the release of 19 October 2026)")
    parser.add_argument("--object", default=DEFAULT_OBJECT, help="the mouse h5ad object")
    parser.add_argument("--manifest", default=DEFAULT_MANIFEST, help="the feature manifest")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    for name, path in (("object", args.object), ("manifest", args.manifest)):
        assert sha256_of(path) == EXPECTED_SHA256[name], f"{path} is not the file the release used"

    obs, feature_names, counts = read_object(args.object)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype={"canonical_id": "string"})
    gene_of_feature, genes = map_features_to_genes(manifest, feature_names)

    pseudobulks, members = describe_pseudobulks(obs, args.design)
    summed = sum_pseudobulks(counts, members, gene_of_feature, len(genes))

    counts_table = pd.DataFrame(summed, columns=genes.ensembl_gene_id)
    counts_table.insert(0, "pseudobulk_id", pseudobulks.pseudobulk_id)
    # The release's TSV format
    tsv = dict(sep="\t", index=False, na_rep="NA", float_format="%.15g", lineterminator="\n")
    counts_table.to_csv(f"{args.out}/xmens_broad_carrier_counts.tsv", **tsv)
    pseudobulks.to_csv(f"{args.out}/xmens_broad_metadata.tsv", **tsv)
    # Empty fields for missing symbols or biotypes; step 2 reads empty fields, and only those, as missing.
    genes.to_csv(f"{args.out}/xmens_broad_genes.tsv", sep="\t", index=False, na_rep="", lineterminator="\n")

    print(pseudobulks[["pseudobulk_id", "n_cells", "n_samples"]].to_string(index=False))
    print(f"{len(genes)} genes from {int((gene_of_feature >= 0).sum())} of {len(feature_names)} features; "
          f"{int(pseudobulks.n_cells.sum())} nuclei; {int(summed.sum())} UMIs", flush=True)


if __name__ == "__main__":
    main()
