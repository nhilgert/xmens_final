"""Build the input object of one MI pseudotime-of-record arm (rule mi_record_input, 11_mi_pseudotime.smk).

Reads a fibroblast/decidual object and the per-nucleus pseudotime table of record (parquet indexed by nucleus name,
config pseudotime_of_record), checks the table's sha256, and writes a copy of the object whose obs[--pt-key] holds
the table's --column. The object and the table must hold the same nuclei. Writes the h5ad (--out) and a JSON record
of the build (--meta).
"""
import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--h5ad", required=True, type=Path, help="fibroblast/decidual object")
parser.add_argument("--table", required=True, type=Path, help="pseudotime table of record (parquet)")
parser.add_argument("--column", required=True, help="table column with the pseudotime of record")
parser.add_argument("--expect-sha256", required=True, help="sha256 of the table pinned in config")
parser.add_argument("--pt-key", required=True, help="obs column that receives the pseudotime of record")
parser.add_argument("--out", required=True, type=Path, help="output h5ad")
parser.add_argument("--meta", required=True, type=Path, help="output JSON record of the build")
args = parser.parse_args()

table_sha256 = sha256_of(args.table)
if table_sha256 != args.expect_sha256:
    raise SystemExit(f"{args.table}: sha256 {table_sha256} differs from the pinned {args.expect_sha256}")

adata = ad.read_h5ad(args.h5ad)
table = pd.read_parquet(args.table)
table.index = table.index.astype(str)

# Every nucleus of the object needs a pseudotime of record, and the record may not describe other nuclei.
missing = adata.obs_names.difference(table.index)
extra = table.index.difference(adata.obs_names)
if len(missing) or len(extra):
    raise SystemExit(
        f"object and table disagree: {len(missing)} nuclei missing from the table, "
        f"{len(extra)} table rows not in the object"
    )

pt_record = table.loc[adata.obs_names, args.column].to_numpy(dtype=float)
pt_object = adata.obs[args.pt_key].to_numpy(dtype=float) if args.pt_key in adata.obs else None
adata.obs[args.pt_key] = pt_record
adata.uns["mi_pseudotime_of_record"] = {
    "table": str(args.table),
    "column": args.column,
    "sha256": table_sha256,
}
adata.write_h5ad(args.out)

n_changed = None if pt_object is None else int((pt_object != pt_record).sum())
meta = {
    "h5ad": str(args.h5ad),
    "table": str(args.table),
    "column": args.column,
    "table_sha256": table_sha256,
    "pt_key": args.pt_key,
    "n_cells": int(adata.n_obs),
    "n_finite_pt": int(np.isfinite(pt_record).sum()),
    "n_changed_vs_object": n_changed,
}
args.meta.write_text(json.dumps(meta, indent=2) + "\n")
print(json.dumps(meta), flush=True)
