"""Write the pseudotime of record into a fibroblast/decidual release object (decision G).

Used by refresh_fibdec_round2.py (04g2 refresh_fibdec_round2) and by its unit test. The run's table (config
pseudotime_of_record.<species>.path, pinned by sha256) holds one row per nucleus, indexed by barcode: the rank
pseudotime (`column`), the raw pseudotime (`raw_column`) and the diffusion components dc0, dc1, ...; run_record is the
run's provenance.json (it names the root nucleus) and script the exact script the run executed. The object takes,
under its release names (the deposit renames them to pseudotime, X_diffmap_pseudotime and uns pseudotime):
  obs  pt_harmonized_rank <- column; pt_harmonized <- raw_column              (float64)
  obsm X_diffmap_harmonized <- dc0 ... dc<n-1>, n as in the object           (float64)
  uns  pt_harmonized <- the run record, laid out as in the deposit copies (prune_release_object.py)
It fails unless the table matches its sha256, holds exactly the object's nuclei and puts the run's root at the lowest
rank.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

RANK_OBS = "pt_harmonized_rank"
RAW_OBS = "pt_harmonized"
DIFFMAP_OBSM = "X_diffmap_harmonized"
RECORD_UNS = "pt_harmonized"
CHUNK_BYTES = 1 << 24


def file_sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def read_run_table(record: dict, obs_names: pd.Index) -> pd.DataFrame:
    """The run's table, checked against its sha256 and the object's nuclei, in the object's nucleus order."""
    table_path = Path(record["path"])
    found = file_sha256(table_path)
    if found != record["sha256"]:
        raise ValueError(f"{table_path}: sha256 {found} is not the configured {record['sha256']}")
    table = pd.read_parquet(table_path)
    if not table.index.is_unique:
        raise ValueError(f"{table_path}: barcodes are not unique")
    missing, extra = obs_names.difference(table.index), table.index.difference(obs_names)
    if len(missing) or len(extra):
        raise ValueError(f"{table_path}: {len(missing)} nuclei of the object are not in the table and {len(extra)} "
                         "table rows are not in the object")
    return table.loc[obs_names]


def apply_pseudotime_of_record(adata, record: dict) -> dict:
    """Replace the object's pseudotime, diffusion map and their uns record by the run's; returns a summary."""
    table = read_run_table(record, adata.obs_names)
    run = json.loads(Path(record["run_record"]).read_text())
    script_path = Path(record["script"])
    n_dcs = np.asarray(adata.obsm[DIFFMAP_OBSM]).shape[1]
    release_record = dict(adata.uns[RECORD_UNS])
    release_rank = adata.obs[RANK_OBS].copy()

    adata.obs[RANK_OBS] = table[record["column"]].to_numpy(dtype="float64")
    adata.obs[RAW_OBS] = table[record["raw_column"]].to_numpy(dtype="float64")
    adata.obsm[DIFFMAP_OBSM] = table[[f"dc{i}" for i in range(n_dcs)]].to_numpy(dtype="float64")

    # Same keys and values as the deposit's substitution (prune_release_object.replace_pseudotime), so the deposit
    # copies the release's record unchanged.
    adata.uns[RECORD_UNS] = {
        "description": f"{release_record['description']}. {record['note']}",
        "n_cells": int(len(table)),
        "n_dcs": int(n_dcs),
        "provenance_json": json.dumps(run),
        "root": run["root"],
        "run": release_record["run"],
        "script": script_path.name,
        "script_sha256": file_sha256(script_path),
        "table": str(record["path"]),
        "table_sha256": record["sha256"],
        **{key: record[key] for key in ("latent_copy", "latent_sha256") if key in record},
    }

    rank = adata.obs[RANK_OBS]
    if rank.loc[run["root"]] != rank.min():
        raise ValueError(f"the run's root {run['root']} is not at the lowest rank of {record['path']}")
    return {
        "n_cells": int(len(table)),
        "n_dcs": int(n_dcs),
        "root": run["root"],
        "spearman_with_release": float(release_rank.corr(rank, method="spearman")),
    }
