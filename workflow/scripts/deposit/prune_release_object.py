"""Deposit copy of one release object: obs, uns and obsm pruned and renamed as configured; nothing else changes.

Called by rule deposit_object (workflow/rules/95_deposit_objects.smk).
Reads:  a release object (data/deliverables/<release>/<object>.h5ad), config/deposit_objects.yaml and, for the mouse
        objects, the zone-module AUCell table of the run of record (results/zone_aucell/<run>/zone_aucell.csv, with the
        run's provenance.json next to it). For an object listed under pseudotime_of_record: that run's pseudotime table,
        run record and script.
Writes: the deposit copy (.h5ad, same compression as the release) and a provenance JSON that lists every column and
        entry removed, renamed or replaced, and the checks below.

Changes, in this order:
  1. obs: the configured columns are removed and the configured renames applied (pt_harmonized_rank -> pseudotime).
  2. obs, mouse objects: every auc_* column of the release (an earlier AUCell run) is replaced by the run of record's
     scores for the configured modules on the nuclei of the configured timepoint; every other nucleus gets NaN.
  3. obsm: the configured embeddings are removed, and the configured renames applied
     (X_diffmap_harmonized -> X_diffmap_pseudotime, the diffusion map of the pseudotime kept).
  4. uns: the configured entries are removed or renamed; a sub-entry naming a removed embedding goes with it, and one
     naming a renamed embedding follows the new name. The mouse objects gain uns['zone_aucell'], the record of the
     scores kept.
  5. Pseudotime of record (objects listed under pseudotime_of_record): obs 'pseudotime', obsm 'X_diffmap_pseudotime'
     and uns 'pseudotime' are replaced by the run's rank pseudotime (pt_rank), its diffusion components and a record
     built from its run record, because the release's values predate the run.
Checks on the written file: same nuclei in the same order, every kept obs column identical to the release, pseudotime
identical to the release's pt_harmonized_rank (for a pseudotime of record: to the run's pt_rank, with the run's root
at the lowest rank), scores present exactly on the timepoint's nuclei, each renamed embedding identical to its release
original (for a pseudotime of record: the diffusion map identical to the run's components), and X, layers, obsp, varm,
var and the other obsm identical in value to the release.

usage: prune_release_object.py --input H5AD --release R --settings YAML --zone-aucell CSV|'' --output H5AD
                               --provenance JSON
"""
import argparse
import datetime
import hashlib
import importlib.metadata
import json
import logging
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml

LOG = logging.getLogger("prune_release_object")
CHUNK = 1 << 24  # rows per read when hashing large arrays


# Hashing ----

def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def array_digests(path: Path, groups=("X", "layers", "obsp", "obsm", "varm")) -> dict:
    """Value digest of every matrix in the given groups. Values are cast to float64 (data) or int64 (sparse indices)
    before hashing, so a copy written with other storage dtypes or chunking still matches when the values match."""
    digests = {}
    with h5py.File(path, "r") as handle:
        for group in groups:
            if group not in handle:
                continue
            node = handle[group]
            members = {group: node} if group == "X" else {f"{group}/{key}": node[key] for key in node}
            for name, member in members.items():
                digest = hashlib.sha256()
                if isinstance(member, h5py.Group):  # sparse matrix: data, indices, indptr
                    parts = [(member["data"], np.float64), (member["indices"], np.int64), (member["indptr"], np.int64)]
                else:
                    parts = [(member, np.float64)]
                for dataset, dtype in parts:
                    for start in range(0, dataset.shape[0], CHUNK):
                        digest.update(np.ascontiguousarray(dataset[start:start + CHUNK], dtype=dtype).tobytes())
                digests[name] = digest.hexdigest()
    return digests


def release_compression(path: Path):
    """Compression of the release's count matrix, so the deposit copy is stored the same way."""
    with h5py.File(path, "r") as handle:
        node = handle["layers/counts"] if "layers/counts" in handle else handle["X"]
        dataset = node["data"] if isinstance(node, h5py.Group) else node
        return dataset.compression


# Changes ----

def prune_obs(obs: pd.DataFrame, drop: list, rename: dict):
    removed = [column for column in drop if column in obs.columns]
    obs = obs.drop(columns=removed)
    renamed = {old: new for old, new in rename.items() if old in obs.columns}
    clashes = [new for new in renamed.values() if new in obs.columns]
    if clashes:
        raise ValueError(f"rename targets already in obs: {clashes}")
    return obs.rename(columns=renamed), removed, renamed


def replace_zone_scores(obs: pd.DataFrame, table_path: Path, modules: list, timepoint: str):
    """Swap the release's auc_* columns for the run of record's scores, at the position the old columns held."""
    old = [column for column in obs.columns if column.startswith("auc_")]
    position = obs.columns.get_loc(old[0]) if old else len(obs.columns)
    obs = obs.drop(columns=old)

    table = pd.read_csv(table_path, index_col="barcode")
    if not table.index.is_unique:
        raise ValueError(f"{table_path}: barcodes are not unique")
    on_timepoint = obs.index[obs["timepoint"].astype(str) == timepoint]
    # The run scored every nucleus of the timepoint in the all-cells object, so every object must find all of its own
    # timepoint nuclei in the table, and no scored nucleus of the object may belong to another timepoint.
    unscored = on_timepoint.difference(table.index)
    if len(unscored):
        raise ValueError(f"{len(unscored)} {timepoint} nuclei of the object have no score in {table_path}")
    off_timepoint = table.index.intersection(obs.index).difference(on_timepoint)
    if len(off_timepoint):
        raise ValueError(f"{len(off_timepoint)} scored nuclei are not {timepoint} nuclei in the object")

    added = []
    for module in modules:
        column = f"auc_{module}"
        values = pd.Series(np.nan, index=obs.index, dtype="float64")
        values.loc[on_timepoint] = table.loc[on_timepoint, column].to_numpy(dtype="float64")
        obs.insert(position + len(added), column, values)
        added.append(column)

    run_provenance = json.loads((table_path.parent / "provenance.json").read_text())
    uns_record = {
        "run": table_path.parent.name,
        "table": str(table_path),
        "table_sha256": file_sha256(table_path),
        "modules": modules,
        "timepoint": timepoint,
        "note": (f"AUCell scores of the zone-module run of record (Fig. S24) on the {timepoint} nuclei; NaN for every "
                 "other nucleus. Run parameters in run_provenance_json."),
        "run_provenance_json": json.dumps(run_provenance),
    }
    record = {"removed": old, "added": added, "n_scored": int(len(on_timepoint)), "table_rows": int(len(table)),
              "table_rows_not_in_object": int(len(table.index.difference(obs.index)))}
    return obs, record, uns_record


def prune_obsm(obsm, drop: list, rename: dict):
    removed = [key for key in drop if key in obsm]
    for key in removed:
        del obsm[key]
    renamed = {old: new for old, new in rename.items() if old in obsm}
    clashes = [new for new in renamed.values() if new in obsm]
    if clashes:
        raise ValueError(f"obsm rename targets already present: {clashes}")
    for old, new in renamed.items():
        obsm[new] = obsm.pop(old)
    return removed, renamed


def prune_uns(uns: dict, drop: list, rename: dict, removed_obsm: list, renamed_obsm: dict):
    removed = [key for key in drop if key in uns]
    for key in removed:
        del uns[key]
    renamed = {}
    for old, new in rename.items():
        if old in uns:
            if new in uns:
                raise ValueError(f"uns rename target already present: {new}")
            uns[new] = uns.pop(old)
            renamed[old] = new
    # A description of a removed embedding goes with it (e.g. embedding_provenance['X_diffmap']); the description of a
    # renamed embedding follows its new name.
    sub_removed, sub_renamed = [], []
    for key, value in uns.items():
        if not isinstance(value, dict):
            continue
        for name in removed_obsm:
            if name in value:
                del value[name]
                sub_removed.append(f"{key}/{name}")
        for old, new in renamed_obsm.items():
            if old in value:
                value[new] = value.pop(old)
                sub_renamed.append(f"{key}/{old} -> {key}/{new}")
    return removed, renamed, sub_removed, sub_renamed


def replace_pseudotime(adata: ad.AnnData, entry: dict):
    """Swap the release's pseudotime, its diffusion map and its uns record for those of a newer run, the pseudotime of
    record. The run's table must hold exactly the object's nuclei and match its configured sha256."""
    table_path, run_path, script_path = Path(entry["table"]), Path(entry["run_record"]), Path(entry["script"])
    table_sha256 = file_sha256(table_path)
    if table_sha256 != entry["table_sha256"]:
        raise ValueError(f"{table_path}: sha256 {table_sha256} is not the configured {entry['table_sha256']}")
    table = pd.read_parquet(table_path)
    if not table.index.is_unique:
        raise ValueError(f"{table_path}: barcodes are not unique")
    missing, extra = adata.obs_names.difference(table.index), table.index.difference(adata.obs_names)
    if len(missing) or len(extra):
        raise ValueError(f"{table_path}: {len(missing)} nuclei of the object are not in the table and {len(extra)} "
                         "table rows are not in the object")
    table = table.loc[adata.obs_names]
    run = json.loads(run_path.read_text())

    # Same layout as the release: the rank pseudotime in obs, one column per diffusion component in obsm, float64.
    n_dcs = adata.obsm["X_diffmap_pseudotime"].shape[1]
    release_pseudotime = adata.obs["pseudotime"].copy()
    release_record = dict(adata.uns["pseudotime"])
    adata.obs["pseudotime"] = table["pt_rank"].to_numpy(dtype="float64")
    adata.obsm["X_diffmap_pseudotime"] = table[[f"dc{i}" for i in range(n_dcs)]].to_numpy(dtype="float64")
    adata.uns["pseudotime"] = {
        "description": f"{release_record['description']}. {entry['note']}",
        "n_cells": int(len(table)),
        "n_dcs": int(n_dcs),
        "provenance_json": json.dumps(run),
        "root": run["root"],
        "run": release_record["run"],
        "script": script_path.name,
        "script_sha256": file_sha256(script_path),
        "table": str(table_path),
        "table_sha256": table_sha256,
        **{key: entry[key] for key in ("latent_copy", "latent_sha256") if key in entry},
    }

    record = {
        "table": str(table_path),
        "table_sha256": table_sha256,
        "run_record": str(run_path),
        "run_record_sha256": file_sha256(run_path),
        "script": str(script_path),
        "script_sha256": adata.uns["pseudotime"]["script_sha256"],
        "root": run["root"],
        "release_root": str(release_record.get("root")),
        "spearman_with_release": float(release_pseudotime.corr(adata.obs["pseudotime"], method="spearman")),
        "replaced": {"obs": ["pseudotime"], "obsm": ["X_diffmap_pseudotime"], "uns": ["pseudotime"]},
    }
    expected = {"obs": adata.obs["pseudotime"].copy(), "obsm": adata.obsm["X_diffmap_pseudotime"].copy(),
                "root": run["root"]}
    return record, expected


# Checks ----

def check_written(output: Path, release_obs: pd.DataFrame, expected_obs: pd.DataFrame, renamed: dict,
                  expected_obsm: list, expected_uns: list, release_var: pd.DataFrame, auc: dict,
                  timepoint: str, replaced: dict = None) -> dict:
    written = ad.read_h5ad(output, backed="r")
    checks = {
        "same_nuclei_same_order": bool(written.obs_names.equals(release_obs.index)),
        "obs_columns_as_expected": list(written.obs.columns) == list(expected_obs.columns),
        "var_identical": bool(written.var.equals(release_var)),
        "obsm_keys_as_expected": sorted(written.obsm.keys()) == sorted(expected_obsm),
        "uns_keys_as_expected": sorted(written.uns.keys()) == sorted(expected_uns),
    }
    unchanged = [c for c in expected_obs.columns if c in release_obs.columns and not c.startswith("auc_")]
    checks["kept_obs_columns_identical"] = all(written.obs[c].equals(release_obs[c]) for c in unchanged)
    # A replaced pseudotime is compared with its run's table, not with the release column it was renamed from.
    checks["renamed_columns_identical"] = all(written.obs[new].equals(release_obs[old].rename(new))
                                              for old, new in renamed.items()
                                              if not (replaced and new == "pseudotime"))
    if replaced:
        pseudotime = written.obs["pseudotime"]
        checks["pseudotime_equals_run_table"] = bool(np.array_equal(pseudotime.to_numpy(), replaced["obs"].to_numpy()))
        checks["run_root_at_lowest_rank"] = bool(pseudotime.loc[replaced["root"]] == pseudotime.min())
    if auc:
        on_timepoint = (written.obs["timepoint"].astype(str) == timepoint).to_numpy()
        present = written.obs[auc["added"]].notna().to_numpy()
        checks["scores_exactly_on_timepoint_nuclei"] = bool((present == on_timepoint[:, None]).all())
    written.file.close()
    return checks


# Run ----

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--release", required=True)
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--zone-aucell", default="", help="zone_aucell.csv of the run of record; '' for human objects")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    settings = yaml.safe_load(args.settings.read_text())["deposit_objects"]

    compression = release_compression(args.input)
    adata = ad.read_h5ad(args.input)
    release_obs, release_var = adata.obs.copy(), adata.var.copy()
    obs_before, uns_before, obsm_before = list(adata.obs.columns), sorted(adata.uns.keys()), sorted(adata.obsm.keys())
    LOG.info("read %s: %d nuclei, %d obs columns, compression %s",
             args.input, adata.n_obs, len(obs_before), compression)

    obs, obs_removed, obs_renamed = prune_obs(adata.obs, settings["obs_drop"], settings["obs_rename"])
    auc, uns_zone = {}, None
    if args.zone_aucell:
        obs, auc, uns_zone = replace_zone_scores(obs, Path(args.zone_aucell), settings["auc"]["modules"],
                                                 settings["auc"]["timepoint"])
    adata.obs = obs

    obsm_removed, obsm_renamed = prune_obsm(adata.obsm, settings["obsm_drop"], settings.get("obsm_rename", {}))
    uns_removed, uns_renamed, uns_sub_removed, uns_sub_renamed = prune_uns(
        adata.uns, settings["uns_drop"], settings["uns_rename"], obsm_removed, obsm_renamed)
    if uns_zone is not None:
        adata.uns["zone_aucell"] = uns_zone
    pt_entry = settings.get("pseudotime_of_record", {}).get(args.input.stem)
    pt_record, pt_expected = replace_pseudotime(adata, pt_entry) if pt_entry else (None, None)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output, compression=compression)
    LOG.info("wrote %s", args.output)
    expected_obs, expected_obsm, expected_uns = adata.obs, list(adata.obsm.keys()), list(adata.uns.keys())
    del adata

    checks = check_written(args.output, release_obs, expected_obs, obs_renamed, expected_obsm, expected_uns,
                           release_var, auc, settings["auc"]["timepoint"], replaced=pt_expected)
    digests_in = array_digests(args.input)
    digests_out = array_digests(args.output)
    # The release's matrices under their deposit names: removed embeddings left out, renamed ones under the new name.
    removed_names = {f"obsm/{key}" for key in obsm_removed}
    new_names = {f"obsm/{old}": f"obsm/{new}" for old, new in obsm_renamed.items()}
    kept = {new_names.get(name, name): digest for name, digest in digests_in.items() if name not in removed_names}
    # A replaced diffusion map is compared with its run's components; every other matrix with the release.
    replaced_names = {"obsm/X_diffmap_pseudotime"} if pt_expected is not None else set()
    checks["matrices_identical_in_value"] = (
        {name: digest for name, digest in kept.items() if name not in replaced_names}
        == {name: digest for name, digest in digests_out.items() if name not in replaced_names})
    renamed_kept = {old: new for old, new in obsm_renamed.items() if f"obsm/{new}" not in replaced_names}
    if renamed_kept:
        checks["renamed_embeddings_identical"] = all(digests_in[f"obsm/{old}"] == digests_out.get(f"obsm/{new}")
                                                     for old, new in renamed_kept.items())
    if pt_expected is not None:
        expected_digest = hashlib.sha256(np.ascontiguousarray(pt_expected["obsm"], dtype=np.float64).tobytes())
        checks["diffmap_equals_run_table"] = digests_out.get("obsm/X_diffmap_pseudotime") == expected_digest.hexdigest()
    LOG.info("checks: %s", checks)

    record = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "script": str(Path(__file__).resolve()),
        "script_sha256": file_sha256(Path(__file__)),
        "settings": str(args.settings),
        "settings_sha256": file_sha256(args.settings),
        "anndata": importlib.metadata.version("anndata"),
        "release": args.release,
        "input": str(args.input),
        "input_sha256": file_sha256(args.input),
        "output": str(args.output),
        "output_sha256": file_sha256(args.output),
        "compression": compression,
        "n_obs": len(release_obs),
        "obs": {"before": obs_before, "after": list(expected_obs.columns), "removed": obs_removed,
                "renamed": obs_renamed, "zone_scores": auc},
        "uns": {"before": uns_before, "after": sorted(expected_uns), "removed": uns_removed, "renamed": uns_renamed,
                "sub_entries_removed": uns_sub_removed, "sub_entries_renamed": uns_sub_renamed,
                "added": ["zone_aucell"] if uns_zone is not None else []},
        "obsm": {"before": obsm_before, "after": sorted(expected_obsm), "removed": obsm_removed,
                 "renamed": obsm_renamed},
        "pseudotime_of_record": pt_record,
        "matrix_digests": digests_out,
        "checks": checks,
    }
    args.provenance.write_text(json.dumps(record, indent=1) + "\n")
    if not all(checks.values()):
        raise SystemExit(f"checks failed: {[name for name, ok in checks.items() if not ok]}")


if __name__ == "__main__":
    main()
