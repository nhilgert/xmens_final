#!/usr/bin/env python
"""Fig. 5E, step 1: Spearman rho of every floor-6 consensus regulon's activity with the pseudotime, both species, and the
quadrant-plot input table (rule programs_regulon_table, 35_fig_programs.smk).

Port of scenic_quadrant_rho_fibdec_latent.py (2026-10-30, sha256 e2ea80fe...), the script that made the table of record
(resources/figures/panels_programs/regulon_quadrants_input_fibdec_latent.csv): same computation, every path an argument.
  activity   the CLI AUCell matrix of the floor-6 consensus (auc_mtx_weighted.csv) per species
  pseudotime mouse: obs[--mouse-pt-col] of the fibroblast/decidual object (apoptotic decidual nuclei are NaN and are
             dropped); human: --human-pt-col of the pseudotime-of-record table joined on obs_names, which must hold exactly
             the barcodes of the human object and of the release's human fib/dec object (nothing is re-ranked)
  orthologs  Ensembl one-to-one human-mouse symbol pairs (rule mi_overlap)
Spearman rho = Pearson correlation of average ranks; activity cast to float32 before ranking. Nuclei are matched to the
activity matrix by obs_names; the script stops unless every nucleus with a pseudotime has an activity score.
Outputs: regulon_rho_pseudotime_{mouse,human}.tsv, regulon_quadrants_input.csv (human_tf, mouse_tf, rho_human, rho_mouse,
label), provenance.json (sha256 of every input, the sanity checks, the human pseudotime side table as human_pt_side_table).
"""
import argparse, datetime, hashlib, json, os
import h5py, numpy as np, pandas as pd
from scipy.stats import rankdata


def sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def strs(ds):
    try:
        return np.asarray(ds.asstr()[:], dtype=object)
    except Exception:
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in ds[:]], dtype=object)


def obs_col(obs, name):
    node = obs[name]
    if isinstance(node, h5py.Group):                       # categorical (anndata >= 0.8)
        cats, codes = strs(node["categories"]), node["codes"][:]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
    v = node[:]
    return v.astype(float) if v.dtype.kind in "fiu" else strs(node)


def read_obs(path, pt_col):
    """obs_names, the pseudotime column and, if present, the decidual subtype and day columns (sanity checks only)."""
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        idx = strs(obs[obs.attrs["_index"]])
        pt = obs_col(obs, pt_col).astype(float)
        cols = {}
        for c in obs.keys():
            if c == obs.attrs["_index"] or c == pt_col:
                continue
            try:
                v = obs_col(obs, c)
            except Exception:
                continue
            if v.dtype == object and np.isin(["decidual_mature", "decidual_boundary"], v).all():
                cols["subtype"] = (c, v)
            if c.lower() in ("day", "days", "day_num", "timepoint_day") and v.dtype.kind in "fiu":
                cols["day"] = (c, v)
    return pd.Index(idx), pt, cols


def spearman(A, pt):
    """Spearman rho of each column of A (nuclei x regulons) with pt."""
    pr = rankdata(pt); pr -= pr.mean()
    RA = rankdata(np.asarray(A, dtype=np.float32), axis=0); RA -= RA.mean(axis=0)
    return (pr @ RA) / (np.sqrt((RA * RA).sum(axis=0)) * np.sqrt((pr * pr).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--activity-mouse", required=True); ap.add_argument("--activity-human", required=True)
    ap.add_argument("--object-mouse", required=True); ap.add_argument("--object-human", required=True)
    ap.add_argument("--mouse-pt-col", default="pt_harmonized_rank")
    ap.add_argument("--human-pt-table", required=True, help="pseudotime of record, parquet indexed by obs_names")
    ap.add_argument("--human-pt-col", default="pt_rank"); ap.add_argument("--human-pt-sha256", required=True)
    ap.add_argument("--human-release-object", required=True, help="the release's human fib/dec object (barcode check)")
    ap.add_argument("--orthologs", required=True, help="ortholog_pairs_one2one.csv with columns human, mouse")
    ap.add_argument("--labelled", required=True, help="comma-separated human TF symbols flagged in the label column")
    ap.add_argument("--expected-n", default="", help="e.g. mouse=42290,human=47533 (recorded, not enforced)")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    labelled = [x for x in a.labelled.split(",") if x]
    expected = {k: int(v) for k, v in (kv.split("=") for kv in a.expected_n.split(",") if kv)}
    activity = {"mouse": a.activity_mouse, "human": a.activity_human}
    objects = {"mouse": a.object_mouse, "human": a.object_human}

    found = sha256(a.human_pt_table)
    if found != a.human_pt_sha256:
        raise SystemExit(f"{a.human_pt_table}: sha256 {found} != pinned {a.human_pt_sha256}")
    prov = {"written": datetime.datetime.now().isoformat(timespec="seconds"),
            "pseudotime_column": {"mouse": a.mouse_pt_col, "human": f"{a.human_pt_col} of the pseudotime-of-record table"},
            "inputs": {}, "species": {}}
    rho = {}
    for sp in ("mouse", "human"):
        idx, pt, cols = read_obs(objects[sp], a.mouse_pt_col)
        if sp == "human":
            side = pd.read_parquet(a.human_pt_table)
            with h5py.File(a.human_release_object, "r") as f:
                rel_idx = pd.Index(strs(f["obs"][f["obs"].attrs["_index"]]))
            for name, other in ((objects[sp], idx), (a.human_release_object, rel_idx)):
                if not (side.index.is_unique and len(side) == len(other) and side.index.isin(other).all()):
                    raise SystemExit(f"{a.human_pt_table} and {name} hold different barcodes")
            pt_old = pt
            pt = side.loc[idx, a.human_pt_col].to_numpy(dtype=float)
            prov["human_pt_side_table"] = {
                "path": a.human_pt_table, "sha256": found, "column": a.human_pt_col,
                "barcodes_match_object": True, "barcodes_match_release_object": a.human_release_object,
                "spearman_with_object_pt": float(pd.Series(pt).rank().corr(pd.Series(pt_old).rank())),
                "median_pt_by_label": {str(k): float(v) for k, v in side.groupby("label", observed=True)[a.human_pt_col].median().items()} if "label" in side else {},
                "median_pt_by_timepoint": {str(k): float(v) for k, v in side.groupby("timepoint", observed=True)[a.human_pt_col].median().items()} if "timepoint" in side else {},
            }
        ok = ~np.isnan(pt)
        chk = {"n_obs": int(len(idx)), "n_with_pseudotime": int(ok.sum()), "expected_n_with_pseudotime": expected.get(sp),
               "pt_min": float(np.nanmin(pt)), "pt_max": float(np.nanmax(pt))}
        if "subtype" in cols:
            c, v = cols["subtype"]
            chk["subtype_column"] = c
            chk["median_pt"] = {s: float(np.median(pt[ok & (v == s)])) for s in ("decidual_boundary", "decidual_proliferative", "decidual_mature") if (ok & (v == s)).any()}
            chk["n_nan_by_subtype"] = {str(s): int(n) for s, n in pd.Series(v[~ok]).value_counts().items()}
        if "day" in cols:
            c, v = cols["day"]
            chk["day_column"] = c
            chk["spearman_pt_day"] = float(pd.Series(pt[ok]).rank().corr(pd.Series(v[ok].astype(float)).rank()))
        auc = pd.read_csv(activity[sp], index_col=0)
        keep = idx[ok]
        miss = keep.difference(auc.index)
        chk["n_with_pseudotime_missing_activity"] = int(len(miss))
        chk["n_activity_nuclei_without_pseudotime"] = int(len(auc.index.difference(keep)))
        if len(miss):
            raise SystemExit(f"{sp}: {len(miss)} nuclei with pseudotime have no activity score, e.g. {list(miss[:3])}")
        ptk = pd.Series(pt[ok], index=keep)
        A = auc.loc[keep]
        r = spearman(A.to_numpy(), ptk.to_numpy())
        t = pd.DataFrame({"regulon": auc.columns, "rho_pseudotime": r})
        t.to_csv(os.path.join(a.out_dir, f"regulon_rho_pseudotime_{sp}.tsv"), sep="\t", index=False, float_format="%.6g")
        rho[sp] = dict(zip(t.regulon.str[:-3], t.rho_pseudotime))
        chk["n_nuclei_correlated"] = int(len(keep)); chk["n_regulons"] = int(auc.shape[1])
        prov["species"][sp] = chk
        prov["inputs"][f"activity_{sp}"] = {"path": activity[sp], "sha256": sha256(activity[sp])}
        prov["inputs"][f"object_{sp}"] = {"path": objects[sp], "sha256": sha256(objects[sp])}
        print(sp, json.dumps(chk), flush=True)

    O = pd.read_csv(a.orthologs)[["human", "mouse"]].drop_duplicates()
    if not (O.human.is_unique and O.mouse.is_unique):
        raise SystemExit(f"{a.orthologs}: not one-to-one")
    prov["inputs"]["orthologs"] = {"path": os.path.abspath(a.orthologs), "sha256": sha256(a.orthologs)}
    d = O[O.human.isin(rho["human"]) & O.mouse.isin(rho["mouse"])].rename(columns={"human": "human_tf", "mouse": "mouse_tf"})
    d = d.sort_values("human_tf").reset_index(drop=True)
    d["rho_human"] = d.human_tf.map(rho["human"])
    d["rho_mouse"] = d.mouse_tf.map(rho["mouse"])
    d["label"] = d.human_tf.isin(labelled).astype(int)
    missing = sorted(set(labelled) - set(d.human_tf))
    d.to_csv(os.path.join(a.out_dir, "regulon_quadrants_input.csv"), index=False, float_format="%.6g")
    prov["n_tf_pairs"] = int(len(d)); prov["labelled"] = sorted(d.human_tf[d.label == 1]); prov["labelled_not_paired"] = missing
    json.dump(prov, open(os.path.join(a.out_dir, "provenance.json"), "w"), indent=1)
    print(f"{len(d)} TF pairs with a regulon in both species; {int(d.label.sum())} labelled; not paired: {missing}", flush=True)


if __name__ == "__main__":
    main()
