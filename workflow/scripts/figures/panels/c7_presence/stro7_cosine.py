"""Contrastive cosine of each human fibroblast/decidual nucleus to the 11 stromal subtype centroids of Ulrich et al.
2024, for the C7-fibroblast figure.

Called by rule c7_stro7_cosine (workflow/rules/41_fig_c7_presence.smk). A port, as a rule, of the human branch (panel 1)
of paper_code/ulrich_mapping/WRITEUP_2026_09_20/02_spatial_human/code/continuous_correspondence_job.py (4 Sep 2026),
which made the stro7_cosine column of the final human object; the statements below are that script's, in its order.
  atlas      Dataset S5 sheet 5B (subtype centroids, log1p CP10k), the 11 stromal subtypes (fibroblast 5, 7, 8, 9, 10,
             11; perivascular and smooth muscle 1, 2, 3, 4, 6)
  markers    Dataset S6 sheet 6C: each subtype's top 50 genes by avg_log2FC, pooled and held out of the panel
  panel      the 3,000 genes with the largest variance across the 11 centroids, after the markers are removed
  nuclei     the human half of the SAMap object's X (species hs; gene names hs_<symbol>), log1p applied when its
             maximum exceeds 50; centred by the mean over the human nuclei, each nucleus scaled to unit length
  centroids  centred by the mean over the 11 centroids, each centroid scaled to unit length
Check:   the Stro-7 column equals obs stro7_cosine of the final human object (--reference) within --tolerance.
Writes:  the 11 cosines per nucleus (obs_name, Stro-5 ... Stro-6) and a record JSON.

usage: stro7_cosine.py --joint H5AD --joint-sha256 S --s5 XLSX --s5-sha256 S --s6 XLSX --s6-sha256 S
                       --reference H5AD --reference-sha256 S [--n-panel 3000 --n-markers 50 --tolerance 1e-5]
                       --cosine CSV.GZ --record JSON
"""
import argparse
import logging
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

import c7_common as common

LOG = logging.getLogger("c7_stro7_cosine")
SPECIES_TAGS = ("hs", "human")


def num_of(s):
    head = str(s).split("_")[0].replace("Stro-", "").strip()
    try:
        return int(head)
    except ValueError:
        return None


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for name in ("joint", "s5", "s6", "reference"):
        p.add_argument(f"--{name}", required=True, type=Path)
        p.add_argument(f"--{name}-sha256", required=True)
    p.add_argument("--fibroblast", default="5,7,8,9,10,11")
    p.add_argument("--controls", default="1,2,3,4,6")
    p.add_argument("--target", default="Stro-7")
    p.add_argument("--n-panel", type=int, default=3000)
    p.add_argument("--n-markers", type=int, default=50)
    p.add_argument("--tolerance", type=float, default=1e-5)
    p.add_argument("--cosine", required=True, type=Path)
    p.add_argument("--record", required=True, type=Path)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    checks = common.Checks(LOG.info)
    fib_num = [int(x) for x in a.fibroblast.split(",")]
    ctrl_num = [int(x) for x in a.controls.split(",")]
    sha = {name: common.verify_sha256(getattr(a, name), getattr(a, f"{name}_sha256"), name)
           for name in ("joint", "s5", "s6", "reference")}

    # atlas centroids and published markers (as the source script)
    A = pd.ExcelFile(a.s5).parse("5B Subtype centroids", header=2, index_col=0)
    A.columns = pd.Index(A.columns).astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
    A.index = A.index.astype(str).str.strip()
    A = A[~A.index.isin(["nan", ""])].apply(pd.to_numeric, errors="coerce")
    A = A[~A.index.duplicated(keep="first")]
    stro = {num_of(c): c for c in A.columns if str(c).startswith("Stro-") and num_of(c) is not None}
    COLS = [stro[n] for n in fib_num + ctrl_num if n in stro]
    SHORT = [f"Stro-{num_of(c)}" for c in COLS]
    checks.equal("stromal subtype columns in S5B", len(COLS), len(fib_num) + len(ctrl_num))
    u6 = pd.ExcelFile(a.s6).parse("6C Stro markers", header=2)
    u6.columns = ["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"]
    u6 = u6.dropna(subset=["gene", "cluster"]); u6["num"] = u6["cluster"].map(num_of)
    markers = {n: list(dict.fromkeys(
                   u6[u6.num == n].sort_values("avg_log2FC", ascending=False)["gene"].astype(str)
               ))[:a.n_markers] for n in fib_num + ctrl_num}

    # the SAMap object: the human half of the species-prefixed gene axis
    with h5py.File(a.joint, "r") as handle:
        jobs = read_elem(handle["obs"]); jobs.index = jobs.index.astype(str)
        vn = pd.Index(read_elem(handle["var"]).index.astype(str))
        X = read_elem(handle["X"])
    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    heads = vn.str.split("_").str[0]
    keep = np.asarray(heads.isin(SPECIES_TAGS))
    sub = vn[keep]
    s2v = pd.Series(np.flatnonzero(keep), index=sub.str.split("_", n=1).str[1].values)
    s2v = s2v[~s2v.index.duplicated(keep="first")]
    rows = np.flatnonzero(jobs["species"].astype(str).isin(SPECIES_TAGS).to_numpy())
    Xh = X[rows]
    names = jobs.index[rows]
    LOG.info("SAMap object %s x %s; human nuclei %d; human genes %d", X.shape[0], X.shape[1], len(rows), len(s2v))
    del X
    x_max = float(Xh.max())
    log1p_applied = x_max > 50
    if log1p_applied:
        Xh = Xh.astype(np.float32) if not np.issubdtype(Xh.dtype, np.floating) else Xh
        np.log1p(Xh.data, out=Xh.data)
        LOG.info("applied log1p (max %.1f)", x_max)

    A_s = A.loc[:, COLS].copy()
    shared = pd.Index(s2v.index).intersection(A_s.index)
    marker_all = sorted({g for v in markers.values() for g in v})
    marker_in = [g for g in marker_all if g in set(shared)]
    var_rank = A_s.loc[shared, COLS].var(axis=1).sort_values(ascending=False)
    excl = set(marker_in)
    panel = [g for g in var_rank.index if g not in excl][:a.n_panel]
    checks.equal("panel size", len(panel), a.n_panel)
    checks.equal("atlas markers in the panel", len(set(panel) & excl), 0)
    checks.equal("C7 held out of the panel", "C7" not in set(panel), True)
    LOG.info("shared genes %d | markers held out %d | panel %d", len(shared), len(marker_in), len(panel))

    Xp = Xh[:, s2v.loc[panel].to_numpy()]
    Xp = Xp.toarray().astype(np.float32)
    Xc = Xp - Xp.mean(axis=0, keepdims=True)
    nr = np.linalg.norm(Xc, axis=1, keepdims=True); nr[nr == 0] = 1.0
    Xn = Xc / nr
    Ap = A_s.loc[panel, COLS].to_numpy(np.float32)
    Ac = np.nan_to_num(Ap - Ap.mean(axis=1, keepdims=True))
    An = Ac / np.linalg.norm(Ac, axis=0, keepdims=True)
    COS = pd.DataFrame(Xn @ An, index=pd.Index(names, name="obs_name"), columns=SHORT)

    ref = common.read_obs(a.reference, ["stro7_cosine"])["stro7_cosine"]
    checks.equal("reference nuclei found in the SAMap object", int(ref.index.isin(COS.index).sum()), len(ref))
    checks.equal("SAMap human nuclei found in the reference", int(COS.index.isin(ref.index).sum()), len(COS))
    diff = (COS[a.target].reindex(ref.index).to_numpy(float) - ref.to_numpy(float))
    checks.within(f"max |recomputed {a.target} - stored stro7_cosine|", float(np.nanmax(np.abs(diff))), a.tolerance)

    a.cosine.parent.mkdir(parents=True, exist_ok=True)
    COS.to_csv(a.cosine)
    inputs = {"joint": a.joint, "s5": a.s5, "s6": a.s6, "reference": a.reference}
    record = common.base_record(Path(__file__), inputs, {"cosine": a.cosine, "record": a.record}, sha)
    record["definition"] = {"subtypes": SHORT, "n_panel": a.n_panel, "n_markers_per_subtype": a.n_markers,
                            "panel_rank": "variance across the 11 centroids", "query_centring": "mean over human nuclei",
                            "atlas_centring": "mean over the 11 centroids", "x_max": x_max,
                            "log1p_applied": bool(log1p_applied)}
    record["counts"] = {"human_nuclei": int(len(COS)), "shared_genes": int(len(shared)),
                        "markers_held_out": int(len(marker_in)), "panel": panel}
    record["median_by_subtype"] = {k: float(v) for k, v in COS.median().items()}
    record["checks"] = checks.rows
    common.write_json(a.record, record)
    checks.fail_if_any("c7_stro7_cosine")


if __name__ == "__main__":
    main()
