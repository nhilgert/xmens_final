#!/usr/bin/env python
"""varD pseudotime from two candidate roots, both at the DC1 MAXIMUM (fibroblast pole).

rootA  global DC1 max, no filter        -> the literal "cell with the maximum DC1"
rootB  uninduced n efib_middle, pct=2   -> the resting-stroma root

Same method as notebook section 12: DPT on the stored diffmap, Palantir on X_scVI.
Orientation gate reports REFUSED; it never flips a pseudotime. Writes per-cell CSV + report;
does NOT modify the varD object.
"""
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, h5py, json, pathlib, warnings
warnings.filterwarnings("ignore")
from anndata.io import read_elem
from scipy.stats import spearmanr
import palantir

X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
P = X / "results/pt_variants/varD_keep_uninduced_no_efibouter/varD_keep_uninduced_no_efibouter_embedded.h5ad"
OUT = X / "results/pt_variants/manual_roots"; OUT.mkdir(parents=True, exist_ok=True)
SUBORDER = ['efib_proliferative','efib_outer','efib_middle','efib_middle_induced',
            'efib_subluminal','decidual_boundary','decidual_proliferative','decidual_mature']
TPORDER = ['uninduced','1_dpi_early','1_dpi_late','2_dpi','3_dpi','4_dpi_bleeding']
RANK = {s: i for i, s in enumerate(SUBORDER)}; TPRANK = {t: i for i, t in enumerate(TPORDER)}
def log(*a): print(*a, flush=True)

with h5py.File(P, "r") as f:
    a = ad.AnnData(obs=read_elem(f["obs"]))
    for k in ("X_diffmap", "X_scVI", "X_umap"): a.obsm[k] = read_elem(f[f"obsm/{k}"])
    for k in ("connectivities", "distances"):
        if k in f["obsp"]: a.obsp[k] = read_elem(f[f"obsp/{k}"])
    for k in ("neighbors", "diffmap_evals"):
        if k in f["uns"]: a.uns[k] = read_elem(f[f"uns/{k}"])
o = a.obs
log(f"varD: {a.n_obs:,} cells | diffmap {a.obsm['X_diffmap'].shape} | neighbors={'neighbors' in a.uns} evals={'diffmap_evals' in a.uns}")

DC1 = a.obsm["X_diffmap"][:, 1]
def resolve(spec):
    m = np.ones(a.n_obs, bool)
    if "subtype" in spec:   m &= (o["cell_subtype"].astype(str) == spec["subtype"]).values
    if "timepoint" in spec: m &= (o["timepoint"].astype(str) == spec["timepoint"]).values
    assert m.any(), spec
    v = DC1[m]
    if spec.get("pct", 0) == 0:
        j = int(np.argmax(v)) if spec["side"] == "high" else int(np.argmin(v))
    else:
        q = np.percentile(v, 100 - spec["pct"] if spec["side"] == "high" else spec["pct"])
        j = int(np.argmin(np.abs(v - q)))
    idx = int(np.flatnonzero(m)[j])
    return o.index[idx], idx

SPECS = {
    "rootA_globalmax": dict(dc=1, side="high", pct=0),
    "rootB_uninduced_efibmiddle": dict(timepoint="uninduced", subtype="efib_middle", dc=1, side="high", pct=2),
}
def gates(pt):
    sub_r = o["cell_subtype"].astype(str).map(RANK).to_numpy(float)
    tp_r = o["timepoint"].astype(str).map(TPRANK).to_numpy(float)
    ok = np.isfinite(pt)
    r_sub = float(spearmanr(pt[ok], sub_r[ok]).statistic)
    m2 = ok & np.isfinite(tp_r); r_tp = float(spearmanr(pt[m2], tp_r[m2]).statistic)
    # decidual-compartment day correlation: the orientation check that is not confounded by composition
    dec = ok & o["cell_subtype"].astype(str).str.startswith("decidual").values & np.isfinite(tp_r)
    r_tp_dec = float(spearmanr(pt[dec], tp_r[dec]).statistic)
    return dict(rho_subtype_rank=round(r_sub, 4), rho_timepoint=round(r_tp, 4),
                rho_timepoint_decidual_only=round(r_tp_dec, 4),
                verdict="OK" if (r_sub > 0 and r_tp_dec > 0) else "REFUSED (orientation gate failed - not flipped)")

report = {}
res = pd.DataFrame(index=a.obs_names); res.index.name = "cell"
res["cell_subtype"] = o["cell_subtype"].astype(str); res["timepoint"] = o["timepoint"].astype(str)
res["sample"] = o["sample"].astype(str); res["DC1"] = DC1

# Palantir diffusion space is root-independent -> build once, reuse for both roots
df = pd.DataFrame(a.obsm["X_scVI"], index=a.obs_names)
dm = palantir.utils.run_diffusion_maps(df, n_components=10, knn=30)
ms = palantir.utils.determine_multiscale_space(dm)
log("palantir multiscale space ready")

for tag, spec in SPECS.items():
    bc, idx = resolve(spec)
    log(f"\n=== {tag}: {bc}  [{o['cell_subtype'].iloc[idx]}, {o['timepoint'].iloc[idx]}]  "
        f"DC1={DC1[idx]:.5g} (global pct {float((DC1 < DC1[idx]).mean()):.3f})")
    ent = {"spec": spec, "root_barcode": bc, "root_idx": idx,
           "root_subtype": str(o["cell_subtype"].iloc[idx]), "root_timepoint": str(o["timepoint"].iloc[idx]),
           "root_DC1": float(DC1[idx]), "root_DC1_global_pct": round(float((DC1 < DC1[idx]).mean()), 4)}
    a.uns["iroot"] = idx
    sc.tl.dpt(a, n_dcs=min(10, a.obsm["X_diffmap"].shape[1]))
    dpt = a.obs["dpt_pseudotime"].to_numpy(float)
    res[f"{tag}_dpt"] = dpt; ent["dpt"] = gates(dpt); log("  DPT      " + json.dumps(ent["dpt"]))
    pr = palantir.core.run_palantir(ms, bc, num_waypoints=1200, n_jobs=8)
    pal = pr.pseudotime.reindex(a.obs_names).to_numpy(float)
    res[f"{tag}_palantir"] = pal; ent["palantir"] = gates(pal); log("  Palantir " + json.dumps(ent["palantir"]))
    for c in (f"{tag}_dpt", f"{tag}_palantir"):
        res[c + "_rank"] = res[c].rank(pct=True)
    ent["spearman_dpt_vs_palantir"] = round(float(spearmanr(dpt, pal).statistic), 4)
    ent["subtype_median_palantir_rank"] = {s: round(float(v), 4) for s, v in
        res.groupby("cell_subtype")[f"{tag}_palantir_rank"].median().reindex(SUBORDER).dropna().items()}
    ent["timepoint_median_palantir_rank"] = {t: round(float(v), 4) for t, v in
        res.groupby("timepoint")[f"{tag}_palantir_rank"].median().reindex(TPORDER).dropna().items()}
    report[tag] = ent

report["agreement_between_roots"] = {
    "dpt": round(float(spearmanr(res.rootA_globalmax_dpt, res.rootB_uninduced_efibmiddle_dpt).statistic), 4),
    "palantir": round(float(spearmanr(res.rootA_globalmax_palantir, res.rootB_uninduced_efibmiddle_palantir).statistic), 4)}
res.to_csv(OUT / "varD_two_roots_pt.csv")
json.dump(report, open("varD_two_roots_report.json", "w"), indent=1, default=str)
import shutil; shutil.copy(OUT / "varD_two_roots_pt.csv", "varD_two_roots_pt.csv")
log("\nagreement between roots: " + json.dumps(report["agreement_between_roots"]))
log("DONE")
