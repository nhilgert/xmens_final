"""Harmonized time-informed diffusion pseudotime for the mouse and human fibroblast/decidual objects.

One procedure for both species:
  1. moscot TemporalProblem between consecutive days on the stromal-lineage scVI latent
     (epsilon 1e-3, tau_a = tau_b = 0.95, scale_cost "mean", growth priors from moscot's species gene sets).
  2. CellRank RealTimeKernel (self_transitions "all", conn_weight 0.2, threshold "auto").
     Same-day kNN graph (15 neighbors): --graph scvi -> on the stromal scVI;
                                        --graph pca  -> scanpy's default 50-PC PCA of .X on the stromal
                                                        lineage's stored (pooled) highly variable genes.
  3. Diffusion map: (T + T^T) / 2, rows renormalized, 15 leading eigenvectors (seed 0).
  4. Root: trajectory component = the DC (1..5) most correlated with the fibroblast -> decidual order;
     root = among earliest-day fibroblasts, the cell at the fibroblast end of that component.
  5. DPT with n_dcs = 15; rank = rank / n over the cells in the run.

    python harmonized_pseudotime.py --species mouse --graph scvi --cells final --out <dir>

Reads the input objects; writes cells.parquet and provenance.json to --out.

Vendored from the 22 Oct 2026 provenance copy (data/deliverables/2026_10_22_flagged_harmonized/provenance/ and
tmp/2026_10_31_netscratch_rescue/; sha256 d4b618ae2b037a833e927948c4095cce19f6028b8123dfda68e60bb16553ffa1).
Changes: lines wrapped at 120 characters, one statement per line, and the optional --object and --latent, which
replace the built-in object and latent paths (rule pt_record_pseudotime, 06c_pseudotime_of_record.smk). Without them
the script runs as on 22 Oct. The human pseudotime of record (30 Oct, harmonized_pseudotime_fibdec.py, sha256
5a6d6af7...) is the 22 Oct script with only the human latent path changed, i.e. this script with --latent.
Variable names are those of the run of record, so the two files can be compared line by line.
"""
import argparse
import datetime
import json
import os
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr

# As in the run of record: moscot, CellRank and scanpy warnings are left out of the log.
warnings.filterwarnings("ignore")

ap = argparse.ArgumentParser()
ap.add_argument("--species", choices=["mouse", "human"], required=True)
ap.add_argument("--graph", choices=["scvi", "pca"], required=True)
ap.add_argument("--cells", choices=["final", "released"], required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--object", default=None, help="input object; default: the species' built-in path for --cells")
ap.add_argument("--latent", default=None, help="stromal scVI latent; default: the species' built-in path")
args = ap.parse_args()

REPO = "/n/lab_storage/eddy_lab/users/nhilgert/xmens_final"
NS = "/n/netscratch/eddy_lab/Everyone/nhilgert"
SPECIES = {
    "mouse": dict(
        objects={
            "final": f"{NS}/2026_10_13_misplaced_removed/Cevrim_XMens_Only_Fib_Dec_Cells_misplaced_removed.h5ad",
            "released": f"{REPO}/data/deliverables/2026_09_02_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad",
        },
        latent=f"{REPO}/data/interim/latent_mouse_stromal.h5ad",
        gene_set="mouse",
        label_col="cell_subtype",
        drop={"cell_subtype": ["decidual_apoptotic"]},
        day={"uninduced": 0.0, "1_dpi_early": 1.0, "1_dpi_late": 1.5, "2_dpi": 2.0, "3_dpi": 3.0,
             "4_dpi_bleeding": 4.0},
        order=["efib_outer", "efib_middle", "efib_middle_induced", "efib_subluminal", "efib_proliferative",
               "decidual_boundary", "decidual_proliferative", "decidual_mature"],
        is_fibroblast=lambda lab: lab.str.startswith("efib")),
    "human": dict(
        objects={
            "final": (f"{NS}/2026_10_14_misplaced_removed_human_v2/"
                      "Cevrim_Human_Only_Fib_Dec_Cells_misplaced_removed.h5ad"),
            "released": f"{REPO}/data/deliverables/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad",
        },
        latent=f"{REPO}/data/interim/latent_human_stromal.h5ad",
        gene_set="human",
        label_col="cell_type",
        drop={},
        day={"secretory": 0.0, "menstrual": 1.0},
        order=["endometrial_fibroblast", "decidual_cell"],
        is_fibroblast=lambda lab: lab == "endometrial_fibroblast"),
}
P = dict(epsilon=1e-3, tau=0.95, scale_cost="mean", self_transitions="all", conn_weight=0.2, threshold="auto",
         n_neighbors=15, n_comps=15, n_dcs=15, seed=0)
cfg = SPECIES[args.species]
os.makedirs(args.out, exist_ok=True)
t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:7.0f}s] {msg}", flush=True)


# Cells, day, latent ----
src = args.object or cfg["objects"][args.cells]
latent_path = args.latent or cfg["latent"]
a = sc.read_h5ad(src)
for col in (cfg["label_col"], "timepoint", "sample"):
    a.obs[col] = a.obs[col].astype(str)
n_in = a.n_obs
for col, values in cfg["drop"].items():
    a = a[~a.obs[col].isin(values)].copy()
log(f"{args.species}/{args.cells}: {n_in} cells in object, {a.n_obs} after exclusions")
assert set(a.obs.timepoint) <= set(cfg["day"]), set(a.obs.timepoint) - set(cfg["day"])
a.obs["day"] = a.obs.timepoint.map(cfg["day"]).astype(float).astype("category")

lat = ad.read_h5ad(latent_path, backed="r")
lat_names = pd.Index(lat.obs_names.astype(str))
idx = lat_names.get_indexer(a.obs_names.astype(str))
assert (idx >= 0).all(), f"{int((idx < 0).sum())} cells missing from the stromal latent"
a.obsm["X_scVI"] = np.asarray(lat.obsm["X_scVI"])[idx]  # replaces whatever X_scVI the object carried
lat_hvg = pd.Series(np.asarray(lat.var["highly_variable"]).astype(bool), index=lat.var_names.astype(str))
lat.file.close()
hvg_genes = lat_hvg.index[lat_hvg.to_numpy()]
a.var["highly_variable"] = a.var_names.isin(hvg_genes)  # pooled stromal HVGs; used only by --graph pca
assert int(a.var.highly_variable.sum()) == len(hvg_genes) == 3000, int(a.var.highly_variable.sum())
log(f"stromal latent attached ({a.obsm['X_scVI'].shape[1]} dims); {len(hvg_genes)} stromal HVGs flagged")

# 1. moscot ----
from moscot.problems.time import TemporalProblem  # noqa: E402 (imported where the run of record imported it)
import jax  # noqa: E402

tp = TemporalProblem(a).score_genes_for_marginals(
    gene_set_proliferation=cfg["gene_set"], gene_set_apoptosis=cfg["gene_set"]
)
tp = tp.prepare(time_key="day", joint_attr="X_scVI")
tp = tp.solve(epsilon=P["epsilon"], tau_a=P["tau"], tau_b=P["tau"], scale_cost=P["scale_cost"])
converged = {str(k): bool(pr.solution.converged) for k, pr in tp.problems.items()}
log(f"moscot on {jax.default_backend()}: {len(converged)} coupling(s), all converged = {all(converged.values())}")

# 2. CellRank real-time kernel ----
import cellrank as cr  # noqa: E402

conn_kwargs = {"n_neighbors": P["n_neighbors"]}
if args.graph == "scvi":
    conn_kwargs["use_rep"] = "X_scVI"
rtk = cr.kernels.RealTimeKernel.from_moscot(tp)
n_dead = 0
try:
    rtk = rtk.compute_transition_matrix(self_transitions=P["self_transitions"], conn_weight=P["conn_weight"],
                                        threshold=P["threshold"], conn_kwargs=conn_kwargs)
except ValueError as e:
    # CellRank rejects a matrix whose rows do not sum to 1; same within-day uniform repair as the released notebooks.
    if "row stochastic" not in str(e):
        raise
    T = csr_matrix(rtk._transition_matrix).tocsr()
    rs = np.asarray(T.sum(1)).ravel()
    dead = np.where(~np.isfinite(rs) | (rs == 0))[0]
    n_dead = len(dead)
    dv = a.obs.day.astype(float).to_numpy()
    L = T.tolil()
    for i in dead:
        same = np.where(dv == dv[i])[0]
        L.rows[i] = list(same)
        L.data[i] = [1.0 / len(same)] * len(same)
    T = L.tocsr()
    T = csr_matrix(T.multiply(1.0 / np.asarray(T.sum(1)).ravel()[:, None]))
    rtk.transition_matrix = T
T = rtk.transition_matrix
dv = a.obs.day.astype(float).to_numpy()
Tc = T.tocoo()
mass = {
    k: float(Tc.data[np.sign(dv[Tc.col] - dv[Tc.row]) == s].sum() / Tc.data.sum())
    for k, s in (("later", 1), ("same", 0), ("earlier", -1))
}
log(f"transition matrix {T.shape}; mass later {mass['later']:.3f} / same {mass['same']:.3f} / "
    f"earlier {mass['earlier']:.3f}; repaired rows {n_dead}")

# 3. diffusion map ----
from scanpy.tools._dpt import DPT  # noqa: E402

S = (T + T.T) / 2
rs = S.sum(axis=1).A1
S.data = S.data / rs[S.nonzero()[0]]
# DPT() requires a graph; it is not used below.
sc.pp.neighbors(a, use_rep="X_scVI", n_neighbors=P["n_neighbors"], random_state=P["seed"])
dpt = DPT(adata=a, neighbors_key="neighbors")
dpt._transitions_sym = S
dpt.compute_eigen(n_comps=P["n_comps"], random_state=P["seed"])
a.obsm["X_diffmap"] = dpt.eigen_basis
a.uns["diffmap_evals"] = dpt.eigen_values
DM = np.asarray(a.obsm["X_diffmap"])

# 4. root ----
lab = a.obs[cfg["label_col"]]
order_rank = lab.map({s: i for i, s in enumerate(cfg["order"])}).to_numpy(float)
assert np.isfinite(order_rank).all(), sorted(set(lab[~np.isfinite(order_rank)]))
rhos = {i: float(spearmanr(DM[:, i], order_rank).statistic) for i in range(1, 6)}
best = max(rhos, key=lambda i: abs(rhos[i]))
fib_end = "max" if rhos[best] < 0 else "min"  # fibroblasts come first in the order, so they sit at the opposite sign
pool = (cfg["is_fibroblast"](lab) & (a.obs.day.astype(float) == a.obs.day.astype(float).min())).to_numpy()
v = np.where(pool, DM[:, best], np.nan)
root = int(np.nanargmax(v) if fib_end == "max" else np.nanargmin(v))
global_extreme = int(np.argmax(DM[:, best]) if fib_end == "max" else np.argmin(DM[:, best]))
log(f"trajectory component DC{best} (rho {rhos[best]:+.3f}); root {a.obs_names[root]} ({lab.iloc[root]}); "
    f"global extreme {'same cell' if root == global_extreme else a.obs_names[global_extreme]}")

# 5. pseudotime ----
a.uns["iroot"] = root
sc.tl.dpt(a, n_dcs=P["n_dcs"])
pt = a.obs["dpt_pseudotime"].to_numpy(float)
assert np.isfinite(pt).all() and pt[root] == 0
cells = pd.DataFrame({"pt_raw": pt, "pt_rank": pd.Series(pt).rank(pct=True).to_numpy(),
                      "label": lab.to_numpy(), "timepoint": a.obs.timepoint.to_numpy(),
                      "sample": a.obs["sample"].to_numpy()},
                     index=a.obs_names)
for i in range(DM.shape[1]):
    cells[f"dc{i}"] = DM[:, i]
cells.to_parquet(os.path.join(args.out, "cells.parquet"))

from importlib.metadata import version as pkg_version  # noqa: E402

st = os.stat(src)
prov = dict(species=args.species, graph=args.graph, cells=args.cells, n_cells=int(a.n_obs), params=P,
            conn_kwargs=conn_kwargs,
            input_object=dict(path=src, bytes=st.st_size,
                              mtime=datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")),
            latent=latent_path, hvg_source="stromal latent var['highly_variable'] (pooled, 3000)", day=cfg["day"],
            moscot_backend=jax.default_backend(), moscot_converged=converged, transition_mass=mass,
            repaired_rows=n_dead, evals=[float(x) for x in a.uns["diffmap_evals"]], trajectory_dc=int(best),
            dc_rho=rhos, fibroblast_end=fib_end, root=str(a.obs_names[root]), root_label=str(lab.iloc[root]),
            root_is_global_extreme=bool(root == global_extreme), root_pool_n=int(pool.sum()),
            versions={p: pkg_version(p) for p in ("scanpy", "anndata", "cellrank", "moscot", "ott-jax", "jax",
                                                  "numpy", "scipy")},
            runtime_s=round(time.time() - t_start), finished=datetime.datetime.now().isoformat(timespec="seconds"))
json.dump(prov, open(os.path.join(args.out, "provenance.json"), "w"), indent=1, default=str)
log("done")
