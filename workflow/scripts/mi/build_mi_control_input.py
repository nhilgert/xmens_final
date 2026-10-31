"""Build the input object for one mouse MI specificity-control arm (rule mi_control_input).

Reads the mouse fibroblast/decidual object and writes an h5ad with raw counts in X and
layers['counts'] and obs pt_norm (pseudotime re-ranked within each sample) for one cell set:
  uninduced        every uninduced nucleus with pseudotime
  induced_matched  a seeded random draw of induced nuclei, size-matched to the uninduced set
Ranking within sample removes between-timepoint differences, so MI in these arms reflects
ordering inside a sample. Also writes build_meta.json next to the h5ad.
"""
import argparse
import json
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def select_cells(timepoint: np.ndarray, has_pt: np.ndarray, cell_set: str, seed: int) -> np.ndarray:
    """Boolean mask of the nuclei in one control arm."""
    is_uninduced = (timepoint == "uninduced") & has_pt
    if cell_set == "uninduced":
        return is_uninduced

    induced = np.flatnonzero((timepoint != "uninduced") & has_pt)
    drawn = np.random.default_rng(seed).choice(induced, size=int(is_uninduced.sum()), replace=False)
    mask = np.zeros(timepoint.shape[0], bool)
    mask[np.sort(drawn)] = True
    return mask


def rank_within_sample(pt: np.ndarray, sample: np.ndarray) -> np.ndarray:
    """Average rank / n of pseudotime, computed separately in each sample."""
    pt_norm = np.full(pt.shape[0], np.nan)
    for s in pd.unique(sample):
        in_sample = sample == s
        pt_norm[in_sample] = pd.Series(pt[in_sample]).rank(method="average").to_numpy() / in_sample.sum()
    return pt_norm


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--h5ad", required=True, help="mouse fibroblast/decidual object")
parser.add_argument("--pt-key", required=True, help="obs column with the global pseudotime")
parser.add_argument("--cell-set", required=True, choices=["uninduced", "induced_matched"])
parser.add_argument("--seed", type=int, required=True, help="seed for the induced_matched draw")
parser.add_argument("--out", required=True, type=pathlib.Path, help="output h5ad")
args = parser.parse_args()

adata = ad.read_h5ad(args.h5ad)
pt = adata.obs[args.pt_key].to_numpy(float)
timepoint = adata.obs["timepoint"].astype(str).to_numpy()
sample = adata.obs["sample"].astype(str).to_numpy()

keep = select_cells(timepoint, np.isfinite(pt), args.cell_set, args.seed)
adata = adata[keep].copy()
pt_norm = rank_within_sample(pt[keep], sample[keep])
print(f"{args.cell_set}: {adata.n_obs:,} nuclei from {len(pd.unique(sample[keep]))} samples", flush=True)

# Raw counts go in X as well as layers['counts'], so every reader of X sees counts.
obs = adata.obs[["timepoint", "sample", "cell_subtype"]].copy()
obs.insert(0, "pt_raw", pt[keep])
obs.insert(0, "pt_norm", pt_norm)
control = ad.AnnData(
    X=adata.layers["counts"].copy(),
    obs=obs,
    var=adata.var.copy(),
    layers={"counts": adata.layers["counts"].copy()},
)
args.out.parent.mkdir(parents=True, exist_ok=True)
control.write_h5ad(args.out)

meta = dict(
    source=args.h5ad,
    cell_set=args.cell_set,
    pt_key_source=args.pt_key,
    pt_key_used="pt_norm",
    normalization="rank within sample",
    seed=args.seed,
    n_cells=int(control.n_obs),
    timepoints=obs["timepoint"].value_counts().to_dict(),
    spearman_pt_norm_vs_pt_raw=round(float(spearmanr(pt_norm, pt[keep]).statistic), 4),
)
(args.out.parent / "build_meta.json").write_text(json.dumps(meta, indent=2))
print(f"wrote {args.out}", flush=True)
