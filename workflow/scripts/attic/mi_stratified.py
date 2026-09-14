import scanpy as sc, numpy as np, pandas as pd, sys, pathlib, json
from sklearn.feature_selection import mutual_info_regression
sp, path, out = sys.argv[1], sys.argv[2], sys.argv[3]
a = sc.read_h5ad(path)
if sp == 'mouse':
    a = a[a.obs['timepoint'].astype(str) != 'uninduced'].copy()
a = a[a.obs['isodepth_consensus'].notna()].copy()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.filter_genes(a, min_cells=int(0.01*a.n_obs))
coord = a.obs['isodepth_consensus'].to_numpy(float)
samp = (a.obs['sample'] if 'sample' in a.obs.columns else a.obs['timepoint']).astype(str)
X = a.X.toarray() if not isinstance(a.X, np.ndarray) else a.X
genes = np.array(a.var_names)
print(f"[{sp}] {X.shape} | {samp.nunique()} samples", flush=True)
rng = np.random.default_rng(0)
per_sample = {}
for s in sorted(samp.unique()):
    m = (samp == s).values
    if m.sum() < 300:
        print(f"  skip {s} (n={m.sum()})", flush=True); continue
    idx = np.where(m)[0]
    if len(idx) > 3000: idx = rng.choice(idx, 3000, replace=False)
    per_sample[s] = mutual_info_regression(X[idx], coord[idx], random_state=0)
    print(f"  {s}: n={len(idx)} done", flush=True)
W = np.array([ (samp == s).sum() for s in per_sample ], float)
M = np.vstack([per_sample[s] for s in per_sample])
strat = (M * (W/W.sum())[:, None]).sum(0)
pooled_idx = rng.choice(len(coord), min(3000, len(coord)), replace=False)
pooled = mutual_info_regression(X[pooled_idx], coord[pooled_idx], random_state=0)
k = int(0.1 * M.shape[1])
top_counts = np.zeros(M.shape[1], int)
for row in M:
    top_counts[np.argsort(-row)[:k]] += 1
df = pd.DataFrame(dict(gene=genes, mi_stratified=strat, mi_pooled=pooled,
                       n_samples_top10=top_counts, n_samples=len(per_sample),
                       sd_across_samples=M.std(0)))
for s, v in per_sample.items():
    df[f'mi_{s}'] = v
df.sort_values('mi_stratified', ascending=False).to_csv(out, index=False)
print(f"[{sp}] wrote {out}: {df.shape}", flush=True)
