#!/usr/bin/env python
"""EXACT code that produced hammoud_gene_drivers_v5.json.

Extracted verbatim (below the header) from the executed job script
jobs/414e5b83-3cb0-4522-a5fb-26cb34b28246/cmd.sh, run 2026-08-31 05:57, env rsc_uv.

Key facts:
- human clusters h0..h9 = 'h' + leiden_1.0 from results/human_fib_annotation/
  fine_cluster_assignments.csv (see hammoud_h_provenance.txt for how leiden_1.0 was made)
- gene panels: human = atlas genes (Ulrich sd05, fib-subtype max > 0.5) present in the
  human deliverable var_names (n=3,099); mouse = same atlas genes mapped through strict
  1:1 Ensembl orthologs (component_id pivot, ortholog_one2one) present in the mouse
  deliverable var_names (n=2,728). Driver gene symbols are reported in HUMAN space.
- centring: contrastive_cosine double-centres over the fib atlas columns only (fib_cols);
  Qz/Az are per-gene z-scores across those columns
- frac_from_top30 = sum of the top-30 per-gene contributions / total cosine. It is a
  RATIO, not a fraction: >1 or <0 when the cosine is near zero (see decidual_apoptotic
  1.519, human h3 -3.856) - do not report it for near-zero matches.
"""

import sys, json, pathlib, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scanpy as sc
HAM = pathlib.Path('/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud')
sys.path.insert(0, str(HAM/'fib_mapping'))
from hammoud_fib_mapping import (load_atlas_centroids, linear_panel_matrix,
                                group_centroids, contrastive_cosine, FIB_SUBTYPES)
X = pathlib.Path('/n/lab_storage/eddy_lab/users/nhilgert/xmens_final')
atlas = load_atlas_centroids(HAM/'pnas.2404775121.sd05.xlsx')
fib_cols = [c for c in atlas.columns if c in FIB_SUBTYPES] or [c for c in atlas.columns if c.startswith('Stro')]
atlas_f = atlas.loc[atlas.loc[:, fib_cols].max(1) > 0.5]
print(f"atlas after floor: {atlas_f.shape[0]:,} genes | fib cols {len(fib_cols)}", flush=True)

o = pd.read_csv('/n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro/orthogroups/ensembl_human_mouse_components.tsv', sep='\t')
o = o[o.homology_type == 'ortholog_one2one']
piv = o.pivot_table(index='component_id', columns='species', values='gene_symbol', aggfunc=lambda s: s.iloc[0] if len(s)==1 else np.nan)
piv = piv.dropna(subset=['human','mouse'])
h2m = piv.drop_duplicates('human').set_index('human')['mouse']
print(f"h2m strict 1:1 pairs: {len(h2m):,} | e.g. {list(h2m.items())[:3]}", flush=True)

res = {}
def drivers(tag, a, key, panel_h, panel_local):
    Xl = linear_panel_matrix(a, panel_local)
    C, ns = group_centroids(Xl, a.obs[key].astype(str), panel_local)
    C.index = panel_h
    A = atlas_f.loc[panel_h, fib_cols]
    M = contrastive_cosine(C, A)
    Qz = C.sub(C.mean(1), axis=0).div(C.std(1).replace(0, np.nan), axis=0).fillna(0)
    Az = A.sub(A.mean(1), axis=0).div(A.std(1).replace(0, np.nan), axis=0).fillna(0)
    print(f"[{tag}] panel {len(panel_h):,} | groups {list(C.columns)}", flush=True)
    for q in C.columns:
        row = M.loc[q] if q in M.index else M[q]
        best = row.idxmax()
        u = Qz[q].to_numpy(float); v = Az[best].to_numpy(float)
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        if nu == 0 or nv == 0: continue
        contrib = (u/nu)*(v/nv); cos = float(contrib.sum())
        df = pd.DataFrame(dict(gene=panel_h, contribution=contrib)).sort_values('contribution', ascending=False)
        df.to_csv(f'drivers_{tag}_{q}.csv', index=False)
        top = df.head(30)
        res[f'{tag}:{q}->{best}'] = dict(cosine=round(cos,4), n_cells=int(ns[q]), n_panel=int(len(panel_h)),
            frac_from_top30=round(float(top.contribution.sum()/cos),3) if cos else None,
            top15=[(r.gene, round(r.contribution,5)) for r in top.itertuples()][:15])
        print(f"[{tag}] {q} -> {best}: cos={cos:.3f} top30 explains {res[f'{tag}:{q}->{best}']['frac_from_top30']}", flush=True)
        print("     ", ", ".join(g for g,_ in res[f'{tag}:{q}->{best}']['top15'][:12]), flush=True)

am = sc.read_h5ad('/n/netscratch/eddy_lab/Everyone/nhilgert/2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad')
pairs = [(h, h2m[h]) for h in atlas_f.index if h in h2m.index and h2m[h] in am.var_names]
print(f"mouse ortholog panel: {len(pairs):,}", flush=True)
assert len(pairs) > 500, f"panel too small ({len(pairs)}) — ortholog join still wrong"
drivers('mouse', am, 'cell_subtype', pd.Index([h for h,_ in pairs]), pd.Index([m for _,m in pairs]))

ah = sc.read_h5ad('/n/netscratch/eddy_lab/Everyone/nhilgert/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad')
fine = pd.read_csv(X/'results/human_fib_annotation/fine_cluster_assignments.csv', index_col=0)
ah.obs['fine'] = ('h' + fine['leiden_1.0'].reindex(ah.obs_names).astype('Int64').astype(str)).values
ah = ah[~ah.obs['fine'].isin(['h<NA>']) & ah.obs['fine'].notna()].copy()
ph = pd.Index([g for g in atlas_f.index if g in ah.var_names])
drivers('human', ah, 'fine', ph, ph)

json.dump(res, open('hammoud_gene_drivers_v5.json','w'), indent=1)
print("DONE", flush=True)