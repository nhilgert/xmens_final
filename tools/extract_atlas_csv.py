#!/usr/bin/env python3
"""Recover the atlas choices CSV from annotation that has ALREADY been done.

Not a re-derivation: it reads the cluster -> lineage mapping straight out of
annotated_{species}.h5ad. EXCLUDE clusters are recovered by set-difference
against latent_{species}.h5ad, since the annotated object has them dropped.

    python extract_atlas_csv.py --species mouse --res leiden_scVI_2.5
"""
import argparse, sys
from pathlib import Path
import scanpy as sc, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--species", required=True)
ap.add_argument("--res", default="leiden_scVI_2.5")
ap.add_argument("--repo", default=".")
ap.add_argument("--label-col", default=None, help="default: lineage, else coarse_refined")
a = ap.parse_args()

R = Path(a.repo)
lat = sc.read_h5ad(R/f"data/interim/latent_{a.species}.h5ad", backed="r")
ann = sc.read_h5ad(R/f"data/interim/annotated_{a.species}.h5ad", backed="r")
if a.res not in lat.obs:
    sys.exit(f"{a.res} not in latent object; have "
             f"{[c for c in lat.obs if c.startswith('leiden')]}")

col = a.label_col or next((c for c in ("lineage", "coarse_refined", "cell_type")
                           if c in ann.obs), None)
if col is None: sys.exit(f"no label column found in annotated obs: {list(ann.obs)[:20]}")
print(f"{a.species}: latent {lat.n_obs:,} cells, annotated {ann.n_obs:,} "
      f"({lat.n_obs - ann.n_obs:,} dropped) | label column '{col}'")

clus = lat.obs[a.res].astype(str)
shared = ann.obs.index.intersection(lat.obs.index)
print(f"shared barcodes: {len(shared):,} / {ann.n_obs:,}")
if len(shared) < 0.9 * ann.n_obs:
    print("!! low overlap -- annotated object may have different cell names")

df = pd.DataFrame({"cluster": clus.loc[shared].values,
                   "label": ann.obs.loc[shared, col].astype(str).values})
vote = df.groupby("cluster")["label"].agg(lambda s: s.value_counts().idxmax())
purity = df.groupby("cluster")["label"].agg(lambda s: s.value_counts(normalize=True).iloc[0])

rows = []
for c in sorted(clus.unique(), key=lambda v: int(v)):
    n = int((clus == c).sum())
    if c in vote.index:
        kept = int((df.cluster == c).sum())
        rows.append({"cluster": c, "label": vote[c], "n": n, "kept": kept,
                     "purity": round(float(purity[c]), 3)})
    else:
        rows.append({"cluster": c, "label": "EXCLUDE", "n": n, "kept": 0, "purity": 1.0})
out = pd.DataFrame(rows)

low = out[(out.label != "EXCLUDE") & (out.purity < 0.9)]
if len(low):
    print(f"\n!! {len(low)} cluster(s) below 0.9 purity -- the original annotation "
          f"may not have been made at {a.res}:")
    print(low.to_string(index=False))

print(f"\n{len(out)} clusters | {int((out.label=='EXCLUDE').sum())} EXCLUDE "
      f"({int(out.loc[out.label=='EXCLUDE','n'].sum()):,} cells)")
print(out.label.value_counts().to_string())

p = R/f"resources/annotations/{a.species}_lineage.csv"
p.parent.mkdir(parents=True, exist_ok=True)
out[["cluster", "label"]].to_csv(p, index=False)
print(f"\nwrote {p}")
