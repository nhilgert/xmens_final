
import h5py, numpy as np, pandas as pd, json
from anndata.io import read_elem
X = "/n/lab_storage/eddy_lab/users/nhilgert/xmens_final"
with h5py.File(f"{X}/data/interim/reannotated/human_stromal.h5ad", "r") as f:
    cols = list(f["obs"].keys())
    qcc = [k for k in cols if k in ("n_genes_by_counts","total_counts","pct_counts_mt","n_genes","doublet_score","predicted_doublet")]
    d = {k: pd.to_numeric(read_elem(f["obs/"+k]), errors="coerce") for k in qcc}
    d["cluster"] = np.asarray(read_elem(f["obs/leiden_scVI_3.0"])).astype(str)
    Q = pd.DataFrame(d, index=read_elem(f["obs/_index"]))
print("QC columns found on the stromal object:", qcc)
if qcc:
    g = Q.groupby("cluster")
    T = pd.concat([g.size().rename("n")] + [g[k].median().rename("med_"+k) for k in qcc], axis=1)
    T["depth_pct"] = T["med_total_counts"].rank(pct=True).round(2) if "med_total_counts" in T else np.nan
    T["genes_pct"] = T["med_n_genes_by_counts"].rank(pct=True).round(2) if "med_n_genes_by_counts" in T else np.nan
    print("\n=== per-cluster depth (sorted by median total_counts) ===")
    print(T.sort_values("med_total_counts" if "med_total_counts" in T else T.columns[1]).round(1).to_string())
    T.to_csv(f"{X}/results/human_boundary/cluster_qc_depth.csv")
    for cid in ("33", "24", "11", "13"):
        if cid in T.index:
            r = T.loc[cid]
            print(f"\ncluster {cid}: n={int(r['n'])}  median total_counts={r.get('med_total_counts', float('nan')):.0f} "
                  f"(pct {r.get('depth_pct')})  median n_genes={r.get('med_n_genes_by_counts', float('nan')):.0f} (pct {r.get('genes_pct')})")
