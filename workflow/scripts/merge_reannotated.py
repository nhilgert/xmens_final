#!/usr/bin/env python3
"""04d -- concatenate the cleaned per-lineage objects into one per-species object.

Each input already has EXCLUDE dropped and carries coarse_refined/fine_refined,
so a throwback lands in its final lineage automatically (the label IS the
placement). Drops the per-lineage embeddings -- separate scVI fits, different
spaces -- and keeps counts + spatial. 04e re-embeds from scratch.

Also writes the guardrail: kept vs published, agreement vs the deposited object,
throwback destinations, composition delta.
"""
from pathlib import Path
import scanpy as sc, anndata as ad, pandas as pd, numpy as np

sm = globals().get("snakemake")
parts_in = list(sm.input.parts); out_h5ad = sm.output.h5ad
report = sm.output.report; species = sm.params.species
ref_path = sm.params.ref; target_n = int(sm.params.target_n or 0)

KEEP_OBS = ["sample", "cell_barcode", "timepoint",
            "coarse_refined", "fine_refined", "origin_lineage"]

parts = []
for p in parts_in:
    a = sc.read_h5ad(p)
    a.obs["origin_lineage"] = Path(p).stem.split("_", 1)[1]
    for k in ("X_scVI", "X_umap_scVI"):
        if k in a.obsm: del a.obsm[k]
    a.obs = a.obs[[c for c in KEEP_OBS if c in a.obs.columns]].copy()
    print(f"  {Path(p).name:34} {a.n_obs:>7,} cells")
    parts.append(a)

merged = ad.concat(parts, axis=0, join="outer", merge="first",
                   uns_merge="first", index_unique=None)
assert "counts" in merged.layers, "lost layers['counts'] -- 04e needs raw counts"
print(f"\nmerged: {merged.n_obs:,} cells x {merged.n_vars:,} genes")

# ---- guardrail vs the deposited object --------------------------------------
def canon(s):
    s = s.astype(str).str.lower().str.replace(r"[^a-z0-9]", "", regex=True)
    return s.str.replace("bleeding", "menstrual", regex=False)
def key(o):
    return canon(o["sample"]) + "||" + \
           o["cell_barcode"].astype(str).str.replace(r"-\d+$", "", regex=True)

rows = {"kept_total": merged.n_obs, "published_n": target_n}
if target_n:
    rows["delta_pct"] = round(100 * (merged.n_obs - target_n) / target_n, 2)
if ref_path and Path(ref_path).exists():
    them = sc.read_h5ad(ref_path, backed="r").obs.copy(); them["k"] = key(them)
    r = key(merged.obs).map(them.set_index("k")["cell_type_annotations"])
    ref = r.astype(object).where(r.notna(), "unlabeled")
    lab = ref.ne("unlabeled").values
    cr = merged.obs.coarse_refined.astype(str)
    rows["labelled_by_them"] = int(lab.sum())
    rows["coarse_agreement_pct"] = round(
        100 * (cr[lab].values == ref.astype(str)[lab].values).mean(), 2) if lab.sum() else np.nan
    rows["kept_but_they_removed"] = int((~lab).sum())
else:
    print("  note: reference object not found -- skipping agreement")

pd.Series(rows).to_frame("value").to_csv(report, sep="\t")
print("\n=== guardrail ===")
print(pd.Series(rows).to_string())
print("\ncomposition:")
print(merged.obs.coarse_refined.value_counts().to_string())

Path(out_h5ad).parent.mkdir(parents=True, exist_ok=True)
merged.write_h5ad(out_h5ad)
print(f"\nwrote {out_h5ad} ({merged.n_obs:,} cells, no embedding -- 04e re-embeds)")
