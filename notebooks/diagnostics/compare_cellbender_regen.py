#!/usr/bin/env python
"""Compare regenerated CellBender outputs vs originals. Robust numeric handling."""
import os, math, numbers, warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import pandas as pd

SAMPLES_TSV = "/n/eddy_lab/users/nhilgert/xmens_final/config/samples.tsv"
ORIG_DIR    = "/n/eddy_lab/users/nhilgert/decidualization/xmens_inputs"
REGEN_DIR   = "/n/eddy_lab/users/nhilgert/decidualization/xmens_inputs_regen"

df = pd.read_csv(SAMPLES_TSV, sep="\t")
scol = "sample_id" if "sample_id" in df.columns else df.columns[0]
df = df.rename(columns={scol: "sample"})

def is_num(v):
    return isinstance(v, numbers.Number) and not (isinstance(v, float) and math.isnan(v))

def count_h5(path):
    if not path or not os.path.exists(path):
        return None
    try:
        return int(sc.read_10x_h5(path).n_obs)
    except Exception as e:
        return f"ERR:{type(e).__name__}"

def find_regen(sample, lane):
    for c in [f"{REGEN_DIR}/{sample}/{sample}_cellbender_filtered.h5",
              f"{REGEN_DIR}/{lane}/{lane}_cellbender_filtered.h5"]:
        if os.path.exists(c):
            return c
    return None

rows = []
for _, r in df.iterrows():
    sample, lane, species = r["sample"], r["coords_id"], r.get("species", "?")
    rows.append(dict(sample=sample, species=species,
                     orig=count_h5(f"{ORIG_DIR}/{sample}.h5"),
                     regen=count_h5(find_regen(sample, lane))))
out = pd.DataFrame(rows)

print("\n=== PER-SAMPLE (regen vs collaborator original) ===")
for _, x in out.iterrows():
    o, rg = x["orig"], x["regen"]
    if is_num(o) and is_num(rg):
        o, rg = int(o), int(rg)
        print(f"  {x['sample']:<24} {x['species']:<6} orig={o:>6}  regen={rg:>6}  {100*(rg-o)/o:+5.1f}%")
    else:
        st = "MISSING" if rg is None else rg
        print(f"  {x['sample']:<24} {x['species']:<6} orig={str(o):>6}  regen={st}")

done = out[out["orig"].apply(is_num) & out["regen"].apply(is_num)].copy()
done["orig"]  = done["orig"].astype(int)
done["regen"] = done["regen"].astype(int)

print("\n=== PER-SPECIES TOTALS (finished samples only) ===")
if len(done):
    for sp, g in done.groupby("species"):
        o, rg = int(g["orig"].sum()), int(g["regen"].sum())
        print(f"  {sp:<6}  n={len(g):>2}  orig={o:>7}  regen={rg:>7}  {100*(rg-o)/o:+.1f}%")
    o, rg = int(done["orig"].sum()), int(done["regen"].sum())
    print(f"  {'ALL':<6}  n={len(done):>2}  orig={o:>7}  regen={rg:>7}  {100*(rg-o)/o:+.1f}%")
else:
    print("  (none finished)")

n_missing = out["regen"].isna().sum() + (out["regen"].astype(str) == "None").sum()
n_pending = sum(1 for v in out["regen"] if v is None)
print(f"\n{len(done)}/{len(out)} regenerated; {n_pending} still missing/pending.")
