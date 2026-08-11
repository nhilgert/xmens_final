#!/usr/bin/env python3
"""Compile a reproducible de-novo annotation notebook from a frozen choices CSV.

Workflow:  interactive 04c notebook  ->  choices CSV  ->  this compiler  ->  de-novo .ipynb

The CSV is the version-controlled source of truth (cluster -> label, hand-editable,
diffable). The emitted notebook is a generated artifact: it loads the object,
scores markers, shows the dotplot, applies the FROZEN map, verifies it against the
marker calls, and writes the annotated object -- a runnable record of the call.

CSV schema (header required):
    cluster,label,note
    0,epithelial cells,
    2,EXCLUDE,doublet
    ...
Optional 4th col `resolution` overrides --res per row (rarely needed).

Usage:
    python compile_denovo.py --species mouse --res leiden_scVI_2.5 \
        --csv results/tables/annot/mouse_lineage.csv \
        --out notebooks/curation/04_annotate_mouse_denovo.ipynb \
        [--level lineage|fine] [--lineage epithelial]   # fine == per-lineage object
"""
import argparse, csv, datetime, json, sys
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell as md, new_code_cell as code

ap = argparse.ArgumentParser()
ap.add_argument("--species", required=True)
ap.add_argument("--res", required=True, help="e.g. leiden_scVI_2.5")
ap.add_argument("--csv", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--level", choices=["lineage", "fine"], default="lineage")
ap.add_argument("--lineage", default=None, help="for --level fine: which lineage object")
ap.add_argument("--repo", default="/n/eddy_lab/users/nhilgert/xmens_final")
ap.add_argument("--spec", default="resources/lineage_spec.yaml",
                help="for --level lineage: source of the lineage->slug map")
a = ap.parse_args()

# ---- read the frozen choices -------------------------------------------------
rows = list(csv.DictReader(open(a.csv)))
if not rows: sys.exit(f"{a.csv} is empty")
need = {"cluster", "label"}
if not need <= set(rows[0]): sys.exit(f"CSV needs columns {need}; has {list(rows[0])}")
ANNOTATION = {r["cluster"].strip(): r["label"].strip() for r in rows}
NOTES = {r["cluster"].strip(): r.get("note", "").strip() for r in rows if r.get("note", "").strip()}
labels = sorted(set(ANNOTATION.values()) - {"EXCLUDE"})
SLUG, ROLLUP = {}, {}
try:
    import yaml
    _spec = yaml.safe_load(open(a.spec))
    if a.level == "fine" and a.lineage:
        _lin = _spec["species"][a.species]["lineages"][a.lineage]
        sys.path.insert(0, "workflow/scripts")
        from apply_annotation import build_rollup
        ROLLUP = build_rollup(_lin.get("fine", {}), _lin["coarse"])
except Exception as e:
    print(f"  note: no rollup from {a.spec} ({e})")
if a.level == "lineage":
    try:
        import yaml
        SLUG = yaml.safe_load(open(a.spec))["atlas"][a.species]["labels"]
    except Exception as e:
        print(f"  note: no slug map from {a.spec} ({e}); lineage subsets will be skipped")
n_excl = sum(v == "EXCLUDE" for v in ANNOTATION.values())
print(f"{len(ANNOTATION)} clusters | {len(labels)} labels | {n_excl} EXCLUDE | from {a.csv}")

# pretty literal for the emitted config cell
ann_lit = "{\n" + "\n".join(f'    "{k}": {json.dumps(v)},' for k, v in
                            sorted(ANNOTATION.items(), key=lambda kv: int(kv[0]))) + "\n}"
notes_lit = json.dumps(NOTES, indent=4) if NOTES else "{}"
obj = (f'latent_{a.species}_{a.lineage}.h5ad' if a.level == "fine"
       else f'latent_{a.species}.h5ad')
label_col = "fine" if a.level == "fine" else "lineage"
stamp = datetime.date.today().isoformat()

C = []
C.append(md(f"""# {'04c fine' if a.level=='fine' else '04'} — {a.species} {label_col} annotation (compiled, de novo)

**Generated** {stamp} by `compile_denovo.py` from `{Path(a.csv).name}` — do not
hand-edit the map here; edit the CSV and recompile. The map is frozen below; this
notebook loads the object, shows the marker dotplot, applies the map, verifies it
against marker-score calls, and writes the annotated object."""))

C.append(md("## 1. Config — the frozen map"))
C.append(code(f'''SPECIES = "{a.species}"
RES     = "{a.res}"
LEVEL   = "{a.level}"            # "lineage" (whole tissue) or "fine" (one lineage object)
LABEL_COL = "{label_col}"
REPO    = "{a.repo}"

# ---- frozen from {Path(a.csv).name} (compiled {stamp}) — edit the CSV, not this ----
ANNOTATION = {ann_lit}
NOTES = {notes_lit}
LABELS = {json.dumps(labels)}
SLUG   = {json.dumps(SLUG, indent=4)}\nROLLUP = {json.dumps(ROLLUP, indent=4)}'''))

C.append(code(f'''import warnings; warnings.filterwarnings("ignore")
import scanpy as sc, pandas as pd, numpy as np
from pathlib import Path
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 60)
REPO   = Path(REPO)
LATENT = REPO / "data/interim" / "{obj}"
PANEL  = REPO / "resources/markers" / f"panel_{{SPECIES}}.tsv"
OUT    = REPO / "data/interim" / ("annotated_{a.species}.h5ad" if LEVEL=="lineage"
                                  else "reannotated/{a.species}_{a.lineage or 'LIN'}.h5ad")
MAPCSV = Path("{a.csv}")
print(f"{{SPECIES}} {{LABEL_COL}} | {{LATENT.name}} | {{RES}} | {{len(ANNOTATION)}} clusters mapped")'''))

C.append(md("## 2. Load + check the map covers the clustering"))
C.append(code('''adata = sc.read_h5ad(LATENT)
assert RES in adata.obs, f"{RES} not in object"
adata.obs[RES] = adata.obs[RES].astype("category").cat.reorder_categories(
    sorted(adata.obs[RES].cat.categories, key=int))
clusters = set(adata.obs[RES].cat.categories)
missing = clusters - set(ANNOTATION); extra = set(ANNOTATION) - clusters
assert not missing, f"clusters with no CSV row: {sorted(missing, key=int)}"
if extra: print("NOTE CSV has rows for absent clusters (ignored):", sorted(extra, key=int))
bad = {c:v for c,v in ANNOTATION.items() if v not in LABELS+["EXCLUDE"]}
assert not bad, f"labels not in LABELS: {bad}"
print(f"{adata.n_obs:,} cells | {adata.obs[RES].nunique()} clusters, all mapped")'''))

C.append(md("## 3. Marker dotplot (evidence for the map)"))
C.append(code('''plot_sets = {}
if PANEL.exists():
    panel = pd.read_csv(PANEL, sep="\\t", comment="#")
    panel["genes"] = panel["genes"].str.split(",")
    plot_sets = {r.celltype:[g for g in r.genes if g in adata.var_names] for _,r in panel.iterrows()}
    plot_sets = {k:v for k,v in plot_sets.items() if v}
    sc.pl.dotplot(adata, plot_sets, groupby=RES, standard_scale="var",
                  figsize=(8, max(6,0.22*adata.obs[RES].nunique())), show=True)
else:
    print("no panel file -- skipping dotplot")'''))

C.append(md("""## 4. Apply the frozen map + verify against marker scores

The apply itself comes from `workflow/scripts/apply_annotation.py` — the *same*
function the Snakemake rules call, so this notebook and the pipeline cannot
drift. Everything else here is evidence (dotplots, marker-score cross-check)."""))
C.append(code('''# marker-score call per cluster (independent check that the manual map is sane)
marker_call = pd.Series("n/a", index=adata.obs[RES].cat.categories)
if plot_sets:
    for ct, genes in plot_sets.items():
        sc.tl.score_genes(adata, genes, score_name=f"__s_{ct}", ctrl_size=50, random_state=0)
    scols = [f"__s_{ct}" for ct in plot_sets]
    cs = adata.obs.groupby(RES, observed=True)[scols].mean()
    cs.columns = [c[4:] for c in cs.columns]
    marker_call = cs.idxmax(axis=1)
    for c in [c for c in adata.obs.columns if c.startswith("__s_")]: del adata.obs[c]

# --- the shared implementation (identical to what the pipeline runs) ---------
import sys; sys.path.insert(0, str(REPO/"workflow/scripts"))
from apply_annotation import apply_map, build_rollup

summary = apply_map(adata, RES, ANNOTATION, rollup=ROLLUP)
final = summary.rename(columns={"label": LABEL_COL, "excluded": "exclude"})
final["marker_call"] = marker_call.reindex(final.index)
final["marker_agrees"] = final.exclude | (final[LABEL_COL] == final.marker_call)
kept = final[~final.exclude]
print(f"{len(final)} clusters | {int(final.exclude.sum())} EXCLUDE ({int(final.loc[final.exclude,'n'].sum()):,} cells)")
print(f"retained {int(kept.n.sum()):,}")
print(f"marker agreement: {int(kept.marker_agrees.sum())}/{len(kept)} clusters\\n")
print(final.to_string())'''))

C.append(md("## 5. Dotplot ordered by assignment (blocks = shared identity)"))
C.append(code('''if plot_sets:
    rank = {l:i for i,l in enumerate(LABELS)}
    order = sorted(final.index, key=lambda c:(len(LABELS) if final.loc[c,"exclude"]
                   else rank.get(final.loc[c,LABEL_COL], len(LABELS)), int(c)))
    def _lbl(c):
        r = final.loc[c]; head = "EXCLUDE" if r.exclude else str(r[LABEL_COL])
        tail = "" if r.marker_agrees else f"  (markers: {r.marker_call})"
        return f"{head} | c{c}{tail}"
    lbl = {c:_lbl(c) for c in final.index}
    adata.obs["cluster_lbl"] = pd.Categorical(adata.obs[RES].map(lbl),
        categories=[lbl[c] for c in order], ordered=True)
    sc.pl.dotplot(adata, plot_sets, groupby="cluster_lbl", standard_scale="var",
                  figsize=(8, max(6,0.22*len(final))), show=True)
adata.obs[LABEL_COL] = adata.obs["coarse_refined"] if LEVEL == "lineage" else adata.obs["label"]
sc.pl.embedding(adata, "X_umap_scVI", color=LABEL_COL, show=True)'''))

C.append(md("## 6. Write the annotated object"))
C.append(code('''ann = adata[~adata.obs.excluded].copy()
if "cluster_lbl" in ann.obs: del ann.obs["cluster_lbl"]
OUT.parent.mkdir(parents=True, exist_ok=True)
ann.write(OUT)
print(f"wrote {OUT}: {ann.n_obs:,} cells")
print(ann.obs[LABEL_COL].value_counts().to_string())

# round-trip the frozen map (audit copy next to the object)
final.assign(source_resolution=RES, derived_from=LATENT.name,
             manual_note=pd.Series(NOTES).reindex(final.index).fillna("")) \\
     .to_csv(OUT.with_suffix(".map.csv"))
print(f"froze map -> {OUT.with_suffix('.map.csv')}")

if LEVEL == "lineage" and SLUG:
    # lineage subsets -- 04b_reintegrate consumes these
    LIN = REPO/"data/interim/lineages"; LIN.mkdir(parents=True, exist_ok=True)
    ann.obs[LABEL_COL] = pd.Categorical(ann.obs[LABEL_COL],
        categories=[l for l in LABELS if l in set(ann.obs[LABEL_COL])])
    rows = []
    for lin in ann.obs[LABEL_COL].cat.categories:
        sub = ann[ann.obs[LABEL_COL] == lin].copy()
        if sub.n_obs == 0: continue
        p = LIN/f"{SPECIES}_{SLUG[lin]}.h5ad"; sub.write(p)
        rows.append({"lineage": lin, "n": sub.n_obs, "file": p.name})
    print(pd.DataFrame(rows).set_index("lineage").to_string())

    # leiden ladder -- read by the R clustree script
    LAD = REPO/"results/tables/clusters"/SPECIES; LAD.mkdir(parents=True, exist_ok=True)
    base = pd.DataFrame({"cell": adata.obs_names,
                         "sample": adata.obs["sample"].astype(str).values,
                         "cell_barcode": adata.obs["cell_barcode"].astype(str).values})
    if "timepoint" in adata.obs: base["timepoint"] = adata.obs["timepoint"].astype(str).values
    res_cols = sorted([c for c in adata.obs if c.startswith("leiden_scVI_")],
                      key=lambda c: float(c.rsplit("_", 1)[1]))
    for col in res_cols:
        base.assign(cluster=adata.obs[col].astype(str).values) \\
            .to_csv(LAD/f"leiden_res{col.rsplit('_',1)[1]}.csv", index=False)
    wide = base.copy()
    for col in res_cols: wide[col] = adata.obs[col].astype(str).values
    wide[LABEL_COL] = adata.obs[LABEL_COL].astype(str).values
    wide.to_csv(LAD/"leiden_all_resolutions.csv", index=False)
    print(f"wrote {len(res_cols)} ladder files -> {LAD}")'''))

nb = new_notebook(cells=C)
nb.metadata["kernelspec"] = {"display_name":"Python 3","language":"python","name":"python3"}
nb.metadata["language_info"] = {"name":"python"}
nbf.validate(nb)
Path(a.out).parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, a.out)
print(f"compiled -> {a.out} ({len(nb.cells)} cells)")
