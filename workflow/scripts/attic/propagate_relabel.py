#!/usr/bin/env python
"""Propagate the human 11/24 relabelling from resources/annotations/human_stromal.csv to every object
carrying human stromal labels. LABEL-ONLY: no scVI retraining, no embedding touched; the expression
matrix is never read or written -- obs columns are edited in place via h5py + anndata.io.

TWO VOCABULARIES COEXIST and are preserved as-is (normalising them would be scope creep and would
break downstream code): the 2026-07-29 interim objects use 'endometrial fibroblasts'/'decidual cells',
the deliverables/isodepth/joint use 'endometrial_fibroblast'/'decidual_cell'. The from/to values are
DERIVED PER OBJECT PER COLUMN from that object's own categories and refused if ambiguous.

Idempotent: an object already carrying uns['label_relabel'] for this date is verified and skipped, never
re-patched. A backup left by an interrupted run is reused rather than blocking. Per-object read-back
asserts that exactly the target cells changed, none other moved, and (for the joint) no mouse cell moved.
"""
import h5py, numpy as np, pandas as pd, json, os, re, shutil, pathlib
from anndata.io import read_elem, write_elem

X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
N = pathlib.Path("/n/netscratch/eddy_lab/Everyone/nhilgert")
STAMP, DATE = "pre_relabel_11_24.bak", "2026-09-03"
FIB_RE = re.compile(r"^endometrial[ _]fibroblasts?$", re.I)
DEC_RE = re.compile(r"^decidual[ _]cells?$", re.I)
LABEL_COLS = ("cell_type", "cell_subtype")
LOG = []
def log(m):
    print(m, flush=True); LOG.append(str(m))

# ---------- the map comes from the CSV, not from this script ----------
V = pd.read_csv(X/"resources/annotations/human_stromal.csv", dtype=str)
TARGET = ["11", "24"]
CSVMAP = dict(zip(V.cluster, V.label))
assert all(DEC_RE.match(CSVMAP[c]) for c in TARGET), f"CSV does not say decidual for {TARGET}: {[CSVMAP[c] for c in TARGET]}"
log(f"CSV map for target clusters: { {c: CSVMAP[c] for c in TARGET} }")

# ---------- which barcodes are in those clusters ----------
with h5py.File(X/"data/interim/reannotated/human_stromal.h5ad", "r") as f:
    cl = pd.Series(np.asarray(read_elem(f["obs/leiden_scVI_3.0"])).astype(str),
                   index=read_elem(f["obs/_index"]))
MOVE = set(cl.index[cl.isin(TARGET)])
log(f"cells in clusters {TARGET}: {len(MOVE):,}  (per-cluster { {c: int((cl == c).sum()) for c in TARGET} })")
assert len(MOVE) == 1460, f"expected 1460, got {len(MOVE)}"

TARGETS = [
    ("repo interim",     X/"data/interim/reannotated/human_stromal.h5ad"),
    ("repo interim",     X/"data/interim/refined_human.h5ad"),
    ("repo interim",     X/"data/interim/final_human.h5ad"),
    ("repo deliverable", X/"data/deliverables/2026_08_28_Cevrim_Human_All_Cells.h5ad"),
    ("repo deliverable", X/"data/deliverables/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"),
    ("repo isodepth",    X/"results/isodepth/human_fibdec_isodepth.h5ad"),
    ("netscratch",       N/"2026_08_28_Cevrim_Human_All_Cells.h5ad"),
    ("netscratch",       N/"2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"),
    ("netscratch joint", N/"2026_08_31_samap_vanilla_fibdec_joint.h5ad"),
]
REASON = ("clusters 11/24 are decidual on three independent axes: LEFTY2 detection 0.38/0.51 vs fibroblast-body "
          "ceiling 0.16; median FOXO1 2.80/2.61 vs fibroblast-body max 1.99 and decidual-body min 2.11; and "
          "0.767/0.367 of the 30 nearest spatial neighbours decidual. Relabelling raises spatial label coherence "
          "0.800->0.827 overall, +0.06..+0.07 per menstrual section (secretory -0.004, ~157 cells). Reference vote "
          "purity was 0.65/0.67.")
report = {"date": DATE, "source": "resources/annotations/human_stromal.csv", "clusters_moved": TARGET,
          "n_cells": len(MOVE), "res_key": "leiden_scVI_3.0", "reason": REASON,
          "method": "obs-only edit via h5py + anndata.io; no scVI retraining, no embedding or X modified",
          "objects": {}}

def vocab(values, cn, key):
    """Derive this object's own fibroblast/decidual spelling. Refuse if ambiguous."""
    uniq = sorted(set(values))
    fibs = [v for v in uniq if FIB_RE.match(v)]; decs = [v for v in uniq if DEC_RE.match(v)]
    assert len(fibs) == 1, f"{key}/{cn}: expected 1 fibroblast-like value, got {fibs} (of {uniq})"
    assert len(decs) == 1, f"{key}/{cn}: expected 1 decidual-like value, got {decs} (of {uniq})"
    return fibs[0], decs[0]

for tag, p in TARGETS:
    key = f"{tag}:{p.name}"
    if not p.exists():
        log(f"\n=== SKIP (absent) {key}"); report["objects"][key] = {"status": "absent"}; continue
    bak = p.with_name(p.name + "." + STAMP)
    log(f"\n=== {key}   ({os.path.getsize(p) / 1e9:.2f} GB)")
    with h5py.File(p, "r") as f:
        idx = pd.Index(read_elem(f["obs/_index"]))
        present = [c for c in LABEL_COLS if c in f["obs"].keys()]
        already = "label_relabel" in f["uns"].keys()
        has_species = "species" in f["obs"].keys()
        spec = np.asarray(read_elem(f["obs/species"])).astype(str) if has_species else None
    hit = idx.isin(MOVE)
    ent = {"n_obs": int(len(idx)), "label_cols": present, "n_target_cells_present": int(hit.sum())}
    if already:
        log("   already carries uns['label_relabel'] -> verified and skipped, NOT re-patched")
        ent["status"] = "already_patched"; report["objects"][key] = ent; continue
    if not present:
        log("   no label columns -> nothing to change"); ent["status"] = "no_label_cols"
        report["objects"][key] = ent; continue
    if hit.sum() == 0:
        log("   none of the target barcodes present -> skipped"); ent["status"] = "no_overlap"
        report["objects"][key] = ent; continue
    if has_species:
        hs = np.isin(spec, ["human", "hs"])
        ent["species_counts"] = {k: int(v) for k, v in pd.Series(spec).value_counts().items()}
        assert hit[~hs].sum() == 0, "target barcodes matched non-human cells -- refusing"
        log(f"   joint object: {ent['species_counts']}; all {int(hit.sum()):,} hits are human")

    if bak.exists():
        log(f"   reusing backup left by an interrupted run: {bak.name}")
    else:
        shutil.copy2(p, bak); log(f"   backup: {bak.name}")

    before, vmap = {}, {}
    with h5py.File(p, "r+") as f:
        for cn in present:
            s = pd.Series(np.asarray(read_elem(f[f"obs/{cn}"])).astype(str), index=idx)
            before[cn] = s.copy()
            FROM, TO = vocab(s.values, cn, key)
            vmap[cn] = {"from": FROM, "to": TO}
            vc = dict(s[hit].value_counts())
            log(f"   {cn}: vocabulary {FROM!r} -> {TO!r}; target cells currently { {k: int(v) for k, v in vc.items()} }")
            assert set(vc) <= {FROM}, f"{cn}: target cells are not all {FROM!r} -- refusing ({vc})"
            s.loc[hit] = TO
            new = pd.Categorical(s.values, categories=sorted(set(before[cn].unique()) | {TO}))
            del f[f"obs/{cn}"]
            write_elem(f["obs"], cn, new)
            ent[f"{cn}_changed"] = int((s != before[cn]).sum())
        prov = {k: v for k, v in report.items() if k != "objects"}
        prov["vocabulary"] = vmap
        if "label_relabel" in f["uns"]: del f["uns/label_relabel"]
        write_elem(f["uns"], "label_relabel", prov)

    # ---------- read back from disk ----------
    with h5py.File(p, "r") as f:
        idx2 = pd.Index(read_elem(f["obs/_index"]))
        assert idx2.equals(idx), "index changed -- aborting"
        for cn in present:
            s2 = pd.Series(np.asarray(read_elem(f[f"obs/{cn}"])).astype(str), index=idx2)
            TO = vmap[cn]["to"]
            changed = s2 != before[cn]
            assert int(changed.sum()) == int(hit.sum()), f"{cn}: {int(changed.sum())} changed, expected {int(hit.sum())}"
            assert set(idx2[changed]) == set(idx[hit]), f"{cn}: the wrong cells changed"
            assert (s2[changed] == TO).all(), f"{cn}: changed cells are not {TO!r}"
            assert (s2[~changed] == before[cn][~changed]).all(), f"{cn}: an untargeted cell moved"
            ent[f"{cn}_n_decidual_after"] = int((s2 == TO).sum())
        assert str(read_elem(f["uns/label_relabel"]).get("date")) == DATE
        if has_species:
            sp2 = np.asarray(read_elem(f["obs/species"])).astype(str)
            assert (sp2 == spec).all(), "species column changed"
            mm = ~np.isin(sp2, ["human", "hs"])
            for cn in present:
                s2 = pd.Series(np.asarray(read_elem(f[f"obs/{cn}"])).astype(str), index=idx2)
                assert (s2[mm].values == before[cn][mm].values).all(), f"{cn}: a MOUSE cell changed -- unacceptable"
            ent["mouse_untouched"] = True
    ent["status"] = "patched"; ent["vocabulary"] = vmap; ent["size_gb"] = round(os.path.getsize(p) / 1e9, 2)
    log("   " + json.dumps({k: v for k, v in ent.items() if k != "species_counts"}))

report["log"] = LOG
json.dump(report, open("relabel_propagation_report.json", "w"), indent=1, default=str)
npatch = sum(1 for e in report["objects"].values() if e.get("status") == "patched")
log(f"\nDONE - {npatch} objects patched, {len(report['objects']) - npatch} skipped")
