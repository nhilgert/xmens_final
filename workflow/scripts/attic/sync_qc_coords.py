#!/usr/bin/env python
"""Sync every current-vintage object to the settled 04f spatial QC.

OVERWRITES: x_um_qc, y_um_qc, spatial_qc_reason, obsm['X_spatial'] (and obsm['spatial'] if present)
ADDS:       tissue_section  (and x_um_dbscan/y_um_dbscan only where absent)
NEVER TOUCHES: an existing x_um_dbscan/y_um_dbscan -- the raw layer is the input to the QC rule
               and must stay reproducible-from-config.
Backup .pre_5mad_coords.bak per file (refuses if one exists); atomic tmp+rename; read-back verify.
"""
import anndata as ad, scanpy as sc, numpy as np, pandas as pd, json, pathlib, os, shutil, warnings
warnings.filterwarnings("ignore")

X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
N = pathlib.Path("/n/netscratch/eddy_lab/Everyone/nhilgert")
COLS = ["x_um_dbscan", "y_um_dbscan", "x_um_qc", "y_um_qc", "spatial_qc_reason", "tissue_section"]
REASONS = ["kept", "far", "degenerate", "boxed", "unmapped"]
NOTE = ("04f run 2026-09-02: far = 5 plain MADs (n_mads 5.0, mad_scale 1.0), 10 exclusion boxes, "
        "degenerate (dup coords / <10um); tissue_section = {sample}__A with uninduced split by y > x - 500. "
        "Human = degenerate-only per user decision; its tissue_section uses the same single-section default. "
        "obsm['X_spatial'] now holds the QC'd coordinates (it previously held raw). "
        "x_um_dbscan/y_um_dbscan left untouched.")
def log(*a): print(*a, flush=True)

# ---------------- barcode-keyed QC map, both species ----------------
qm = sc.read_h5ad(X / "data/interim/qc_spatial_mouse.h5ad", backed="r").obs
assert {"x_um_qc", "y_um_qc", "spatial_qc_reason", "tissue_section"} <= set(qm.columns)
M = qm[["x_um_dbscan", "y_um_dbscan", "x_um_qc", "y_um_qc"]].apply(pd.to_numeric)
M["spatial_qc_reason"] = qm.spatial_qc_reason.astype(str)
M["tissue_section"] = qm.tissue_section.astype(object)

qh = sc.read_h5ad(X / "data/interim/qc_spatial_human.h5ad", backed="r").obs
assert {"x_um_qc", "y_um_qc", "spatial_qc_reason", "sample"} <= set(qh.columns), sorted(qh.columns)[:40]
H = qh[["x_um_dbscan", "y_um_dbscan", "x_um_qc", "y_um_qc"]].apply(pd.to_numeric)
H["spatial_qc_reason"] = qh.spatial_qc_reason.astype(str)
# human has no 04f sections entry (mouse-only species list); apply 04f's own single-section default
H["tissue_section"] = np.where(H.spatial_qc_reason.values == "kept",
                               qh["sample"].astype(str).values + "__A", None)

QC = pd.concat([M[COLS], H[COLS]])
assert QC.index.is_unique, "barcode collision between species"
log(f"QC map: {len(QC):,} barcodes (mouse {len(M):,} + human {len(H):,})")
log(f"  mouse kept {int((M.spatial_qc_reason=='kept').sum()):,} | sections {M.tissue_section.nunique()}")
log(f"  human kept {int((H.spatial_qc_reason=='kept').sum()):,} | sections {pd.Series(H.tissue_section).nunique()}")

TARGETS = [
    ("lab  mouse All   2026_09_02", X / "data/deliverables/2026_09_02_Cevrim_XMens_All_Cells.h5ad"),
    ("lab  mouse FibDec 2026_09_02", X / "data/deliverables/2026_09_02_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"),
    ("lab  mouse All   2026_08_28", X / "data/deliverables/2026_08_28_Cevrim_XMens_All_Cells.h5ad"),
    ("lab  mouse FibDec 2026_08_28", X / "data/deliverables/2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"),
    ("lab  human All   2026_08_28", X / "data/deliverables/2026_08_28_Cevrim_Human_All_Cells.h5ad"),
    ("lab  human FibDec 2026_08_28", X / "data/deliverables/2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"),
    ("lab  mouse isodepth 2026_09_02", X / "results/isodepth/mouse_fibdec_isodepth_2026_09_02.h5ad"),
    ("lab  mouse isodepth", X / "results/isodepth/mouse_fibdec_isodepth.h5ad"),
    ("lab  human isodepth", X / "results/isodepth/human_fibdec_isodepth.h5ad"),
    ("net  mouse All   2026_09_02", N / "2026_09_02_Cevrim_XMens_All_Cells.h5ad"),
    ("net  mouse FibDec 2026_09_02", N / "2026_09_02_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"),
    ("net  mouse All   2026_08_28", N / "2026_08_28_Cevrim_XMens_All_Cells.h5ad"),
    ("net  mouse FibDec 2026_08_28", N / "2026_08_28_Cevrim_XMens_Only_Fib_Dec_Cells.h5ad"),
    ("net  human All   2026_08_28", N / "2026_08_28_Cevrim_Human_All_Cells.h5ad"),
    ("net  human FibDec 2026_08_28", N / "2026_08_28_Cevrim_Human_Only_Fib_Dec_Cells.h5ad"),
    ("net  samap vanilla joint", N / "2026_08_31_samap_vanilla_fibdec_joint.h5ad"),
    ("net  samap raw all", N / "2026_09_02_samap_raw_all.h5ad"),
]
report = {"note": NOTE, "objects": {}, "skipped": {
    "20251008_* (4 objects, Aug 21 vintage)": "bare 16-mer barcodes with no -1-<sample> suffix: 0 of them "
    "appear in the 04f index, and cell counts differ (57,242 vs 57,033). Cannot be joined by barcode; "
    "would need sample-aware reconstruction.",
    "donttouch/, heca/, *.bak": "left alone"}}

for tag, p in TARGETS:
    if not p.exists():
        report["objects"][("lab:" if str(p).startswith(str(X)) else "net:") + p.name] = {"tag": tag, "status": "MISSING"}; log(f"!! MISSING {p}"); continue
    bak = p.with_name(p.name + ".pre_5mad_coords.bak")
    log(f"\n=== {tag}: {p.name}")
    # key MUST include the location: six filenames exist in BOTH lab_storage and netscratch, and a
    # filename-only key silently drops one copy from the report.
    key = ("lab:" if str(p).startswith(str(X)) else "net:") + p.name
    ent = {"tag": tag, "path": str(p), "backup": bak.name}
    # resume: an object already carrying the 2026-09-02 push is verified and reported, not rewritten
    if bak.exists():
        done = sc.read_h5ad(p, backed="r")
        ent["n_obs"] = int(done.n_obs)
        if str(done.uns.get("spatial_qc_push", {}).get("date", "")) == "2026-09-02":
            o = done.obs
            kept = (o.spatial_qc_reason.astype(str) == "kept").to_numpy()
            bb = np.isfinite(pd.to_numeric(o.x_um_dbscan, errors="coerce").to_numpy(float))
            qq = np.isfinite(pd.to_numeric(o.x_um_qc, errors="coerce").to_numpy(float))
            Msp = np.asarray(done.obsm["X_spatial"], dtype=float)
            ok = (o.tissue_section.notna().to_numpy().sum() == kept.sum()
                  and np.array_equal(Msp[:, 0], pd.to_numeric(o.x_um_qc, errors="coerce").to_numpy(float), equal_nan=True))
            sec = o.tissue_section.value_counts()
            ent.update({"status": "already synced (verified, not rewritten)", "joined_frac": 1.0,
                        "n_qc_coords_after": int(qq.sum()), "raw_finite": int(bb.sum()), "qc_finite": int(qq.sum()),
                        "qc_finite_but_raw_nan": int((qq & ~bb).sum()), "readback_ok": bool(ok),
                        "obsm_updated": ["X_spatial"], "size_gb": round(os.path.getsize(p) / 1e9, 2),
                        "reasons": {str(k): int(v) for k, v in o.spatial_qc_reason.value_counts().items() if v},
                        "sections": {str(k): int(v) for k, v in sec.items() if v}})
            if ent["qc_finite_but_raw_nan"]:
                ent["raw_layer_warning"] = ("x_um_dbscan has fewer placed cells than x_um_qc: this object's "
                                            "'raw' column is itself filtered and is NOT the true raw layer")
            report["objects"][key] = ent
            log(f"   already synced; verified readback_ok={ok} qc={ent['qc_finite']:,} sections={len(ent['sections'])}")
            del done
            continue
        del done
        raise AssertionError(f"backup exists but object lacks the 2026-09-02 push, refusing: {p}")
    shutil.copy2(p, bak)
    a = ad.read_h5ad(p)
    ent["n_obs"] = int(a.n_obs)

    hit = a.obs_names.isin(QC.index)
    ent["joined_frac"] = round(float(hit.mean()), 6)
    assert hit.all(), f"{p.name}: only {hit.mean():.4f} of barcodes in the QC map"
    J = QC.reindex(a.obs_names)

    # raw layer: never overwrite; add only if absent; report any divergence from the QC source
    for c, src in (("x_um_dbscan", "x_um_dbscan"), ("y_um_dbscan", "y_um_dbscan")):
        if c in a.obs.columns:
            have = pd.to_numeric(a.obs[c], errors="coerce").to_numpy(float)
            want = pd.to_numeric(J[src]).to_numpy(float)
            ent[f"{c}_preserved"] = True
            ent[f"{c}_differs_from_04f_source"] = int((~np.isclose(have, want, equal_nan=True)).sum())
            a.obs[c] = have                      # unchanged values, dtype normalised
        else:
            a.obs[c] = pd.to_numeric(J[src]).to_numpy(float)
            ent[f"{c}_added"] = True

    xq = pd.to_numeric(J.x_um_qc).to_numpy(float)
    yq = pd.to_numeric(J.y_um_qc).to_numpy(float)
    ent["n_qc_coords_before"] = int(pd.to_numeric(a.obs["x_um_qc"], errors="coerce").notna().sum()) if "x_um_qc" in a.obs.columns else None
    a.obs["x_um_qc"], a.obs["y_um_qc"] = xq, yq
    a.obs["spatial_qc_reason"] = pd.Categorical(J.spatial_qc_reason.values, categories=REASONS)
    a.obs["tissue_section"] = pd.Categorical(J.tissue_section.where(J.tissue_section.notna(), None).values)
    ent["n_qc_coords_after"] = int(np.isfinite(xq).sum())
    ent["reasons"] = {k: int(v) for k, v in pd.Series(J.spatial_qc_reason.values).value_counts().items()}
    sec = pd.Series(J.tissue_section.values).value_counts()
    ent["sections"] = {str(k): int(v) for k, v in sec.items()}

    ent["obsm_updated"] = []
    for k in ("X_spatial", "spatial"):
        if k in a.obsm:
            a.obsm[k] = np.c_[xq, yq]; ent["obsm_updated"].append(k)
    if not ent["obsm_updated"]:
        a.obsm["X_spatial"] = np.c_[xq, yq]; ent["obsm_updated"].append("X_spatial (added)")

    a.uns["spatial_qc_push"] = {"date": "2026-09-02", "note": NOTE,
                                "source_mouse": "data/interim/qc_spatial_mouse.h5ad",
                                "source_human": "data/interim/qc_spatial_human.h5ad",
                                "n_joined": str(int(hit.sum())), "n_with_qc_coords": str(ent["n_qc_coords_after"]),
                                "n_sections": str(int((sec > 0).sum())), "backup": bak.name}
    if "tissue_section_provenance" in a.uns:
        a.uns["tissue_section_provenance_superseded"] = {
            "note": "2026-09-01 DBSCAN {sample}__s{k} definition retired; replaced 2026-09-02 by the 04f "
                    "line-split definition recorded in spatial_qc_push"}
        del a.uns["tissue_section_provenance"]
        ent["retired_old_section_provenance"] = True

    tmp = p.with_suffix(".h5ad.tmp"); a.write_h5ad(tmp); os.replace(tmp, p)
    del a

    b = sc.read_h5ad(p, backed="r")
    o = b.obs
    kept = (o.spatial_qc_reason.astype(str) == "kept").to_numpy()
    assert int(pd.to_numeric(o.x_um_qc, errors="coerce").notna().sum()) == ent["n_qc_coords_after"]
    assert o.tissue_section.notna().to_numpy().sum() == kept.sum(), "sections must exist exactly on kept cells"
    Msp = np.asarray(b.obsm["X_spatial"], dtype=float)
    assert np.array_equal(Msp[:, 0], pd.to_numeric(o.x_um_qc, errors="coerce").to_numpy(float), equal_nan=True)
    # Raw coverage is REPORTED, not asserted. On the deliverables raw is a superset of QC'd, but the
    # SAMap joint's mouse-half x_um_dbscan turned out to hold already-QC'd values under a raw name,
    # so "raw >= qc" is not an invariant of every object. Recorded so the anomaly stays visible.
    bb = np.isfinite(pd.to_numeric(o.x_um_dbscan, errors="coerce").to_numpy(float))
    qq = np.isfinite(pd.to_numeric(o.x_um_qc, errors="coerce").to_numpy(float))
    ent["raw_finite"], ent["qc_finite"] = int(bb.sum()), int(qq.sum())
    ent["qc_finite_but_raw_nan"] = int((qq & ~bb).sum())
    if ent["qc_finite_but_raw_nan"]:
        ent["raw_layer_warning"] = ("x_um_dbscan has fewer placed cells than x_um_qc: this object's 'raw' "
                                    "column is itself filtered and is NOT the true raw layer")
    ent["readback_ok"] = True
    ent["size_gb"] = round(os.path.getsize(p) / 1e9, 2)
    report["objects"][key] = ent
    log("   " + json.dumps({k: ent[k] for k in ("joined_frac", "n_qc_coords_before", "n_qc_coords_after",
                                                "obsm_updated", "size_gb", "readback_ok") if k in ent}))
    del b

json.dump(report, open("coord_sync_report.json", "w"), indent=1, default=str)
log(f"\nDONE - {sum(1 for e in report['objects'].values() if e.get('readback_ok'))} objects synced")
