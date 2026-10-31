#!/usr/bin/env python
"""Second pass of facts for the methods-number registry (xmens / Cevrim et al.).

Usage: python collect_facts2.py OUT_DIR    (OUT_DIR must not exist; nothing is written elsewhere)
Part A reads small text/CSV/JSON files. Part B opens the joint SAMap object read-only and computes the
cross-species graph statistics quoted in the methods (edge totals, alignment scores, label transfer), with the cell
types of the release config.yaml names.
"""
import os, sys, re, json, glob, time, subprocess, traceback, hashlib, shutil
import numpy as np, pandas as pd

from release_objects import read_release

T0 = time.time()
OUT = sys.argv[1]; os.makedirs(OUT); os.makedirs(f"{OUT}/files")
L = "/n/lab_storage/eddy_lab/users/nhilgert"; X = f"{L}/xmens_final"; N = "/n/netscratch/eddy_lab/Everyone/nhilgert"
RELEASE = read_release(X)
# Supplementary tables of 28 Oct, built on the X-Mens fit with the uninduced fibroblasts (read-only; Drive mirror).
# The build's directory lost its date prefix (README: built 2026-10-28; Wang provenance built 2026-10-27T22:20:53 by
# supp_tables_2026_10_28); the CSVs are byte-identical to the copy in 2026_10_29_supplementary_tables_v7.
SUPP_CSV = f"{L}/xmens_figures/supplementary_tables/csv"
F = {"_meta": {"started": time.strftime("%Y-%m-%d %H:%M:%S"), "host": os.uname().nodename}}


def rec(name, fn):
    t = time.time()
    try:
        F[name] = fn()
    except Exception as e:
        F[name] = {"_error": f"{type(e).__name__}: {e}"[:400], "_trace": traceback.format_exc()[-800:]}
    print(f"[{time.time() - T0:6.0f}s] {name} ({time.time() - t:.1f}s){' ERROR' if isinstance(F[name], dict) and '_error' in F[name] else ''}", flush=True)


def sh(c, t=120):
    try:
        return subprocess.run(["bash", "-c", c], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=t).stdout
    except subprocess.TimeoutExpired:
        return f"TIMEOUT: {c[:120]}"


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def keep(p, name=None):
    """Copy a small source file into OUT/files and record its md5 and mtime."""
    name = name or os.path.basename(p)
    shutil.copy2(p, f"{OUT}/files/{name}")
    return {"path": p, "copied_as": name, "md5": md5(p), "mtime": time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p))), "bytes": os.path.getsize(p)}


# ================================================================== Part A: small files
def samap_code():
    out = {"rule": keep(f"{X}/workflow/rules/10b_samap_raw.smk")}
    txt = open(f"{X}/workflow/rules/10b_samap_raw.smk").read()
    scripts = sorted(set(re.findall(r"workflow/scripts/[\w/.\-]+\.py", txt)))
    out["scripts"] = {}
    for s in scripts:
        p = f"{X}/{s}"
        if os.path.exists(p): out["scripts"][s] = keep(p, s.replace("/", "__"))
    out["results_tree"] = sh(f"find {X}/results/samap_raw {X}/results/samap_transfer -maxdepth 3 2>/dev/null | head -80")
    out["label_files"] = sh(f"find {X}/results {N} -maxdepth 4 -iname '*label*transfer*' -o -maxdepth 4 -iname '*transfer*label*' 2>/dev/null | head -20").split()
    return out
rec("samap_code", samap_code)


def misc_small():
    out = {}
    cb = glob.glob(f"{L}/.conda/envs/cellbender_022/lib/python3*/site-packages/cellbender/remove_background/argparse.py")
    out["cellbender_argparse"] = cb[0] if cb else None
    out["cellbender_epochs"] = sh(f"grep -n -A2 -- '--epochs' {cb[0]} | head -6") if cb else "argparse.py not found"
    out["cellbender_rule"] = keep(f"{X}/workflow/rules/00_cellbender.smk")
    out["proteomes"] = sh(f"timeout 30 ls {L}/data/proteomes 2>&1 | head -30; for f in {L}/data/proteomes/README* {L}/data/proteomes/*.md {L}/data/proteomes/*.txt; do test -f \"$f\" && (echo \"## $f\"; timeout 10 head -40 \"$f\"); done", t=90)
    out["proteome_release_grep"] = sh(f"timeout 40 grep -rhoE 'release-[0-9]+|release_[0-9]+|Ensembl [0-9]+' {L}/data/proteomes/*.md {L}/data/proteomes/*.txt {L}/data/proteomes/*.sh {L}/data/proteomes/*.json {X}/workflow/rules/09_homology.smk {X}/resources/orthology/README.md 2>/dev/null | sort | uniq -c | head -12", t=90)
    out["mi_null_pct"] = sh(f"grep -rnE '97\\.5|2\\.5' {X}/workflow/scripts/mi {X}/workflow/scripts/overlap {X}/workflow/scripts 2>/dev/null --include=*.py | grep -iE 'percentile|quantile|null' | cut -c1-220 | head -12")
    out["overlap_scripts"] = sh(f"grep -rlE 'overlap_cutoff_sweep|overlap_permutation_null' {X}/workflow 2>/dev/null | head -8").split()
    for p in out["overlap_scripts"][:4]:
        out.setdefault("overlap_kept", {})[p] = keep(p, "overlap__" + os.path.basename(p))
    out["heatmap_code_find"] = sh(f"grep -rlE '27,?064|evenly spaced|spline|dispersion' {N}/figures/pt_gene_heatmaps {N}/figures/*.py {X}/workflow/scripts 2>/dev/null --include=*.py --include=*.R | head -10").split()
    hm = f"{N}/figures/pt_gene_heatmaps/fig_pt_gene_heatmaps.py"
    if os.path.exists(hm): out["heatmap_fig"] = keep(hm)
    out["heatmap_prov"] = {os.path.basename(p): json.load(open(p)) for p in sorted(glob.glob(f"{N}/figures/pt_gene_heatmaps/*.provenance.json"))}
    return out
rec("misc_small", misc_small)


def mouse_reference():
    import yaml
    cfg = yaml.safe_load(open(f"{X}/config/config.yaml"))
    g = cfg["genome_references"]["mouse"]
    out = {"config": g}
    cands = [v for v in g.values() if isinstance(v, str) and v.startswith("/")]
    for c in cands:
        for p in (c, f"{c}/reference.json", os.path.join(os.path.dirname(c), "reference.json")):
            if os.path.isfile(p) and p.endswith(".json"):
                out["reference_json"] = {"path": p, "content": json.load(open(p))}
                break
    gtf = [f"{c}/genes/genes.gtf" for c in cands if os.path.isfile(f"{c}/genes/genes.gtf")] + [p for c in cands for p in glob.glob(f"{os.path.dirname(c)}/*/genes/genes.gtf")]
    if not gtf:
        rj = [v for v in [out.get("reference_json", {}).get("path")] if v]
        gtf = [os.path.join(os.path.dirname(p), "genes", "genes.gtf") for p in rj if os.path.isfile(os.path.join(os.path.dirname(p), "genes", "genes.gtf"))]
    if gtf:
        with open(gtf[0]) as fh: out["gtf_header"] = {"path": gtf[0], "lines": [next(fh).rstrip("\n") for _ in range(5)]}
    if "reference_json" not in out:
        out["find"] = sh(f"find {' '.join(cands)} -maxdepth 2 -name 'reference.json' 2>/dev/null | head -3").split()
        for p in out["find"][:1]: out["reference_json"] = {"path": p, "content": json.load(open(p))}
    if "gtf_header" not in out and "reference_json" in out:
        g_ = os.path.join(os.path.dirname(out["reference_json"]["path"]), "genes", "genes.gtf")
        if os.path.isfile(g_):
            with open(g_) as fh: out["gtf_header"] = {"path": g_, "lines": [next(fh).rstrip("\n") for _ in range(5)]}
    return out
rec("mouse_reference", mouse_reference)


def scenic_mouse_manifest():
    R = f"{L}/decidualization/revisions/scenic_v2_2026_08_28/scenic_seeds/mouse"
    out = {}
    for p in [f"{R}/mouse_dsc_pooled.manifest.json"] + sorted(glob.glob(f"{R}/looms_sketch15k/*.manifest.json"))[:2]:
        if os.path.exists(p): out[p] = json.load(open(p))
    return out
rec("scenic_mouse_manifest", scenic_mouse_manifest)


def ivd_release():
    out = {}
    R = f"{N}/cagri_ivd_ulrich_deseq2_runs"
    out["listing"] = {d: sorted(os.listdir(f"{R}/{d}")) for d in sorted(os.listdir(R))}
    mt = {}
    for d in sorted(os.listdir(R)):
        p = f"{R}/{d}/xmens_broad_all_genes.csv"
        if os.path.exists(p): mt[d] = md5(p)
    out["xmens_table_md5"] = mt
    # The Wang in vitro support rows (gene mapping to X-Mens identifiers) are read from the last run directory;
    # their md5 in every directory shows whether that choice matters.
    out["wang_ivd_support_md5"] = {d: md5(f"{R}/{d}/wang_ivd_support_rows.csv") for d in sorted(os.listdir(R))
                                   if os.path.exists(f"{R}/{d}/wang_ivd_support_rows.csv")}
    rel = sorted(mt)[-1] if mt else None
    out["release_used"] = rel
    if rel:
        for f in sorted(os.listdir(f"{R}/{rel}")):
            p = f"{R}/{rel}/{f}"
            if os.path.isfile(p) and f.endswith((".csv", ".tsv", ".json", ".md", ".yaml", ".txt")) and os.path.getsize(p) < 40e6:
                out.setdefault("kept", {})[f] = keep(p, "ivdrel__" + f)
    out["corrected_find"] = sh(f"find {L}/cagri_ivd_end_to_end_20260815 {L}/human_ivd_repro {N} -maxdepth 4 -name '*stadtmauer_authors_criteria*' 2>/dev/null | head -12").split()
    for p in out["corrected_find"]:
        if os.path.isfile(p) and os.path.getsize(p) < 20e6: out.setdefault("corrected_kept", {})[p] = keep(p, "corr__" + os.path.basename(p))
    C = f"{L}/cagri_ivd_end_to_end_20260815"
    for f in ("resources/wang_pregnancy.csv", "resources/wang_oil.csv", "resources/wang_ivd.csv", "config.yaml", "README.md"):
        if os.path.exists(f"{C}/{f}"): out.setdefault("upstream", {})[f] = keep(f"{C}/{f}", "ivdup__" + f.replace("/", "__"))
    out["pb_floor_grep"] = sh(f"grep -nE 'min_cells|MIN_CELLS|floor|>= ?25|25' {C}/run_pipeline.py {C}/config.yaml 2>/dev/null | cut -c1-200 | head -20")
    return out
rec("ivd_release", ivd_release)


def supp_tables():
    """Wang comparison summary and the methods-text checks of the 28 Oct supplementary tables."""
    out = {}
    for name in ("SuppTable_Wang_comparison__Summary.csv", "SuppTable_IVD_methods_checks.csv"):
        out[name] = keep(f"{SUPP_CSV}/{name}", "supp__" + name)
    return out
rec("supp_tables", supp_tables)


# ================================================================== Part B: joint SAMap graph
def samap_graph():
    import anndata as ad, scipy.sparse as ss, h5py
    from anndata.io import read_elem
    JOINT = f"{N}/2026_09_02_samap_raw_all.h5ad"
    MOUSE = RELEASE["paths"]["mouse_fibdec"]
    HUMAN = RELEASE["paths"]["human_fibdec"]
    out = {"joint": JOINT, "joint_md5_head": md5(JOINT) if os.path.getsize(JOINT) < 2e9 else "skipped (large)", "mtime": time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(JOINT)))}
    joint = ad.read_h5ad(JOINT, backed="r")
    out["obs_cols"] = list(joint.obs.columns); out["uns_keys"] = list(joint.uns.keys()); out["obsp_keys"] = list(joint.obsp.keys())
    sp_ = joint.obs["species"].astype(str).values; names = joint.obs_names.astype(str).values
    K = int(np.asarray(joint.uns["mapping_K"]).ravel()[0]); out["K"] = K
    conn = joint.obsp["connectivities"].tocsr()
    hs_m, mm_m = sp_ == "hs", sp_ == "mm"
    out["n_hs"], out["n_mm"] = int(hs_m.sum()), int(mm_m.sum())
    Chm = conn[hs_m][:, mm_m].tocsr(); Cmh = conn[mm_m][:, hs_m].tocsr()
    with h5py.File(HUMAN, "r") as f:
        hlab = pd.Series(np.asarray(read_elem(f["obs/cell_type"])).astype(str), index=pd.Index(read_elem(f["obs/_index"])).astype(str))
    with h5py.File(MOUSE, "r") as f:
        mlab = pd.Series(np.asarray(read_elem(f["obs/cell_type"])).astype(str), index=pd.Index(read_elem(f["obs/_index"])).astype(str))
    lh = hlab.reindex(names[hs_m]).values.astype(str); lm = mlab.reindex(names[mm_m]).values.astype(str)
    out["unmatched_labels"] = {"human": int((lh == "nan").sum()), "mouse": int((lm == "nan").sum())}
    # per-cell cross-species edge totals
    th = np.asarray(Chm.sum(1)).ravel(); tm = np.asarray(Cmh.sum(1)).ravel()
    out["edge_total_mean"] = {f"hs:{a}": float(th[lh == a].mean()) for a in sorted(set(lh))}
    out["edge_total_mean"].update({f"mm:{b}": float(tm[lm == b].mean()) for b in sorted(set(lm))})
    out["edge_total_max_possible"] = K
    out["cells_without_cross_edges"] = {"human": int((th == 0).sum()), "mouse": int((tm == 0).sum())}
    # per-direction alignment (mean weight / K) and max-symmetrised score
    def onehot(lab):
        cats = pd.Index(sorted(pd.unique(lab)))
        return cats, ss.csr_matrix((np.ones(len(lab)), (np.arange(len(lab)), cats.get_indexer(lab))), shape=(len(lab), len(cats)))
    ch, Mh = onehot(lh); cm, Mm = onehot(lm)
    Sh = np.asarray((Chm @ Mm).todense()); Sm = np.asarray((Cmh @ Mh).todense())
    HM = pd.DataFrame(np.vstack([Sh[lh == a].mean(0) for a in ch]) / K, index=ch, columns=cm)
    MH = pd.DataFrame(np.vstack([Sm[lm == b].mean(0) for b in cm]).T / K, index=ch, columns=cm)
    out["align_h2m"] = HM.round(6).to_dict(); out["align_m2h"] = MH.round(6).to_dict()
    out["align_max_sym"] = np.maximum(HM, MH).round(6).to_dict()
    # label transfer: argmax of connectivity-weighted vote of mouse neighbours' coarse label, for every human cell
    labs = [b for b in cm if b != "nan"]
    V = Sh[:, [list(cm).index(b) for b in labs]]
    has = V.sum(1) > 0
    lab_t = np.where(has, np.array(labs)[V.argmax(1)], "")
    out["label_vote_cols"] = labs
    out["label_covered"] = int(has.sum()); out["label_total"] = int(hs_m.sum())
    for a in ("endometrial_fibroblast", "decidual_cell"):
        m = has & (lh == a); out[f"label_agree_{a}"] = float((lab_t[m] == a).mean()) if m.any() else None
    m = has & np.isin(lh, labs); out["label_agree_all"] = float((lab_t[m] == lh[m]).mean())
    out["label_denominators"] = {"covered_with_native_label": int(m.sum()), "native_label_counts": {a: int((lh == a).sum()) for a in sorted(set(lh))}}
    # any stored label-transfer column
    cand = [c for c in joint.obs.columns if re.search(r"label|transfer|mapped|samap", c, re.I)]
    out["stored_label_cols"] = {c: joint.obs[c].astype(str).value_counts().head(6).to_dict() for c in cand[:8]}
    pd.DataFrame({"cell": names[hs_m], "native": lh, "transferred": lab_t, "covered": has, "edge_total": th}).to_csv(f"{OUT}/samap_hs_label_transfer.csv.gz", index=False, compression="gzip")
    return out
rec("samap_graph", samap_graph)

F["_meta"]["finished"] = time.strftime("%Y-%m-%d %H:%M:%S"); F["_meta"]["runtime_s"] = round(time.time() - T0)
json.dump(F, open(f"{OUT}/facts2.json", "w"), indent=1, default=str)
print("done", round(time.time() - T0), "s |", sum(1 for v in F.values() if isinstance(v, dict) and "_error" in v), "items failed", flush=True)
