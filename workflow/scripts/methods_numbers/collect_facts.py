#!/usr/bin/env python
"""Collect the raw facts behind the numbers in the computational methods (xmens / Cevrim et al.).

Usage: python collect_facts.py OUT_DIR        (OUT_DIR must not exist)
Reads pipeline outputs, logs, configs, run records and installed environments; never writes outside OUT_DIR.
The four objects are those of the release config.yaml names (deliverables.release; paths from its MANIFEST.json).
Writes OUT_DIR/facts.json plus per-nucleus pseudotime tables (pt_cells_{mouse,human}.csv.gz) and per-file texts.
Every item is collected independently: a failure is recorded as {"_error": ...} and the rest still runs.
"""
import os, sys, re, json, glob, time, subprocess, traceback
import gzip
import shutil
import numpy as np, pandas as pd

from release_objects import read_release

T0 = time.time()
OUT = sys.argv[1]; os.makedirs(OUT)
L = "/n/lab_storage/eddy_lab/users/nhilgert"; X = f"{L}/xmens_final"; N = "/n/netscratch/eddy_lab/Everyone/nhilgert"
O22 = f"{N}/2026_10_22_flagged_harmonized"  # provenance/ of the pseudotime of record
RELEASE = read_release(X)
OBJ = RELEASE["paths"]  # registry tag -> released object

# Executed runs of the epithelium notebooks (04c_*); their outputs are not yet written under results/.
EPI_HUMAN = f"{N}/2026_10_29b_human_epithelial_layers"
# v9 run (2026-10-30): v8 decisions (cluster 28 EXCLUDE), cell types around each nucleus from release 2026_10_30;
# its CSV outputs are byte-identical to the v8 run (2026_10_29_v8_mouse_epithelial_compartments).
EPI_MOUSE = f"{N}/2026_10_30_v9_mouse_epithelial_compartments"
# Release 2026_10_29 reinstated human epithelial clusters 11, 22 and 44 (04g reinstate_release_cells) and wrote the
# record; 2026_10_29b and 2026_10_30 were rebuilt by other rules, which write no reinstatement record, and keep the same
# epithelial nuclei (checked 30 Oct 2026: per cluster, the nuclei in each release equal n_in_release + n_added).
REINSTATED_RELEASE = "2026_10_29"
SCENIC_REV = f"{L}/decidualization/revisions/scenic_v2_2026_08_28"
IVD = f"{L}/cagri_ivd_end_to_end_20260815"; HIR = f"{L}/human_ivd_repro"
SAMAP_RUN = f"{L}/code/samap-scvi/full_fibdec"
F = {"_meta": {"started": time.strftime("%Y-%m-%d %H:%M:%S"), "host": os.uname().nodename, "python": sys.version.split()[0]}}


def rec(name, fn):
    t = time.time()
    try:
        F[name] = fn()
    except Exception as e:
        F[name] = {"_error": f"{type(e).__name__}: {e}"[:400], "_trace": traceback.format_exc()[-600:]}
    print(f"[{time.time() - T0:6.0f}s] {name} ({time.time() - t:.1f}s){' ERROR' if isinstance(F[name], dict) and '_error' in F[name] else ''}", flush=True)


def sh(cmd):
    return subprocess.run(["bash", "-c", cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True).stdout


def jload(p):
    return json.load(open(p))


def text(p, n=None):
    s = open(p, errors="replace").read()
    return s if n is None else s[:n]


def mtime(p):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p))) if os.path.exists(p) else None


# ------------------------------------------------------------------ configs and code texts
import yaml
rec("config", lambda: yaml.safe_load(open(f"{X}/config/config.yaml")))
rec("lineage_spec", lambda: yaml.safe_load(open(f"{X}/resources/lineage_spec.yaml")))
rec("git_head", lambda: sh(f"cd {X} && git log -1 --format='%h %ci %s'").strip())
rec("zone_genesets", lambda: {k: list(v) for k, v in yaml.safe_load(open(f"{X}/resources/signatures/aucell_modules_mouse.yaml")).items()})
CODE = {"harmonized_pseudotime.py": f"{O22}/provenance/harmonized_pseudotime.py",
        "integrate.py": f"{X}/workflow/scripts/atlas/integrate.py", "reintegrate.py": f"{X}/workflow/scripts/atlas/reintegrate.py",
        "qc.py": f"{X}/workflow/scripts/ingest_qc/qc.py", "09_homology.smk": f"{X}/workflow/rules/09_homology.smk",
        "scenic_make_loom.py": f"{X}/workflow/scripts/scenic/scenic_make_loom.py", "ctx_only.py": f"{SCENIC_REV}/scripts/ctx_only.py",
        "build_consensus_aucell.py": f"{SCENIC_REV}/scripts/build_consensus_aucell.py", "scenic_README.md": f"{SCENIC_REV}/README.md",
        "crossspecies_gene_selection.py": f"{X}/results/mi_dpt/overlap_harmonized_nall/crossspecies_heatmap_selection/crossspecies_gene_selection.py",
        "samap_FULL_REPORT.md": f"{SAMAP_RUN}/FULL_REPORT.md", "samap_INTEGRATION_SUMMARY.md": f"{SAMAP_RUN}/INTEGRATION_SUMMARY.md",
        "samap_HANDOFF.md": f"{SAMAP_RUN}/HANDOFF.md", "label_stain_full.py": f"{SAMAP_RUN}/label_stain_full.py",
        "samap_supp_README.md": f"{N}/figures/samap_supp/README.md", "SAMAP_FIGURES_HANDOFF.md": f"{N}/figures/samap_supp/SAMAP_FIGURES_HANDOFF.md",
        "scenic_PROVENANCE.md": f"{N}/figures/_data/scenic/PROVENANCE.md", "ivd_README.md": f"{IVD}/README.md",
        "build_xmens_v2_tables.py": f"{IVD}/build_xmens_v2_tables.py", "hir_build_consensus.py": f"{HIR}/build_consensus.py",
        "hir_methods.md": f"{HIR}/methods.md", "orthology_README.md": f"{X}/resources/orthology/README.md"}
os.makedirs(f"{OUT}/code"); code_index = {}
for k, p in CODE.items():
    if os.path.exists(p):
        open(f"{OUT}/code/{k}", "w").write(text(p)); code_index[k] = {"path": p, "mtime": mtime(p), "bytes": os.path.getsize(p)}
    else:
        code_index[k] = {"path": p, "missing": True}
F["code_index"] = code_index
rec("label_stain_log", lambda: {"path": sorted(glob.glob(f"{SAMAP_RUN}/logs/label_stain_*.out"))[-1],
                                "tail": text(sorted(glob.glob(f"{SAMAP_RUN}/logs/label_stain_*.out"))[-1])[-4000:]})


# ------------------------------------------------------------------ QC logs
def qc_logs():
    out = {}
    for p in sorted(glob.glob(f"{X}/logs/01_qc_*.log")):
        s = text(p); smp = os.path.basename(p)[len("01_qc_"):-4]
        def g(pat, cast=int):
            m = re.search(pat, s)
            return cast(m.group(1)) if m else None
        out[smp] = dict(path=p, species=g(r"species=(\w+)", str), loaded=g(r"loaded (\d+) cells x"), n_genes=g(r"cells x (\d+) genes"),
                        min_counts=g(r"min (\d+) counts: removed"), umi_removed=g(r"min \d+ counts: removed (\d+)"),
                        mito_removed=g(r"filter pct_counts_mt: removed (\d+)"), counts_removed=g(r"filter total_counts: removed (\d+)"),
                        genes_removed=g(r"filter n_genes_by_counts: removed (\d+)"), doublets=g(r"flagged (\d+) doublets"),
                        final=g(r"final: (\d+) cells"), dd_line=g(r"(DoubletDetection [^\n]+)", str), mito_line=g(r"(mito > [^\n]+)", str),
                        mad_line=g(r"(\d+-MAD total_counts[^\n]*)", str))
    return out
rec("qc_logs", qc_logs)


# ------------------------------------------------------------------ annotation decision tables and atlas clusters
def decision_table_paths():
    """Path of each annotation decision table, keyed by file name, as release v9 uses it: the round-2 tables of
    config annotation_round2 for the atlas ({species}_lineage.csv) and the round-2 lineages, the July tables under
    resources/annotations/ for every other step. small_reads.py applies the same rule to the table texts, from which
    the registry takes the atlas EXCLUDE clusters."""
    paths = {os.path.basename(p): p for p in sorted(glob.glob(f"{X}/resources/annotations/*.csv"))}
    round2 = yaml.safe_load(open(f"{X}/config/config.yaml"))["annotation_round2"]
    for sp, table in round2["atlas_decisions"].items():
        paths[f"{sp}_lineage.csv"] = f"{X}/{table}"
    for sp, lineages in round2["lineages"].items():
        for lineage in lineages:
            paths[f"{sp}_{lineage}.csv"] = f"{X}/{round2['decisions_dir']}/{sp}_{lineage}.csv"
    return paths


def annot_tables():
    out = {}
    for b, p in decision_table_paths().items():
        if b.startswith("_"): continue
        d = pd.read_csv(p, dtype=str)
        lab = d["label"].str.strip() if "label" in d else pd.Series(dtype=str)
        out[b] = dict(path=p, cols=list(d.columns), n_rows=len(d), n_exclude=int(lab.str.upper().str.startswith("EXCLU").sum()),
                      labels=sorted(set(lab[~lab.str.upper().str.startswith("EXCLU")])), mtime=mtime(p),
                      res=sorted(set(d["res"].dropna())) if "res" in d else None)
    return out
rec("annot_tables", annot_tables)


def atlas_clusters():
    out = {}
    for sp in ("mouse", "human"):
        D = f"{X}/results/tables/clusters/{sp}"
        allr = pd.read_csv(f"{D}/leiden_all_resolutions.csv", nrows=3)
        c = pd.read_csv(f"{D}/leiden_res2.5.csv")
        cc = [x for x in c.columns if x not in (c.columns[0],)][0]
        out[sp] = dict(file=f"{D}/leiden_res2.5.csv", cols=list(c.columns), n_cells=len(c), sizes={str(k): int(v) for k, v in c[cc].value_counts().items()},
                       all_res_cols=list(allr.columns), all_res_file=f"{D}/leiden_all_resolutions.csv", mtime=mtime(f"{D}/leiden_res2.5.csv"))
        c.to_csv(f"{OUT}/atlas_res2.5_{sp}.csv.gz", index=False, compression="gzip")
    return out
rec("atlas_clusters", atlas_clusters)


# ------------------------------------------------------------------ h5ad readers (obs/uns/obsm/var only, plus count detection)
import h5py
from anndata.io import read_elem


def obs_cols(h, cols):
    idx = pd.Index(np.asarray(read_elem(h["obs/_index"])).astype(str)) if "_index" in h["obs"] else pd.Index(np.asarray(read_elem(h["obs"][h["obs"].attrs["_index"]])).astype(str))
    d = {}
    for c in cols:
        if c in h["obs"]: d[c] = read_elem(h[f"obs/{c}"])
    return pd.DataFrame(d, index=idx)


def to_bool(v):
    v = pd.Series(v)
    return v.astype(bool).to_numpy() if v.dtype == bool else v.astype(str).isin(["True", "true", "1", "1.0"]).to_numpy()


def jsonable(x, depth=0):
    if depth > 6: return str(x)[:200]
    if isinstance(x, dict): return {str(k): jsonable(v, depth + 1) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [jsonable(v, depth + 1) for v in x][:200]
    if isinstance(x, np.ndarray): return jsonable(x.tolist(), depth + 1) if x.size <= 200 else f"array{list(x.shape)}"
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.floating,)): return float(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    if isinstance(x, bytes): return x.decode(errors="replace")
    if isinstance(x, pd.DataFrame): return {"_df_shape": list(x.shape), "head": x.head(5).astype(str).to_dict()}
    return x if isinstance(x, (str, int, float, bool)) or x is None else str(x)[:300]


OBS_WANT = ["cell_type", "cell_subtype", "timepoint", "sample", "spatial_qc_reason", "x_um_qc", "y_um_qc", "embedding_outlier",
            "pt_harmonized", "pt_harmonized_rank", "dpt_pseudotime", "dpt_pseudotime_no_apoptotic", "dpt_pseudotime_rank", "lineage", "doublet_score"]


def objects():
    out = {}
    fibdec_names = {}  # species -> obs names of its fibroblast/decidual object
    for tag, p in OBJ.items():
        sp = tag.split("_")[0]
        with h5py.File(p, "r") as h:
            obs_keys = list(h["obs"].keys())
            O = obs_cols(h, OBS_WANT)
            r = dict(path=p, mtime=mtime(p), n_obs=len(O), obs_keys=obs_keys, n_vars=int(h["var/_index"].shape[0]) if "_index" in h["var"] else None,
                     obsm={k: list(h[f"obsm/{k}"].shape) if isinstance(h[f"obsm/{k}"], h5py.Dataset) else "group" for k in h["obsm"].keys()},
                     layers=list(h["layers"].keys()) if "layers" in h else [], uns_keys=list(h["uns"].keys()))
            for c in ("cell_type", "cell_subtype", "timepoint", "spatial_qc_reason", "lineage", "sample"):
                if c in O: r[f"vc_{c}"] = {str(k): int(v) for k, v in O[c].astype(str).value_counts().items()}
            if "spatial_qc_reason" in O and "timepoint" in O:
                r["spatial_by_tp"] = {f"{a}|{b}": int(v) for (a, b), v in O.groupby([O.spatial_qc_reason.astype(str), O.timepoint.astype(str)]).size().items()}
            if "sample" in O and "timepoint" in O:
                # 10x lanes (obs['sample']) of each timepoint, i.e. of each array or section
                samples = O.groupby(O.timepoint.astype(str))["sample"].unique()
                r["samples_by_tp"] = {tp: sorted(str(s) for s in names) for tp, names in samples.items()}
            if "x_um_qc" in O:
                r["n_positioned"] = int(np.isfinite(O.x_um_qc.astype(float)).sum())
            if "embedding_outlier" in O:
                fl = to_bool(O.embedding_outlier); r["n_embedding_outlier"] = int(fl.sum())
                r["embedding_outlier_names"] = sorted(O.index[fl])
            for u in ("embedding_outlier_rule", "pt_harmonized", "epithelial_reinstatement"):
                if u in h["uns"]:
                    try: r[f"uns_{u}"] = jsonable(read_elem(h[f"uns/{u}"]))
                    except Exception as e: r[f"uns_{u}"] = f"unreadable: {e}"
            if "var" in h and "highly_variable" in h["var"]:
                r["n_highly_variable"] = int(np.asarray(read_elem(h["var/highly_variable"])).astype(bool).sum())
            if tag.endswith("fibdec"):
                keep = [c for c in ("cell_type", "cell_subtype", "timepoint", "sample", "pt_harmonized", "pt_harmonized_rank", "dpt_pseudotime",
                                    "dpt_pseudotime_no_apoptotic", "dpt_pseudotime_rank", "embedding_outlier") if c in O]
                O[keep].to_csv(f"{OUT}/pt_cells_{sp}.csv.gz", compression="gzip")
                fibdec_names[sp] = set(O.index)

                # Nuclei per cell type and timepoint: the pseudobulks of the decidual-versus-fibroblast comparison.
                celltype_by_tp = O.groupby([O.cell_type.astype(str), O.timepoint.astype(str)]).size()
                r["celltype_by_tp"] = {f"{ct}|{tp}": int(n) for (ct, tp), n in celltype_by_tp.items() if n > 0}

                # gene detection on raw counts (SCENIC loom rule and n_vars)
                lay = "counts" if "layers" in h and "counts" in h["layers"] else None
                M = read_elem(h[f"layers/{lay}"]) if lay else read_elem(h["X"])
                M = M.tocsr() if hasattr(M, "tocsr") else M
                nz = np.bincount(M.indices[M.data > 0], minlength=M.shape[1]) if hasattr(M, "indices") else (np.asarray(M) > 0).sum(0)
                genes = np.asarray(read_elem(h["var/_index"])).astype(str)
                pd.DataFrame({"gene": genes, "n_cells": nz}).to_csv(f"{OUT}/detect_{sp}.csv.gz", index=False, compression="gzip")
                r["detect_layer"] = lay or "X"; r["detect_shape"] = list(M.shape)
                r["is_integer_counts"] = bool(np.allclose(M.data[:100000], np.round(M.data[:100000]))) if hasattr(M, "data") else None
                del M
        out[tag] = r

    # The all-cells and fibroblast/decidual objects carry their own embedding-outlier flags (the human Fib/Dec object
    # keeps the flags of the atlas fitted before the epithelial reinstatement), so their overlap is recorded here.
    for sp, names in fibdec_names.items():
        all_flags = set(out[f"{sp}_all"]["embedding_outlier_names"])
        fibdec_flags = set(out[f"{sp}_fibdec"]["embedding_outlier_names"])
        out[f"{sp}_all"]["n_embedding_outlier_on_fibdec"] = len(all_flags & names)
        out[f"{sp}_fibdec"]["n_embedding_outlier_shared_with_all"] = len(fibdec_flags & all_flags)
    return out
rec("release", lambda: {"release": RELEASE["release"], "dir": RELEASE["dir"], "manifest": RELEASE["manifest"]})
rec("objects", objects)


def doublets():
    """Per-nucleus DoubletDetection scores of the post-QC nuclei, split by atlas EXCLUDE clusters."""
    out = {}
    for sp in ("mouse", "human"):
        cands = [f"{X}/data/interim/latent_{sp}.h5ad", f"{X}/data/interim/annotated_{sp}.h5ad"]
        src = None
        for p in cands:
            if not os.path.exists(p): continue
            with h5py.File(p, "r") as h:
                ks = [k for k in h["obs"].keys() if "doublet" in k.lower()]
                if ks:
                    O = obs_cols(h, ks); src = (p, ks); break
        if src is None:
            parts = []
            for p in sorted(glob.glob(f"{X}/data/interim/qc/*.h5ad")):
                with h5py.File(p, "r") as h:
                    ks = [k for k in h["obs"].keys() if "doublet" in k.lower()]
                    if ks: parts.append(obs_cols(h, ks))
            O = pd.concat(parts) if parts else None; src = ("data/interim/qc/*.h5ad", list(O.columns) if O is not None else [])
        out[sp] = {"source": src[0], "columns": src[1], "n": None if O is None else len(O)}
        if O is not None:
            O.to_csv(f"{OUT}/doublet_scores_{sp}.csv.gz", compression="gzip")
    return out
rec("doublets", doublets)

# ------------------------------------------------------------------ pseudotime run records
rec("pt_provenance", lambda: {sp: jloads for sp, jloads in ((sp, jload(f"{O22}/provenance/{sp}_pca_released_provenance.json")) for sp in ("mouse", "human"))})
rec("o22_manifest", lambda: jload(f"{O22}/MANIFEST.json"))


# ------------------------------------------------------------------ mutual information and overlap
def mi():
    out = {}
    for run in ("pt_harmonized_mouse_nall", "pt_harmonized_human_nall", "dpt_noapop_mouse_nall", "dpt_full_human_nall"):
        D = f"{X}/results/mi_dpt/{run}"
        out[run] = {"manifest": jload(f"{D}/manifest.json"), "mi_all": f"{D}/mi_all.csv", "mtime": mtime(f"{D}/mi_all.csv")}
        pd.read_csv(f"{D}/mi_all.csv").to_csv(f"{OUT}/mi_{run}.csv.gz", index=False, compression="gzip")
    O = f"{X}/results/mi_dpt/overlap_harmonized_nall"
    out["overlap_report"] = jload(f"{O}/overlap_report.json")
    sw = pd.read_csv(f"{O}/overlap_cutoff_sweep.csv"); out["sweep"] = {"file": f"{O}/overlap_cutoff_sweep.csv", "cols": list(sw.columns),
                                                                       "frac_min": float(sw.frac.min()), "frac_max": float(sw.frac.max()), "n": int(sw.frac.nunique())}
    S = f"{O}/crossspecies_heatmap_selection"
    out["heatmap_selection"] = jload(f"{S}/selection_summary.json")
    g = pd.read_csv(f"{S}/conserved_profiles_selected_genes.csv"); out["heatmap_selected"] = {"cols": list(g.columns), "n": len(g), "rows": g.head(40).astype(str).to_dict("records")}
    orth = pd.read_csv(f"{X}/results/mi_dpt/overlap/ortholog_pairs_one2one.csv")
    out["orth"] = {"file": f"{X}/results/mi_dpt/overlap/ortholog_pairs_one2one.csv", "n_pairs": len(orth), "cols": list(orth.columns), "mtime": mtime(f"{X}/results/mi_dpt/overlap/ortholog_pairs_one2one.csv")}
    orth.to_csv(f"{OUT}/orth.csv.gz", index=False, compression="gzip")
    return out
rec("mi", mi)

# ------------------------------------------------------------------ SAMap
def samap():
    D = f"{N}/figures/samap_supp/data"
    out = {"edge_share": pd.read_csv(f"{D}/samap_compartment_edge_share.csv").to_dict("records"),
           "transfer_summary": jload(f"{X}/results/samap_raw/pt_transfer_noapop_joinfix/summary.json"),
           "sfig_provenance_rho": jload(f"{N}/figures/samap_supp/sfig_samap.provenance.json").get("spearman_rho")}
    hs = pd.read_csv(f"{D}/hs_transfer_spatial.csv")
    out["hs_transfer_cols"] = list(hs.columns); out["hs_transfer_n"] = len(hs)
    hs.to_csv(f"{OUT}/hs_transfer_spatial.csv.gz", index=False, compression="gzip")
    out["run_dir_listing"] = sorted(os.listdir(SAMAP_RUN))[:200]
    out["meta_files"] = sh(f"find {SAMAP_RUN} {L}/code/samap-scvi -maxdepth 3 \\( -name 'meta*.json' -o -name 'eval*.json' -o -name 'mapping*.csv' -o -name 'degree_check*.json' -o -name '*alignment*' \\) 2>/dev/null | head -40").split()
    for p in out["meta_files"][:12]:
        if p.endswith(".json"):
            try: out.setdefault("meta_json", {})[p] = jsonable(jload(p))
            except Exception as e: out.setdefault("meta_json", {})[p] = str(e)
        elif p.endswith(".csv"):
            try: out.setdefault("meta_csv", {})[p] = pd.read_csv(p, nrows=12).astype(str).to_dict("records")
            except Exception as e: out.setdefault("meta_csv", {})[p] = str(e)
    out["grep_numbers"] = sh(f"grep -rnE '0\\.58|12\\.6|13\\.5|47,?509|84\\.2|81\\.3|87\\.0' {SAMAP_RUN}/*.md {N}/figures/samap_supp/*.md {N}/figures/samap_supp/*.json 2>/dev/null | cut -c1-260 | head -30")
    return out
rec("samap", samap)


# ------------------------------------------------------------------ SCENIC
def scenic():
    D = f"{N}/figures/_data/scenic"; out = {}
    for sp in ("mouse", "human"):
        r = pd.read_csv(f"{D}/{sp}_recurrence_report.tsv", sep="\t")
        out[sp] = {"report_cols": list(r.columns), "n_rows": len(r), "summary": text(f"{D}/{sp}_recurrence_summary.txt", 3000),
                   "report_file": f"{D}/{sp}_recurrence_report.tsv", "mtime": mtime(f"{D}/{sp}_recurrence_report.tsv")}
        r.to_csv(f"{OUT}/scenic_recurrence_{sp}.csv.gz", index=False, compression="gzip")
    out["ensembl116_one2one"] = f"{D}/ensembl116_one2one.csv"
    pd.read_csv(f"{D}/ensembl116_one2one.csv").to_csv(f"{OUT}/scenic_ensembl116_one2one.csv.gz", index=False, compression="gzip")
    out["databases"] = sh(f"ls {SCENIC_REV}/databases | head -40").split()
    for t in ("allTFs_mm.txt", "allTFs_hg38.txt"):
        p = sh(f"find {SCENIC_REV}/databases -name '{t}' | head -1").strip()
        if p: out[t] = {"path": p, "tfs": [l.strip() for l in open(p) if l.strip()]}
    out["aucell_log_lines"] = sh(f"grep -rhE 'AUCell \\(|auc_threshold' {SCENIC_REV}/slurm_logs 2>/dev/null | cut -c1-200 | sort | uniq -c | head -10")
    site = glob.glob(f"{L}/.conda/envs/pyscenic/lib/python3*/site-packages/pyscenic")
    out["pyscenic_defaults"] = sh(f"grep -rnE '(rank_threshold|auc_threshold|nes_threshold)\\s*[:=]\\s*[0-9.]+' {site[0]}/*.py 2>/dev/null | cut -c1-200 | head -20") if site else "no pyscenic site-packages"
    out["floor_grep"] = sh(f"grep -rnE 'floor|min_occurrence|occ6|> ?5' {SCENIC_REV}/README.md {N}/figures/_data/scenic/PROVENANCE.md 2>/dev/null | cut -c1-220 | head -20")
    return out
rec("scenic", scenic)


# ------------------------------------------------------------------ in vitro comparison
def ivd():
    out = {}; V2 = f"{IVD}/results/xmens_v2"
    for f in ("consensus_ortholog_xmens_broad_support_rows_v2.csv", "wang_ivd_support_rows_v2.csv", "xmens_broad_all_genes_v2.csv"):
        d = pd.read_csv(f"{V2}/{f}", low_memory=False); out[f] = {"cols": list(d.columns), "n": len(d), "mtime": mtime(f"{V2}/{f}")}
        d.to_csv(f"{OUT}/ivd_{f}.gz", index=False, compression="gzip")
    for f in ("fit_inputs/metadata.tsv", "verify_inputs.tsv", "verify_tables.tsv", "xmens_v2_eligibility.tsv"):
        d = pd.read_csv(f"{V2}/{f}", sep="\t"); out[f] = {"cols": list(d.columns), "n": len(d)}
        d.to_csv(f"{OUT}/ivd_{os.path.basename(f)}.gz", index=False, compression="gzip", sep="\t")
    out["resources"] = sh(f"ls {IVD}/resources | head -60").split()
    for f in out["resources"]:
        if f.endswith((".csv", ".tsv")) and os.path.getsize(f"{IVD}/resources/{f}") < 30e6:
            d = pd.read_csv(f"{IVD}/resources/{f}", sep="\t" if f.endswith(".tsv") else ",", low_memory=False)
            out.setdefault("resource_tables", {})[f] = {"cols": list(d.columns)[:40], "n": len(d)}
            d.to_csv(f"{OUT}/ivdres_{f}.gz", index=False, compression="gzip", sep="\t" if f.endswith(".tsv") else ",")
    out["summary_release"] = pd.read_csv(f"{IVD}/results/cagri_ivd_end_to_end_20260815_ortholog114/ivd_support_summary.csv").to_dict("records")
    out["hir_signatures"] = {}
    for p in sorted(glob.glob(f"{HIR}/signatures/*.csv")) + sorted(glob.glob(f"{HIR}/papers/stadtmauer/*.csv")):
        d = pd.read_csv(p, low_memory=False); out["hir_signatures"][p] = {"cols": list(d.columns)[:30], "n": len(d)}
        d.to_csv(f"{OUT}/hir_{os.path.basename(p)}.gz", index=False, compression="gzip")
    out["hir_results"] = sh(f"ls -R {HIR}/results | head -80")
    return out
rec("ivd", ivd)


# ------------------------------------------------------------------ installed versions and references
def versions():
    want = re.compile(r"^(scanpy|anndata|scvi_tools|torch|scikit_learn|statsmodels|scipy|pydeseq2|moscot|cellrank|pyscenic|ctxcore|arboreto|geosketch|"
                      r"doubletdetection|cellbender|sc_samap|samap|igraph|leidenalg|umap_learn|numpy|pandas|ott_jax|jax)-([^-]+)\.dist-info$", re.I)
    envs = {"rsc_uv": f"{L}/.uv/envs/rsc_uv", "samap_uv": f"{L}/.uv/envs/samap_uv", "pyscenic": f"{L}/.conda/envs/pyscenic",
            "cellrank": f"{L}/.conda/envs/cellrank", "cellbender_022": f"{L}/.conda/envs/cellbender_022", "SAMap": f"{L}/.conda/envs/SAMap",
            "cuda11.8_env": f"{L}/.conda/envs/cuda11.8_env", "mi_env": f"{L}/.conda/envs/mi_env"}
    out = {}
    for e, p in envs.items():
        d = {}
        for f in glob.glob(f"{p}/lib/python3*/site-packages/*.dist-info"):
            m = want.match(os.path.basename(f))
            if m: d[m.group(1).lower().replace("_", "-")] = m.group(2)
        out[e] = {"path": p, "packages": d}
    cb = glob.glob(f"{L}/.conda/envs/cellbender_022/lib/python3*/site-packages/cellbender/remove_background/*.py")
    out["cellbender_epochs_default"] = sh(f"grep -nE 'epochs' {' '.join(cb)} 2>/dev/null | grep -iE 'default|EPOCH' | cut -c1-200 | head -8") if cb else ""
    return out
rec("versions", versions)


def umap_defaults():
    import inspect, scanpy as sc
    sig = inspect.signature(sc.tl.umap); pca = inspect.signature(sc.pp.pca); nb = inspect.signature(sc.pp.neighbors)
    return {"scanpy": sc.__version__, "umap_min_dist": sig.parameters["min_dist"].default, "umap_spread": sig.parameters["spread"].default,
            "pca_n_comps_arg_default": str(pca.parameters["n_comps"].default), "settings_N_PCS": int(sc.settings.N_PCS),
            "neighbors_n_neighbors": nb.parameters["n_neighbors"].default}
rec("scanpy_defaults", umap_defaults)

# ------------------------------------------------------------------ spatial QC table and misc
rec("spatial_qc_table", lambda: pd.read_csv(f"{X}/results/tables/spatial_qc_mouse.csv").to_dict("records"))
rec("misc_find", lambda: {"gpu_provenance": sh(f"find {X} {L}/xmens_final_repro* -maxdepth 4 -name 'gpu_provenance*.json' 2>/dev/null | head -3").split(),
                          "heatmap_code": sh(f"grep -rlE 'plot_pseudotime_heatmap|vgam|VGAM|ns\\(|dispersion_trend' {N}/figures {X}/workflow {X}/notebooks --include=*.py --include=*.R --include=*.ipynb 2>/dev/null | head -10").split(),
                          "yang_files": sh(f"find {N} {X}/resources -maxdepth 4 -iname '*yang*' 2>/dev/null | head -10").split()})


# ------------------------------------------------------------------ v7 epithelium: human layers, mouse compartments
EPI_CLUSTER_KEY = "leiden_scVI_3.0"  # resolution both notebooks annotate at
# Markers behind the reinstatement of the human clusters excluded for PECAM1 co-detection (methods, epithelium).
REINSTATEMENT_GENES = ["PECAM1", "FLT1", "VWF", "CDH5", "PAX8", "PAEP"]


def copy_gz(src, name):
    """Copy a text file into OUT as name.gz, byte for byte."""
    with open(src, "rb") as fh_in, gzip.open(f"{OUT}/{name}.gz", "wb") as fh_out:
        shutil.copyfileobj(fh_in, fh_out)
    return {"path": src, "mtime": mtime(src), "bytes": os.path.getsize(src)}


def notebook_code(path):
    """Code cells of an executed notebook, joined; the parameters it ran with are in its injected cell."""
    cells = json.load(open(path))["cells"]
    return "\n\n".join("".join(cell["source"]) for cell in cells if cell["cell_type"] == "code")


def detection_by_cluster(path, genes, clusters):
    """Share of nuclei with a raw count above zero for each gene (columns), in each listed cluster (rows)."""
    with h5py.File(path, "r") as h:
        obs = obs_cols(h, [EPI_CLUSTER_KEY])
        var_names = pd.Index(np.asarray(read_elem(h["var"][h["var"].attrs["_index"]])).astype(str))
        counts = read_elem(h["layers/counts"]).tocsr()

    gene_columns = var_names.get_indexer(genes)
    if (gene_columns < 0).any():
        raise KeyError(f"genes missing from {path}: {[g for g, i in zip(genes, gene_columns) if i < 0]}")
    is_detected = pd.DataFrame((counts[:, gene_columns] > 0).toarray(), index=obs.index, columns=genes)

    cluster = obs[EPI_CLUSTER_KEY].astype(str)
    in_listed = cluster.isin(clusters)
    return is_detected[in_listed].groupby(cluster[in_listed]).mean()


def epithelial_human():
    """Outputs of the executed 04c_human_epithelial_layers notebook, and marker detection in the reinstated clusters."""
    run_out = f"{EPI_HUMAN}/human_epithelial_layers"
    annotated = f"{run_out}/annotated_human_epithelial.h5ad"
    out = {"run_dir": EPI_HUMAN, "annotated": annotated, "copied": {}}

    out["copied"]["per_cluster"] = copy_gz(f"{run_out}/human_epithelial_layers_per_cluster.csv", "epi_human_per_cluster.csv")
    out["copied"]["ranking"] = copy_gz(f"{run_out}/de_inner_vs_outer_ranking.csv", "epi_human_ranking.csv")
    out["copied"]["by_donor"] = copy_gz(f"{run_out}/de_inner_vs_outer_by_donor.csv", "epi_human_by_donor.csv")

    # Clusters reinstated (n_added > 0), from the record of the release that made the reinstatement.
    reinstated_path = f"{X}/results/tables/reinstated_cells_{REINSTATED_RELEASE}.tsv"
    out["copied"]["reinstated"] = copy_gz(reinstated_path, "epi_human_reinstated.tsv")
    reinstated = pd.read_csv(reinstated_path, sep="\t")
    out["reinstated_clusters"] = reinstated.loc[reinstated.n_added > 0, "cluster"].astype(str).tolist()

    detection = detection_by_cluster(annotated, REINSTATEMENT_GENES, out["reinstated_clusters"])
    detection.to_csv(f"{OUT}/epi_human_reinstated_detection.csv")

    with h5py.File(annotated, "r") as h:
        obs = obs_cols(h, [EPI_CLUSTER_KEY, "cell_subtype", "timepoint", "layer_of_nucleus"])
    obs.to_csv(f"{OUT}/epi_human_obs.csv.gz", compression="gzip")

    notebook = f"{EPI_HUMAN}/04c_human_epithelial_layers.executed.ipynb"
    open(f"{OUT}/code/epi_human_notebook.py", "w").write(notebook_code(notebook))
    out["notebook"] = {"path": notebook, "mtime": mtime(notebook)}
    return out
rec("epithelial_human", epithelial_human)


def epithelial_mouse():
    """Outputs of the executed 04c_mouse_epithelial_compartments notebook."""
    run_out = f"{EPI_MOUSE}/mouse_epithelial_compartments"
    annotated = f"{run_out}/annotated_mouse_epithelial.h5ad"
    out = {"run_dir": EPI_MOUSE, "annotated": annotated, "copied": {}}

    out["copied"]["per_cluster"] = copy_gz(f"{run_out}/mouse_epithelial_compartments_per_cluster.csv", "epi_mouse_per_cluster.csv")
    out["copied"]["ranking"] = copy_gz(f"{run_out}/de_luminal_vs_glandular_ranking.csv", "epi_mouse_ranking.csv")
    out["copied"]["proliferative"] = copy_gz(f"{run_out}/de_proliferative_vs_other_induced.csv", "epi_mouse_proliferative.csv")

    # The lineage object keeps the pre-v5 timepoint names; "uninduced" is the same in both vocabularies.
    with h5py.File(annotated, "r") as h:
        obs = obs_cols(h, [EPI_CLUSTER_KEY, "cell_subtype", "timepoint", "position_vote"])
    obs.to_csv(f"{OUT}/epi_mouse_obs.csv.gz", compression="gzip")

    notebook = f"{EPI_MOUSE}/04c_mouse_epithelial_compartments.executed.ipynb"
    open(f"{OUT}/code/epi_mouse_notebook.py", "w").write(notebook_code(notebook))
    out["notebook"] = {"path": notebook, "mtime": mtime(notebook)}
    return out
rec("epithelial_mouse", epithelial_mouse)

F["_meta"]["finished"] = time.strftime("%Y-%m-%d %H:%M:%S"); F["_meta"]["runtime_s"] = round(time.time() - T0)
json.dump(F, open(f"{OUT}/facts.json", "w"), indent=1, default=str)
print("done", round(time.time() - T0), "s |", sum(1 for v in F.values() if isinstance(v, dict) and "_error" in v), "items failed", flush=True)
