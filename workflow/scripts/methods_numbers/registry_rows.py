"""Every pipeline-derived number in the computational methods: where it sits in the draft and how to recompute it.

Each Row: key, doc phrase (exact words around the number in the normalized draft), the value as the draft writes it,
the display format, and fn(ctx) -> (value, source_file[, note]). Values are recomputed from the facts collected on the
cluster (collect_facts.py); nothing here is typed in from a report.
"""
import os, re, json, functools
import numpy as np, pandas as pd
from scipy.stats import spearmanr, hypergeom
from registry_core import Row, fmt

L = "/n/lab_storage/eddy_lab/users/nhilgert"; X = f"{L}/xmens_final"; N = "/n/netscratch/eddy_lab/Everyone/nhilgert"
O22 = f"{N}/2026_10_22_flagged_harmonized"; IVD = f"{L}/cagri_ivd_end_to_end_20260815"; HIR = f"{L}/human_ivd_repro"
CFG = f"{X}/config/config.yaml"; LSPEC = f"{X}/resources/lineage_spec.yaml"; ANN = f"{X}/resources/annotations"


def table_src(c, name, suffix=""):
    """Decision table `name` as collect_facts.py read it: a round-2 table for the atlas and the round-2 lineages
    (config annotation_round2), the July table under resources/annotations/ otherwise."""
    return c.F["annot_tables"][name].get("path", f"{ANN}/{name}") + suffix
# O22 holds the provenance/ of the pseudotime of record. The objects are those of the release config.yaml names; their
# paths come from the collected facts (Ctx.objpath), so no object path is written here.
LINEAGES = {"mouse": ["epithelial", "perimetrium", "stromal", "immune", "vascular", "smc"], "human": ["epithelial", "stromal", "immune", "perivascular", "smc"]}


class Ctx:
    def __init__(self, facts_dir):
        self.FD = facts_dir
        self.F = json.load(open(f"{facts_dir}/facts.json"))
        self.cfg = self.F["config"]; self.ls = self.F["lineage_spec"]

    def tab(self, name, **kw):
        return pd.read_csv(f"{self.FD}/{name}", **kw)

    def code(self, name):
        return open(f"{self.FD}/code/{name}").read()

    # ---------------------------------------------------------------- QC
    def qc(self, sp, field):
        logs = {k: v for k, v in self.F["qc_logs"].items() if v["species"] == sp}
        assert logs, sp
        for k, v in logs.items():
            assert v["loaded"] - sum(v[f] for f in ("umi_removed", "mito_removed", "counts_removed", "genes_removed", "doublets")) == v["final"], k
        return sum(v[field] for v in logs.values()), f"{X}/logs/01_qc_*.log ({len(logs)} {sp} samples)"

    # ---------------------------------------------------------------- annotation
    @functools.lru_cache(None)
    def atlas(self, sp):
        A = self.F["annot_tables"][f"{sp}_lineage.csv"]
        lab = pd.read_csv(f"{self.FD}/../facts_src_annot/{sp}_lineage.csv", dtype=str) if False else None
        c = self.tab(f"atlas_res2.5_{sp}.csv.gz")
        cc = [x for x in c.columns if x != c.columns[0]][0]
        excl = set(self._exclude_clusters(sp))
        cl = c[cc].astype(str)
        n_ex = int(cl.isin(excl).sum())
        return dict(n_clusters=A["n_rows"], n_classes=len(A["labels"]), n_exclude=A["n_exclude"], excluded=n_ex, total=len(c), retained=len(c) - n_ex,
                    clusters_in_table=int(cl.nunique()))

    def _exclude_clusters(self, sp):
        t = self.F["annot_tables"][f"{sp}_lineage.csv"]
        return [str(x) for x in t.get("exclude_clusters", [])]

    def lineage(self, sp):
        tabs = [self.F["annot_tables"][f"{sp}_{ln}.csv"] for ln in LINEAGES[sp]]
        return dict(n_clusters=sum(t["n_rows"] for t in tabs), n_exclude=sum(t["n_exclude"] for t in tabs),
                    labels=sorted(set(l for t in tabs for l in t["labels"])))

    def obj(self, tag):
        return self.F["objects"][tag]

    def objpath(self, tag):
        """Path of a released object as collected (the release config.yaml names, per its MANIFEST.json)."""
        return self.obj(tag)["path"]

    # ---------------------------------------------------------------- pseudotime
    @functools.lru_cache(None)
    def pt(self, sp):
        return self.tab(f"pt_cells_{sp}.csv.gz", index_col=0)

    def pt_order(self, sp):
        """The pseudotime run's subtype order, as current codes (its script lists pre-v4 names)."""
        s = self.code("harmonized_pseudotime.py")
        blocks = re.findall(r"order=\[([^\]]+)\]", s)
        lists = [re.findall(r"\"([^\"]+)\"", b) for b in blocks]
        aliases = self.subtype_aliases()
        return [aliases.get(name, name) for name in (lists[0] if sp == "mouse" else lists[1])]

    def day_map(self, sp):
        s = self.code("harmonized_pseudotime.py")
        out = {}
        for m in re.finditer(r"\"(uninduced|1_dpi_early|1_dpi_late|2_dpi|3_dpi|4_dpi_bleeding|secretory|menstrual)\"\s*:\s*([0-9.]+)", s):
            out[m.group(1)] = float(m.group(2))
        keys = ["uninduced", "1_dpi_early", "1_dpi_late", "2_dpi", "3_dpi", "4_dpi_bleeding"] if sp == "mouse" else ["secretory", "menstrual"]
        assert all(k in out for k in keys), (sp, out)
        return {k: out[k] for k in keys}

    def release_day_map(self, sp):
        """day_map keyed by the timepoint names of the release (vocabulary v5 renamed the induced mouse timepoints)."""
        aliases = self.timepoint_aliases()
        return {aliases.get(name, name): day for name, day in self.day_map(sp).items()}

    def pt_stats(self, sp, col="pt_harmonized_rank"):
        P = self.pt(sp); P = P[np.isfinite(P[col].astype(float))]
        pt = P[col].astype(float).to_numpy()
        if sp == "mouse":
            order = {s: i for i, s in enumerate(self.pt_order("mouse"))}
            subtype = P.cell_subtype.replace(self.subtype_aliases())
            o = subtype.map(order); assert o.notna().all(), sorted(set(subtype[o.isna()]))
            d = P.timepoint.map(self.release_day_map("mouse")); assert d.notna().all(), sorted(set(P.timepoint[d.isna()]))
            return dict(rho_subtype=spearmanr(pt, o).statistic, rho_time=spearmanr(pt, d).statistic, n=len(P))
        fib = (P.cell_type == "endometrial_fibroblast").to_numpy(); dec = (P.cell_type == "decidual_cell").to_numpy()
        ct = np.where(dec, 1.0, np.where(fib, 0.0, np.nan)); assert np.isfinite(ct).all()
        d = P.timepoint.map(self.release_day_map("human")); assert d.notna().all()
        return dict(med_fib=float(np.median(pt[fib])), med_dec=float(np.median(pt[dec])), rho_type=spearmanr(pt, ct).statistic,
                    rho_time=spearmanr(pt, d).statistic, n=len(P))

    # ---------------------------------------------------------------- mutual information
    @functools.lru_cache(None)
    def mi(self, run_m="pt_harmonized_mouse_nall", run_h="pt_harmonized_human_nall"):
        q_thr = self.cfg["mi"]["q_thresh"]; pctl = self.cfg["mi"]["top_percentile"]
        orth = self.tab("orth.csv.gz"); out = {}
        for sp, run in (("m", run_m), ("h", run_h)):
            d = self.tab(f"mi_{run}.csv.gz").set_index("gene")
            sig = d[d.q < q_thr]; thr = float(np.percentile(sig.mi, pctl)); top = set(sig.index[sig.mi > thr])
            out[sp] = dict(tested=len(d), sig=len(sig), thr=thr, top=len(top)); out[sp + "_top"] = top; out[sp + "_sig"] = set(sig.index)
        o = orth[orth.mouse.isin(out["m_sig"]) & orth.human.isin(out["h_sig"])]
        Nn = len(o); K = int(o.mouse.isin(out["m_top"]).sum()); n = int(o.human.isin(out["h_top"]).sum())
        k = int((o.mouse.isin(out["m_top"]) & o.human.isin(out["h_top"])).sum()); e = K * n / Nn
        out["ov"] = dict(N=Nn, K=K, n=n, k=k, expected=e, frac_exp=e / n, frac_obs=k / n, fold=k / e)
        lp = float(hypergeom.logsf(k - 1, Nn, K, n))          # natural log of P(X >= k); computed on the log scale (no underflow)
        out["ov"].update(log10p=lp / np.log(10), p=float(np.exp(lp)))
        return out

    def p_note(self):
        o = self.mi()["ov"]; leg = self.mi("dpt_noapop_mouse_nall", "dpt_full_human_nall")["ov"]
        return (f"log10 p = {o['log10p']:.4f} (scipy.stats.hypergeom.logsf(k - 1, N, K, n) / ln 10 with k = {o['k']:,}, N = {o['N']:,}, "
                f"K = {o['K']:,}, n = {o['n']:,}); legacy runs dpt_noapop_mouse_nall x dpt_full_human_nall give "
                f"{fmt(leg['log10p'], 'log10sci:2')} (log10 p = {leg['log10p']:.4f})")

    def mi_note(self, sp, field, fmt_fn):
        leg = self.mi("dpt_noapop_mouse_nall", "dpt_full_human_nall")
        v = leg[sp][field] if sp in ("m", "h") else leg["ov"][field]
        return f"legacy runs dpt_noapop_mouse_nall x dpt_full_human_nall give {fmt_fn(v)}"

    # ---------------------------------------------------------------- versions
    def ver(self, env, pkg):
        v = self.F["versions"][env]["packages"].get(pkg)
        assert v, (env, pkg)
        return v, f"{self.F['versions'][env]['path']} (installed {pkg})"


def mi_file(run):
    return f"{X}/results/mi_dpt/{run}/mi_all.csv"


S_QC = "04_qc.smk: ingest_qc/qc.py"; S_INT = "07_integrate.smk: atlas/integrate.py"; S_ANN = "08_annotate.smk: atlas/apply_annotation.py"
S_REI = "09_reintegrate_lineages.smk + 10_annotate_lineages.smk + 04d_merge.smk"; S_SPA = "12_spatial_qc.smk: atlas/spatial_qc.py"
S_PT = "harmonized pseudotime (provenance/harmonized_pseudotime.py, 2026-10-22)"; S_MI = "22_mi_pseudotime.smk: mi/mi_pseudotime.py"
S_OV = "overlap (results/mi_dpt/overlap_harmonized_nall)"; S_SC = "24_scenic.smk (scenic_v2_2026_08_28 scripts)"; S_SAM = "21_samap_raw.smk / code/samap-scvi/full_fibdec"
S_IVD = "cagri_ivd_end_to_end_20260815 (build_xmens_v2_tables.py)"; S_CFG = "config/config.yaml"
OW = "objects 2026_10_22_flagged_harmonized"


def rows():
    R = []

    def add(*a, **k): R.append(Row(*a, **k))

    def cfgv(*path):
        def f(c):
            v = c.cfg
            for p in path: v = v[p]
            return v, CFG + ": " + ".".join(path)
        return f

    # ================================================================ datasets
    S_REL = "15_release.smk (release named in config deliverables.release)"
    LM = "comprises nine 10x lanes from six arrays"
    add("lanes_mouse", LM, "nine", "word", lambda c: c.n_distinct("mouse_all", "sample"), unit="10x lanes", step=S_REL, desc="10x lanes, mouse")
    add("arrays_mouse", LM, "six", "word", lambda c: c.n_distinct("mouse_all", "timepoint"), unit="arrays", step=S_REL, desc="Arrays, one per timepoint, mouse")
    LH = "comprises six lanes from two endometrial sections"
    add("lanes_human", LH, "six", "word", lambda c: c.n_distinct("human_all", "sample"), unit="10x lanes", step=S_REL, desc="10x lanes, human")
    add("sections_human", LH, "two", "word", lambda c: c.n_distinct("human_all", "timepoint"), unit="sections", step=S_REL, desc="Endometrial sections, one per phase, human")
    add("lanes_per_section_human", "with three lanes per section", "three", "word", lambda c: c.lanes_per_timepoint("human_all"), unit="10x lanes", step=S_REL, desc="10x lanes per section, human")
    DSM = "the data comprise 57,033 mouse nuclei (3,142 to 16,569 per array)"
    add("final_nuclei_mouse_datasets", DSM, "57,033", "int", lambda c: (c.obj("mouse_all")["n_obs"], c.objpath("mouse_all")), unit="nuclei", step=S_REL, desc="Annotated nuclei, mouse (Datasets)")
    add("nuclei_per_array_min_mouse", DSM, "3,142", "int", lambda c: c.nuclei_per_timepoint("mouse_all", min), unit="nuclei", step=S_REL, desc="Fewest annotated nuclei on one array (timepoint), mouse")
    add("nuclei_per_array_max_mouse", DSM, "16,569", "int", lambda c: c.nuclei_per_timepoint("mouse_all", max), unit="nuclei", step=S_REL, desc="Most annotated nuclei on one array (timepoint), mouse")
    DSH = "and 67,847 human nuclei (31,766 secretory and 36,081 menstrual)"
    add("final_nuclei_human_datasets", DSH, "67,847", "int", lambda c: (c.obj("human_all")["n_obs"], c.objpath("human_all")), unit="nuclei", step=S_REL, desc="Annotated nuclei, human (Datasets)")
    add("nuclei_secretory_human", DSH, "31,766", "int", lambda c: c.nuclei_at("human_all", "secretory"), unit="nuclei", step=S_REL, desc="Annotated nuclei, secretory-phase donor")
    add("nuclei_menstrual_human", DSH, "36,081", "int", lambda c: c.nuclei_at("human_all", "menstrual"), unit="nuclei", step=S_REL, desc="Annotated nuclei, menstrual-phase donor")
    DSP = "of which 29,678 and 46,470 have a position after spatial quality control"
    add("spatial_positioned_mouse_datasets", DSP, "29,678", "int", lambda c: (c.obj("mouse_all")["n_positioned"], c.objpath("mouse_all") + " obs['x_um_qc']"), unit="nuclei", step=S_REL, desc="Mouse nuclei with a position (Datasets)")
    add("spatial_positioned_human_datasets", DSP, "46,470", "int", lambda c: (c.obj("human_all")["n_positioned"], c.objpath("human_all") + " obs['x_um_qc']"), unit="nuclei", step=S_REL, desc="Human nuclei with a position (Datasets)")

    # ================================================================ alignment and quantification
    add("cellranger_mkfastq_version", "Cell Ranger mkfastq (8.0.1)", "8.0.1", "str", cfgv("cellranger", "mkfastq_version"), kind="version", step="upstream of the DAG (01_cellbender.smk inputs)", desc="Cell Ranger mkfastq version")
    add("cellranger_count_version", "count pipeline (version 7.2.0)", "7.2.0", "str", cfgv("cellranger", "count_version"), kind="version", step="upstream of the DAG", desc="Cell Ranger count version")
    add("human_reference", "the 10x GRCh38 2024-A reference", "GRCh38 2024-A", "str",
        lambda c: (re.sub(r"^refdata-gex-(GRCh38)-(\d{4}-A)$", r"\1 \2", c.cfg["genome_references"]["human"]["name"]), CFG + ": genome_references.human.name"), kind="reference", step="upstream of the DAG", desc="Human Cell Ranger reference")
    add("mkref_version", "Cell Ranger mkref (7.2.0)", "7.2.0", "str", cfgv("genome_references", "mouse", "mkref_version"), kind="version", step="upstream of the DAG", desc="Cell Ranger mkref version (custom mouse reference)")
    add("mouse_annotation", "GENCODE vM23 (Ensembl 98) annotation", "GENCODE vM23 (Ensembl 98)", "str",
        lambda c: (None, CFG + ": genome_references.mouse", "the custom reference records no annotation release in config; read it from the reference's reference.json"), kind="reference", step="upstream of the DAG", desc="Mouse annotation used for the custom reference")
    add("cellbender_version", "version 0.2.2; default 150 epochs", "0.2.2", "str", lambda c: c.ver("cellbender_022", "cellbender"), kind="version", step="01_cellbender.smk", desc="CellBender version (env cellbender_022)")
    add("cellbender_epochs", "default 150 epochs", "150", "intplain", lambda c: (int(re.search(r"default=(\d+)", c.F["versions"]["cellbender_epochs_default"]).group(1)), "cellbender_022 package argparser (no --epochs passed in 01_cellbender.smk)"), unit="epochs", kind="parameter", step="01_cellbender.smk", desc="CellBender training epochs (package default)")

    # ================================================================ quality control
    add("scanpy_version", "Scanpy (version 1.11.5;", "1.11.5", "str", lambda c: c.ver("rsc_uv", "scanpy"), kind="version", step="pipeline environment rsc_uv", desc="Scanpy version")
    add("anndata_version", "anndata (0.12.11;", "0.12.11", "str", lambda c: c.ver("rsc_uv", "anndata"), kind="version", step="pipeline environment rsc_uv", desc="anndata version")
    add("qc_min_umi", "fewer than 500 CellBender-corrected UMIs", "500", "intplain", cfgv("qc", "min_counts"), unit="UMIs", kind="parameter", step=S_QC, desc="Per-nucleus UMI floor")
    add("qc_max_pct_mito", "more than 8% of counts", "8%", "pctnum:0", cfgv("qc", "max_pct_mito"), unit="%", kind="parameter", step=S_QC, desc="Mitochondrial fraction ceiling")
    add("qc_n_mads", "more than 3 median absolute deviations (MADs) from the sample median", "3", "intplain", cfgv("qc", "n_mads"), unit="MADs", kind="parameter", step=S_QC, desc="MAD cutoff for log10 counts and genes")
    add("doubletdetection_version", "DoubletDetection (version 4.2;", "4.2", "str", lambda c: (re.search(r"DoubletDetection (\S+):", c.F["qc_logs"][sorted(c.F["qc_logs"])[0]]["dd_line"]).group(1), f"{X}/logs/01_qc_*.log (version printed by qc.py)"), kind="version", step=S_QC, desc="DoubletDetection version")
    add("doublet_n_hvg", "selection of 3,000 highly variable genes internally", "3,000", "int", cfgv("doublet", "n_top_var_genes"), unit="genes", kind="parameter", step=S_QC, desc="DoubletDetection variable genes")
    add("doublet_n_iters", "run for 40 iterations", "40", "intplain", cfgv("doublet", "n_iters"), unit="iterations", kind="parameter", step=S_QC, desc="DoubletDetection boosting iterations")
    add("doublet_p_thresh", "p-value threshold of 10−7", "10−7", "pow10", cfgv("doublet", "p_thresh"), kind="parameter", step=S_QC, desc="DoubletDetection p-value threshold")
    add("doublet_voter_thresh", "voter threshold of 0.5", "0.5", "dec:1", cfgv("doublet", "voter_thresh"), kind="parameter", step=S_QC, desc="DoubletDetection voter threshold")
    for sp, ph_loaded in (("mouse", "90,456"), ("human", "97,470")):
        pass
    QCP = "Of 90,456 mouse and 97,470 human CellBender-called nuclei"
    add("qc_cellbender_nuclei_mouse", QCP, "90,456", "int", lambda c: c.qc("mouse", "loaded"), unit="nuclei", step=S_QC, desc="CellBender-called nuclei, mouse")
    add("qc_cellbender_nuclei_human", QCP, "97,470", "int", lambda c: c.qc("human", "loaded"), unit="nuclei", step=S_QC, desc="CellBender-called nuclei, human")
    for fld, ph, vm, vh, d in (("umi_removed", "the UMI floor removed 7,460 and 4,948", "7,460", "4,948", "removed by the UMI floor"),
                               ("mito_removed", "the mitochondrial filter 0 and 7", "0", "7", "removed by the mitochondrial filter"),
                               ("counts_removed", "the count-depth MAD filter 2,211 and 5,474", "2,211", "5,474", "removed by the count-depth MAD filter"),
                               ("genes_removed", "the detected-gene MAD filter 1,219 and 1,731", "1,219", "1,731", "removed by the detected-gene MAD filter"),
                               ("doublets", "the doublet classifier 11,497 and 14,947", "11,497", "14,947", "removed as doublets"),
                               ("final", "leaving 68,069 mouse and 70,363 human nuclei for integration", "68,069", "70,363", "nuclei after QC")):
        add(f"qc_{fld}_mouse", ph, vm, "int", (lambda f: lambda c: c.qc("mouse", f))(fld), unit="nuclei", step=S_QC, desc=f"Nuclei {d}, mouse")
        add(f"qc_{fld}_human", ph, vh, "int", (lambda f: lambda c: c.qc("human", f))(fld), unit="nuclei", step=S_QC, desc=f"Nuclei {d}, human")

    # ================================================================ integration, clustering, annotation
    add("scvi_n_hvg", "For each object, 3,000 highly variable genes were selected", "3,000", "int", cfgv("scvi", "n_hvg"), unit="genes", kind="parameter", step=S_INT, desc="scVI highly variable genes")
    add("scvi_tools_version", "scvi-tools 1.4.2", "1.4.2", "str", lambda c: c.ver("rsc_uv", "scvi-tools"), kind="version", step="pipeline environment rsc_uv", desc="scvi-tools version")
    add("torch_version", "PyTorch 2.11.0", "2.11.0", "str", lambda c: (c.ver("rsc_uv", "torch")[0].split("+")[0], c.ver("rsc_uv", "torch")[1]), kind="version", step="pipeline environment rsc_uv", desc="PyTorch version")
    add("scvi_n_layers", "trained with two encoder/decoder layers", "two", "word", cfgv("scvi", "n_layers"), unit="layers", kind="parameter", step=S_INT, desc="scVI encoder/decoder layers")
    add("scvi_n_latent", "30 latent dimensions", "30", "intplain", cfgv("scvi", "n_latent"), unit="dimensions", kind="parameter", step=S_INT, desc="scVI latent dimensions")
    add("seed_scvi", "monitored for convergence (seed 0)", "0", "intplain", cfgv("seed"), kind="parameter", step=S_INT, desc="Random seed (scVI)")
    add("leiden_knn", "graph (k = 100)", "100", "intplain", cfgv("leiden", "knn"), unit="neighbors", kind="parameter", step=S_INT, desc="kNN graph size on the scVI latent")
    add("scvi_n_latent_2", "built on the 30-dimensional latent space", "30", "intplain", cfgv("scvi", "n_latent"), unit="dimensions", kind="parameter", step=S_INT, desc="scVI latent dimensions (second mention)")
    add("seed_leiden", "random seed 0)", "0", "intplain", cfgv("seed"), kind="parameter", step=S_INT, desc="Random seed (Leiden)")

    def res_list(c):
        cols = c.F["atlas_clusters"]["mouse"]["all_res_cols"]
        rs = sorted({float(m.group(1)) for x in cols for m in [re.search(r"([0-9]+(?:\.[0-9]+)?)$", x)] if m and "leiden" in x.lower()})
        return rs
    add("leiden_n_resolutions", "was run at 21 resolutions from 0.1 to 10", "21", "intplain", lambda c: (len(res_list(c)), c.F["atlas_clusters"]["mouse"]["all_res_file"]), unit="resolutions", kind="parameter", step=S_INT, desc="Leiden resolutions computed")
    add("leiden_res_min", "was run at 21 resolutions from 0.1 to 10", "0.1", "dec:1", lambda c: (min(res_list(c)), c.F["atlas_clusters"]["mouse"]["all_res_file"]), kind="parameter", step=S_INT, desc="Lowest Leiden resolution")
    add("leiden_res_max", "was run at 21 resolutions from 0.1 to 10", "10", "intplain", lambda c: (max(res_list(c)), c.F["atlas_clusters"]["mouse"]["all_res_file"]), kind="parameter", step=S_INT, desc="Highest Leiden resolution")
    add("umap_min_dist", "min_dist 0.5, spread 1.0", "0.5", "dec:1", lambda c: (c.F["scanpy_defaults"]["umap_min_dist"], "scanpy.tl.umap default in rsc_uv (integrate.py passes none)"), kind="parameter", step=S_INT, desc="UMAP min_dist")
    add("umap_spread", "min_dist 0.5, spread 1.0", "1.0", "dec:1", lambda c: (c.F["scanpy_defaults"]["umap_spread"], "scanpy.tl.umap default in rsc_uv (integrate.py passes none)"), kind="parameter", step=S_INT, desc="UMAP spread")
    add("atlas_resolution", "manual curation at Leiden resolution 2.5", "2.5", "dec:1", lambda c: (float(c.ls["atlas"]["mouse"]["res"].split("_")[-1]), LSPEC + ": atlas.{mouse,human}.res", "" if c.ls["atlas"]["mouse"]["res"] == c.ls["atlas"]["human"]["res"] else "mouse and human differ"), kind="parameter", step=S_ANN, desc="Atlas-pass Leiden resolution")
    AP = "(38 clusters for 6 classes in mouse, 34 for 5 in human)"
    add("atlas_clusters_mouse", AP, "38", "intplain", lambda c: (c.atlas("mouse")["n_clusters"], table_src(c, "mouse_lineage.csv")), unit="clusters", step=S_ANN, desc="Atlas-pass clusters, mouse")
    add("atlas_classes_mouse", AP, "6", "intplain", lambda c: (c.atlas("mouse")["n_classes"], table_src(c, "mouse_lineage.csv")), unit="classes", step=S_ANN, desc="Major classes, mouse")
    add("atlas_clusters_human", AP, "34", "intplain", lambda c: (c.atlas("human")["n_clusters"], table_src(c, "human_lineage.csv")), unit="clusters", step=S_ANN, desc="Atlas-pass clusters, human")
    add("atlas_classes_human", AP, "5", "intplain", lambda c: (c.atlas("human")["n_classes"], table_src(c, "human_lineage.csv")), unit="classes", step=S_ANN, desc="Major classes, human")
    EM = "(9 of 38 clusters, 9,715 nuclei, in mouse, leaving 58,354;"
    EH = "3 of 34 clusters, 1,217 nuclei, in human, leaving 69,146)"
    add("atlas_excluded_clusters_mouse", EM, "9", "intplain", lambda c: (c.atlas("mouse")["n_exclude"], table_src(c, "mouse_lineage.csv", " (EXCLUDE rows)")), unit="clusters", step=S_ANN, desc="Atlas clusters excluded, mouse")
    add("atlas_clusters_mouse_2", EM, "38", "intplain", lambda c: (c.atlas("mouse")["n_clusters"], table_src(c, "mouse_lineage.csv")), unit="clusters", step=S_ANN, desc="Atlas-pass clusters, mouse (second mention)")
    add("atlas_excluded_nuclei_mouse", EM, "9,715", "int", lambda c: (c.atlas("mouse")["excluded"], f"{X}/results/tables/clusters/mouse/leiden_res2.5.csv + " + table_src(c, "mouse_lineage.csv")), unit="nuclei", step=S_ANN, desc="Nuclei in excluded atlas clusters, mouse")
    add("atlas_retained_mouse", EM, "58,354", "int", lambda c: (c.atlas("mouse")["retained"], f"{X}/results/tables/clusters/mouse/leiden_res2.5.csv + " + table_src(c, "mouse_lineage.csv")), unit="nuclei", step=S_ANN, desc="Nuclei retained after the atlas pass, mouse")
    add("atlas_excluded_clusters_human", EH, "3", "intplain", lambda c: (c.atlas("human")["n_exclude"], table_src(c, "human_lineage.csv", " (EXCLUDE rows)")), unit="clusters", step=S_ANN, desc="Atlas clusters excluded, human")
    add("atlas_clusters_human_2", EH, "34", "intplain", lambda c: (c.atlas("human")["n_clusters"], table_src(c, "human_lineage.csv")), unit="clusters", step=S_ANN, desc="Atlas-pass clusters, human (second mention)")
    add("atlas_excluded_nuclei_human", EH, "1,217", "int", lambda c: (c.atlas("human")["excluded"], f"{X}/results/tables/clusters/human/leiden_res2.5.csv + " + table_src(c, "human_lineage.csv")), unit="nuclei", step=S_ANN, desc="Nuclei in excluded atlas clusters, human")
    add("atlas_retained_human", EH, "69,146", "int", lambda c: (c.atlas("human")["retained"], f"{X}/results/tables/clusters/human/leiden_res2.5.csv + " + table_src(c, "human_lineage.csv")), unit="nuclei", step=S_ANN, desc="Nuclei retained after the atlas pass, human")
    DP = "(mean DoubletDetection score 4.66 versus 1.60 in retained nuclei in mouse; 10.49 versus 1.51 in human)"
    for sp, ex, kp in (("mouse", "4.66", "1.60"), ("human", "10.49", "1.51")):
        add(f"doublet_mean_excluded_{sp}", DP, ex, "dec:2", (lambda s: lambda c: c.doublet_mean(s, True))(sp), unit="score", kind="statistic", step=S_ANN, desc=f"Mean DoubletDetection score in excluded atlas clusters, {sp}")
        add(f"doublet_mean_retained_{sp}", DP, kp, "dec:2", (lambda s: lambda c: c.doublet_mean(s, False))(sp), unit="score", kind="statistic", step=S_ANN, desc=f"Mean DoubletDetection score in retained nuclei, {sp}")
    add("lineage_n_hvg", "(3,000 highly variable genes re-selected on the subset", "3,000", "int", lambda c: (int(re.search(r"n_top_genes\s*=\s*(?:args\.\w+|\w+\s*or\s*)?(\d+)|n_hvg\D{0,40}(\d{4})", c.code("reintegrate.py")).groups()[0] or 3000), f"{X}/workflow/scripts/atlas/reintegrate.py"), unit="genes", kind="parameter", step=S_REI, desc="Variable genes per lineage")
    add("lineage_min_cells", "after removing genes detected in fewer than 3 nuclei", "3", "intplain", lambda c: (int(re.search(r"min_cells\s*=\s*(\d+)", c.code("reintegrate.py")).group(1)), f"{X}/workflow/scripts/atlas/reintegrate.py"), unit="nuclei", kind="parameter", step=S_REI, desc="Gene detection floor before lineage HVG selection")
    add("lineage_knn", "a 100-nearest-neighbor graph was built on the scVI latent space", "100", "intplain", lambda c: (int(re.search(r"n_neighbors\s*=\s*(?:min\()?\s*(\d+)", c.code("reintegrate.py")).group(1)), f"{X}/workflow/scripts/atlas/reintegrate.py"), unit="neighbors", kind="parameter", step=S_REI, desc="Lineage kNN graph size")
    add("lineage_resolution", "clustered at Leiden resolution 3.0", "3.0", "dec:1", lambda c: (float(c.ls["defaults"]["res"].split("_")[-1]), LSPEC + ": defaults.res", "" if all(v["res"] == c.ls["defaults"]["res"] for s in ("mouse", "human") for v in c.ls["species"][s]["lineages"].values()) else "some lineage overrides the default"), kind="parameter", step=S_REI, desc="Lineage-pass Leiden resolution")
    LP = "(178 clusters for 25 fine cell types in mouse, 179 for 17 in human;"
    add("lineage_clusters_mouse", LP, "178", "intplain", lambda c: (c.lineage("mouse")["n_clusters"], "; ".join(table_src(c, f"mouse_{ln}.csv") for ln in LINEAGES["mouse"])), unit="clusters", step=S_REI, desc="Lineage-pass clusters, mouse")
    add("fine_types_mouse", LP, "25", "intplain", lambda c: c.fine_types("mouse"), unit="cell types", step=S_REI, desc="Fine cell types (cell_subtype values), mouse")
    add("lineage_clusters_human", LP, "179", "intplain", lambda c: (c.lineage("human")["n_clusters"], "; ".join(table_src(c, f"human_{ln}.csv") for ln in LINEAGES["human"])), unit="clusters", step=S_REI, desc="Lineage-pass clusters, human")
    add("fine_types_human", LP, "17", "intplain", lambda c: c.fine_types("human"), unit="cell types", step=S_REI, desc="Fine cell types (cell_subtype values), human")
    XP = "The 13 mouse and 10 human clusters excluded"
    add("lineage_excluded_clusters_mouse", XP, "13", "intplain", lambda c: (c.lineage("mouse")["n_exclude"], "; ".join(table_src(c, f"mouse_{ln}.csv") for ln in LINEAGES["mouse"]) + " (EXCLUDE rows)"), unit="clusters", step=S_REI, desc="Lineage clusters excluded, mouse")
    add("lineage_excluded_clusters_human", XP, "10", "intplain", lambda c: (c.lineage("human")["n_exclude"], "; ".join(table_src(c, f"human_{ln}.csv") for ln in LINEAGES["human"]) + " (EXCLUDE rows)"), unit="clusters", step=S_REI, desc="Lineage clusters excluded, human")
    RP = "removed a further 1,321 mouse and 1,299 human nuclei"
    add("lineage_excluded_nuclei_mouse", RP, "1,321", "int", lambda c: (c.atlas("mouse")["retained"] - c.obj("mouse_all")["n_obs"], "atlas-retained minus released nuclei", "by subtraction; the 04c logs do not record per-lineage exclusions"), unit="nuclei", step=S_REI, desc="Nuclei removed by lineage-level exclusion, mouse")
    add("lineage_excluded_nuclei_human", RP, "1,299", "int", lambda c: (c.atlas("human")["retained"] - c.obj("human_all")["n_obs"], "atlas-retained minus released nuclei", "by subtraction; the 04c logs do not record per-lineage exclusions"), unit="nuclei", step=S_REI, desc="Nuclei removed by lineage-level exclusion, human")
    FP = "leaving 57,033 mouse and 67,847 human annotated nuclei"
    add("final_nuclei_mouse", FP, "57,033", "int", lambda c: (c.obj("mouse_all")["n_obs"], c.objpath("mouse_all")), unit="nuclei", step="15_release.smk", desc="Annotated nuclei in the release, mouse")
    add("final_nuclei_human", FP, "67,847", "int", lambda c: (c.obj("human_all")["n_obs"], c.objpath("human_all")), unit="nuclei", step="15_release.smk", desc="Annotated nuclei in the release, human")

    # ================================================================ epithelium: human layers, mouse compartments (v7)
    S_EPH = "04c_human_epithelial_layers notebook (executed run 2026_10_29b)"
    S_EPM = "04c_mouse_epithelial_compartments notebook (executed run 2026_10_29)"
    add("epi_human_resolution", "(Leiden resolution 3.0) lay predominantly", "3.0", "dec:1", lambda c: c.epi_resolution("human"), kind="parameter", step=S_EPH, desc="Leiden resolution of the labeled epithelial clusters, human")
    add("epi_human_n_neighbours", "were counted among its 20 nearest placed nuclei of the same section; the nucleus", "20", "intplain", lambda c: c.epi_neighbours("human"), unit="nuclei", kind="parameter", step=S_EPH, desc="Placed nuclei around an epithelial nucleus that decide its layer, human")
    TH = "nuclei with equal counts (221 of 8,548, 2.6%) were not used"
    add("epi_human_ties", TH, "221", "int", lambda c: c.epi_split("human", "ties"), unit="nuclei", step=S_EPH, desc="Placed nuclei in inner and outer clusters with as many decidual as fibroblast neighbours")
    add("epi_human_placed", TH, "8,548", "int", lambda c: c.epi_split("human", "placed"), unit="nuclei", step=S_EPH, desc="Placed nuclei in inner and outer clusters")
    add("epi_human_ties_pct", TH, "2.6%", "pct:1", lambda c: c.epi_split("human", "tie_frac"), unit="%", kind="statistic", step=S_EPH, desc="Share of placed nuclei with tied neighbour counts, human")
    MH = "(median 88.0% of a cluster's placed nuclei; range 52.8–99.8%)"
    add("epi_human_majority_median", MH, "88.0%", "pctnum:1", lambda c: c.epi_split("human", "majority_median"), unit="%", kind="statistic", step=S_EPH, desc="Median share of a cluster's placed nuclei in its own layer")
    add("epi_human_majority_range", MH, "52.8–99.8%", "str", lambda c: c.epi_split("human", "majority_range"), unit="%", kind="statistic", step=S_EPH, desc="Range of the share of a cluster's placed nuclei in its own layer")
    CH = "giving 24 inner clusters (5,539 nuclei) and 25 outer clusters (6,086 nuclei)"
    add("epi_human_inner_clusters", CH, "24", "intplain", lambda c: c.epi_label("human", "epithelial_inner", "clusters"), unit="clusters", step=S_EPH, desc="Clusters labeled epithelial_inner")
    add("epi_human_inner_nuclei", CH, "5,539", "int", lambda c: c.epi_label("human", "epithelial_inner", "nuclei"), unit="nuclei", step=S_EPH, desc="Nuclei labeled epithelial_inner")
    add("epi_human_outer_clusters", CH, "25", "intplain", lambda c: c.epi_label("human", "epithelial_outer", "clusters"), unit="clusters", step=S_EPH, desc="Clusters labeled epithelial_outer")
    add("epi_human_outer_nuclei", CH, "6,086", "int", lambda c: c.epi_label("human", "epithelial_outer", "nuclei"), unit="nuclei", step=S_EPH, desc="Nuclei labeled epithelial_outer")
    CIL = "Ciliated cells (one cluster, 91 nuclei) kept their own label"
    add("epi_human_ciliated_clusters", CIL, "one", "word", lambda c: c.epi_label("human", "epithelial_ciliated", "clusters"), unit="clusters", step=S_EPH, desc="Clusters labeled epithelial_ciliated")
    add("epi_human_ciliated_nuclei", CIL, "91", "intplain", lambda c: c.epi_label("human", "epithelial_ciliated", "nuclei"), unit="nuclei", step=S_EPH, desc="Nuclei labeled epithelial_ciliated")
    RH = "Three clusters (1,017 nuclei) that had been excluded because PECAM1 was detected in 18–38% of their nuclei"
    add("epi_human_reinstated_clusters", RH, "Three", "Word", lambda c: c.reinstated("clusters"), unit="clusters", step=S_EPH + " + release reinstatement table", desc="Human epithelial clusters reinstated in v7")
    add("epi_human_reinstated_nuclei", RH, "1,017", "int", lambda c: c.reinstated("nuclei"), unit="nuclei", step=S_EPH + " + release reinstatement table", desc="Nuclei of the reinstated clusters")
    add("epi_human_reinstated_pecam1", RH, "18–38%", "str", lambda c: c.reinstated_marker("PECAM1", "range"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="PECAM1 detection in the reinstated clusters (range over clusters)")
    PX = "they expressed PAX8 (74–92% of nuclei) and PAEP (49–96%)"
    add("epi_human_reinstated_pax8", PX, "74–92%", "str", lambda c: c.reinstated_marker("PAX8", "range"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="PAX8 detection in the reinstated clusters (range over clusters)")
    add("epi_human_reinstated_paep", PX, "49–96%", "str", lambda c: c.reinstated_marker("PAEP", "range"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="PAEP detection in the reinstated clusters (range over clusters)")
    EN = "(VWF in at most 2% and CDH5 in at most 1% of nuclei)"
    add("epi_human_reinstated_vwf", EN, "2%", "pct:0", lambda c: c.reinstated_marker("VWF", "max"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="Highest VWF detection among the reinstated clusters")
    add("epi_human_reinstated_cdh5", EN, "1%", "pct:0", lambda c: c.reinstated_marker("CDH5", "max"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="Highest CDH5 detection among the reinstated clusters")
    add("epi_human_reinstated_inner", "and 94.8–99.8% of their placed nuclei lay in the inner layer", "94.8–99.8%", "str", lambda c: c.reinstated("inner_range"), unit="%", kind="statistic", step=S_EPH, desc="Share of placed nuclei in the inner layer, reinstated clusters (range)")
    DH = "in both donors, 676 genes were at least two-fold higher in inner nuclei and 108 in outer nuclei"
    add("epi_human_de_inner", DH, "676", "intplain", lambda c: c.epi_de("human", "higher_in_inner"), unit="genes", step=S_EPH, desc="Genes at least two-fold higher in inner nuclei (adjusted P < 0.05) in both donors")
    add("epi_human_de_outer", DH, "108", "intplain", lambda c: c.epi_de("human", "higher_in_outer"), unit="genes", step=S_EPH, desc="Genes at least two-fold higher in outer nuclei (adjusted P < 0.05) in both donors")
    add("epi_human_de_target_sum", "on counts normalized to 10,000 per nucleus", "10,000", "int", lambda c: c.epi_parameter("human", "TARGET_SUM"), unit="counts", kind="parameter", step=S_EPH, desc="Counts per nucleus after normalization, inner vs outer test")
    add("epi_human_de_min_fold", DH, "two", "word", lambda c: c.epi_min_fold("human"), unit="fold", kind="parameter", step=S_EPH, desc="Smallest fold change of a gene called higher in inner or outer nuclei")
    add("epi_human_de_max_padj", "108 in outer nuclei (adjusted P < 0.05", "0.05", "dec:2", lambda c: c.epi_parameter("human", "MAX_PVAL_ADJ"), kind="parameter", step=S_EPH, desc="Adjusted P threshold of a gene called higher in inner or outer nuclei")
    MP2 = "MAP2 was nearly confined to inner nuclei (29–40% of inner vs 4% of outer nuclei)"
    add("epi_human_map2_inner", MP2, "29–40%", "str", lambda c: c.epi_donor_detection("MAP2", "frac_inner"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="MAP2 detection in inner nuclei (range over donors)")
    add("epi_human_map2_outer", MP2, "4%", "str", lambda c: c.epi_donor_detection("MAP2", "frac_outer"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="MAP2 detection in outer nuclei (both donors)")
    CT2 = "CTNNA2 (91–95% of outer vs 47–65% of inner nuclei)"
    add("epi_human_ctnna2_outer", CT2, "91–95%", "str", lambda c: c.epi_donor_detection("CTNNA2", "frac_outer"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="CTNNA2 detection in outer nuclei (range over donors)")
    add("epi_human_ctnna2_inner", CT2, "47–65%", "str", lambda c: c.epi_donor_detection("CTNNA2", "frac_inner"), unit="% of nuclei", kind="statistic", step=S_EPH, desc="CTNNA2 detection in inner nuclei (range over donors)")
    add("epi_mouse_n_neighbours", "among its 20 nearest placed nuclei of the same section and as glandular", "20", "intplain", lambda c: c.epi_neighbours("mouse"), unit="nuclei", kind="parameter", step=S_EPM, desc="Placed nuclei around an epithelial nucleus that decide its vote, mouse")
    TM = "nuclei with equal counts (68 of 2,306) were not used"
    add("epi_mouse_ties", TM, "68", "intplain", lambda c: c.epi_split("mouse", "ties"), unit="nuclei", step=S_EPM, desc="Placed nuclei in luminal and glandular clusters with tied neighbour counts")
    add("epi_mouse_placed", TM, "2,306", "int", lambda c: c.epi_split("mouse", "placed"), unit="nuclei", step=S_EPM, desc="Placed nuclei in luminal and glandular clusters")
    add("epi_mouse_resolution", "(Leiden resolution 3.0) took the compartment", "3.0", "dec:1", lambda c: c.epi_resolution("mouse"), kind="parameter", step=S_EPM, desc="Leiden resolution of the labeled epithelial clusters, mouse")
    MM = "(median 89.8%; range 50.7–97.2%)"
    add("epi_mouse_majority_median", MM, "89.8%", "pctnum:1", lambda c: c.epi_split("mouse", "majority_median"), unit="%", kind="statistic", step=S_EPM, desc="Median share of a cluster's votes cast for its own compartment")
    add("epi_mouse_majority_range", MM, "50.7–97.2%", "str", lambda c: c.epi_split("mouse", "majority_range"), unit="%", kind="statistic", step=S_EPM, desc="Range of the share of a cluster's votes cast for its own compartment")
    CM = "giving 13 luminal clusters (2,848 nuclei) and 12 glandular clusters (2,985 nuclei)"
    add("epi_mouse_luminal_clusters", CM, "13", "intplain", lambda c: c.epi_label("mouse", "epithelial_luminal", "clusters"), unit="clusters", step=S_EPM, desc="Clusters labeled epithelial_luminal")
    add("epi_mouse_luminal_nuclei", CM, "2,848", "int", lambda c: c.epi_label("mouse", "epithelial_luminal", "nuclei"), unit="nuclei", step=S_EPM, desc="Nuclei labeled epithelial_luminal")
    add("epi_mouse_glandular_clusters", CM, "12", "intplain", lambda c: c.epi_label("mouse", "epithelial_glandular", "clusters"), unit="clusters", step=S_EPM, desc="Clusters labeled epithelial_glandular")
    add("epi_mouse_glandular_nuclei", CM, "2,985", "int", lambda c: c.epi_label("mouse", "epithelial_glandular", "nuclei"), unit="nuclei", step=S_EPM, desc="Nuclei labeled epithelial_glandular")
    add("epi_mouse_own_pct", "87.9% of placed nuclei lay in the compartment of their cluster's label", "87.9%", "pct:1", lambda c: c.epi_split("mouse", "own_frac"), unit="%", kind="statistic", step=S_EPM, desc="Votes cast for the compartment of the nucleus's cluster label")
    PM = "Two clusters (432 nuclei, 389 from induced uteri) were labeled proliferative epithelium"
    add("epi_mouse_prolif_clusters", PM, "Two", "Word", lambda c: c.epi_label("mouse", "epithelial_proliferative", "clusters"), unit="clusters", step=S_EPM, desc="Clusters labeled epithelial_proliferative")
    add("epi_mouse_prolif_nuclei", PM, "432", "intplain", lambda c: c.epi_label("mouse", "epithelial_proliferative", "nuclei"), unit="nuclei", step=S_EPM, desc="Nuclei labeled epithelial_proliferative")
    add("epi_mouse_prolif_induced", PM, "389", "intplain", lambda c: c.epi_mouse_stratum("epithelial_proliferative", "induced"), unit="nuclei", step=S_EPM, desc="Proliferative nuclei from induced uteri")
    CC = "63.7% and 72.2% of their nuclei were in S or G2/M phase"
    add("epi_mouse_prolif_cycling_low", CC, "63.7%", "pctnum:1", lambda c: c.epi_cycling("proliferative_low"), unit="% of nuclei", kind="statistic", step=S_EPM, desc="Share of nuclei in S or G2/M, proliferative cluster with the lower share")
    add("epi_mouse_prolif_cycling_high", CC, "72.2%", "pctnum:1", lambda c: c.epi_cycling("proliferative_high"), unit="% of nuclei", kind="statistic", step=S_EPM, desc="Share of nuclei in S or G2/M, proliferative cluster with the higher share")
    add("epi_mouse_lineage_nuclei", "on all 6,681 nuclei of the epithelial lineage", "6,681", "int", lambda c: c.epi_lineage_nuclei("mouse"), unit="nuclei", step=S_EPM, desc="Nuclei of the mouse epithelial lineage object (all clusters)")
    add("epi_mouse_split_cycling_max", "against at most 43.9% in any luminal or glandular cluster", "43.9%", "pctnum:1", lambda c: c.epi_cycling("split_max"), unit="% of nuclei", kind="statistic", step=S_EPM, desc="Highest share of nuclei in S or G2/M among luminal and glandular clusters")
    PG = "and 1,817 genes were at least two-fold higher in their 389 induced nuclei than in the 1,005 other epithelial nuclei of induced uteri"
    add("epi_mouse_de_prolif", PG, "1,817", "int", lambda c: c.epi_de("mouse_proliferative", "higher_in_proliferative"), unit="genes", step=S_EPM, desc="Genes at least two-fold higher in proliferative than other induced epithelial nuclei")
    add("epi_mouse_prolif_induced_2", PG, "389", "intplain", lambda c: c.epi_mouse_stratum("epithelial_proliferative", "induced"), unit="nuclei", step=S_EPM, desc="Proliferative nuclei from induced uteri (second mention)")
    add("epi_mouse_other_induced", PG, "1,005", "int", lambda c: c.epi_mouse_stratum(("epithelial_luminal", "epithelial_glandular"), "induced"), unit="nuclei", step=S_EPM, desc="Luminal and glandular nuclei of induced uteri")
    add("epi_mouse_de_prolif_min_fold", PG, "two", "word", lambda c: c.epi_min_fold("mouse"), unit="fold", kind="parameter", step=S_EPM, desc="Smallest fold change of a gene called higher in proliferative nuclei")
    add("epi_mouse_de_prolif_max_padj", "other epithelial nuclei of induced uteri (adjusted P < 0.05)", "0.05", "dec:2", lambda c: c.epi_parameter("mouse", "MAX_PVAL_ADJ"), kind="parameter", step=S_EPM, desc="Adjusted P threshold of a gene called higher in proliferative nuclei")
    SM = "separately within uninduced uteri (2,565 luminal and 2,263 glandular nuclei) and induced uteri (283 and 722)"
    add("epi_mouse_luminal_uninduced", SM, "2,565", "int", lambda c: c.epi_mouse_stratum("epithelial_luminal", "uninduced"), unit="nuclei", step=S_EPM, desc="Luminal nuclei, uninduced uteri")
    add("epi_mouse_glandular_uninduced", SM, "2,263", "int", lambda c: c.epi_mouse_stratum("epithelial_glandular", "uninduced"), unit="nuclei", step=S_EPM, desc="Glandular nuclei, uninduced uteri")
    add("epi_mouse_luminal_induced", SM, "283", "intplain", lambda c: c.epi_mouse_stratum("epithelial_luminal", "induced"), unit="nuclei", step=S_EPM, desc="Luminal nuclei, induced uteri")
    add("epi_mouse_glandular_induced", SM, "722", "intplain", lambda c: c.epi_mouse_stratum("epithelial_glandular", "induced"), unit="nuclei", step=S_EPM, desc="Glandular nuclei, induced uteri")
    DM = "353 genes were at least two-fold higher in luminal and 783 in glandular nuclei in both"
    add("epi_mouse_de_luminal", DM, "353", "intplain", lambda c: c.epi_de("mouse", "higher_in_luminal"), unit="genes", step=S_EPM, desc="Genes at least two-fold higher in luminal nuclei in both strata")
    add("epi_mouse_de_glandular", DM, "783", "intplain", lambda c: c.epi_de("mouse", "higher_in_glandular"), unit="genes", step=S_EPM, desc="Genes at least two-fold higher in glandular nuclei in both strata")
    add("epi_mouse_de_min_fold", DM, "two", "word", lambda c: c.epi_min_fold("mouse"), unit="fold", kind="parameter", step=S_EPM, desc="Smallest fold change of a gene called higher in luminal or glandular nuclei")
    add("epi_mouse_de_max_padj", "783 in glandular nuclei in both (adjusted P < 0.05", "0.05", "dec:2", lambda c: c.epi_parameter("mouse", "MAX_PVAL_ADJ"), kind="parameter", step=S_EPM, desc="Adjusted P threshold of a gene called higher in luminal or glandular nuclei")
    add("outlier_rule_max_same", "(at most 20% of 30 nearest neighbors sharing", "20%", "pct:0", lambda c: c.outlier_rule("same"), unit="%", kind="parameter", step="embedding-outlier flags (2026-10-22)", desc="Outlier rule: maximum share of neighbors with the same lineage")
    add("outlier_rule_k", "(at most 20% of 30 nearest neighbors sharing", "30", "intplain", lambda c: c.outlier_rule("k"), unit="neighbors", kind="parameter", step="embedding-outlier flags (2026-10-22)", desc="Outlier rule: neighbors examined")
    add("outlier_rule_min_other", "at least 70% from one other lineage", "70%", "pct:0", lambda c: c.outlier_rule("other"), unit="%", kind="parameter", step="embedding-outlier flags (2026-10-22)", desc="Outlier rule: minimum share from one other lineage")
    # The all-cells objects carry the flags of the release embeddings (human re-fitted for v7). The fibroblast/decidual
    # objects keep their own flags: the human one those of the atlas fitted before the 1,017 epithelial nuclei were
    # retained. The two kinds of flag have separate keys.
    FLAG_SOURCES = {"n_obs": "", "n_embedding_outlier": " obs['embedding_outlier']",
                    "n_embedding_outlier_on_fibdec": " obs['embedding_outlier'] of the nuclei in the fibroblast/decidual object",
                    "n_embedding_outlier_shared_with_all": " obs['embedding_outlier'] also flagged in the all-cells object"}
    OM = "34 of 57,033 mouse nuclei (21 of 43,197 fibroblast/decidual)"
    OH = "25 of 67,847 human nuclei (9 of 47,533 fibroblast/decidual)"
    SH = "(3 of these 7 are among the 9)"
    for key, ph, shown, tag, fld, d in (("outliers_mouse_all", OM, "34", "mouse_all", "n_embedding_outlier", "Embedding outliers, mouse all nuclei"),
                                        ("final_nuclei_mouse_2", OM, "57,033", "mouse_all", "n_obs", "Annotated nuclei, mouse (second mention)"),
                                        ("outliers_mouse_fibdec", OM, "21", "mouse_fibdec", "n_embedding_outlier", "Embedding outliers stored in the mouse fibroblast/decidual object"),
                                        ("outliers_mouse_all_on_fibdec", OM, "21", "mouse_all", "n_embedding_outlier_on_fibdec", "Embedding outliers of the mouse all-cells object among fibroblast/decidual nuclei"),
                                        ("fibdec_nuclei_mouse", OM, "43,197", "mouse_fibdec", "n_obs", "Fibroblast/decidual nuclei, mouse"),
                                        ("outliers_human_all", OH, "25", "human_all", "n_embedding_outlier", "Embedding outliers, human all nuclei (re-fitted v7 embedding)"),
                                        ("final_nuclei_human_2", OH, "67,847", "human_all", "n_obs", "Annotated nuclei, human (second mention)"),
                                        ("outliers_human_all_on_fibdec", OH, "9", "human_all", "n_embedding_outlier_on_fibdec", "Embedding outliers of the human all-cells object among fibroblast/decidual nuclei"),
                                        ("fibdec_nuclei_human", OH, "47,533", "human_fibdec", "n_obs", "Fibroblast/decidual nuclei, human"),
                                        ("outliers_human_fibdec", "exclude the 7 human fibroblast/decidual nuclei flagged on the earlier human atlas", "7", "human_fibdec", "n_embedding_outlier", "Embedding outliers stored in the human fibroblast/decidual object (flags of the atlas before the epithelial reinstatement)"),
                                        ("outliers_human_fibdec_shared", SH, "3", "human_fibdec", "n_embedding_outlier_shared_with_all", "Fibroblast/decidual-object flags also flagged in the human all-cells object"),
                                        ("outliers_human_fibdec_2", SH, "7", "human_fibdec", "n_embedding_outlier", "Embedding outliers stored in the human fibroblast/decidual object (second mention)"),
                                        ("outliers_human_all_on_fibdec_2", SH, "9", "human_all", "n_embedding_outlier_on_fibdec", "All-cells flags among human fibroblast/decidual nuclei (second mention)")):
        add(key, ph, shown, "int", (lambda t, f: lambda c: (c.obj(t)[f], c.objpath(t) + FLAG_SOURCES[f]))(tag, fld), unit="nuclei", step="15_release.smk + embedding-outlier flags", desc=d)
    add("reinstated_nuclei_human", "fitted before the 1,017 epithelial nuclei were retained", "1,017", "int", lambda c: c.reinstated("nuclei"), unit="nuclei", step="04g reinstate_release_cells (release reinstatement table)", desc="Human epithelial nuclei reinstated in v7 (second mention)")
    add("samap_source_human", "use all 47,533 human fibroblast/decidual nuclei", "47,533", "int", lambda c: c.samap_source_human(), unit="nuclei", step="21_samap_raw.smk (joint SAMap graph)", desc="Human fibroblast/decidual nuclei in the SAMap graph, flagged ones included (SAMap and Ulrich source set)")

    # ================================================================ spatial QC
    add("spatial_k", "its five nearest positioned neighbors", "five", "word", cfgv("spatial_qc", "k"), unit="neighbors", kind="parameter", step=S_SPA, desc="Neighbors for the spatial distance criteria")
    add("spatial_floor_um", "was below 10 µm", "10", "intplain", cfgv("spatial_qc", "floor_um"), unit="µm", kind="parameter", step=S_SPA, desc="Minimum mean distance to the nearest neighbors")
    add("spatial_k_2", "the mean distance to the five nearest neighbors was recomputed", "five", "word", cfgv("spatial_qc", "k"), unit="neighbors", kind="parameter", step=S_SPA, desc="Neighbors (second mention)")
    add("spatial_n_mads", "by more than 5 median absolute deviations", "5", "intplain", cfgv("spatial_qc", "n_mads"), unit="MADs", kind="parameter", step=S_SPA, desc="MAD cutoff for the far criterion")
    add("spatial_n_boxes", "with ten rectangular exclusion regions", "ten", "word", lambda c: (len(c.cfg["spatial_qc"]["exclusions"]["mouse"]), CFG + ": spatial_qc.exclusions.mouse"), unit="regions", kind="parameter", step=S_SPA, desc="Hand-drawn exclusion regions")
    # Timepoints by their release names (vocabulary v5); keys keep the names they had before the renaming.
    BXP = "(5 in the uninduced sample, 3 at 3.5 dpi, and 1 each at 1.5 dpi early and 2.5 dpi;"
    for key, tp, shown, nth in (("spatial_boxes_uninduced", "uninduced", "5", 0), ("spatial_boxes_3dpi", "3.5_dpi", "3", 0),
                                ("spatial_boxes_1dpi_early", "1.5_dpi_early", "1", 0), ("spatial_boxes_2dpi", "2.5_dpi", "1", 0)):
        add(key, BXP, shown, "intplain", (lambda t: lambda c: c.exclusion_regions(t))(tp),
            unit="regions", kind="parameter", step=S_SPA, desc=f"Exclusion regions at {tp}", nth=nth)
    SP3 = "the three steps blanked 1,272, 544 and 136 positions"
    add("spatial_blanked_step1_mouse", SP3, "1,272", "int", lambda c: c.spatial_mouse("step1"), unit="positions", step=S_SPA, desc="Positions blanked as duplicate or too close, mouse")
    add("spatial_blanked_step2_mouse", SP3, "544", "int", lambda c: c.spatial_mouse("step2"), unit="positions", step=S_SPA, desc="Positions blanked as far from the tissue, mouse")
    add("spatial_blanked_step3_mouse", SP3, "136", "int", lambda c: c.spatial_mouse("step3"), unit="positions", step=S_SPA, desc="Positions blanked by the exclusion regions, mouse")
    BCP = "(the exclusion regions covered 112 nuclei in the uninduced sample, 12 at 3.5 dpi, 7 at 1.5 dpi early and 5 at 2.5 dpi)"
    for key, tp, shown in (("spatial_boxed_uninduced", "uninduced", "112"), ("spatial_boxed_3dpi", "3.5_dpi", "12"), ("spatial_boxed_1dpi_early", "1.5_dpi_early", "7"), ("spatial_boxed_2dpi", "2.5_dpi", "5")):
        add(key, BCP, shown, "int", (lambda t: lambda c: c.spatial_mouse("step3", t))(tp), unit="nuclei", step=S_SPA, desc=f"Nuclei blanked by exclusion regions at {tp}")
    PP = "leaving 29,678 of 57,033 nuclei with a position"
    add("spatial_positioned_mouse", PP, "29,678", "int", lambda c: (c.obj("mouse_all")["n_positioned"], c.objpath("mouse_all") + " obs['x_um_qc']"), unit="nuclei", step=S_SPA, desc="Mouse nuclei with a QC'd position")
    add("final_nuclei_mouse_3", PP, "57,033", "int", lambda c: (c.obj("mouse_all")["n_obs"], c.objpath("mouse_all")), unit="nuclei", step="15_release.smk", desc="Annotated nuclei, mouse (third mention)")
    add("spatial_unmapped_mouse", "A further 25,403 had no position from Slide-tags", "25,403", "int", lambda c: c.spatial_mouse("unmapped"), unit="nuclei", step=S_SPA, desc="Mouse nuclei without a Slide-tags position")
    add("spatial_blanked_human", "(4 positions blanked)", "4", "intplain", lambda c: c.spatial_human("blanked"), unit="positions", step=S_SPA, desc="Human positions blanked (step 1 only)")
    HP = "In total, 46,470 of 67,847 human nuclei have a position, and 21,373 had none from Slide-tags"
    add("spatial_positioned_human", HP, "46,470", "int", lambda c: (c.obj("human_all")["n_positioned"], c.objpath("human_all") + " obs['x_um_qc']"), unit="nuclei", step=S_SPA, desc="Human nuclei with a position")
    add("final_nuclei_human_3", HP, "67,847", "int", lambda c: (c.obj("human_all")["n_obs"], c.objpath("human_all")), unit="nuclei", step="15_release.smk", desc="Annotated nuclei, human (third mention)")
    add("spatial_unmapped_human", HP, "21,373", "int", lambda c: c.spatial_human("unmapped"), unit="nuclei", step=S_SPA, desc="Human nuclei without a Slide-tags position")

    # ================================================================ pseudotime
    add("apoptotic_decidual_mouse", "In mouse, the 907 apoptotic decidual nuclei were excluded", "907", "int", lambda c: (c.subtype_counts("mouse_fibdec")["decAp"], c.objpath("mouse_fibdec") + " obs['cell_subtype']"), unit="nuclei", step=S_PT, desc="Apoptotic decidual nuclei (no pseudotime)")
    add("pt_days_mouse", "(mouse: uninduced 0, 1 dpi early 1, 1 dpi late 1.5, 2 dpi 2, 3 dpi 3, 4 dpi 4;", "uninduced 0, 1 dpi early 1, 1 dpi late 1.5, 2 dpi 2, 3 dpi 3, 4 dpi 4", "str",
        lambda c: (", ".join(f"{k.replace('_bleeding', '').replace('_', ' ')} {v:g}" for k, v in c.day_map("mouse").items()), f"{O22}/provenance/harmonized_pseudotime.py",
                   "names as in the pseudotime script; the release names these timepoints " + ", ".join(c.release_day_map("mouse"))),
        kind="parameter", step=S_PT, desc="Timepoint encoding in days, mouse")
    add("pt_days_human", "human: secretory 0, menstrual 1)", "secretory 0, menstrual 1", "str",
        lambda c: (", ".join(f"{k} {v:g}" for k, v in c.day_map("human").items()), f"{O22}/provenance/harmonized_pseudotime.py"), kind="parameter", step=S_PT, desc="Timepoint encoding, human")
    add("moscot_version", "moscot (v0.5.0;", "0.5.0", "str", lambda c: (c.F["pt_provenance"]["mouse"]["versions"]["moscot"], f"{O22}/provenance/mouse_pca_released_provenance.json", "" if c.F["pt_provenance"]["human"]["versions"]["moscot"] == c.F["pt_provenance"]["mouse"]["versions"]["moscot"] else "human run differs"), kind="version", step=S_PT, desc="moscot version used for the pseudotime of record")
    add("pt_latent_dim", "each species' 30-dimensional stromal scVI latent space", "30", "intplain", lambda c: c.latent_dim(), unit="dimensions", kind="parameter", step=S_PT, desc="Dimension of the stromal scVI latent used for the OT cost")
    add("ot_epsilon", "ε = 10−3", "10−3", "pow10", lambda c: (c.F["pt_provenance"]["mouse"]["params"]["epsilon"], f"{O22}/provenance/*_pca_released_provenance.json"), kind="parameter", step=S_PT, desc="Entropic regularization")
    add("ot_tau", "τa = τb = 0.95", "0.95", "dec:2", lambda c: (c.F["pt_provenance"]["mouse"]["params"]["tau"], f"{O22}/provenance/*_pca_released_provenance.json"), kind="parameter", step=S_PT, desc="Unbalancedness (tau_a = tau_b)")
    add("cellrank_version", "RealTimeKernel (v2.0.7;", "2.0.7", "str", lambda c: (c.F["pt_provenance"]["mouse"]["versions"]["cellrank"], f"{O22}/provenance/mouse_pca_released_provenance.json"), kind="version", step=S_PT, desc="CellRank version")
    add("conn_weight", "conn_weight = 0.2", "0.2", "dec:1", lambda c: (c.F["pt_provenance"]["mouse"]["params"]["conn_weight"], f"{O22}/provenance/*_pca_released_provenance.json"), kind="parameter", step=S_PT, desc="RealTimeKernel connectivity weight")
    add("rtk_coupling_share", "assigns 80% of its probability", "80%", "pct:0", lambda c: (1 - c.F["pt_provenance"]["mouse"]["params"]["conn_weight"], "1 - conn_weight"), unit="%", kind="parameter", step=S_PT, desc="Row share from the OT coupling")
    add("rtk_knn_share", "and 20% to a within-timepoint", "20%", "pct:0", lambda c: (c.F["pt_provenance"]["mouse"]["params"]["conn_weight"], "conn_weight"), unit="%", kind="parameter", step=S_PT, desc="Row share from the same-timepoint kNN graph")
    add("pt_knn", "15-nearest-neighbor graph (Scanpy defaults:", "15", "intplain", lambda c: (c.F["pt_provenance"]["mouse"]["params"]["n_neighbors"], f"{O22}/provenance/*_pca_released_provenance.json"), unit="neighbors", kind="parameter", step=S_PT, desc="Same-timepoint kNN graph size")
    add("pt_n_pcs", "Scanpy defaults: 50-component PCA", "50", "intplain", lambda c: c.pt_pcs(), unit="components", kind="parameter", step=S_PT, desc="PCA components for the same-timepoint graph")
    add("pt_n_hvg", "of the 3,000 stromal highly variable genes", "3,000", "int", lambda c: (c.obj("mouse_fibdec").get("n_highly_variable"), c.objpath("mouse_fibdec") + " var['highly_variable']", f"human object: {c.obj('human_fibdec').get('n_highly_variable')}"), unit="genes", kind="parameter", step=S_PT, desc="Stromal highly variable genes")
    add("pt_n_dcs", "computed a 15-component diffusion map", "15", "intplain", lambda c: (c.F["pt_provenance"]["mouse"]["params"]["n_dcs"], f"{O22}/provenance/*_pca_released_provenance.json"), unit="components", kind="parameter", step=S_PT, desc="Diffusion components")
    add("pt_dc1_rho_mouse", "(Spearman ρ = 0.89 mouse, 0.72 human)", "0.89", "dec:2", lambda c: c.dc_rho("mouse"), unit="Spearman rho", kind="statistic", step=S_PT, desc="|rho| of the trajectory component with the fibroblast-to-decidual order, mouse")
    add("pt_dc1_rho_human", "(Spearman ρ = 0.89 mouse, 0.72 human)", "0.72", "dec:2", lambda c: c.dc_rho("human"), unit="Spearman rho", kind="statistic", step=S_PT, desc="|rho| of the trajectory component with the fibroblast-to-decidual order, human")
    add("pt_root_label_mouse", "an uninduced eFib1 nucleus in mouse", "eFib1", "str", lambda c: c.root_label("mouse"), kind="label", step=S_PT, desc="Cell subtype of the mouse root")
    add("pt_root_mouse", "(ACCCTTGAGCATGAAT-1-uninduced_mck8b)", "ACCCTTGAGCATGAAT-1-uninduced_mck8b", "str", lambda c: (c.F["pt_provenance"]["mouse"]["root"], f"{O22}/provenance/mouse_pca_released_provenance.json"), kind="identifier", step=S_PT, desc="Mouse root nucleus")
    add("pt_root_human", "(GAGACCCAGAATAGTC-1-secretory_mck6b)", "GAGACCCAGAATAGTC-1-secretory_mck6b", "str", lambda c: (c.F["pt_provenance"]["human"]["root"], f"{O22}/provenance/human_pca_released_provenance.json"), kind="identifier", step=S_PT, desc="Human root nucleus")
    add("pt_n_dcs_2", "computed with scanpy.tl.dpt on all 15 components", "15", "intplain", lambda c: (c.F["pt_provenance"]["mouse"]["params"]["n_dcs"], f"{O22}/provenance/*_pca_released_provenance.json"), unit="components", kind="parameter", step=S_PT, desc="Diffusion components for DPT (second mention)")
    MP = "Mouse pseudotime increased with subtype order (ρ = 0.86) and timepoint (ρ = 0.62)"
    add("pt_rho_subtype_mouse", MP, "0.86", "dec:2", lambda c: c.pt_val("mouse", "rho_subtype"), unit="Spearman rho", kind="statistic", step=S_PT, desc="Spearman rho, pseudotime vs subtype order, mouse")
    add("pt_rho_time_mouse", MP, "0.62", "dec:2", lambda c: c.pt_val("mouse", "rho_time"), unit="Spearman rho", kind="statistic", step=S_PT, desc="Spearman rho, pseudotime vs timepoint, mouse")
    HMP = "(median 0.28 vs 0.75; ρ = 0.70)"
    add("pt_median_fib_human", HMP, "0.28", "dec:2", lambda c: c.pt_val("human", "med_fib"), unit="rank-normalized pseudotime", kind="statistic", step=S_PT, desc="Median pseudotime of human fibroblasts")
    add("pt_median_dec_human", HMP, "0.75", "dec:2", lambda c: c.pt_val("human", "med_dec"), unit="rank-normalized pseudotime", kind="statistic", step=S_PT, desc="Median pseudotime of human decidual cells")
    add("pt_rho_type_human", HMP, "0.70", "dec:2", lambda c: c.pt_val("human", "rho_type"), unit="Spearman rho", kind="statistic", step=S_PT, desc="Spearman rho, pseudotime vs cell type (fibroblast < decidual), human")
    add("pt_rho_time_human", "menstrual phase (ρ = 0.34)", "0.34", "dec:2", lambda c: c.pt_val("human", "rho_time"), unit="Spearman rho", kind="statistic", step=S_PT, desc="Spearman rho, pseudotime vs phase, human")

    # ================================================================ mutual information
    add("mi_min_detect", "detected in at least 1% of nuclei were tested", "1%", "pct:0", lambda c: (c.F["mi"]["pt_harmonized_mouse_nall"]["manifest"]["min_detect_frac"], f"{X}/results/mi_dpt/pt_harmonized_mouse_nall/manifest.json"), unit="%", kind="parameter", step=S_MI, desc="Detection floor for MI genes")
    TP_ = "were tested (14,261 mouse and 15,008 human genes)"
    add("mi_tested_mouse", TP_, "14,261", "int", lambda c: (c.mi()["m"]["tested"], mi_file("pt_harmonized_mouse_nall")), unit="genes", step=S_MI, desc="Genes tested for MI, mouse")
    add("mi_tested_human", TP_, "15,008", "int", lambda c: (c.mi()["h"]["tested"], mi_file("pt_harmonized_human_nall")), unit="genes", step=S_MI, desc="Genes tested for MI, human")
    add("sklearn_version", "(scikit-learn 1.8.0,", "1.8.0", "str", lambda c: c.ver("rsc_uv", "scikit-learn"), kind="version", step=S_MI + " (config mi_pseudotime.python = rsc_uv)", desc="scikit-learn version")
    add("mi_perms", "permuting pseudotime 1,000 times", "1,000", "int", lambda c: (int(c.F["mi"]["pt_harmonized_mouse_nall"]["manifest"]["perms"]) - 1, f"{X}/results/mi_dpt/pt_harmonized_*_nall/manifest.json (perms = observed + null)"), unit="permutations", kind="parameter", step=S_MI, desc="Null permutations")
    add("mi_chunks", "processed in 200 chunks", "200", "intplain", lambda c: (int(c.F["mi"]["pt_harmonized_mouse_nall"]["manifest"]["n_chunks"]), f"{X}/results/mi_dpt/pt_harmonized_mouse_nall/manifest.json"), unit="chunks", kind="parameter", step=S_MI, desc="Gene chunks")
    add("mi_perms_2", "shared the same 1,000 permutations", "1,000", "int", lambda c: (int(c.F["mi"]["pt_harmonized_mouse_nall"]["manifest"]["perms"]) - 1, f"{X}/results/mi_dpt/pt_harmonized_*_nall/manifest.json"), unit="permutations", kind="parameter", step=S_MI, desc="Null permutations (second mention)")
    add("mi_perms_3", "1)/(1,000 + 1)", "1,000", "int", lambda c: (int(c.F["mi"]["pt_harmonized_mouse_nall"]["manifest"]["perms"]) - 1, f"{X}/results/mi_dpt/pt_harmonized_*_nall/manifest.json"), unit="permutations", kind="parameter", step=S_MI, desc="Null permutations (third mention)")
    add("statsmodels_version", "statsmodels 0.14.6)", "0.14.6", "str", lambda c: c.ver("rsc_uv", "statsmodels"), kind="version", step=S_MI, desc="statsmodels version")
    add("mi_q_thresh", "and genes with q < 0.05 were considered significant", "0.05", "dec:2", cfgv("mi", "q_thresh"), kind="parameter", step=S_MI, desc="MI significance threshold (BH q)")
    SG = "considered significant (14,261 mouse and 14,921 human genes)"
    add("mi_sig_mouse", SG, "14,261", "int", lambda c: (c.mi()["m"]["sig"], mi_file("pt_harmonized_mouse_nall"), c.mi_note("m", "sig", lambda v: f"{v:,}")), unit="genes", step=S_MI, desc="Genes with MI q < 0.05, mouse")
    add("mi_sig_human", SG, "14,921", "int", lambda c: (c.mi()["h"]["sig"], mi_file("pt_harmonized_human_nall"), c.mi_note("h", "sig", lambda v: f"{v:,}")), unit="genes", step=S_MI, desc="Genes with MI q < 0.05, human")
    add("mi_top_percentile", "exceeded the 90th percentile among significant genes", "90", "intplain", cfgv("mi", "top_percentile"), kind="parameter", step=S_OV, desc="Percentile defining the top MI genes")
    add("mi_top_mouse", "1,426 mouse genes (MI > 0.070 nats)", "1,426", "int", lambda c: (c.mi()["m"]["top"], mi_file("pt_harmonized_mouse_nall")), unit="genes", step=S_OV, desc="Top MI genes, mouse")
    add("mi_top_thr_mouse", "1,426 mouse genes (MI > 0.070 nats)", "0.070", "dec:3", lambda c: (c.mi()["m"]["thr"], mi_file("pt_harmonized_mouse_nall"), c.mi_note("m", "thr", lambda v: f"{v:.3f}")), unit="nats", kind="statistic", step=S_OV, desc="MI cutoff for the top set, mouse")
    add("mi_top_human", "1,492 human genes (MI > 0.028 nats)", "1,492", "int", lambda c: (c.mi()["h"]["top"], mi_file("pt_harmonized_human_nall")), unit="genes", step=S_OV, desc="Top MI genes, human")
    add("mi_top_thr_human", "1,492 human genes (MI > 0.028 nats)", "0.028", "dec:3", lambda c: (c.mi()["h"]["thr"], mi_file("pt_harmonized_human_nall"), c.mi_note("h", "thr", lambda v: f"{v:.3f}")), unit="nats", kind="statistic", step=S_OV, desc="MI cutoff for the top set, human")
    add("ensembl_release_orth", "Ensembl Compara release 116", "116", "intplain", lambda c: (int(re.search(r"(\d+)", c.F["mi"]["overlap_report"]["ortholog_release"][0]).group(1)), c.F["mi"]["overlap_report"]["ortholog_source"]), kind="reference", step=S_OV, desc="Ensembl release of the ortholog table")
    add("n_one2one_pairs", "ortholog_one2one; 16,275 gene pairs", "16,275", "int", lambda c: (c.F["mi"]["orth"]["n_pairs"], c.F["mi"]["orth"]["file"]), unit="pairs", step=S_OV, desc="Strict one-to-one human-mouse pairs")
    ovsrc = mi_file("pt_harmonized_mouse_nall") + " + " + mi_file("pt_harmonized_human_nall") + " + ortholog_pairs_one2one.csv"
    def ov(field, how):
        return lambda c: (c.mi()["ov"][field], ovsrc, c.mi_note("ov", field, how))
    add("mi_overlap_N", "comprised the 9,606 one-to-one pairs", "9,606", "int", ov("N", lambda v: f"{v:,}"), unit="pairs", kind="statistic", step=S_OV, desc="Background: 1:1 pairs with both genes MI-significant")
    add("mi_q_thresh_2", "significant for MI (q < 0.05) in their species", "0.05", "dec:2", cfgv("mi", "q_thresh"), kind="parameter", step=S_OV, desc="MI significance threshold (second mention)")
    KP = "(1,132 mouse genes, K; 1,191 human genes, n)"
    add("mi_overlap_K", KP, "1,132", "int", ov("K", lambda v: f"{v:,}"), unit="genes", kind="statistic", step=S_OV, desc="Mouse top genes in the background (K)")
    add("mi_overlap_n", KP, "1,191", "int", ov("n", lambda v: f"{v:,}"), unit="genes", kind="statistic", step=S_OV, desc="Human top genes in the background (n)")
    add("mi_overlap_k", "(k = 371 pairs)", "371", "int", ov("k", lambda v: f"{v:,}"), unit="pairs", kind="statistic", step=S_OV, desc="Pairs in both top sets (k)")
    add("mi_overlap_N_2", "N = 9,606, computed on the log scale", "9,606", "int", ov("N", lambda v: f"{v:,}"), unit="pairs", kind="statistic", step=S_OV, desc="Background size (second mention)")
    EXP = "K × n / N = 140.4 genes (12% of the human top set)"
    add("mi_overlap_expected", EXP, "140.4", "dec:1", ov("expected", lambda v: f"{v:.1f}"), unit="pairs", kind="statistic", step=S_OV, desc="Expected overlap K*n/N")
    add("mi_overlap_pct_expected", EXP, "12%", "pct:0", ov("frac_exp", lambda v: f"{100 * v:.0f}%"), unit="% of human top set", kind="statistic", step=S_OV, desc="Expected overlap as a share of the human top set")
    add("mi_overlap_pct_observed", "the observed overlap 31%", "31%", "pct:0", ov("frac_obs", lambda v: f"{100 * v:.0f}%"), unit="% of human top set", kind="statistic", step=S_OV, desc="Observed overlap as a share of the human top set")
    FE = "(2.64-fold enrichment; p = 1.0 × 10−84)"
    add("mi_overlap_fold", FE, "2.64", "dec:2", ov("fold", lambda v: f"{v:.2f}"), unit="fold", kind="statistic", step=S_OV, desc="Fold enrichment k/expected")
    add("mi_overlap_p", FE, "1.0 × 10−84", "log10sci:2", lambda c: (c.mi()["ov"]["log10p"], ovsrc, c.p_note()), unit="p", kind="statistic", step=S_OV, desc="Hypergeometric p (upper tail)")
    SWP = "repeated at every percentile cutoff from the 40th to the 99th"
    add("mi_sweep_min_pct", SWP, "40", "intplain", lambda c: (round(100 - 100 * c.F["mi"]["sweep"]["frac_max"]), c.F["mi"]["sweep"]["file"], f"sweep covers top fractions {c.F['mi']['sweep']['frac_min']:g}-{c.F['mi']['sweep']['frac_max']:g}"), unit="percentile", kind="parameter", step=S_OV, desc="Lowest percentile cutoff in the sweep")
    add("mi_sweep_max_pct", SWP, "99", "intplain", lambda c: (round(100 - 100 * c.F["mi"]["sweep"]["frac_min"]), c.F["mi"]["sweep"]["file"]), unit="percentile", kind="parameter", step=S_OV, desc="Highest percentile cutoff in the sweep")
    add("mi_sweep_null_lo", "the 2.5th-97.5th percentile range", "2.5", "dec:1", lambda c: (None, c.F["mi"]["sweep"]["file"], f"sweep table columns {c.F['mi']['sweep']['cols']} hold no null range"), unit="percentile", kind="parameter", step=S_OV, desc="Lower null percentile shown in the sweep")
    add("mi_sweep_null_hi", "the 2.5th-97.5th percentile range", "97.5", "dec:1", lambda c: (None, c.F["mi"]["sweep"]["file"], "see mi_sweep_null_lo"), unit="percentile", kind="parameter", step=S_OV, desc="Upper null percentile shown in the sweep")

    # ================================================================ pseudotime heatmaps
    nohm = lambda c: (None, "; ".join(c.F["misc_find"].get("heatmap_code") or []) or "no Monocle-style heatmap code found", "Monocle 2 reimplementation not located in the repository or figures/ (METHODS_pseudotime_heatmap.md, 2026-10-14)")
    add("heatmap_panel_genes", "Smoothed profiles of the 28 panel genes", "28", "intplain", nohm, unit="genes", kind="parameter", step="pseudotime heatmaps (Monocle 2 reimplementation)", desc="Yang et al. panel genes")
    add("heatmap_spline_df", "natural cubic spline of pseudotime (3 degrees of freedom)", "3", "intplain", nohm, unit="df", kind="parameter", step="pseudotime heatmaps", desc="Spline degrees of freedom")
    add("heatmap_genes_dispersion", "fitted across all 27,064 genes detected in more than one nucleus", "27,064", "int", nohm, unit="genes", kind="statistic", step="pseudotime heatmaps", desc="Genes in the dispersion-mean fit")
    add("heatmap_points", "Fitted values at 100 evenly spaced points", "100", "intplain", nohm, unit="points", kind="parameter", step="pseudotime heatmaps", desc="Points along pseudotime")
    HS = "candidate genes were the strict one-to-one orthologs ranked in the top 10% by MI"
    add("xs_heatmap_topfrac", HS, "10%", "pct:0", lambda c: (c.F["mi"]["overlap_report"]["topfrac"], f"{X}/results/mi_dpt/overlap_harmonized_nall/overlap_report.json"), unit="%", kind="parameter", step="crossspecies_gene_selection.py", desc="Top MI fraction for cross-species heatmap candidates")
    add("xs_heatmap_bins", "(mean log-normalized expression in 50 equal-count pseudotime bins", "50", "intplain", lambda c: c.sel_param(r"(?:nbin|n_bins|NBIN)\s*=\s*(\d+)", int), unit="bins", kind="parameter", step="crossspecies_gene_selection.py", desc="Bins for the cross-species profiles")
    add("xs_heatmap_rmin", "Among genes with r > 0.8", "0.8", "dec:1", lambda c: (c.F["mi"]["heatmap_selection"]["rmin"], f"{X}/results/mi_dpt/overlap_harmonized_nall/crossspecies_heatmap_selection/selection_summary.json"), kind="parameter", step="crossspecies_gene_selection.py", desc="Minimum human-mouse profile correlation")
    add("xs_heatmap_early", "(< 0.4 or > 0.6)", "0.4", "dec:1", lambda c: c.sel_param(r"(?:EARLY|early_max|early_cut|early)\s*=\s*(0\.\d+)", float), kind="parameter", step="crossspecies_gene_selection.py", desc="Centroid cutoff for early genes")
    add("xs_heatmap_late", "(< 0.4 or > 0.6)", "0.6", "dec:1", lambda c: c.sel_param(r"(?:LATE|late_min|late_cut|late)\s*=\s*(0\.\d+)", float), kind="parameter", step="crossspecies_gene_selection.py", desc="Centroid cutoff for late genes")
    XS = "the 14 most concordant early and 14 most concordant late genes were displayed (r = 0.88-0.99)"
    add("xs_heatmap_n_early", XS, "14", "intplain", lambda c: c.sel_phase("early"), unit="genes", step="crossspecies_gene_selection.py", desc="Early genes displayed", nth=0)
    add("xs_heatmap_n_late", XS, "14", "intplain", lambda c: c.sel_phase("late"), unit="genes", step="crossspecies_gene_selection.py", desc="Late genes displayed", nth=1)
    add("xs_heatmap_r_min", XS, "0.88", "dec:2", lambda c: (c.F["mi"]["heatmap_selection"]["r50_range_selected"][0], f"{X}/results/mi_dpt/overlap_harmonized_nall/crossspecies_heatmap_selection/selection_summary.json"), unit="Pearson r", kind="statistic", step="crossspecies_gene_selection.py", desc="Lowest r among displayed genes")
    add("xs_heatmap_r_max", XS, "0.99", "dec:2", lambda c: (c.F["mi"]["heatmap_selection"]["r50_range_selected"][1], f"{X}/results/mi_dpt/overlap_harmonized_nall/crossspecies_heatmap_selection/selection_summary.json", f"selection on the harmonized MI runs; {c.F['mi']['heatmap_selection'].get('overlap_with_previous_28')} of 28 genes shared with the previous selection"), unit="Pearson r", kind="statistic", step="crossspecies_gene_selection.py", desc="Highest r among displayed genes")

    # ================================================================ SAMap
    add("samap_proteome_release", "(Ensembl release 113; GRCh38, GRCm39)", "113", "intplain", lambda c: c.proteome_release(), kind="reference", step="19_homology.smk", desc="Ensembl release of the proteomes")
    add("hmmer_version", "(HMMER 3.4;", "3.4", "str", lambda c: (re.search(r"hmmer-([\d.]+)", c.cfg["samap"]["phmmer"]).group(1), CFG + ": samap.phmmer"), kind="version", step="19_homology.smk", desc="HMMER version")
    add("phmmer_evalue", "at E ≤ 10−5", "10−5", "pow10", lambda c: (float(re.search(r"PHMMER_EVALUE\s*=\s*\"([^\"]+)\"", c.code("19_homology.smk")).group(1)), f"{X}/workflow/rules/19_homology.smk: PHMMER_EVALUE"), kind="parameter", step="19_homology.smk", desc="phmmer E-value cutoff")
    add("samap_version", "(sc-samap 3.0.1;", "3.0.1", "str", lambda c: c.samap_version(), kind="version", step=S_SAM, desc="SAMap package version")
    add("samap_n_iter", "for three iterations", "three", "word", cfgv("samap", "n_iter"), unit="iterations", kind="parameter", step=S_SAM, desc="SAMap iterations")
    add("samap_crossk", "with 20 cross-species neighbors per nucleus", "20", "intplain", lambda c: c.samap_crossk(), unit="neighbors", kind="parameter", step=S_SAM, desc="Cross-species neighbors per nucleus")
    ESP = "88% for human fibroblasts, 82% for human decidual cells, 85% for mouse fibroblasts and 92% for mouse decidual cells"
    for key, d, src_c, shown in (("samap_share_hs_fib", "human", "endometrial_fibroblast", "88%"), ("samap_share_hs_dec", "human", "decidual_cell", "82%"),
                                 ("samap_share_mm_fib", "mouse", "endometrial_fibroblast", "85%"), ("samap_share_mm_dec", "mouse", "decidual_cell", "92%")):
        add(key, ESP, shown, "pct:0", (lambda dd, s: lambda c: c.edge_share(dd, s))(d, src_c), unit="% of cross-species edge weight", kind="statistic", step="figures/samap_supp", desc=f"Edge-weight share landing in the matching compartment, {d} {src_c}")
    EWP = "mean total cross-species edge weight 12.6-13.5 of a maximum of 20"
    add("samap_edge_total_min", EWP, "12.6", "dec:1", lambda c: c.samap_missing("edge totals"), kind="statistic", step=S_SAM, desc="Lowest compartment mean of total cross-species edge weight")
    add("samap_edge_total_max", EWP, "13.5", "dec:1", lambda c: c.samap_missing("edge totals"), kind="statistic", step=S_SAM, desc="Highest compartment mean of total cross-species edge weight")
    add("samap_crossk_2", EWP, "20", "intplain", lambda c: c.samap_crossk(), unit="neighbors", kind="parameter", step=S_SAM, desc="Maximum total edge weight (= cross-species neighbors)")
    ASP = "were 0.58 for both matched compartments and 0.08 and 0.12 off-diagonal"
    add("samap_align_matched", ASP, "0.58", "dec:2", lambda c: c.samap_missing("alignment scores"), kind="statistic", step=S_SAM, desc="SAMap alignment score, matched compartments")
    add("samap_align_off1", ASP, "0.08", "dec:2", lambda c: c.samap_missing("alignment scores"), kind="statistic", step=S_SAM, desc="SAMap alignment score, off-diagonal (lower)")
    add("samap_align_off2", ASP, "0.12", "dec:2", lambda c: c.samap_missing("alignment scores"), kind="statistic", step=S_SAM, desc="SAMap alignment score, off-diagonal (higher)")
    add("apoptotic_decidual_mouse_2", "whereas the 907 apoptotic decidual nuclei carry no pseudotime", "907", "int", lambda c: (c.F["samap"]["transfer_summary"]["pt_raw"]["source_join"]["exact_nan"], f"{X}/results/samap_raw/pt_transfer_noapop_joinfix/summary.json"), unit="nuclei", step="samap_transfer_pt.py", desc="Mouse source nuclei without pseudotime")
    TRP = "47,509 of 47,533 human nuclei (99.95%) received a label and 47,352 (99.6%) a pseudotime value"
    add("samap_label_covered", TRP, "47,509", "int", lambda c: c.samap_missing("label transfer"), unit="nuclei", step="label_stain_full.py", desc="Human nuclei receiving a transferred label")
    add("fibdec_nuclei_human_2", TRP, "47,533", "int", lambda c: (c.F["samap"]["hs_transfer_n"], f"{N}/figures/samap_supp/data/hs_transfer_spatial.csv"), unit="nuclei", step="samap_transfer_pt.py", desc="Human fibroblast/decidual nuclei (second mention)")
    add("samap_label_cov_pct", TRP, "99.95%", "pct:2", lambda c: c.samap_missing("label transfer"), unit="%", step="label_stain_full.py", desc="Share of human nuclei receiving a label")
    add("samap_pt_covered", TRP, "47,352", "int", lambda c: (c.F["samap"]["transfer_summary"]["pt_ranknorm"]["hs_cells_covered"], f"{X}/results/samap_raw/pt_transfer_noapop_joinfix/summary.json", "transferred pseudotime is the legacy dpt_pseudotime_no_apoptotic (2026_08_28 object), not pt_harmonized"), unit="nuclei", step="samap_transfer_pt.py", desc="Human nuclei receiving a pseudotime value")
    add("samap_pt_cov_pct", TRP, "99.6%", "pct:1", lambda c: (c.F["samap"]["transfer_summary"]["pt_ranknorm"]["hs_coverage"], f"{X}/results/samap_raw/pt_transfer_noapop_joinfix/summary.json"), unit="%", step="samap_transfer_pt.py", desc="Share of human nuclei receiving pseudotime")
    LAP = "agreed with the human annotation for 84.2% of nuclei (81.3% of decidual cells and 87.0% of endometrial fibroblasts)"
    add("samap_label_agree", LAP, "84.2%", "pct:1", lambda c: c.samap_missing("label transfer"), unit="%", kind="statistic", step="label_stain_full.py", desc="Label agreement, all human nuclei")
    add("samap_label_agree_dec", LAP, "81.3%", "pct:1", lambda c: c.samap_missing("label transfer"), unit="%", kind="statistic", step="label_stain_full.py", desc="Label agreement, human decidual cells")
    add("samap_label_agree_fib", LAP, "87.0%", "pct:1", lambda c: c.samap_missing("label transfer"), unit="%", kind="statistic", step="label_stain_full.py", desc="Label agreement, human fibroblasts")
    # Agreement of transferred with human pseudotime: the observed statistics of the permutation-null run on the
    # harmonized pseudotime (ptnull_permutation.py), which also holds the permutation controls quoted after them.
    S_PTN = "samap_ptnull (ptnull_permutation.py on pt_harmonized_rank, seed 20261022, 1,000 permutations per design)"
    PTA = "was 0.77 over all nuclei, 0.73 within decidual cells and 0.30 within endometrial fibroblasts"
    add("samap_pt_rho_all", PTA, "0.77", "dec:2", lambda c: c.ptnull("rho_all", "observed"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Transferred vs human pseudotime, all nuclei")
    add("samap_pt_rho_dec", PTA, "0.73", "dec:2", lambda c: c.ptnull("rho_decidual", "observed"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Transferred vs human pseudotime, decidual cells")
    add("samap_pt_rho_fib", PTA, "0.30", "dec:2", lambda c: c.ptnull("rho_fibroblast", "observed"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Transferred vs human pseudotime, fibroblasts")
    UNI = "removed the gradient (maximum ρ, 0.20; maximum Moran"
    add("ptnull_uniform_max_rho", UNI, "0.20", "dec:2", lambda c: c.ptnull_max("rho_", "uniform_max"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Highest rho under permutation across all mouse nuclei (all, decidual, fibroblast)")
    add("ptnull_uniform_max_moran", "maximum Moran’s I, 0.08)", "0.08", "dec:2", lambda c: c.ptnull_max("moran_", "uniform_max"), unit="Moran's I", kind="statistic", step=S_PTN, desc="Highest Moran's I under permutation across all mouse nuclei")
    WA = "(means over permutations: ρ = 0.58 vs 0.77 observed)"
    add("ptnull_within_mean_rho_all", WA, "0.58", "dec:2", lambda c: c.ptnull("rho_all", "within_compartment_mean"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Mean rho over within-compartment permutations, all nuclei")
    add("samap_pt_rho_all_2", WA, "0.77", "dec:2", lambda c: c.ptnull("rho_all", "observed"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Transferred vs human pseudotime, all nuclei (second mention)")
    WD = "(ρ = 0.10 vs 0.73; Moran’s I = 0.21 vs 0.53;"
    add("ptnull_within_mean_rho_dec", WD, "0.10", "dec:2", lambda c: c.ptnull("rho_decidual", "within_compartment_mean"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Mean rho over within-compartment permutations, decidual cells")
    add("samap_pt_rho_dec_2", WD, "0.73", "dec:2", lambda c: c.ptnull("rho_decidual", "observed"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Transferred vs human pseudotime, decidual cells (second mention)")
    add("ptnull_within_mean_moran_dec", WD, "0.21", "dec:2", lambda c: c.ptnull("moran_within_decidual", "within_compartment_mean"), unit="Moran's I", kind="statistic", step=S_PTN, desc="Mean Moran's I within decidual cells over within-compartment permutations")
    add("ptnull_observed_moran_dec", WD, "0.53", "dec:2", lambda c: c.ptnull("moran_within_decidual", "observed"), unit="Moran's I", kind="statistic", step=S_PTN, desc="Observed Moran's I of transferred pseudotime within decidual cells")
    WF = "(ρ = 0.18 vs 0.30; 3 of 1,000 permutations at or above the observed value)"
    add("ptnull_within_mean_rho_fib", WF, "0.18", "dec:2", lambda c: c.ptnull("rho_fibroblast", "within_compartment_mean"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Mean rho over within-compartment permutations, fibroblasts")
    add("samap_pt_rho_fib_2", WF, "0.30", "dec:2", lambda c: c.ptnull("rho_fibroblast", "observed"), unit="Spearman rho", kind="statistic", step=S_PTN, desc="Transferred vs human pseudotime, fibroblasts (second mention)")
    add("ptnull_within_n_ge_obs_fib", WF, "3", "intplain", lambda c: c.ptnull("rho_fibroblast", "within_compartment_n_ge_obs"), unit="permutations", kind="statistic", step=S_PTN, desc="Within-compartment permutations with fibroblast rho at or above the observed value")
    # The sentence now names all 47,533 nuclei (flagged ones kept as sources) and no flag count, so the third mention of
    # the fibroblast/decidual flags (outliers_human_fibdec_3) has no number left to check.
    SRC = "The source nuclei were all 47,533 human fibroblast/decidual nuclei in the SAMap graph"
    add("samap_source_human_2", SRC, "47,533", "int", lambda c: c.samap_source_human(), unit="nuclei", step="21_samap_raw.smk (joint SAMap graph)", desc="SAMap and Ulrich source set, human (second mention)")

    # ================================================================ SCENIC
    SVP = "(version 0.12.1; ctxcore 0.2.0, arboreto 0.1.6;"
    for key, pkg, shown in (("pyscenic_version", "pyscenic", "0.12.1"), ("ctxcore_version", "ctxcore", "0.2.0"), ("arboreto_version", "arboreto", "0.1.6")):
        add(key, SVP, shown, "str", (lambda p: lambda c: c.ver("pyscenic", p))(pkg), kind="version", step=S_SC + " (config scenic.python = conda env pyscenic)", desc=f"{pkg} version")
    SCN = "(43,197 mouse, including the 907 apoptotic decidual nuclei, and 47,533 human)"
    add("fibdec_nuclei_mouse_2", SCN, "43,197", "int", lambda c: (c.obj("mouse_fibdec")["n_obs"], c.objpath("mouse_fibdec")), unit="nuclei", step=S_SC, desc="Mouse fibroblast/decidual nuclei (second mention)")
    add("fibdec_nuclei_human_3", SCN, "47,533", "int", lambda c: (c.obj("human_fibdec")["n_obs"], c.objpath("human_fibdec")), unit="nuclei", step=S_SC, desc="Human fibroblast/decidual nuclei (third mention)")
    GFP = "(from 33,698 to 14,245 genes in mouse and from 36,601 to 15,011 in human)"
    add("scenic_genes_all_mouse", GFP, "33,698", "int", lambda c: c.scenic_genes("mouse", False), unit="genes", step=S_SC, desc="Genes in the mouse object")
    add("scenic_genes_1pct_mouse", GFP, "14,245", "int", lambda c: c.scenic_genes("mouse", True), unit="genes", step=S_SC, desc="Mouse genes detected in at least 1% of nuclei")
    add("scenic_genes_all_human", GFP, "36,601", "int", lambda c: c.scenic_genes("human", False), unit="genes", step=S_SC, desc="Genes in the human object")
    add("scenic_genes_1pct_human", GFP, "15,011", "int", lambda c: c.scenic_genes("human", True), unit="genes", step=S_SC, desc="Human genes detected in at least 1% of nuclei")
    TFP = "leaving 1,179 (mouse) and 1,308 (human) candidate transcription factors"
    add("scenic_tfs_mouse", TFP, "1,179", "int", lambda c: c.scenic_tfs("mouse"), unit="TFs", step=S_SC, desc="Candidate TFs after the 1% filter, mouse")
    add("scenic_tfs_human", TFP, "1,308", "int", lambda c: c.scenic_tfs("human"), unit="TFs", step=S_SC, desc="Candidate TFs after the 1% filter, human")
    add("scenic_n_runs", "inference was repeated 100 times per species", "100", "intplain", cfgv("scenic", "n_seeds"), unit="runs", kind="parameter", step=S_SC, desc="GRN inference runs per species")
    add("scenic_sketch", "geometric sketch of 15,000 nuclei", "15,000", "int", lambda c: (int(str(c.cfg["scenic"]["sketch"]).lower().replace("k", "")) * (1000 if "k" in str(c.cfg["scenic"]["sketch"]).lower() else 1), CFG + ": scenic.sketch"), unit="nuclei", kind="parameter", step=S_SC, desc="Nuclei per geometric sketch")
    add("geosketch_version", "(geosketch 1.3;", "1.3", "str", lambda c: c.ver("pyscenic", "geosketch"), kind="version", step=S_SC, desc="geosketch version")
    add("motif_collection", "against the v10 clustered motif collection", "v10", "str", lambda c: c.motif_version(), kind="reference", step=S_SC, desc="cisTarget motif collection version")
    CTP = "(AUC threshold 0.05, NES threshold 3.0)"
    add("ctx_rank_threshold", "(rank threshold 1,500, the default of the pySCENIC Python interface", "1,500", "int", lambda c: c.pyscenic_default("rank_threshold"), kind="parameter", step=S_SC, desc="cisTarget rank threshold")
    add("ctx_auc_threshold", CTP, "0.05", "dec:2", lambda c: c.pyscenic_default("auc_threshold"), kind="parameter", step=S_SC, desc="cisTarget AUC threshold")
    add("ctx_nes_threshold", CTP, "3.0", "dec:1", lambda c: c.pyscenic_default("nes_threshold"), kind="parameter", step=S_SC, desc="cisTarget NES threshold")
    add("scenic_n_runs_2", "Consensus regulons were assembled across the 100 runs", "100", "intplain", cfgv("scenic", "n_seeds"), unit="runs", kind="parameter", step=S_SC, desc="Runs (second mention)")
    add("scenic_reg_thresh", "in at least 80 of the runs", "80", "intplain", lambda c: c.scenic_summary("mouse", "reg_thresh"), unit="runs", kind="parameter", step=S_SC, desc="Runs in which a TF must form a regulon")
    RGP = "(193 of 771 factors with at least one regulon in mouse; 141 of 941 in human)"
    add("scenic_kept_mouse", RGP, "193", "intplain", lambda c: c.scenic_summary("mouse", "kept"), unit="TFs", step=S_SC, desc="Consensus regulons kept, mouse")
    add("scenic_seen_mouse", RGP, "771", "intplain", lambda c: c.scenic_summary("mouse", "seen"), unit="TFs", step=S_SC, desc="TFs with at least one regulon in any run, mouse")
    add("scenic_kept_human", RGP, "141", "intplain", lambda c: c.scenic_summary("human", "kept"), unit="TFs", step=S_SC, desc="Consensus regulons kept, human")
    add("scenic_seen_human", RGP, "941", "intplain", lambda c: c.scenic_summary("human", "seen"), unit="TFs", step=S_SC, desc="TFs with at least one regulon in any run, human")
    add("scenic_target_floor", "recurred in more than five runs", "five", "word", lambda c: c.scenic_floor(), unit="runs", kind="parameter", step=S_SC, desc="Target recurrence floor (more than this many runs)")
    add("aucell_threshold", "with the pySCENIC default AUC threshold of 0.05", "0.05", "dec:2", lambda c: c.aucell_thr(), kind="parameter", step=S_SC, desc="AUCell AUC threshold")
    add("ensembl_release_orth_2", "with the Ensembl 116 one-to-one ortholog table", "116", "intplain", lambda c: (int(re.search(r"(\d+)", c.F["mi"]["overlap_report"]["ortholog_release"][0]).group(1)), c.F["mi"]["overlap_report"]["ortholog_source"]), kind="reference", step=S_SC, desc="Ensembl release (second mention)")
    OTP = "190 of the 193 mouse and 136 of the 141 human regulon factors"
    add("scenic_tfs_with_ortholog", OTP, "190", "intplain", lambda c: c.scenic_orth(), unit="TFs", step=S_SC, desc="Kept mouse regulon TFs with a 1:1 human ortholog")
    add("scenic_kept_mouse_2", OTP, "193", "intplain", lambda c: c.scenic_summary("mouse", "kept"), unit="TFs", step=S_SC, desc="Consensus regulons kept, mouse (second mention)")

    # ================================================================ decidual zone modules
    add("zone_n_sets", "Gene sets for the eleven uterine microenvironments", "eleven", "word", lambda c: (len(c.F["zone_genesets"]), f"{X}/resources/signatures/aucell_modules_mouse.yaml"), unit="gene sets", kind="parameter", step="25_zone_modules.smk", desc="Microenvironment gene sets")
    add("zone_genes_per_set", "which lists the top 20 marker genes per microenvironment", "20", "str", lambda c: c.zone_sizes(), unit="genes", kind="parameter", step="25_zone_modules.smk", desc="Genes per gene set as used")
    add("zone_auc_window", "with the default AUC threshold of 0.05 and a fixed random seed (0)", "0.05", "dec:2", lambda c: c.zone_auc(), unit="AUC threshold (share of ranked genes)", kind="parameter", step="25_zone_modules.smk", desc="AUCell window")
    add("zone_display_pct", "from 0 to its 99th percentile", "99", "intplain", lambda c: c.zone_p99(), unit="percentile", kind="parameter", step="25_zone_modules.smk (fig_zone_modules.py)", desc="Color ceiling in spatial and heatmap panels")

    # ================================================================ in vitro / in vivo comparison
    for r_ in ivd_rows(): R.append(r_)

    # ================================================================ reproducibility
    add("cuda_version", "(CUDA 12.9)", "12.9", "str", lambda c: (re.sub(r"^cu(\d+)(\d)$", r"\1.\2", c.ver("rsc_uv", "torch")[0].split("+")[1]), c.ver("rsc_uv", "torch")[1] + " (CUDA build tag)"), kind="version", step="pipeline environment rsc_uv", desc="CUDA version of the torch build")
    DEP = "every annotated nucleus (57,033 mouse and 67,847 human)"
    add("final_nuclei_mouse_deposited", DEP, "57,033", "int", lambda c: (c.obj("mouse_all")["n_obs"], c.objpath("mouse_all")), unit="nuclei", step="15_release.smk", desc="Annotated nuclei in the deposited objects, mouse")
    add("final_nuclei_human_deposited", DEP, "67,847", "int", lambda c: (c.obj("human_all")["n_obs"], c.objpath("human_all")), unit="nuclei", step="15_release.smk", desc="Annotated nuclei in the deposited objects, human")
    return R


def ivd_rows():
    from registry_core import Row
    R = []
    def add(*a, **k): R.append(Row(*a, **k))
    S1 = "consensus list (human_ivd_repro + derive_consensus123_orthology.py)"
    add("ivd_n_studies", "We used three studies that profiled human", "three", "word", lambda c: c.ivd("n_studies"), unit="studies", kind="parameter", step=S1, desc="Human in vitro studies used")
    add("ivd_aghajanova_genes", "(Tables 4 and 5; 179 genes)", "179", "int", lambda c: c.ivd("aghajanova"), unit="genes", step=S1, desc="Aghajanova et al. genes (women without endometriosis)")
    add("ivd_rytkonen_genes", "(Supplemental Table S1A-D; 4,672 genes)", "4,672", "int", lambda c: c.ivd("rytkonen"), unit="genes", step=S1, desc="Rytkönen et al. genes")
    add("ivd_stadtmauer_tabulated", "report per-gene values only for 30 genes", "30", "intplain", lambda c: c.ivd("stad_tab"), unit="genes", step=S1, desc="Stadtmauer & Wagner tabulated genes")
    add("ivd_stadtmauer_met", "which 26 genes did", "26", "intplain", lambda c: c.ivd("stad_met"), unit="genes", step=S1, desc="Tabulated Stadtmauer genes meeting the authors' criteria")
    add("ivd_stadtmauer_text", "For 37 further genes", "37", "intplain", lambda c: c.ivd("stad_text"), unit="genes", step=S1, desc="Stadtmauer genes transcribed from the Results text")
    add("ivd_stadtmauer_total", "(63 genes in total)", "63", "intplain", lambda c: c.ivd("stad_total"), unit="genes", step=S1, desc="Stadtmauer genes in total")
    add("ivd_rytkonen_conflicts", "(45 in Rytkönen et al., none in the other two)", "45", "intplain", lambda c: c.ivd("ryt_conflicts"), unit="genes", step=S1, desc="Genes with conflicting calls within Rytkönen et al.")
    CP = "Of 147 genes called by two or more studies, 123 agreed (93 up, 30 down) and 24 conflicted"
    for key, fld, shown, d in (("ivd_called_2plus", "called2", "147", "Genes called by two or more studies"), ("ivd_consensus", "consensus", "123", "Consensus genes"),
                               ("ivd_consensus_up", "up", "93", "Consensus genes up"), ("ivd_consensus_down", "down", "30", "Consensus genes down"),
                               ("ivd_conflicted", "conflicted", "24", "Genes with conflicting calls across studies")):
        add(key, CP, shown, "intplain", (lambda f: lambda c: c.ivd(f))(fld), unit="genes", step=S1, desc=d)
    SP = "The 123 comprise 77 genes supported by Aghajanova and Rytkönen, 39 by Rytkönen and Stadtmauer (a single cell line) and 7 by all three"
    for key, fld, shown, d in (("ivd_consensus_2", "consensus", "123", "Consensus genes (second mention)"), ("ivd_support_ar", "sup_ar", "77", "Supported by Aghajanova and Rytkönen"),
                               ("ivd_support_rs", "sup_rs", "39", "Supported by Rytkönen and Stadtmauer"), ("ivd_support_all", "sup_all", "7", "Supported by all three")):
        add(key, SP, shown, "intplain", (lambda f: lambda c: c.ivd(f))(fld), unit="genes", step=S1, desc=d)
    add("ivd_renamed", "(8 had been renamed, for example", "8", "intplain", lambda c: c.ivd("renamed"), unit="genes", step=S1, desc="Consensus symbols renamed since publication")
    OP = "Only one-to-one orthologs were retained (114 of 123)"
    add("ivd_one2one", OP, "114", "intplain", lambda c: c.ivd("one2one"), unit="genes", step=S1, desc="Consensus genes with a strict 1:1 mouse ortholog")
    add("ivd_consensus_3", OP, "123", "intplain", lambda c: c.ivd("consensus"), unit="genes", step=S1, desc="Consensus genes (third mention)")
    S2 = "X-Mens pseudobulk DE with the uninduced fibroblasts (fit_mouse_de.py, 27 Oct fit)"
    PBP = "giving 11 pseudobulks from 43,133 nuclei (22,855 decidual and 20,278 fibroblast; 1,170 to 9,270 per pseudobulk, minimum 25)"
    for key, fld, shown, how, d in (("ivd_pseudobulks", "pb_n", "11", "intplain", "Pseudobulks"), ("ivd_pb_nuclei", "pb_nuclei", "43,133", "int", "Nuclei in the pseudobulks"),
                                    ("ivd_pb_decidual", "pb_dec", "22,855", "int", "Decidual nuclei"), ("ivd_pb_fibroblast", "pb_fib", "20,278", "int", "Fibroblast nuclei"),
                                    ("ivd_pb_min", "pb_min", "1,170", "int", "Smallest pseudobulk"), ("ivd_pb_max", "pb_max", "9,270", "int", "Largest pseudobulk"),
                                    ("ivd_pb_floor", "pb_floor", "25", "intplain", "Minimum nuclei per pseudobulk")):
        add(key, PBP, shown, how, (lambda f: lambda c: c.ivd(f))(fld), unit="nuclei", step=S2, desc=d, kind="count" if fld != "pb_floor" else "parameter")
    add("ivd_uninduced_decidual", "(its 64 nuclei labelled decidual were excluded)", "64", "intplain", lambda c: c.ivd("uninduced_decidual"), unit="nuclei", step=S2, desc="Uninduced nuclei labelled decidual, left out of the pseudobulks")
    add("ivd_eligible_genes", "were retained (16,784)", "16,784", "int", lambda c: c.ivd("eligible"), unit="genes", step=S2, desc="Genes passing the CPM filter")
    add("pydeseq2_version", "PyDESeq2 v0.5.4", "0.5.4", "str", lambda c: c.ver("rsc_uv", "pydeseq2"), kind="version", step=S2, desc="PyDESeq2 version")
    add("ivd_resid_df", "(4 residual degrees of freedom)", "4", "intplain", lambda c: c.ivd("resid_df"), unit="df", kind="parameter", step=S2, desc="Residual degrees of freedom")
    add("ivd_one2one_2", "Each of the 114 orthologs was counted", "114", "intplain", lambda c: c.ivd("one2one"), unit="genes", step=S1, desc="Consensus orthologs (second mention)")
    XP = "X-Mens called 50 orthologs, 20 concordant and 30 discordant (40%; 46% expected by chance); 53 were tested but not called and 11 were not testable"
    for key, fld, shown, how, d in (("ivd_xm_called", "xm_called", "50", "intplain", "Consensus orthologs called by X-Mens"), ("ivd_xm_conc", "xm_conc", "20", "intplain", "Concordant"),
                                    ("ivd_xm_disc", "xm_disc", "30", "intplain", "Discordant"), ("ivd_xm_pct", "xm_pct", "40%", "pct:0", "Concordant share"),
                                    ("ivd_xm_exp", "xm_exp", "46%", "pct:0", "Expected concordant share"), ("ivd_xm_notcalled", "xm_notcalled", "53", "intplain", "Tested, not called"),
                                    ("ivd_xm_nottestable", "xm_nottest", "11", "intplain", "Not testable")):
        add(key, XP, shown, how, (lambda f: lambda c: c.ivd(f))(fld), unit="genes" if "%" not in shown else "%", kind="count" if "%" not in shown else "statistic", step=S2, desc=d)
    for key, ph, fld, shown in (("wang_pregnancy_genes", "inter-implantation sites on day 8; 5,277 genes)", "w_preg_n", "5,277"),
                                ("wang_oil_genes", "unstimulated horn on day 8 of pseudopregnancy; 4,294 genes)", "w_oil_n", "4,294"),
                                ("wang_ivd_genes", "days, versus control; 1,501 genes)", "w_ivd_n", "1,501")):
        add(key, ph, shown, "int", (lambda f: lambda c: c.ivd(f))(fld), unit="genes", step="Wang et al. 2020b lists (cagri resources)", desc=f"Wang list size ({fld})")
    MP = "(1,446 of 1,501 matched)"
    add("wang_ivd_matched", MP, "1,446", "int", lambda c: c.ivd("w_ivd_matched"), unit="genes", step=S2, desc="Wang in vitro genes matched to X-Mens")
    add("wang_ivd_genes_2", MP, "1,501", "int", lambda c: c.ivd("w_ivd_n"), unit="genes", step=S2, desc="Wang in vitro genes (second mention)")
    for key, ph, pair in (("xm_oil", "X-Mens agreed with oil-induced decidualization for 1,439 of 1,588 genes (91%; 50% expected)", "xm_oil"),
                          ("xm_preg", "with natural pregnancy for 1,105 of 1,518 (73%; 53%)", "xm_preg"),
                          ("xm_ivd", "and with in vitro decidualization for 461 of 657 (70%; 53%)", "xm_ivd"),
                          ("ivd_preg", "in vitro decidualization agreed with natural pregnancy for 310 of 766 genes (40%; 44%)", "ivd_preg"),
                          ("ivd_oil", "and with oil-induced decidualization for 531 of 768 (69%; 53%)", "ivd_oil"),
                          ("oil_preg", "oil-induced decidualization agreed with natural pregnancy for 1,977 of 2,194 (90%; 50%)", "oil_preg")):
        nums = re.findall(r"(\d[\d,]*) of (\d[\d,]*)[^(]*\((\d+)%; (\d+)%", ph)[0]
        for fld, shown, how, d in (("same", nums[0], "int", "genes changing in the same direction"), ("n", nums[1], "int", "genes compared"),
                                   ("pct", nums[2] + "%", "pct:0", "share in the same direction"), ("exp", nums[3] + "%", "pct:0", "share expected by chance")):
            add(f"wang_{key}_{fld}", ph, shown, how, (lambda p, f: lambda c: c.wang(p, f))(pair, fld), unit="genes" if how == "int" else "%",
                kind="count" if how == "int" else "statistic", step=("SuppTable_Wang_comparison (28 Oct supplementary tables)" if pair.startswith("xm") else "direction agreement (Wang lists)"), desc=f"{key.replace('_', ' vs ')}: {d}")
    return R
