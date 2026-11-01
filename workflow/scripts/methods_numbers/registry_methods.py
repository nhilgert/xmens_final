"""Recompute methods for the methods-number registry (extends registry_rows.Ctx).

Every method returns (value, source_file[, note]). Inputs: the facts directories written by collect_facts.py (FD) and
collect_facts2.py (FD2), small_reads.json (SR), and FIX = the 2026-10-22 Stadtmauer-correction outputs of the in vitro
consensus (job_outputs of consensus123_stadtmauer_fix), which hold the corrected list the draft describes. EI is the
directory given with --extra-inputs: project artifacts staged outside the repository (MANIFEST.json gives each file's
artifact version and sha256), namely the X-Mens fit with the uninduced fibroblasts (27 Oct). The permutation-null statistics of the SAMap
pseudotime transfer are a pipeline output (PTNULL_TABLE, rule samap_ptnull).
"""
import ast
import io
import os, re, json, functools
import zipfile
import numpy as np, pandas as pd
import registry_rows as rr

L = rr.L; X = rr.X; N = rr.N; O22 = rr.O22
HIR = f"{L}/human_ivd_repro"; CAG = f"{L}/cagri_ivd_end_to_end_20260815"; SCR = f"{L}/decidualization/revisions/scenic_v2_2026_08_28"

# Permutation-null statistics of the SAMap pseudotime transfer (rule samap_ptnull, workflow/rules/29_fig_samap.smk).
PTNULL_TABLE = f"{X}/results/samap_figures/ptnull/ptnull_statistics_table.csv"

# Inputs staged with --extra-inputs.
FIT_BUNDLE = "xmens_with_uninduced_20261027.zip"  # fit_mouse_de.py with the uninduced fibroblasts: diagnostics, pseudobulks, DE

# Epithelial labels that follow tissue position, and the per-cluster vote columns of the notebooks, in the same order.
EPI_SPLIT_LABELS = {"human": ("epithelial_inner", "epithelial_outer"), "mouse": ("epithelial_luminal", "epithelial_glandular")}
EPI_VOTE_COLUMNS = {"human": ("n_inner", "n_outer"), "mouse": ("n_luminal", "n_glandular")}

# Rows of the Summary sheet of SuppTable_Wang_comparison, by registry pair.
WANG_SUMMARY_ROWS = {"xm_oil": "X-Mens vs Wang oil-induced", "xm_preg": "X-Mens vs Wang natural pregnancy",
                     "xm_ivd": "X-Mens vs Wang in vitro", "ivd_preg": "Wang in vitro vs Wang natural pregnancy",
                     "ivd_oil": "Wang in vitro vs Wang oil-induced", "oil_preg": "Wang oil-induced vs Wang natural pregnancy"}


def pct_range(percents, decimals):
    """A range of percentages as the methods write it: '18–38%', or '4%' when both ends round alike."""
    low = f"{min(percents):.{decimals}f}"
    high = f"{max(percents):.{decimals}f}"
    return f"{low}%" if low == high else f"{low}–{high}%"


class Ctx2(rr.Ctx):
    def __init__(self, facts_dir, facts2_dir=None, small_reads=None, fix_dir=None, fix_label="", artifact_inputs=None, proteome_check=None,
                 extra_inputs=None):
        super().__init__(facts_dir)
        self.AI = artifact_inputs
        self.AIM = json.load(open(f"{artifact_inputs}/MANIFEST.json")) if artifact_inputs else {}
        self.EI = extra_inputs
        self.EIM = json.load(open(f"{extra_inputs}/MANIFEST.json")) if extra_inputs else {}
        self.PCHK = proteome_check
        self.FD2 = facts2_dir
        self.F2 = json.load(open(f"{facts2_dir}/facts2.json")) if facts2_dir and os.path.exists(f"{facts2_dir}/facts2.json") else {}
        self.SR = json.load(open(small_reads)) if small_reads else {}
        self.FIX = fix_dir; self.FIX_LABEL = fix_label

    # ------------------------------------------------------------------ helpers
    def f2(self, *path):
        v = self.F2
        for p in path:
            if not isinstance(v, dict) or p not in v: raise KeyError("/".join(map(str, path)) + " not collected")
            v = v[p]
        if isinstance(v, dict) and "_error" in v: raise RuntimeError(v["_error"])
        return v

    def fix(self, name, **kw):
        if str(self.FIX).endswith(".zip"):
            import zipfile, io
            z = zipfile.ZipFile(self.FIX); m = [n for n in z.namelist() if n.endswith("/job_outputs/" + name)]
            assert len(m) == 1, (name, m)
            return pd.read_csv(io.BytesIO(z.read(m[0])), **kw)
        return pd.read_csv(f"{self.FIX}/{name}", **kw)

    def fixsrc(self, name):
        return f"{self.FIX_LABEL}/{name}"

    # ------------------------------------------------------------------ atlas (fixes: cluster column, EXCLUDE list)
    def _exclude_clusters(self, sp):
        t = self.SR["annot"][f"{sp}_lineage.csv"]
        d = pd.read_csv(pd.io.common.StringIO(t), dtype=str)
        return d.loc[d.label.str.strip().str.upper().str.startswith("EXCLU"), "cluster"].str.strip().tolist()

    @functools.lru_cache(None)
    def atlas(self, sp):
        A = self.F["annot_tables"][f"{sp}_lineage.csv"]
        c = self.tab(f"atlas_res2.5_{sp}.csv.gz")
        cl = c["cluster"].astype(str)
        excl = set(self._exclude_clusters(sp)); n_ex = int(cl.isin(excl).sum())
        assert set(cl.unique()) <= set(pd.read_csv(pd.io.common.StringIO(self.SR["annot"][f"{sp}_lineage.csv"]), dtype=str).cluster.str.strip()), "cluster ids not in table"
        return dict(n_clusters=A["n_rows"], n_classes=len(A["labels"]), n_exclude=A["n_exclude"], excluded=n_ex, total=len(c), retained=len(c) - n_ex,
                    clusters_in_table=int(cl.nunique()), excluded_ids=sorted(excl, key=int))

    def doublet_mean(self, sp, excluded):
        c = self.tab(f"atlas_res2.5_{sp}.csv.gz").set_index("cell")
        d = self.tab(f"doublet_scores_{sp}.csv.gz", index_col=0)
        assert c.index.equals(d.index) or set(c.index) == set(d.index), "atlas table and doublet scores cover different nuclei"
        d = d.reindex(c.index)
        ex = c.cluster.astype(str).isin(set(self._exclude_clusters(sp)))
        v = d.doubletdetection_doublet_score.astype(float)[ex if excluded else ~ex].mean()
        # the atlas table small_reads.py read (round 2 for v9), from which the EXCLUDE clusters above come
        table = self.SR.get("annot_paths", {}).get(f"{sp}_lineage.csv", f"{rr.ANN}/{sp}_lineage.csv")
        return float(v), f"{self.F['doublets'][sp]['source']} obs['doubletdetection_doublet_score'] + {X}/results/tables/clusters/{sp}/leiden_res2.5.csv + {table}"

    # ------------------------------------------------------------------ spatial QC
    def spatial_mouse(self, step, tp=None):
        o = self.obj("mouse_all"); src = f"{self.objpath('mouse_all')} obs['spatial_qc_reason']"
        key = {"step1": "degenerate", "step2": "far", "step3": "boxed", "unmapped": "unmapped"}[step]
        if tp is None: return o["vc_spatial_qc_reason"].get(key, 0), src
        if tp not in o["vc_timepoint"]:
            raise KeyError(f"timepoint {tp} not in the release (it has {sorted(o['vc_timepoint'])})")
        return o["spatial_by_tp"].get(f"{key}|{tp}", 0), src + " by timepoint"

    def exclusion_regions(self, timepoint):
        """Hand-drawn exclusion regions at one release timepoint; the config lists them under the sample-sheet names."""
        aliases = self.timepoint_aliases()
        regions = self.cfg["spatial_qc"]["exclusions"]["mouse"]
        n_regions = sum(1 for region in regions if aliases.get(region.get("timepoint"), region.get("timepoint")) == timepoint)
        return n_regions, f"{rr.CFG}: spatial_qc.exclusions.mouse (timepoints renamed by the vocabulary's legacy aliases)"

    def spatial_human(self, what):
        o = self.obj("human_all"); src = f"{self.objpath('human_all')} obs['spatial_qc_reason']"
        key = {"blanked": "degenerate", "unmapped": "unmapped"}[what]
        return o["vc_spatial_qc_reason"].get(key, 0), src

    # ------------------------------------------------------------------ pseudotime
    def pt_val(self, sp, field):
        s = self.pt_stats(sp)
        return float(s[field]), f"{self.objpath(sp + '_fibdec')} obs['pt_harmonized_rank'] (+ subtype order and day map from provenance/harmonized_pseudotime.py)", f"n = {s['n']:,} nuclei with a pseudotime"

    def dc_rho(self, sp):
        p = self.F["pt_provenance"][sp]; k = str(p["trajectory_dc"])
        return abs(float(p["dc_rho"][k])), f"{O22}/provenance/{sp}_pca_released_provenance.json: dc_rho[{k}]"

    def latent_dim(self):
        m = re.findall(r"stromal latent attached \((\d+) dims\)", self.SR["pt_runs_log"])
        assert m and len(set(m)) == 1, m
        return int(m[0]), f"{N}/2026_10_22_harmonized_pt/runs/*/run_*.log ('stromal latent attached')"

    def pt_pcs(self):
        assert "Falling back to preprocessing with `sc.pp.pca` and default params" in self.SR["pt_runs_log"]
        return int(self.F["scanpy_defaults"]["settings_N_PCS"]), f"run logs (sc.pp.pca default) + scanpy {self.F['scanpy_defaults']['scanpy']} settings.N_PCS"

    def fine_types(self, sp):
        return len(self.obj(f"{sp}_all")["vc_cell_subtype"]), f"{self.objpath(sp + '_all')} obs['cell_subtype'] (distinct values)"

    def outlier_rule(self, which):
        t = self.obj("mouse_all")["uns_embedding_outlier_rule"]
        t = t if isinstance(t, str) else json.dumps(t)
        pat = {"k": r"(\d+) nearest neighbou?rs", "other": r">= ?(\d+)% to one other lineage", "same": r"<= ?(\d+)% to the cell's own lineage"}[which]
        v = int(re.search(pat, t).group(1))
        return (v / 100 if which != "k" else v), f"{self.objpath('mouse_all')} uns['embedding_outlier_rule']"

    # ------------------------------------------------------------------ cross-species heatmap selection
    def sel_phase(self, which):
        rows = self.F["mi"]["heatmap_selected"]["rows"]; assert len(rows) == self.F["mi"]["heatmap_selected"]["n"]
        selected = self.overlap_file("crossspecies_heatmap_selection/conserved_profiles_selected_genes.csv")
        return sum(r["phase"] == which for r in rows), selected

    def sel_code(self, pat, cast):
        m = re.search(pat, self.code("crossspecies_gene_selection.py")); assert m, pat
        return cast(m.group(1)), f"{X}/workflow/scripts/mi/crossspecies_gene_selection.py"

    # ------------------------------------------------------------------ SAMap
    def edge_share(self, species, src):
        direction = {"human": "human_to_mouse", "mouse": "mouse_to_human"}.get(species, species)
        r = [x for x in self.F["samap"]["edge_share"] if x["direction"] == direction and x["source_compartment"] == src]
        assert len(r) == 1
        return float(r[0][src]), f"{N}/figures/samap_supp/data/samap_compartment_edge_share.csv"

    def samap_version(self):
        v = sorted(set(re.findall(r"Version: (\S+)", self.SR["samap_env"])))
        assert len(v) == 1, v
        return v[0], f"{X}/envs/samap/.venv (sc_samap dist-info)"

    def samap_crossk(self):
        try:
            return int(self.f2("samap_graph", "K")), f"{N}/2026_09_02_samap_raw_all.h5ad uns['mapping_K']"
        except Exception:
            m = re.search(r"k=(\d+)", self.SR["samap_settings"]); return int(m.group(1)), f"{L}/code/samap-scvi/full_fibdec/run_arm.py RUN['k']"

    def samap_graph(self, what):
        g = self.f2("samap_graph"); src = f"{g['joint']} obsp['connectivities'] (cross-species block, K = {g['K']})"
        if what in ("edge_min", "edge_max"):
            v = g["edge_total_mean"]; v = {k: x for k, x in v.items() if not k.endswith(":nan")}
            return (min(v.values()) if what == "edge_min" else max(v.values())), src, "compartment means: " + ", ".join(f"{k} {x:.2f}" for k, x in v.items())
        A = pd.DataFrame(g["align_max_sym"])           # rows human compartment, columns mouse compartment
        A = A.loc[["endometrial_fibroblast", "decidual_cell"], ["endometrial_fibroblast", "decidual_cell"]]
        if what == "align_matched":
            a, b = A.iloc[0, 0], A.iloc[1, 1]
            note = f"fibroblast-fibroblast {a:.3f}, decidual-decidual {b:.3f}"
            return (float(a) if round(a, 2) == round(b, 2) else None), src, note + ("" if round(a, 2) == round(b, 2) else "; the two differ at 2 decimals")
        off = sorted([A.iloc[0, 1], A.iloc[1, 0]])
        return float(off[0] if what == "align_off1" else off[1]), src, f"human fibroblast-mouse decidual {A.iloc[0, 1]:.3f}, human decidual-mouse fibroblast {A.iloc[1, 0]:.3f}"

    def samap_label(self, what):
        g = self.f2("samap_graph"); src = f"{g['joint']} (argmax of connectivity-weighted vote of mouse neighbors' cell_type)"
        if what in ("covered", "total", "cov_pct"):
            v = {"covered": g["label_covered"], "total": g["label_total"], "cov_pct": g["label_covered"] / g["label_total"]}[what]
            st = g["stored_label_cols"].get("transferred_mm_coarse", {}); unm = st.get("unmapped")
            return v, src, ("" if unm is None or unm == g["label_total"] - g["label_covered"] else f"stored column transferred_mm_coarse has {unm} unmapped")
        t = pd.read_csv(f"{self.FD2}/samap_hs_label_transfer.csv.gz")
        if what != "agree": t = t[t.native == {"agree_dec": "decidual_cell", "agree_fib": "endometrial_fibroblast"}[what]]
        ok = (t.transferred.fillna("") == t.native)
        cov = t.covered.astype(bool)
        return float(ok.mean()), src, f"denominator = all {len(t):,} human nuclei (unlabeled count as disagreement); over labeled nuclei only {100 * ok[cov].mean():.2f}%"

    # ------------------------------------------------------------------ SCENIC
    def _scsum(self, sp):
        """The floor-6 consensus the methods describe (consensus_sketch15k_weighted_soft_occ6)."""
        return [v for k, v in self.SR["occ6"][sp].items() if k.endswith("recurrence_summary.txt")][0]

    def _scnote(self, sp):
        f = self.F["scenic"][sp]["summary"]; m = re.search(r"min_occurrence=(\d+)", f).group(1)
        if m == "6": return ""
        fig = self.tab(f"scenic_recurrence_{sp}.csv.gz"); a = set(fig.loc[fig.kept.astype(bool), "name"].astype(str))
        b = set(self.SR["occ6_kept"][sp]["kept"])
        return (f"figures/_data/scenic/{sp}_recurrence_report.tsv is the min_occurrence={m} run; its kept TF names "
                + ("are identical to the floor-6 run's" if a == b else f"differ from the floor-6 run's ({len(a - b)} only there, {len(b - a)} only in floor-6)"))

    def scenic_summary(self, sp, field):
        s = self._scsum(sp); src = f"{SCR}/scenic_seeds/{sp}/consensus_sketch15k_weighted_soft_occ6/recurrence_summary.txt"
        if field == "reg_thresh":
            return int(round(float(re.search(r"reg_thresh=([0-9.]+)", s).group(1)) * int(re.search(r"seeds: (\d+)", s).group(1)))), src, "reg_thresh x seeds"
        pat = {"kept": r"kept regulons: (\d+)", "seen": r"regulon names seen: (\d+)", "seeds": r"seeds: (\d+)"}[field]
        return int(re.search(pat, s).group(1)), src

    def scenic_floor(self):
        s = self._scsum("mouse"); m = int(re.search(r"min_occurrence=(\d+)", s).group(1)); assert m == int(re.search(r"min_occurrence=(\d+)", self._scsum("human")).group(1))
        return m - 1, f"{SCR}/scenic_seeds/*/consensus_sketch15k_weighted_soft_occ6/recurrence_summary.txt (min_occurrence={m})", self._scnote("mouse")

    def aucell_thr(self):
        v = sorted(set(re.findall(r"auc_threshold=([0-9.]+)", self.SR["occ6_logs"])))
        assert len(v) == 1, v
        return float(v[0]), f"{SCR}/consensus_variant_logs/*occ6* (AUCell line)"

    def pyscenic_default(self, name):
        s = self.F["scenic"]["pyscenic_defaults"]
        v = sorted(set(re.findall(rf"transform\.py:\d+:\s+{name}=([0-9.]+)", s)))
        assert len(v) == 1, (name, v)
        assert name not in self.code("ctx_only.py").split("prune2df(")[1].split(")")[0], f"ctx_only.py passes {name} explicitly"
        return float(v[0]), f"pyscenic env transform.py defaults (ctx_only.py passes no {name})"

    def motif_version(self):
        v = sorted(set(re.findall(r"_(v\d+)_clust|motifs-(v\d+)nr", " ".join(self.F["scenic"]["databases"]))))
        v = sorted(set(a or b for a, b in v)); assert len(v) == 1, v
        return v[0], f"{SCR}/databases (ranking and motif-table file names)"

    def _loom_manifest(self, sp):
        if sp == "human":
            p = [k for k in self.SR["loom_json"] if k.endswith("human_dsc_pooled.manifest.json")][0]
            return self.SR["loom_json"][p], p
        man = self.f2("scenic_mouse_manifest"); p = [k for k in man if k.endswith("mouse_dsc_pooled.manifest.json")][0]
        return man[p], p

    def scenic_genes(self, sp, filtered):
        m, p = self._loom_manifest(sp)
        v = m["n_genes_grn"] if filtered else m["n_genes_full"]
        # independent check from the per-gene detection counts of the released object
        d = self.tab(f"detect_{sp}.csv.gz"); n_obs = self.obj(f"{sp}_fibdec")["n_obs"]
        mc = max(1, int(round(0.01 * n_obs))); chk = int((d.n_cells >= mc).sum()) if filtered else len(d)
        return int(v), p, ("" if chk == v else f"detection counts of the Oct 22 object give {chk:,} (>= {mc} of {n_obs:,} nuclei)")

    def scenic_tfs(self, sp):
        m, p = self._loom_manifest(sp)
        return int(m["n_tfs_retained"]), p

    def scenic_orth(self):
        r = self.tab("scenic_recurrence_mouse.csv.gz"); k = r[r.kept.astype(bool)]
        t = k.name.astype(str).str.replace(r"\(\+\)$", "", regex=True)
        o = self.tab("scenic_ensembl116_one2one.csv.gz")
        m1 = o.mouse.map(o.mouse.value_counts()).eq(1); h1 = o.human.map(o.human.value_counts()).eq(1)
        return int(t.isin(set(o[m1 & h1].mouse)).sum()), f"{N}/figures/_data/scenic/mouse_recurrence_report.tsv (kept) + ensembl116_one2one.csv", f"of {len(t)} kept"

    # ------------------------------------------------------------------ zone modules
    def zone_sizes(self):
        s = sorted(set(len(v) for v in self.F["zone_genesets"].values()))
        return (str(s[0]) if len(s) == 1 else f"{s[0]}-{s[-1]}"), f"{X}/resources/signatures/aucell_modules_mouse.yaml"

    def zone_auc(self):
        runs = self.cfg["zone_modules"]["runs"]; v = sorted(set(r["auc_threshold"] for r in runs.values()))
        assert len(v) == 1, v
        return float(v[0]), f"{rr.CFG}: zone_modules.runs.*.auc_threshold"

    def zone_p99(self):
        m = re.findall(r"(\d+)(?:th|st|nd|rd) percentile", self.SR["zone_script"]); m = [x for x in m if x != "2"]
        assert m, "no percentile in figure script"
        return int(m[0]), f"{N}/figures/zone_modules/fig_zone_modules.py"

    # ------------------------------------------------------------------ release vocabulary and datasets
    def timepoint_aliases(self):
        """Pre-v5 timepoint name -> release name (legacy_aliases of the annotation vocabulary)."""
        return self.vocab()[0]["legacy_aliases"]["timepoint"]

    def n_distinct(self, tag, column):
        """Distinct values of an obs column of a released object (a 10x lane is one obs['sample'])."""
        return len(self.obj(tag)[f"vc_{column}"]), f"{self.objpath(tag)} obs['{column}'] (distinct values)"

    def lanes_per_timepoint(self, tag):
        per_timepoint = {tp: len(samples) for tp, samples in self.obj(tag)["samples_by_tp"].items()}
        assert len(set(per_timepoint.values())) == 1, per_timepoint
        note = ", ".join(f"{tp} {n}" for tp, n in per_timepoint.items())
        return next(iter(per_timepoint.values())), f"{self.objpath(tag)} obs['sample'] by obs['timepoint']", note

    def nuclei_per_timepoint(self, tag, pick):
        """Fewest or most nuclei of one timepoint (pick = min or max); a mouse timepoint is one array."""
        counts = self.obj(tag)["vc_timepoint"]
        note = ", ".join(f"{tp} {n:,}" for tp, n in counts.items())
        return int(pick(counts.values())), f"{self.objpath(tag)} obs['timepoint']", note

    def nuclei_at(self, tag, timepoint):
        return int(self.obj(tag)["vc_timepoint"][timepoint]), f"{self.objpath(tag)} obs['timepoint']"

    def samap_source_human(self):
        """Human fibroblast/decidual nuclei in the joint SAMap graph: the source set of the human-to-mouse transfer and of
        the Ulrich comparison. The graph was fitted before any nucleus was flagged as an embedding outlier and holds every
        nucleus of the fibroblast/decidual object, flagged ones included (checked for release 2026_10_30 on 30 Oct 2026,
        when the Ulrich per-nucleus cosines also covered all of these nuclei)."""
        g = self.f2("samap_graph")
        fibdec = self.obj("human_fibdec")
        note = (f"{g['n_hs']:,} human nuclei in the SAMap graph; the fibroblast/decidual object has {fibdec['n_obs']:,}, "
                f"{fibdec['n_embedding_outlier']} of them flagged and kept as sources")
        return int(g["n_hs"]), f"{g['joint']} obs['species'] == 'hs'", note

    # ------------------------------------------------------------------ v7 epithelium (executed 04c notebooks)
    def _epi_src(self, sp, item):
        return self.F[f"epithelial_{sp}"]["copied"][item]["path"]

    @functools.lru_cache(None)
    def _epi_clusters(self, sp):
        """Per-cluster table the epithelium notebook wrote (labels, votes, shares)."""
        return self.tab(f"epi_{sp}_per_cluster.csv.gz")

    def epi_neighbours(self, sp):
        """Placed neighbours that decide a nucleus's layer or compartment, as the executed notebook ran it."""
        values = re.findall(r"^N_NEIGHBOURS = (\d+)", self.code(f"epi_{sp}_notebook.py"), re.M)
        assert values, "N_NEIGHBOURS is not set in the notebook"
        # A parameters cell injected at run time follows the defaults, so the last assignment is the one that ran.
        return int(values[-1]), self.F[f"epithelial_{sp}"]["notebook"]["path"] + " (N_NEIGHBOURS)"

    def epi_parameter(self, sp, name):
        """A constant the executed epithelium notebook ran with (its last assignment, as for N_NEIGHBOURS)."""
        values = re.findall(rf"^{name} = ([^#\n]+)", self.code(f"epi_{sp}_notebook.py"), re.M)
        assert values, f"{name} is not set in the notebook"
        return ast.literal_eval(values[-1].strip()), self.F[f"epithelial_{sp}"]["notebook"]["path"] + f" ({name})"

    def epi_resolution(self, sp):
        """Leiden resolution of the clusters the notebook labels; its cluster key is leiden_scVI_<resolution>."""
        cluster_key, source = self.epi_parameter(sp, "CLUSTER_KEY")
        return float(cluster_key.rsplit("_", 1)[1]), source

    def epi_min_fold(self, sp):
        """Smallest fold change a DE gene of the notebook must reach (2 ** MIN_LOG2FC)."""
        min_log2fc, source = self.epi_parameter(sp, "MIN_LOG2FC")
        return 2 ** min_log2fc, source

    def epi_label(self, sp, label, what):
        """Clusters or nuclei with one label in the notebook's per-cluster table."""
        clusters = self._epi_clusters(sp)
        of_label = clusters[clusters.label == label]
        value = len(of_label) if what == "clusters" else int(of_label.n_nuclei.sum())
        return value, self._epi_src(sp, "per_cluster")

    def epi_split(self, sp, field):
        """Ties and majority shares over the clusters whose label follows tissue position (inner/outer, luminal/glandular)."""
        first, second = EPI_SPLIT_LABELS[sp]
        votes_first, votes_second = EPI_VOTE_COLUMNS[sp]
        clusters = self._epi_clusters(sp)
        split = clusters[clusters.label.isin([first, second])]

        # A cluster's own votes are those for the layer or compartment its label names.
        own_votes = np.where(split.label == first, split[votes_first], split[votes_second])
        all_votes = (split[votes_first] + split[votes_second]).to_numpy()
        n_placed = int((all_votes + split.n_tie).sum())
        majority_pct = 100 * own_votes / all_votes
        src = self._epi_src(sp, "per_cluster")

        if field == "ties": return int(split.n_tie.sum()), src
        if field == "placed": return n_placed, src
        if field == "tie_frac": return split.n_tie.sum() / n_placed, src
        if field == "majority_median": return float(np.median(majority_pct)), src
        if field == "majority_range": return pct_range(majority_pct, 1), src
        if field == "own_frac": return own_votes.sum() / all_votes.sum(), src
        raise KeyError(field)

    def reinstated(self, field):
        """The human epithelial clusters the release reinstated (results/tables/reinstated_cells_<release>.tsv, n_added > 0)."""
        table = self.tab("epi_human_reinstated.tsv.gz", sep="\t")
        added = table[table.n_added > 0]
        cluster_list = "clusters " + ", ".join(map(str, added.cluster))
        src = self._epi_src("human", "reinstated")
        if field == "clusters": return len(added), src, cluster_list
        if field == "nuclei": return int(added.n_added.sum()), src, cluster_list
        if field == "inner_range":
            pct_inner = self._epi_clusters("human").set_index("cluster").loc[added.cluster, "pct_inner_layer"]
            return pct_range(pct_inner, 1), self._epi_src("human", "per_cluster"), cluster_list
        raise KeyError(field)

    def reinstated_marker(self, gene, how):
        """Share of nuclei detecting a marker in each reinstated cluster: the range over clusters, or the maximum."""
        detection = self.tab("epi_human_reinstated_detection.csv", index_col=0)[gene]
        src = f"{self.F['epithelial_human']['annotated']} layers['counts']"
        note = ", ".join(f"cluster {cluster} {100 * share:.1f}%" for cluster, share in detection.items())
        if how == "range": return pct_range(100 * detection, 0), src, note
        if how == "max": return float(detection.max()), src, note
        raise KeyError(how)

    def epi_de(self, table, column):
        """Genes passing the two-fold, adjusted P < 0.05 rule in an epithelium DE table (count of a boolean column)."""
        file_name, item = {"human": ("epi_human_ranking.csv.gz", "ranking"), "mouse": ("epi_mouse_ranking.csv.gz", "ranking"),
                           "mouse_proliferative": ("epi_mouse_proliferative.csv.gz", "proliferative")}[table]
        values = self.tab(file_name)[column]
        assert values.dtype == bool, (file_name, column, values.dtype)
        return int(values.sum()), self._epi_src(table.split("_")[0], item)

    def epi_donor_detection(self, gene, column):
        """Share of inner (frac_inner) or outer (frac_outer) nuclei detecting a gene in each donor, as a range."""
        by_donor = self.tab("epi_human_by_donor.csv.gz", header=[0, 1], index_col=0)
        per_donor = by_donor.xs(column, axis=1, level=1).loc[gene]
        note = ", ".join(f"{donor} {100 * share:.1f}%" for donor, share in per_donor.items())
        return pct_range(100 * per_donor, 0), self._epi_src("human", "by_donor"), note

    def epi_mouse_stratum(self, labels, stratum):
        """Mouse epithelial nuclei with the given label(s) in uninduced or in induced uteri (annotated lineage object)."""
        labels = [labels] if isinstance(labels, str) else list(labels)
        obs = self.tab("epi_mouse_obs.csv.gz", index_col=0)
        is_uninduced = obs.timepoint.astype(str) == "uninduced"
        in_stratum = is_uninduced if stratum == "uninduced" else ~is_uninduced
        n_nuclei = int((obs.cell_subtype.isin(labels) & in_stratum).sum())
        return n_nuclei, f"{self.F['epithelial_mouse']['annotated']} obs['cell_subtype'] by obs['timepoint']"

    def epi_cycling(self, field):
        """Share of nuclei in S or G2/M per cluster (cell-cycle scoring of the mouse epithelial lineage)."""
        clusters = self._epi_clusters("mouse")
        proliferative = clusters.loc[clusters.label == "epithelial_proliferative", "pct_S_or_G2M"].sort_values()
        split = clusters.loc[clusters.label.isin(EPI_SPLIT_LABELS["mouse"]), "pct_S_or_G2M"]
        assert len(proliferative) == 2, proliferative  # the methods quote one share per proliferative cluster
        src = self._epi_src("mouse", "per_cluster")
        if field == "proliferative_low": return float(proliferative.iloc[0]), src
        if field == "proliferative_high": return float(proliferative.iloc[1]), src
        if field == "split_max": return float(split.max()), src
        raise KeyError(field)

    def epi_lineage_nuclei(self, sp):
        """Nuclei of the epithelial lineage object, every cluster included."""
        n_table = int(self._epi_clusters(sp).n_nuclei.sum())
        n_obs = len(self.tab(f"epi_{sp}_obs.csv.gz", index_col=0))
        assert n_table == n_obs, (n_table, n_obs)
        return n_table, self._epi_src(sp, "per_cluster")

    # ------------------------------------------------------------------ in vitro comparison
    @functools.lru_cache(None)
    def _ivd_tables(self):
        T = {"vote": self.fix("consensus_full_vote_table_stadtmauer_authors_criteria.csv"),
             "hits": self.fix("stadtmauer2022_hits_authors_criteria.csv"),
             "tables": self.fix("stadtmauer_tables_S2_S3_parsed_and_called.csv"),
             "cons": self.fix("consensus123_stadtmauer_authors_criteria.csv"),
             "orth": self.fix("consensus123_stadtmauer_authors_criteria_mouse_orthology.csv"),
             "ag": self.tab("hir_aghajanova2011_invitro_decidualization_nonendo.csv.gz"),
             "ry": self.tab("hir_rytkonen2019_complete_signature.csv.gz")}
        return T

    # ------------------------------------------------------------------ inputs staged with --extra-inputs
    def extra(self, name):
        """Path of an input staged with --extra-inputs."""
        if self.EI is None:
            raise FileNotFoundError(f"{name} is held as a project artifact; stage it and pass --extra-inputs (README)")
        return f"{self.EI}/{name}"

    def extra_src(self, name, member=None):
        entry = self.EIM[name]
        label = f"artifact {name} (version {entry['latest_version_id']}, sha256 {entry['sha256'][:12]}..., staged in {self.EI})"
        return label if member is None else f"{label}: {member}"

    @functools.lru_cache(None)
    def _bundle_member(self, suffix):
        """(member name, bytes) of the one file in the 27 Oct fit bundle whose name ends with suffix."""
        bundle = zipfile.ZipFile(self.extra(FIT_BUNDLE))
        names = [name for name in bundle.namelist() if name.endswith(suffix)]
        assert len(names) == 1, (suffix, names)
        return names[0], bundle.read(names[0])

    @functools.lru_cache(None)
    def _xmens(self):
        """X-Mens decidual-versus-fibroblast DE table of the fit with the uninduced fibroblasts (27 Oct)."""
        name, data = self._bundle_member("xmens_broad_all_genes_with_uninduced.csv")
        return pd.read_csv(io.BytesIO(data), index_col=0, low_memory=False), self.extra_src(FIT_BUNDLE, name)

    @functools.lru_cache(None)
    def _ptnull_table(self):
        return pd.read_csv(PTNULL_TABLE, index_col="statistic")

    def ptnull(self, statistic, column):
        """One cell of the permutation-null statistics: the observed value or a summary over 1,000 permutations."""
        value = self._ptnull_table().loc[statistic, column]
        return float(value), PTNULL_TABLE + f" [{statistic}, {column}]"

    def ptnull_max(self, prefix, column):
        """Largest value of a column over the statistics named prefix* (every rho, or every Moran's I)."""
        table = self._ptnull_table()
        of_prefix = table[table.index.str.startswith(prefix)]
        return float(of_prefix[column].max()), PTNULL_TABLE + f" [max over {', '.join(of_prefix.index)}, {column}]"

    def ivd(self, field):
        T = self._ivd_tables(); v = T["vote"]; S = self.fixsrc
        if field == "n_studies":
            st = sorted({s.strip() for x in T["cons"].mature_papers for s in x.split(";")}); return len(st), S("consensus123_stadtmauer_authors_criteria.csv"), ", ".join(st)
        if field == "aghajanova": return T["ag"].gene.nunique(), f"{HIR}/signatures/aghajanova2011_invitro_decidualization_nonendo.csv"
        if field == "rytkonen": return T["ry"].gene.nunique(), f"{HIR}/signatures/rytkonen2019_complete_signature.csv"
        if field == "stad_tab": return len(T["tables"]), S("stadtmauer_tables_S2_S3_parsed_and_called.csv")
        if field == "stad_met": return int(T["tables"].new_study_call.isin(["up", "down"]).sum()), S("stadtmauer_tables_S2_S3_parsed_and_called.csv")
        if field == "stad_text": return int(T["hits"].source.str.contains("Results text").sum()), S("stadtmauer2022_hits_authors_criteria.csv")
        if field == "stad_total": return T["hits"].gene.nunique(), S("stadtmauer2022_hits_authors_criteria.csv")
        if field == "ryt_conflicts":
            d = T["ry"].groupby("gene").direction.agg(lambda s: s.str.lower().nunique()); return int((d > 1).sum()), f"{HIR}/signatures/rytkonen2019_complete_signature.csv"
        c2 = v[(v.n_up + v.n_down) >= 2]; agree = c2[(c2.n_up == 0) | (c2.n_down == 0)]
        if field == "called2": return len(c2), S("consensus_full_vote_table_stadtmauer_authors_criteria.csv")
        if field == "conflicted": return len(c2) - len(agree), S("consensus_full_vote_table_stadtmauer_authors_criteria.csv")
        cons = T["cons"]; assert len(cons) == len(agree) and set(cons.gene) == set(agree.gene)
        if field == "consensus": return len(cons), S("consensus123_stadtmauer_authors_criteria.csv")
        if field in ("up", "down"): return int((cons.consensus_direction == field).sum()), S("consensus123_stadtmauer_authors_criteria.csv")
        if field in ("sup_ar", "sup_rs", "sup_all"):
            want = {"sup_ar": "Aghajanova; Rytkonen", "sup_rs": "Rytkonen; Stadtmauer", "sup_all": "Aghajanova; Rytkonen; Stadtmauer"}[field]
            return int((cons.mature_papers == want).sum()), S("consensus123_stadtmauer_authors_criteria.csv")
        o = T["orth"]
        if field == "renamed": return int((o.current_human_symbol != o.human_gene).sum()), S("consensus123_stadtmauer_authors_criteria_mouse_orthology.csv")
        if field == "one2one": return int(o.strict_one_to_one.astype(bool).sum()), S("consensus123_stadtmauer_authors_criteria_mouse_orthology.csv")
        if field.startswith("pb_") or field in ("eligible", "resid_df", "uninduced_decidual"):
            return self._ivd_fit(field)
        if field.startswith("xm_"):
            return self._ivd_score(field)
        if field.startswith("w_"):
            return self._wang_sizes(field)
        raise KeyError(field)

    def _pseudobulk_check(self, pseudobulks):
        """Compare each pseudobulk's nuclei with the release's nuclei of that cell type and timepoint."""
        aliases = self.timepoint_aliases()
        cell_types = {"decidual": "decidual_cell", "fibroblast": "endometrial_fibroblast"}
        in_release = self.obj("mouse_fibdec")["celltype_by_tp"]
        differ = []
        for pb in pseudobulks.itertuples():
            key = f"{cell_types[pb.group]}|{aliases.get(pb.timepoint, pb.timepoint)}"
            if in_release.get(key) != pb.n_cells:
                differ.append(f"{key}: {pb.n_cells} vs {in_release.get(key)}")
        if differ:
            return "pseudobulk sizes differ from the release: " + "; ".join(differ)
        return f"every pseudobulk equals the release's nuclei of its cell type and timepoint ({self.objpath('mouse_fibdec')})"

    def _ivd_fit(self, field):
        """The X-Mens fit with the uninduced fibroblasts (27 Oct): fit diagnostics and pseudobulk table of the bundle."""
        diag_name, diag_bytes = self._bundle_member("fit_diagnostics_with_uninduced.json")
        diag = json.loads(diag_bytes)
        dsrc = self.extra_src(FIT_BUNDLE, diag_name)
        pb_name, pb_bytes = self._bundle_member("xmens_broad_metadata_with_uninduced.tsv")
        pseudobulks = pd.read_csv(io.BytesIO(pb_bytes), sep="\t")
        pbsrc = self.extra_src(FIT_BUNDLE, pb_name)

        # The diagnostics and the pseudobulk table describe the same pseudobulks.
        n_nuclei = sum(diag["nuclei_per_arm"].values())
        assert len(pseudobulks) == diag["n_pseudobulks"], (len(pseudobulks), diag["n_pseudobulks"])
        assert int(pseudobulks.n_cells.sum()) == n_nuclei, (int(pseudobulks.n_cells.sum()), n_nuclei)

        if field == "pb_n": return int(diag["n_pseudobulks"]), dsrc
        if field == "pb_nuclei": return int(n_nuclei), dsrc
        if field == "pb_dec": return int(diag["nuclei_per_arm"]["decidual"]), dsrc
        if field == "pb_fib": return int(diag["nuclei_per_arm"]["fibroblast"]), dsrc
        if field == "pb_min": return int(pseudobulks.n_cells.min()), pbsrc, self._pseudobulk_check(pseudobulks)
        if field == "pb_max": return int(pseudobulks.n_cells.max()), pbsrc, self._pseudobulk_check(pseudobulks)
        if field == "resid_df":
            assert diag["residual_df"] == diag["n_pseudobulks"] - diag["design_rank"]
            return int(diag["residual_df"]), dsrc, f"design {diag['design']}, rank {diag['design_rank']}"
        if field == "pb_floor":
            import yaml
            y = yaml.safe_load(open(f"{self.FD2}/files/ivdup__config.yaml"))
            return int(y["xmens_broad"]["cell_floor"]), f"{CAG}/config.yaml: xmens_broad.cell_floor"
        if field == "eligible":
            X_, xsrc = self._xmens()
            n_eligible = int(X_.production_eligible.astype(bool).sum())
            assert n_eligible == diag["eligible_genes"], (n_eligible, diag["eligible_genes"])
            return int(diag["eligible_genes"]), dsrc, "equals production_eligible in the fit's DE table"
        if field == "uninduced_decidual":
            in_release = self.obj("mouse_fibdec")["celltype_by_tp"]
            return int(in_release.get("decidual_cell|uninduced", 0)), f"{self.objpath('mouse_fibdec')} obs['cell_type'] by obs['timepoint']"
        raise KeyError(field)

    @functools.lru_cache(None)
    def _ivd_scored(self):
        """Score the corrected consensus orthologs against the PyDESeq2 X-Mens release table (same rule as run_pipeline)."""
        X_, src = self._xmens(); T = self._ivd_tables(); o = T["orth"]
        if "consensus_direction" not in o:
            o = o.merge(T["cons"][["gene", "consensus_direction"]].rename(columns={"gene": "human_gene"}), on="human_gene", how="left", validate="one_to_one")
        assert o.consensus_direction.notna().all()
        o = o[o.strict_one_to_one.astype(bool)].copy()
        r = X_.reindex(o.mouse_ensembl_gene_id)
        elig = r.production_eligible.fillna(False).astype(bool).to_numpy() & r.estimable.fillna(False).astype(bool).to_numpy()
        call = elig & r.significant_fdr05_lfc1.fillna(False).astype(bool).to_numpy()
        sgn = np.sign(r.log2FoldChange_unshrunk.to_numpy())
        want = o.consensus_direction.map({"up": 1, "down": -1}).to_numpy()
        status = np.where(~elig, "ineligible_or_unavailable", np.where(~call, "tested_not_called", np.where(sgn == want, "supported", "reversed")))
        o["status"] = status; o["xm_dir"] = sgn
        return o, src

    def _ivd_score(self, field):
        o, src = self._ivd_scored(); s = o.status
        cs = self.fix("consensus_ortholog_xmens_broad_support_rows_stadtmauer_authors_criteria.csv").support_status.value_counts().to_dict()
        mine = s.value_counts().to_dict(); note = "" if mine == cs else f"bundle support rows give {cs}"
        src = f"{src} + {self.fixsrc('consensus123_stadtmauer_authors_criteria_mouse_orthology.csv')}"
        called = o[s.isin(["supported", "reversed"])]
        if field == "xm_called": return len(called), src, note
        if field == "xm_conc": return int((s == "supported").sum()), src, note
        if field == "xm_disc": return int((s == "reversed").sum()), src, note
        if field == "xm_notcalled": return int((s == "tested_not_called").sum()), src, note
        if field == "xm_nottest": return int((s == "ineligible_or_unavailable").sum()), src, note
        if field == "xm_pct": return (s == "supported").sum() / len(called), src, note
        if field == "xm_exp":
            a = called.consensus_direction.eq("up").mean(); b = (called.xm_dir > 0).mean()
            return a * b + (1 - a) * (1 - b), src, "pup(A)pup(B) + pdown(A)pdown(B) over the called genes"
        raise KeyError(field)

    def _wang_list(self, k):
        w = self.tab(f"ivdres_wang_{k}.csv.gz")
        return w.drop_duplicates("gene", keep=False).set_index("gene").direction.str.lower().map({"up": 1, "down": -1}), f"{CAG}/resources/wang_{k}.csv"

    def _wang_sizes(self, field):
        k = {"w_preg_n": "pregnancy", "w_oil_n": "oil", "w_ivd_n": "ivd"}.get(field)
        if k:
            w, src = self._wang_list(k); return len(w), src
        if field == "w_ivd_matched":
            rel = self.f2("ivd_release")["kept"]; key = [x for x in rel if x.startswith("wang_ivd_support_rows")][0]
            ws = pd.read_csv(f"{self.FD2}/files/{rel[key]['copied_as']}")
            support_md5 = self.f2("ivd_release", "wang_ivd_support_md5")
            if len(set(support_md5.values())) == 1:
                note = f"wang_ivd_support_rows.csv is identical in all {len(support_md5)} run directories"
            else:
                note = "wang_ivd_support_rows.csv differs between run directories: " + ", ".join(f"{d} {m[:8]}" for d, m in support_md5.items())
            return int(ws.xmens_selected_ensembl_gene_id.notna().sum()), rel[key]["path"], note
        raise KeyError(field)

    @functools.lru_cache(None)
    def _wang_pairs(self):
        X_, xsrc = self._xmens()
        xe = X_[X_.estimable.astype(bool)]
        xe = xe[~xe.gene_symbol.astype(str).duplicated(keep=False)].set_index(xe.gene_symbol.astype(str)[~xe.gene_symbol.astype(str).duplicated(keep=False)])
        xc = xe[xe.significant_fdr05_lfc1.astype(bool)]
        xm = np.sign(xc.log2FoldChange_unshrunk)
        W = {k: self._wang_list(k)[0] for k in ("oil", "pregnancy", "ivd")}
        rel = self.f2("ivd_release")["kept"]; key = [x for x in rel if x.startswith("wang_ivd_support_rows")][0]
        ws = pd.read_csv(f"{self.FD2}/files/{rel[key]['copied_as']}").dropna(subset=["xmens_selected_ensembl_gene_id"])
        g = ws[ws.xmens_selected_ensembl_gene_id.isin(X_.index[X_.significant_fdr05_lfc1.astype(bool) & X_.estimable.astype(bool)])]
        xi_w = g.wang_ivd_direction.str.lower().map({"up": 1, "down": -1}).to_numpy(); xi_x = np.sign(X_.log2FoldChange_unshrunk.reindex(g.xmens_selected_ensembl_gene_id).to_numpy())
        def pair(a, b):
            g_ = a.index.intersection(b.index); return a[g_].to_numpy(), b[g_].to_numpy()
        P = {"xm_oil": pair(xm, W["oil"]), "xm_preg": pair(xm, W["pregnancy"]), "xm_ivd": (xi_x, xi_w),
             "ivd_preg": pair(W["ivd"], W["pregnancy"]), "ivd_oil": pair(W["ivd"], W["oil"]), "oil_preg": pair(W["oil"], W["pregnancy"])}
        return P, xsrc

    def _wang_summary(self, pair):
        """One comparison of the Summary sheet of SuppTable_Wang_comparison (28 Oct supplementary tables)."""
        kept = self.f2("supp_tables")["SuppTable_Wang_comparison__Summary.csv"]
        table = pd.read_csv(f"{self.FD2}/files/{kept['copied_as']}").set_index("Comparison")
        row = table.loc[WANG_SUMMARY_ROWS[pair]]
        n_compared = int(row["Genes compared"])
        up_a = row["List A up among genes compared"] / n_compared
        up_b = row["List B up among genes compared"] / n_compared
        stats = {"same": int(row["Same direction"]), "n": n_compared, "pct": row["Same direction"] / n_compared,
                 "exp": up_a * up_b + (1 - up_a) * (1 - up_b)}
        return stats, f"{kept['path']} (row '{WANG_SUMMARY_ROWS[pair]}')"

    def wang(self, pair, field):
        P, xsrc = self._wang_pairs(); a, b = P[pair]
        same = int((a == b).sum()); n = len(a)
        src = (xsrc + " + " if pair.startswith("xm") else "") + f"{CAG}/resources/wang_*.csv" + (" + wang_ivd_support_rows.csv" if pair == "xm_ivd" else "")
        pa, pb = (a > 0).mean(), (b > 0).mean()
        recomputed = {"same": same, "n": n, "pct": same / n, "exp": pa * pb + (1 - pa) * (1 - pb)}
        summary, ssrc = self._wang_summary(pair)
        agree = (summary["same"], summary["n"]) == (same, n)

        # X-Mens comparisons are read from the supplementary table built on the 27 Oct fit; the registry's own
        # recomputation from that fit's DE table is the check. Comparisons among the Wang lists are recomputed.
        if pair.startswith("xm"):
            note = f"recomputed from the fit's DE table: {same:,} of {n:,}" + (" (same)" if agree else " (DIFFERS)")
            return summary[field], ssrc, note
        note = "" if agree else f"SuppTable_Wang_comparison gives {summary['same']:,} of {summary['n']:,}"
        return recomputed[field], src, note


# ---------------------------------------------------------------------- row functions replaced after construction
def overrides():
    O = {}
    O["samap_edge_total_min"] = lambda c: c.samap_graph("edge_min")
    O["samap_edge_total_max"] = lambda c: c.samap_graph("edge_max")
    O["samap_align_matched"] = lambda c: c.samap_graph("align_matched")
    O["samap_align_off1"] = lambda c: c.samap_graph("align_off1")
    O["samap_align_off2"] = lambda c: c.samap_graph("align_off2")
    O["samap_label_covered"] = lambda c: c.samap_label("covered")
    O["samap_label_cov_pct"] = lambda c: c.samap_label("cov_pct")
    O["samap_label_agree"] = lambda c: c.samap_label("agree")
    O["samap_label_agree_dec"] = lambda c: c.samap_label("agree_dec")
    O["samap_label_agree_fib"] = lambda c: c.samap_label("agree_fib")
    O["xs_heatmap_early"] = lambda c: c.sel_code(r"pt_centroid < (0\.\d+), \"early\"", float)
    O["xs_heatmap_late"] = lambda c: c.sel_code(r"pt_centroid > (0\.\d+), \"late\"", float)
    O["xs_heatmap_bins"] = lambda c: c.sel_code(r"Mh(\d+), Zh\d+ = profiles\(", int)
    O["lineage_knn"] = lambda c: (int(c.cfg["leiden"]["knn"]), f"{rr.CFG}: leiden.knn (read by workflow/scripts/atlas/reintegrate.py)")
    O["cellbender_epochs"] = lambda c: (int(re.search(r"default=(\d+)", c.f2("misc_small", "cellbender_epochs")).group(1)),
                                        f"{L}/.conda/envs/cellbender_022 cellbender/remove_background/argparser.py (--epochs default; 01_cellbender.smk passes none)")
    O["mouse_annotation"] = lambda c: c.mouse_annotation()
    O["samap_proteome_release"] = lambda c: c.proteome_release()
    O["mi_sweep_null_lo"] = lambda c: c.mi_null_pct(0)
    O["mi_sweep_null_hi"] = lambda c: c.mi_null_pct(1)
    return O


def _mouse_annotation(self):
    r = self.f2("mouse_reference")
    j = r.get("reference_json", {}).get("content", {}); txt = json.dumps(j)
    gv = re.findall(r"gencode\.?v?(M\d+)|(vM\d+)", txt, re.I); ens = re.findall(r"(?:release[-_ ]|Ensembl ?)(\d{2,3})", txt, re.I)
    gv = sorted(set((a or b).lstrip("v") for a, b in gv)); ens = sorted(set(ens))
    if not gv: return None, r.get("reference_json", {}).get("path", "reference.json not found"), f"no GENCODE version in reference.json: {txt[:200]}"
    return f"GENCODE v{gv[0]}" + (f" (Ensembl {ens[0]})" if ens else ""), r["reference_json"]["path"], ("" if len(gv) == 1 else f"several versions: {gv}")


def _proteome_release(self):
    m = self.f2("misc_small"); txt = m.get("proteome_release_grep", "") + "\n" + m.get("proteomes", "")
    rel = sorted(set(re.findall(r"release[-_ ](\d{2,3})", txt)))
    if len(rel) != 1: return None, f"{L}/data/proteomes", f"release not recorded unambiguously: {rel or 'none found'}"
    return int(rel[0]), f"{L}/data/proteomes (README / download record)"


def _mi_null_pct(self, i):
    m = self.f2("misc_small"); g = m.get("mi_null_pct", "")
    v = re.findall(r"(?:percentile|quantile)\([^)]*?\[?\s*([0-9.]+)\s*,\s*([0-9.]+)", g)
    if not v: return None, "; ".join(m.get("overlap_scripts", [])) or "overlap scripts not found", "no 2.5/97.5 percentile call found in the overlap code"
    lo, hi = float(v[0][0]), float(v[0][1])
    if hi <= 1: lo, hi = 100 * lo, 100 * hi
    return (lo, hi)[i], "; ".join(m.get("overlap_scripts", [])[:2])


Ctx2.mouse_annotation = _mouse_annotation
Ctx2.proteome_release = _proteome_release
Ctx2.mi_null_pct = _mi_null_pct


# ---------------------------------------------------------------------- inputs held as project artifacts (not in the repository)
def _ai(self, name):
    return f"{self.AI}/{name}"


def _ai_src(self, name):
    m = self.AIM[name]
    return f"artifact {name} (version {m['latest_version_id']}; outside the repository)"


def _heat(self, what):
    if what == "panel_genes":
        d = pd.read_csv(self.ai("gene_peak_order.tsv"), sep="\t"); return len(d), self.ai_src("gene_peak_order.tsv")
    if what == "dispersion_genes":
        d = pd.read_csv(self.ai("disp_table_mouse.tsv.gz"), sep="\t"); assert (d.n_detected > 1).all()
        return len(d), self.ai_src("disp_table_mouse.tsv.gz"), "genes detected in more than one nucleus (n_detected > 1 for every row)"
    R = open(self.ai("monocle_heatmap_repro.R")).read(); src = self.ai_src("monocle_heatmap_repro.R")
    pat = {"spline_df": r"sm\.ns\(Pseudotime, df ?= ?(\d+)\)", "points": r"length\.out ?= ?(\d+)"}[what]
    m = sorted(set(re.findall(pat, R))); assert len(m) == 1, (what, m)
    return int(m[0]), src


def _s16(self, what):
    from scipy.stats import hypergeom
    S = pd.read_csv(self.ai("figS16_revision_adopted_cutoff_sweep.csv")); src = self.ai_src("figS16_revision_adopted_cutoff_sweep.csv")
    ov = self.mi()["ov"]; r = S[S.percentile == 90].iloc[0]
    tie = (int(r.N), int(r.K_mouse), int(r.n_human), int(r.overlap)) == (ov["N"], ov["K"], ov["n"], ov["k"])
    note = "its 90th-percentile row equals the recomputation from results/mi_dpt/pt_harmonized_*_nall/mi_all.csv" if tie else "WARNING: its 90th-percentile row differs from the current MI runs"
    if what == "min": return int(S.percentile.min()), src, note
    if what == "max": return int(S.percentile.max()), src, note
    q, col = {"null_lo": (0.025, "null_fold_lo"), "null_hi": (0.975, "null_fold_hi")}[what]
    v = hypergeom.ppf(q, S.N, S.K_mouse, S.n_human) / S.expected
    assert np.allclose(v, S[col], atol=1e-6), "null interval is not the hypergeometric 2.5th/97.5th percentile"
    return 100 * q, src, f"{col} equals hypergeom.ppf({q}) / expected at all {len(S)} cutoffs"


def _cellbender_epochs(self):
    m = self.f2("misc_small"); g = m["cellbender_epochs"]
    v = int(re.search(r"default=(\d+)", g).group(1))
    rule = open(f"{self.FD2}/files/{m['cellbender_rule']['copied_as']}").read()
    assert not re.search(r"--epochs\b", rule.replace("# epochs", "")), "01_cellbender.smk passes --epochs"
    return v, f"{m['cellbender_argparse']} (--epochs default; workflow/rules/01_cellbender.smk passes no --epochs)"


def _mouse_annotation2(self):
    r = self.f2("mouse_reference"); h = r.get("gtf_header")
    if not h: return None, r.get("reference_json", {}).get("path", ""), "genes.gtf header not collected"
    txt = "\n".join(h["lines"]); v = re.search(r"version (M\d+) \(Ensembl (\d+)\)", txt); prov = re.search(r"##provider: (\S+)", txt)
    return f"{prov.group(1)} v{v.group(1)} (Ensembl {v.group(2)})", h["path"] + " (header)"


def _proteome_release2(self):
    t = open(self.PCHK).read()
    loc = dict((b, a) for a, b in re.findall(r"^(\d{4}-\d\d-\d\d \d\d:\d\d) (\S+\.pep\.all\.fa)$", t, re.M))
    rel = re.findall(r"^## release (\d+) (\w+): (\d{4}-\d\d-\d\d \d\d:\d\d) UTC \((\S+)\.gz", t, re.M)
    hits = sorted({r for r, sp, ts, fn in rel if loc.get(fn) == ts})
    full = [r for r in hits if all(loc.get(fn) == ts for r2, sp, ts, fn in rel if r2 == r)]
    if len(full) != 1: return None, "data/proteomes/fasta", f"no single release matches the file times: {hits}"
    return int(full[0]), "/n/lab_storage/eddy_lab/users/nhilgert/data/proteomes/fasta/*.pep.all.fa (file times vs Ensembl FTP listings; check_proteome_release.sh)", \
        "inferred from timestamps: both local files carry the release-%s FTP times to the minute; no download record exists" % full[0]


Ctx2.ai = _ai; Ctx2.ai_src = _ai_src; Ctx2.heat = _heat; Ctx2.s16 = _s16
Ctx2.cellbender_epochs = _cellbender_epochs; Ctx2.mouse_annotation = _mouse_annotation2; Ctx2.proteome_release = _proteome_release2
_ov0 = overrides


def overrides():
    O = _ov0()
    O.update({"heatmap_panel_genes": lambda c: c.heat("panel_genes"), "heatmap_spline_df": lambda c: c.heat("spline_df"),
              "heatmap_genes_dispersion": lambda c: c.heat("dispersion_genes"), "heatmap_points": lambda c: c.heat("points"),
              "mi_sweep_min_pct": lambda c: c.s16("min"), "mi_sweep_max_pct": lambda c: c.s16("max"),
              "mi_sweep_null_lo": lambda c: c.s16("null_lo"), "mi_sweep_null_hi": lambda c: c.s16("null_hi"),
              "cellbender_epochs": lambda c: c.cellbender_epochs(), "mouse_annotation": lambda c: c.mouse_annotation(),
              "samap_proteome_release": lambda c: c.proteome_release()})
    return O


# ---------------------------------------------------------------------- annotation vocabulary (current subtype codes)
def _vocab(self):
    import yaml
    v = self.SR.get("vocab")
    if not isinstance(v, dict): raise KeyError("annotation vocabulary not collected (small_reads.py: vocab)")
    return yaml.safe_load(v["text"]), v["path"]


def _subtype_aliases(self):
    """Superseded cell_subtype name -> current code, from the collected annotation vocabulary.

    The Oct 22 objects and the pseudotime provenance store pre-v4 names (efib_outer, decidual_apoptotic, ...), while
    the methods use the codes (eFib1, decAp, ...), so subtype names read from them are mapped to codes first.
    """
    return self.vocab()[0]["legacy_aliases"]["cell_subtype"]


def _subtype_counts(self, tag):
    """obs['cell_subtype'] counts of one object, keyed by current code."""
    aliases, counts = self.subtype_aliases(), {}
    for name, n in self.obj(tag)["vc_cell_subtype"].items():
        code = aliases.get(name, name)
        counts[code] = counts.get(code, 0) + n
    return counts


def _root_label(self, sp):
    """Cell subtype of the pseudotime root, as its current code."""
    stored = self.F["pt_provenance"][sp]["root_label"]
    code = self.subtype_aliases().get(stored, stored)
    note = f"the provenance stores the pre-v4 name {stored}" if code != stored else ""
    return code, f"{O22}/provenance/{sp}_pca_released_provenance.json", note


Ctx2.vocab = _vocab; Ctx2.subtype_aliases = _subtype_aliases; Ctx2.subtype_counts = _subtype_counts; Ctx2.root_label = _root_label
