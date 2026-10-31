"""Read-only collection of small text/CSV/JSON sources for the methods-number registry. Usage: python small_reads.py OUT.json"""
import os, sys, glob, json, subprocess, time


def sh(c, t=120):
    try:
        return subprocess.run(["bash", "-c", c], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=t).stdout
    except subprocess.TimeoutExpired:
        return f"TIMEOUT: {c[:120]}"


L = "/n/lab_storage/eddy_lab/users/nhilgert"; X = f"{L}/xmens_final"; N = "/n/netscratch/eddy_lab/Everyone/nhilgert"
R = f"{L}/decidualization/revisions/scenic_v2_2026_08_28"; S = f"{N}/figures/samap_supp"
out = {}


def safe(k, fn):
    try:
        out[k] = fn()
    except Exception as e:
        out[k] = f"ERROR {type(e).__name__}: {e}"
    print(k, "done", flush=True)


safe("scenic_seed_dirs", lambda: {sp: sorted(os.listdir(f"{R}/scenic_seeds/{sp}"))[:80] for sp in ("mouse", "human")})
safe("loom_json_find", lambda: sh(f"find -L {R}/scenic_seeds {R}/inputs -maxdepth 3 -name '*.json' -size -300k 2>/dev/null | xargs grep -l 'n_tfs_retained' 2>/dev/null | head -10").split())


def loom_json():
    d = {}
    for p in out.get("loom_json_find", [])[:8]:
        try: d[p] = json.load(open(p))
        except Exception as e: d[p] = str(e)
    return d
safe("loom_json", loom_json)
safe("loom_logs", lambda: sh(f"grep -rhE '^genes: |TF list|species=.* frac=' {R}/slurm_logs {R}/logs 2>/dev/null | sort | uniq -c | head -12"))


def occ6():
    d = {}
    for sp in ("mouse", "human"):
        dd = f"{R}/scenic_seeds/{sp}/consensus_sketch15k_weighted_soft_occ6"
        d[sp] = {"listing": sorted(os.listdir(dd)) if os.path.isdir(dd) else "missing"}
        for f in sorted(glob.glob(dd + "/*summary*"))[:3]:
            d[sp][f] = open(f, errors="replace").read()[:2000]
    return d
safe("occ6", occ6)
def occ6_kept():
    d = {}
    for sp in ("mouse", "human"):
        p = f"{R}/scenic_seeds/{sp}/consensus_sketch15k_weighted_soft_occ6/recurrence_report.tsv"
        rows = [l.rstrip("\n").split("\t") for l in open(p)]
        h = rows[0]; i_name, i_kept = h.index("name"), h.index("kept")
        d[sp] = {"path": p, "kept": sorted(r[i_name] for r in rows[1:] if r[i_kept] in ("True", "true", "1"))}
    return d
safe("occ6_kept", occ6_kept)
def vocab():
    p = f"{X}/resources/annotation_vocabulary.yaml"
    return {"path": p, "mtime": time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p))), "text": open(p).read()}
safe("vocab", vocab)
safe("occ6_logs", lambda: sh(f"grep -hE 'reg_thresh|kept regulons|regulon names seen|seeds:|AUCell' {R}/consensus_variant_logs/*occ6* 2>/dev/null | head -16"))
safe("pt_runs_log", lambda: sh(f"grep -rhE 'stromal latent attached|PCA|pca' {N}/2026_10_22_harmonized_pt/runs --include=*.log --include=*.out --include=*.txt 2>/dev/null | sort | uniq -c | head -12"))
safe("pt_runs_listing", lambda: sh(f"ls -R {N}/2026_10_22_harmonized_pt/runs 2>/dev/null | head -40"))
safe("pt_pca_code", lambda: sh(f"grep -nE 'sc\\.pp\\.pca|n_comps|use_rep|neighbors|graph ==' {N}/2026_10_22_flagged_harmonized/provenance/harmonized_pseudotime.py"))
safe("samap_env", lambda: sh(f"for m in {L}/.uv/envs/samap_harmony_uv/lib/python3*/site-packages/*samap*.dist-info/METADATA {X}/envs/samap/.venv/lib/python3*/site-packages/*samap*.dist-info/METADATA; do echo \"$m\"; grep -E '^(Name|Version):' $m; done 2>/dev/null"))
safe("samap_settings", lambda: sh(f"grep -nE 'crossK|cross_k|NH1|NH2|n_iter|NUMITERS|neigh|k=|K=' {L}/code/samap-scvi/full_fibdec/run_arm.py {L}/code/samap-scvi/full_fibdec/samap_scvi.py {X}/workflow/rules/10b_samap_raw.smk 2>/dev/null | cut -c1-220 | head -30"))
safe("samap_fitlog", lambda: sh(f"grep -hE 'sc-samap|samap.__version__|SAMap version|version' {L}/code/samap-scvi/full_fibdec/logs/*.out 2>/dev/null | sort | uniq -c | head -10"))
safe("samap_supp_files", lambda: sorted(os.listdir(S)))
safe("samap_supp_exporters", lambda: {os.path.basename(p): open(p).read() for p in sorted(glob.glob(f"{S}/*.py")) if os.path.getsize(p) < 60000})
safe("proteomes", lambda: sh(f"ls -la {L}/data/proteomes {L}/data/proteomes/fasta 2>/dev/null | head -24; for f in {L}/data/proteomes/README* {L}/data/proteomes/*.md {L}/data/proteomes/*.txt {L}/data/proteomes/*.json; do test -f $f && (echo \"## $f\"; head -40 $f); done; grep -rhoE 'release-[0-9]+' {L}/data/proteomes --include=*.md --include=*.txt --include=*.json --include=*.sh --include=*.log 2>/dev/null | sort | uniq -c | head"))
safe("heatmap_dir", lambda: sh(f"ls -la {N}/figures/pt_gene_heatmaps | head -60"))
def decision_tables():
    """Text of each annotation decision table, keyed by file name as before, as release v9 uses it: the round-2 tables
    of config annotation_round2 for the atlas ({species}_lineage.csv) and the round-2 lineages, the July tables under
    resources/annotations/ for every other step (the rule of collect_facts.py). The registry takes the atlas EXCLUDE
    clusters, and so the atlas nuclei counts and doublet means, from these texts."""
    import yaml
    paths = {os.path.basename(p): p for p in sorted(glob.glob(f"{X}/resources/annotations/*.csv"))}
    round2 = yaml.safe_load(open(f"{X}/config/config.yaml"))["annotation_round2"]
    for sp, table in round2["atlas_decisions"].items():
        paths[f"{sp}_lineage.csv"] = f"{X}/{table}"
    for sp, lineages in round2["lineages"].items():
        for lineage in lineages:
            paths[f"{sp}_{lineage}.csv"] = f"{X}/{round2['decisions_dir']}/{sp}_{lineage}.csv"
    out["annot_paths"] = paths
    return {name: open(p).read() for name, p in paths.items() if os.path.getsize(p) < 200000}


safe("annot", decision_tables)
safe("tf_lists", lambda: {t: [l.strip() for l in open(f"{R}/databases/{t}") if l.strip()] for t in ("allTFs_mm.txt", "allTFs_hg38.txt")})
safe("mi_code", lambda: {os.path.basename(p): open(p).read() for p in sorted(glob.glob(f"{X}/workflow/scripts/mi/*.py")) if os.path.getsize(p) < 80000})
safe("mi_rule", lambda: sh(f"ls {X}/workflow/rules | grep -iE 'mi|pseudo'; cat {X}/workflow/rules/11_mi_pseudotime.smk 2>/dev/null | head -80"))
safe("ivd_code", lambda: open(f"{L}/cagri_ivd_end_to_end_20260815/build_xmens_v2_tables.py").read())
safe("hir_code", lambda: open(f"{L}/human_ivd_repro/build_consensus.py").read())
safe("hir_derive", lambda: open(f"{L}/cagri_ivd_end_to_end_20260815/derive_consensus123_orthology.py").read())
safe("zone_cfg", lambda: sh(f"grep -nE 'zone|auc' {X}/config/config.yaml | head -30; ls {X}/workflow/rules | grep -iE 'zone|16_'; ls {N}/figures/zone_modules | head -30"))
safe("zone_script", lambda: sh(f"grep -nE 'auc_threshold|aucell|percentile|p99|99' {N}/figures/zone_modules/fig_zone_modules.py {N}/figures/zone_modules/score_zone_aucell.py 2>/dev/null | cut -c1-200 | head -30"))
json.dump(out, open(sys.argv[1], "w"), default=str)
print("wrote", sys.argv[1], flush=True)
