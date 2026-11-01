"""Review kit for annotation decisions: per-cluster QC, marker statistics and per-nucleus checks.

One decision step = one clustering plus the table that labels its clusters: the two atlas tables
(whole-tissue Leiden 2.5) and one table per lineage (Leiden 3.0). For each step the kit computes one row
per cluster and the marker statistics of the review dot plot (review_dotplot.py), so every decision -
kept, relabelled or removed - is shown with the same evidence.

Definitions (as in the 29 Oct per-nucleus check):
  detected            a gene with >= 1 count in a nucleus (layers['counts'])
  expression          mean over a cluster's nuclei of log(1 + 10,000 x counts / total counts of the nucleus)
  expresses a panel   >= 2 of the lineage panel's genes detected
  two lineages        the nucleus expresses two or more lineage panels (the doublet signature)

Fixed flag rules for a cluster, against the kept clusters of the same step:
  mixed          >= 15% of nuclei express two lineages and >= 3x the kept median
  other lineage  (lineage steps only) a lineage panel that kept clusters rarely express (kept median < 15%)
                 is expressed by >= 40% of the cluster's nuclei and >= 4x the kept median, or one of its genes
                 is detected in >= 50% of the cluster's nuclei against < 10% (and 5x less) in kept clusters
  doublet score  mean DoubletDetection score >= 5 and >= 5x the kept clusters' median
  doublets       >= 20% of nuclei called doublets (rarely fires: called doublets were removed at QC)
  low RNA        median detected genes < half the kept median
  small          fewer than 20 nuclei
  label          another label's mean profile fits better than the cluster's own (r >= 0.85, ahead by >= 0.05)
An EXCLUDE cluster is "supported" when it is mixed, of another lineage, doublet-like or low-RNA; otherwise it
is a reinstatement candidate (profile r >= 0.85 with a kept label) or needs a look.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

# Lineage panels for the mixing measure (26 Oct review): fibroblast and decidual are one stromal continuum,
# smooth muscle and mural cells one panel.
LINEAGE_PANELS = {
    "mouse": {"epithelial": ["Cdh1", "Pax8", "Epcam", "Krt8", "Krt18"],
              "stromal": ["Pdgfra", "Col1a1", "Col3a1", "Dcn", "Lum", "Prl8a2", "Wnt4", "Bmp2", "Hand2"],
              "endothelial": ["Pecam1", "Cdh5", "Flt1", "Emcn", "Kdr", "Cldn5"],
              "smc_mural": ["Acta2", "Myh11", "Des", "Cnn1", "Pdgfrb", "Rgs5", "Kcnj8", "Notch3"],
              "immune": ["Ptprc", "Cd74", "Lyz2", "Csf1r", "Cd52"],
              "mesothelial": ["Msln", "Upk3b", "Wt1", "Lrrn4", "Upk1b"]},
    "human": {"epithelial": ["CDH1", "PAX8", "EPCAM", "KRT8", "KRT18"],
              "stromal": ["PDGFRA", "DCN", "LUM", "COL1A1", "COL1A2", "IGFBP1", "PRL", "LEFTY2"],
              "endothelial": ["PECAM1", "VWF", "CDH5", "FLT1", "CLDN5"],
              "smc_mural": ["ACTA2", "MYH11", "CNN1", "DES", "PDGFRB", "RGS5", "NOTCH3"],
              "immune": ["PTPRC", "CD74", "LYZ", "CD52"]},
}
# The lineage panel(s) a step's own cells express; never "other lineage" in that step.
OWN_PANELS = {"immune": {"immune"}, "vascular": {"endothelial"}, "perivascular": {"endothelial"}, "smc": {"smc_mural"},
              "epithelial": {"epithelial"}, "stromal": {"stromal"}, "perimetrium": {"mesothelial"}}
UNLABELLED = ("?", "to annotate")   # provisional rows: no label yet, so not a reference
MITO_PREFIX = {"mouse": "mt-", "human": "MT-"}
RIBO_PREFIX = {"mouse": ("Rps", "Rpl"), "human": ("RPS", "RPL")}
TARGET_SUM = 1e4
K_NEIGHBOURS = 15
N_PROFILE_GENES = 2000
SMALL_CLUSTER = 20


def cluster_order(names):
    return sorted(names, key=lambda s: (0, int(s)) if str(s).isdigit() else (1, str(s)))


def counts_of(adata):
    """Raw counts of the counts layer. X is log-normalized in the lineage objects, so it is never used instead."""
    if "counts" not in adata.layers:
        raise KeyError("object has no layers['counts']; detection, UMIs and CP10k need raw counts")
    return sp.csr_matrix(adata.layers["counts"]).astype(np.float32)


def nucleus_qc(adata, species):
    """Per-nucleus QC and lineage-panel calls, from the raw counts."""
    counts = counts_of(adata)
    genes = pd.Index(adata.var_names.astype(str))
    detected = (counts > 0).tocsr()
    umis = np.asarray(counts.sum(axis=1)).ravel()
    qc = pd.DataFrame(index=adata.obs_names)
    qc["n_genes"] = np.asarray(detected.sum(axis=1)).ravel()
    qc["umis"] = umis
    qc["pct_mito"] = np.asarray(counts[:, genes.str.startswith(MITO_PREFIX[species])].sum(axis=1)).ravel() / np.maximum(umis, 1) * 100
    qc["pct_ribo"] = np.asarray(counts[:, genes.str.startswith(RIBO_PREFIX[species])].sum(axis=1)).ravel() / np.maximum(umis, 1) * 100
    qc["doublet_score"] = adata.obs["doubletdetection_doublet_score"].astype(float).to_numpy()
    qc["predicted_doublet"] = adata.obs["doubletdetection_predicted_doublet"].astype(float).to_numpy()
    qc["sample"] = adata.obs["sample"].astype(str).to_numpy()
    qc["timepoint"] = adata.obs["timepoint"].astype(str).to_numpy() if "timepoint" in adata.obs else ""
    for panel, panel_genes in LINEAGE_PANELS[species].items():
        columns = [genes.get_loc(g) for g in panel_genes if g in genes]
        need = min(2, len(columns))
        qc[f"expr_{panel}"] = (np.asarray(detected[:, columns].sum(axis=1)).ravel() >= need) if columns else False
    qc["n_lineages"] = qc[[c for c in qc if c.startswith("expr_")]].sum(axis=1)
    return qc, counts, detected, genes


def cluster_table(adata, cluster_col, labels, species, latent_key=None, lineage_step=True):
    """One row per cluster: QC, lineage mixing, and identity (neighbours, best-fitting label, top genes).

    labels: dict cluster -> label ("EXCLUDE" for removed clusters); clusters missing from it get "?".
    lineage_step: the lineage's name (e.g. "immune"), or False for the atlas, where kept clusters span every
    lineage (no "other lineage" rule).
    """
    qc, counts, detected, genes = nucleus_qc(adata, species)
    cluster = adata.obs[cluster_col].astype(str).to_numpy()
    names = cluster_order(set(cluster))
    label_of = {n: labels.get(n, "?") for n in names}

    # Mean log1p CP10k profile and detection per cluster.
    umis = qc["umis"].to_numpy()
    norm = counts.multiply(1.0 / np.maximum(umis, 1)[:, None]).tocsr() * TARGET_SUM
    norm.data = np.log1p(norm.data)
    member = sp.csr_matrix((np.ones(len(cluster)), (pd.Index(names).get_indexer(cluster), np.arange(len(cluster)))),
                           shape=(len(names), len(cluster)))
    sizes = np.asarray(member.sum(axis=1)).ravel()
    mean_profile = np.asarray((member @ norm).todense()) / sizes[:, None]
    pct_detected = np.asarray((member @ detected.astype(np.float32)).todense()) / sizes[:, None] * 100
    kept = np.array([label_of[n] not in ("EXCLUDE",) + UNLABELLED for n in names])
    informative = pct_detected.max(axis=0) >= 10
    variance = np.where(informative, mean_profile[kept].var(axis=0) if kept.any() else 0, -1)
    profile_genes = np.argsort(variance)[::-1][:N_PROFILE_GENES]

    latent_key = latent_key or next(k for k in adata.obsm if "scvi" in k.lower() and "umap" not in k.lower())
    latent = np.asarray(adata.obsm[latent_key])
    _, neighbours = NearestNeighbors(n_neighbors=K_NEIGHBOURS + 1).fit(latent).kneighbors(latent)
    neighbour_cluster = cluster[neighbours[:, 1:]]

    rows = []
    for i, name in enumerate(names):
        m = cluster == name
        q = qc[m]
        timepoints = q["timepoint"].value_counts(normalize=True)
        samples = q["sample"].value_counts(normalize=True)
        row = dict(cluster=name, label=label_of[name], excluded=label_of[name] == "EXCLUDE", n_nuclei=int(m.sum()),
                   median_genes=float(q.n_genes.median()), median_umis=float(q.umis.median()),
                   median_pct_mito=float(q.pct_mito.median()), mean_pct_mito=float(q.pct_mito.mean()),
                   median_pct_ribo=float(q.pct_ribo.median()),
                   mean_doublet=float(q.doublet_score.mean()), pct_predicted_doublet=100 * float(q.predicted_doublet.mean()),
                   pct_two_lineages=100 * float((q.n_lineages >= 2).mean()),
                   top_sample=samples.index[0], top_sample_share=float(samples.iloc[0]), n_samples=int(q["sample"].nunique()),
                   timepoints="; ".join(f"{t} {v:.0%}" for t, v in timepoints.head(4).items()))
        for panel, panel_genes in LINEAGE_PANELS[species].items():
            row[f"pct_{panel}"] = 100 * float(q[f"expr_{panel}"].mean())
            present = [genes.get_loc(g) for g in panel_genes if g in genes]
            best = max(present, key=lambda j: pct_detected[i, j]) if present else None
            row[f"top_gene_{panel}"] = genes[best] if best is not None else ""
            row[f"top_gene_pct_{panel}"] = float(pct_detected[i, best]) if best is not None else 0.0

        # Neighbours in the step's latent: own-cluster share and the labels of the others.
        nbr = neighbour_cluster[m].ravel()
        others = pd.Series([label_of.get(c, "?") for c in nbr[nbr != name]]).value_counts() / max(nbr.size, 1)
        row["frac_self"] = float((nbr == name).mean())
        row["neighbour_labels"] = "; ".join(f"{k} {v:.2f}" for k, v in others.head(3).items())

        # Best-fitting label: correlation with each label's centroid, this cluster left out.
        r_by_label = {}
        for lab in sorted({label_of[n] for j, n in enumerate(names) if kept[j] and n != name}):
            idx = [j for j, n in enumerate(names) if n != name and label_of[n] == lab]
            centroid = np.average(mean_profile[idx][:, profile_genes], axis=0, weights=sizes[idx])
            r_by_label[lab] = float(np.corrcoef(mean_profile[i, profile_genes], centroid)[0, 1])
        ranked = sorted(r_by_label.items(), key=lambda kv: -kv[1])
        row["best_label"], row["best_r"] = ranked[0] if ranked else ("", np.nan)
        row["second_label"], row["second_r"] = ranked[1] if len(ranked) > 1 else ("", np.nan)
        row["own_label_r"] = r_by_label.get(label_of[name], np.nan)

        # Top genes by detection difference against the rest of the step's nuclei.
        rest = ((pct_detected * sizes[:, None]).sum(axis=0) - pct_detected[i] * sizes[i]) / max(sizes.sum() - sizes[i], 1)
        diff = np.where(pct_detected[i] >= 25, pct_detected[i] - rest, -np.inf)
        top = [j for j in np.argsort(diff)[::-1][:10] if np.isfinite(diff[j])]
        row["top_genes"] = ", ".join(f"{genes[j]} {pct_detected[i, j]:.0f}/{rest[j]:.0f}%" for j in top)
        rows.append(row)

    table = pd.DataFrame(rows)
    panel_genes = sorted({g for genes_ in LINEAGE_PANELS[species].values() for g in genes_ if g in genes})
    detection = pd.DataFrame(pct_detected[:, [genes.get_loc(g) for g in panel_genes]], columns=panel_genes, index=table.index)
    return add_flags(table, qc, cluster, species, lineage_step, detection), qc


def add_flags(table, qc, cluster, species, lineage_step=True, detection=None):
    """Fixed flag rules (module docstring), against the kept clusters of this step.

    detection: % of each cluster's nuclei detecting each lineage-panel gene (rows aligned with table).
    """
    kept = table[~table.excluded & ~table.label.isin(UNLABELLED)]
    own = OWN_PANELS.get(lineage_step, set()) if lineage_step else set()
    kept_nuclei = np.isin(cluster, kept.cluster)
    ref_genes = float(qc.n_genes[kept_nuclei].median()) if kept_nuclei.any() else np.nan
    ref_two = float(kept.pct_two_lineages.median()) if len(kept) else np.nan
    ref_doublet = float(kept.mean_doublet.median()) if len(kept) else np.nan
    ref_panel = {p: float(kept[f"pct_{p}"].median()) if len(kept) else np.nan for p in LINEAGE_PANELS[species]}
    kept_detection = {}   # kept-cluster median detection of each panel gene, for the single-marker rule
    for p in LINEAGE_PANELS[species]:
        for gene in set(table[f"top_gene_{p}"]) - {""}:
            kept_detection[gene] = float(detection[gene][kept.index].median()) if len(kept) and detection is not None else np.nan
    flags, verdicts = [], []
    for _, row in table.iterrows():
        mixed = row.pct_two_lineages >= 15 and row.pct_two_lineages >= 3 * ref_two
        foreign = [f"{p} {row[f'pct_{p}']:.0f}% (kept {ref:.0f}%)" for p, ref in ref_panel.items()
                   if lineage_step and p not in own and ref < 15 and row[f"pct_{p}"] >= 40 and row[f"pct_{p}"] >= 4 * ref]
        for p in set(LINEAGE_PANELS[species]) - own:
            gene, pct = row[f"top_gene_{p}"], row[f"top_gene_pct_{p}"]
            kept_pct = kept_detection.get(gene, np.nan)
            if (lineage_step and gene and pct >= 50 and kept_pct < 10 and pct >= 5 * max(kept_pct, 1e-9)
                    and not any(f.startswith(p + " ") for f in foreign)):
                foreign.append(f"{p} ({gene} {pct:.0f}% vs kept {kept_pct:.0f}%)")
        # Listed by lineage name: the loop above walks a set, whose order depends on PYTHONHASHSEED.
        foreign.sort(key=lambda entry: entry.split(" ", 1)[0])
        doublet_score = row.mean_doublet >= 5 and row.mean_doublet >= 5 * ref_doublet
        doublets = row.pct_predicted_doublet >= 20
        low_rna = row.median_genes < 0.5 * ref_genes
        small = row.n_nuclei < SMALL_CLUSTER
        other_label = (not row.excluded and row.label not in UNLABELLED and row.best_r >= 0.85 and row.best_label != row.label
                       and (np.isnan(row.own_label_r) or row.best_r - row.own_label_r >= 0.05))
        evidence = [("mixed", mixed), ("other lineage: " + ", ".join(foreign), bool(foreign)),
                    (f"doublet score {row.mean_doublet:.1f} (kept {ref_doublet:.1f})", doublet_score),
                    ("doublets", doublets), ("low RNA", low_rna)]
        found = [name for name, hit in evidence + [("small", small), (f"fits {row.best_label}", other_label)] if hit]
        flags.append("; ".join(found))
        if row.excluded:
            supporting = [name for name, hit in evidence if hit]
            if supporting:
                verdicts.append("exclusion supported: " + "; ".join(supporting))
            elif row.best_r >= 0.85:
                verdicts.append(f"reinstatement candidate: {row.best_label} (r {row.best_r:.2f})")
            else:
                verdicts.append(f"needs a look: no kept label fits (best {row.best_label} r {row.best_r:.2f})")
        else:
            verdicts.append("flag: " + "; ".join(found) if found else "fine")
    table["flags"], table["verdict"] = flags, verdicts
    table["ref_kept_median_genes"], table["ref_kept_median_two_lineages"] = ref_genes, ref_two
    table["ref_kept_median_doublet"] = ref_doublet
    return table


def marker_stats(adata, cluster_col, marker_genes):
    """Per (cluster, gene): % of nuclei detecting the gene and mean log1p CP10k (library over all genes)."""
    counts = counts_of(adata)
    genes = pd.Index(adata.var_names.astype(str))
    present = [g for g in dict.fromkeys(marker_genes) if g in genes]
    library = np.asarray(counts.sum(axis=1)).ravel()
    selected = counts[:, genes.get_indexer(present)].toarray()
    log_cp10k = np.log1p(selected / np.maximum(library, 1)[:, None] * TARGET_SUM)
    cluster = adata.obs[cluster_col].astype(str).to_numpy()
    rows = []
    for name in cluster_order(set(cluster)):
        m = cluster == name
        pct = 100.0 * (selected[m] > 0).mean(axis=0)
        mean = log_cp10k[m].mean(axis=0)
        rows += [dict(cluster=name, gene=g, pct_detected=float(p), mean_log1p_cp10k=float(e)) for g, p, e in zip(present, pct, mean)]
    missing = [g for g in marker_genes if g not in genes]
    return pd.DataFrame(rows), missing
