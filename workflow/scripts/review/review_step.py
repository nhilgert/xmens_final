"""Review one annotation decision step: cluster table, marker statistics and dot plot, saved together.

    from review_step import review_step
    result = review_step(species="mouse", step="smc", object_path="data/interim/latent_mouse_smc.h5ad",
                         cluster_col="leiden_scVI_3.0", decisions="resources/annotations/mouse_smc.csv",
                         out_dir="results/annotation_review", vocab_path="resources/annotation_vocabulary.yaml",
                         panels_path="resources/review/panels_<tag>.json")

Writes {species}_{step}_clusters.csv (one row per cluster: QC, mixing, identity, flags, verdict),
{species}_{step}_marker_stats.csv (every plotted value) and {species}_{step}_dotplot.{png,svg,pdf}.
Panels come from resources/review/panels_<tag>.json (the decision set's version); display names from the
annotation vocabulary.
Called by the review notebooks (notebooks/review/, rule render_review_notebook in 99b_review_notebooks.smk).
"""
import json
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from annotation_review import cluster_table, marker_stats
from review_dotplot import EXCLUDE_COLOUR, FIGURE_RC, GROUP_PALETTE, render

FLAG_COLOUR = "#9A5B00"     # dark amber: flags on kept clusters
ATLAS_ORDER = {"human": ["stromal", "epithelial", "perivascular", "smc", "immune"],
               "mouse": ["stromal", "epithelial", "vascular", "perimetrium", "smc", "immune"]}
ATLAS_DISPLAY = {"stromal": "stromal (fibroblast/decidual)", "epithelial": "epithelial", "vascular": "vascular (endothelial)",
                 "perivascular": "perivascular (endothelial)", "perimetrium": "perimetrium", "smc": "smooth muscle",
                 "immune": "immune"}


def load_panels(species, step, path):
    spec = json.loads(Path(path).read_text())
    for figure in spec["figures"]:
        if figure["species"] == species and figure["step"] == step:
            return figure["panels"]
    raise KeyError(f"no panels for {species} {step} in {path}")


def read_decisions(decisions):
    """cluster -> label from a decision CSV path, a DataFrame with cluster/label columns, or a dict."""
    if isinstance(decisions, dict):
        return {str(k): str(v) for k, v in decisions.items()}
    table = decisions if isinstance(decisions, pd.DataFrame) else pd.read_csv(decisions, dtype=str)
    return dict(zip(table["cluster"].astype(str).str.strip(), table["label"].astype(str).str.strip()))


def label_order_and_display(species, step, labels, vocab):
    """Kept labels in figure order and the text shown for each (vocabulary display names)."""
    kept = [l for l in dict.fromkeys(labels.values()) if l not in ("EXCLUDE", "?")]
    if step == "atlas":
        order = ATLAS_ORDER[species]
        return sorted(kept, key=lambda l: order.index(l) if l in order else len(order)), ATLAS_DISPLAY
    order = []
    for cell_type, subtypes in vocab["allowed"].items():
        order += list(subtypes or []) + [cell_type]
    names = {**vocab["display_names"]["cell_type"], **vocab["display_names"]["cell_subtype"]}
    return sorted(kept, key=lambda l: order.index(l) if l in order else len(order)), names


def review_step(species, step, object_path, cluster_col, decisions, out_dir, vocab_path, panels_path, title=None,
                subtitle=None, notes=(), reasons=None, reasons_header="exclusion reason [record]", extra_text=None,
                step_name=None, adata=None):
    """Compute and save the review of one step. Returns dict(table, qc, stats, missing_genes, figure, paths).

    step_name: the step as named in the title (default: step); qc: the per-nucleus table of cluster_table.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = adata if adata is not None else ad.read_h5ad(object_path)
    labels = read_decisions(decisions)
    vocab = yaml.safe_load(Path(vocab_path).read_text())
    panels = load_panels(species, step, panels_path)

    table, qc = cluster_table(adata, cluster_col, labels, species, lineage_step=step if step != "atlas" else False)
    genes = [g for panel in panels for group in panel["groups"] for g in group["genes"]]
    stats, missing = marker_stats(adata, cluster_col, genes)
    order, display = label_order_and_display(species, step, labels, vocab)

    text_columns = []
    if extra_text:
        text_columns.append(extra_text)
    checks = {r.cluster: r.verdict for r in table.itertuples() if r.excluded or r.verdict != "fine"}
    text_columns.append(("checks (fixed rules)", checks, FLAG_COLOUR))
    if reasons:
        text_columns.append((reasons_header, {k: v for k, v in reasons.items() if labels.get(k) == "EXCLUDE"},
                             EXCLUDE_COLOUR))

    kept, removed = table[~table.excluded], table[table.excluded]
    level = "Leiden 2.5" if step == "atlas" else "Leiden 3.0"
    title = title or (f"{species.capitalize()} {step_name or step} ({level}): {len(kept)} clusters kept ({kept.n_nuclei.sum():,} nuclei) "
                      f"as {kept.label.nunique()} label{'s' if kept.label.nunique() != 1 else ''}; "
                      f"{len(removed)} removed ({removed.n_nuclei.sum():,} nuclei)")
    subtitle = subtitle or f"Clusters: {object_path} ({cluster_col}); decisions: {decisions if isinstance(decisions, (str, Path)) else 'draft table'}"

    # Fixed SVG ids and no dates or software stamps in the files, so a re-render of the same inputs is byte-identical.
    file_metadata = {"png": {"Software": None}, "svg": {"Date": None}, "pdf": {"CreationDate": None}}
    with mpl.rc_context({**FIGURE_RC, "svg.hashsalt": "review_step"}):
        fig, plotted = render(table, stats, panels, order, display, title, subtitle, notes, text_columns)
        stem = out_dir / f"{species}_{step}_dotplot"
        for suffix in ("png", "svg", "pdf"):
            fig.savefig(f"{stem}.{suffix}", dpi=300 if suffix == "png" else None, metadata=file_metadata[suffix])
    table.to_csv(out_dir / f"{species}_{step}_clusters.csv", index=False)
    plotted.to_csv(out_dir / f"{species}_{step}_marker_stats.csv", index=False)
    return dict(table=table, qc=qc, stats=plotted, missing_genes=missing, figure=fig,
                paths=dict(png=f"{stem}.png", svg=f"{stem}.svg", pdf=f"{stem}.pdf",
                           clusters=str(out_dir / f"{species}_{step}_clusters.csv"),
                           marker_stats=str(out_dir / f"{species}_{step}_marker_stats.csv")))


def draw_decision_umap(adata, cluster_col, labels, label_order, display, embedding_key="X_umap_scVI"):
    """UMAP of the step's own embedding: every nucleus by decision (a) and the removed clusters alone (b).

    label_order and display are those of the dot plot (label_order_and_display), so the label colours match its
    row-group blocks. Returns the figure.
    """
    xy = np.asarray(adata.obsm[embedding_key])[:, :2]
    cluster = adata.obs[cluster_col].astype(str).to_numpy()
    decision = np.array([labels[c] for c in cluster])
    removed = decision == "EXCLUDE"

    # Kept labels take the dot plot's group colours; red stays reserved for EXCLUDE.
    colours = {label: GROUP_PALETTE[i % len(GROUP_PALETTE)] for i, label in enumerate(label_order)}
    colours["EXCLUDE"] = EXCLUDE_COLOUR
    point_size = 1.0 if len(cluster) > 20_000 else 4.0   # smaller points where nuclei overlap heavily

    with mpl.rc_context(FIGURE_RC):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")

        # Panel a: nuclei by decision, cluster numbers at the cluster medians.
        for label in list(label_order) + ["EXCLUDE"]:
            member = decision == label
            if member.any():
                axes[0].scatter(*xy[member].T, s=point_size, color=colours[label], lw=0, rasterized=True,
                                label=display.get(label, label))
        for name in np.unique(cluster):
            x, y = np.median(xy[cluster == name], axis=0)
            axes[0].text(x, y, name, fontsize=6, ha="center", va="center")
        axes[0].legend(markerscale=6, loc="center left", bbox_to_anchor=(1.0, 0.5))
        axes[0].set_title(f"a  decision per nucleus ({cluster_col})", loc="left")

        # Panel b: removed clusters in red over the kept nuclei in grey, so their positions stand out.
        axes[1].scatter(*xy[~removed].T, s=point_size, color="#D0D0D0", lw=0, rasterized=True)
        axes[1].scatter(*xy[removed].T, s=2 * point_size, color=EXCLUDE_COLOUR, lw=0, rasterized=True)
        for name in np.unique(cluster[removed]):
            x, y = np.median(xy[cluster == name], axis=0)
            axes[1].text(x, y, name, fontsize=7, ha="center", va="center", color=EXCLUDE_COLOUR, fontweight="bold")
        axes[1].set_title(f"b  removed clusters ({int(removed.sum()):,} nuclei)", loc="left")

        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(f"{embedding_key} 1")
            ax.set_ylabel(f"{embedding_key} 2")
            ax.set_aspect("equal", adjustable="datalim")
    return fig
