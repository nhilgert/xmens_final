"""Write 04c_human_epithelial_layers.ipynb, the notebook that annotates the human epithelial lineage by tissue layer.

The notebook belongs in notebooks/public/ of the repository, where the render_notebook rule
(workflow/rules/99_reports.smk) executes it with papermill from the repository root. The builder writes it to the
current directory.
"""
import nbformat

NOTEBOOK_NAME = "04c_human_epithelial_layers.ipynb"

INTRO = """\
# Human epithelium: inner and outer layers

Reads the unannotated human epithelial lineage object (Leiden clusters at resolution 3.0) and writes it back with an
annotation per nucleus: `cell_subtype` is the label of its cluster in the decisions table, and `cell_type` is
`epithelial_cell` except in excluded clusters.

The decisions table records the label of every cluster: `epithelial_inner`, `epithelial_outer`, `epithelial_ciliated`
or `EXCLUDE`. This notebook shows how the inner and outer labels follow from tissue position, and checks them.

- A spatially placed nucleus is in the **inner** layer when decidual cells outnumber endometrial fibroblasts among its
  20 nearest placed nuclei in the same section, and in the **outer** layer when fibroblasts outnumber decidual cells;
  nuclei with equal counts are not used. A cluster's layer is the layer that holds most of its placed nuclei.
- Every `epithelial_inner` and `epithelial_outer` cluster must lie in the layer its label names; otherwise the
  notebook stops and lists the clusters that disagree.
- `epithelial_ciliated` and `EXCLUDE` clusters keep their decision.

Clusters 11, 22 and 44 co-detect the endothelial marker PECAM1 but express the epithelial markers PAX8 and PAEP. They
were reinstated as inner epithelium; their notes in the decisions table give the marker fractions.

Two descriptions follow that change no label: agreement with the two main expression groups of the same embedding at
resolution 0.5, and the genes that distinguish the layers (Wilcoxon rank-sum test within each donor).

Writes the annotated object, a per-cluster table, the differential-expression tables and three figures."""

PARAMETERS = """\
# Paths are relative to REPO, the repository root; an absolute path is used as given.
REPO = "."
EPITHELIAL_PATH = "data/interim/latent_human_epithelial.h5ad"  # unannotated lineage object: counts, UMAP, clusters
SPATIAL_PATH = "data/interim/qc_spatial_human.h5ad"  # QC-placed coordinates of every human nucleus
# Cell types of the neighbours; the stromal labels are the same in later releases.
TISSUE_PATH = "data/deliverables/2026_10_28/Cevrim_Human_All_Cells.h5ad"
DECISIONS_PATH = "resources/annotations/human_epithelial.csv"  # label of every cluster
OUT_DIR = "results/annotation/human_epithelial_layers"
OUT_OBJECT_PATH = "results/annotation/human_epithelial_layers/annotated_human_epithelial.h5ad"

CLUSTER_KEY = "leiden_scVI_3.0"  # resolution the annotation is made at
CHECK_KEY = "leiden_scVI_0.5"  # two-group split of the same embedding, used only as a check
N_NEIGHBOURS = 20  # placed nuclei around an epithelial nucleus that decide its layer"""

IMPORTS = """\
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

# Joining onto REPO leaves an absolute path unchanged, so any path parameter may point outside the repository.
repo = Path(REPO)
out_dir = repo / OUT_DIR
out_dir.mkdir(parents=True, exist_ok=True)"""

LOAD = """\
# Labels of the decisions table. Inner and outer are checked against tissue position; the other two are kept as decided.
LAYER_LABELS = ["epithelial_inner", "epithelial_outer"]
LABELS = LAYER_LABELS + ["epithelial_ciliated", "EXCLUDE"]


def read_obs(path: Path, columns: list[str]) -> pd.DataFrame:
    \"\"\"obs columns of an h5ad file, without loading its matrices.\"\"\"
    adata_backed = ad.read_h5ad(path, backed="r")
    obs = adata_backed.obs[columns].copy()
    adata_backed.file.close()
    return obs


adata = ad.read_h5ad(repo / EPITHELIAL_PATH)
cluster = adata.obs[CLUSTER_KEY].astype(str).astype(int)

# One decision per cluster, the table's cluster sizes must be this object's, and every label must be a known one.
decisions = pd.read_csv(repo / DECISIONS_PATH, index_col="cluster")
duplicated = sorted(decisions.index[decisions.index.duplicated()])
assert not duplicated, f"clusters with more than one decision: {duplicated}"
unmatched = sorted(set(cluster) ^ set(decisions.index))
assert not unmatched, f"clusters in only one of the object and the decisions table: {unmatched}"
size_in_object = cluster.value_counts().reindex(decisions.index)
assert (decisions.n_cells == size_in_object).all(), "decision table cluster sizes differ from the object"
unknown_labels = sorted(set(decisions.label) - set(LABELS))
assert not unknown_labels, f"decision table labels outside {LABELS}: {unknown_labels}"

spatial = read_obs(repo / SPATIAL_PATH, ["sample", "timepoint", "x_um_qc", "y_um_qc"])
tissue_cell_type = read_obs(repo / TISSUE_PATH, ["cell_type"]).cell_type.astype(str)

print(f"{adata.n_obs:,} epithelial lineage nuclei in {cluster.nunique()} clusters")
print(decisions.label.value_counts().to_string())"""

NUCLEUS_LAYER_MD = """\
## Layer of each placed nucleus

Neighbours are drawn from every QC-placed nucleus of the same section, whatever its cell type. Nuclei outside the
released annotation (for example the excluded stromal-mix cluster) count as neither decidual nor fibroblast."""

NUCLEUS_LAYER = """\
placed = spatial[spatial.x_um_qc.notna()].join(tissue_cell_type)

nuclei = pd.DataFrame({"cluster": cluster, "timepoint": adata.obs.timepoint.astype(str)})
nuclei = nuclei.join(placed[["sample", "x_um_qc", "y_um_qc"]])

neighbour_counts = []
for sample, section in placed.groupby("sample", observed=True):
    query = nuclei[nuclei["sample"] == sample]
    tree = cKDTree(section[["x_um_qc", "y_um_qc"]].to_numpy())
    _, idx = tree.query(query[["x_um_qc", "y_um_qc"]].to_numpy(), k=N_NEIGHBOURS + 1)

    # The nearest hit is the nucleus itself; the next N_NEIGHBOURS are its neighbours.
    assert (section.index.to_numpy()[idx[:, 0]] == query.index.to_numpy()).all()
    neighbour_types = section.cell_type.to_numpy()[idx[:, 1:]]

    neighbour_counts.append(pd.DataFrame(
        {
            "n_decidual": (neighbour_types == "decidual_cell").sum(axis=1),
            "n_fibroblast": (neighbour_types == "endometrial_fibroblast").sum(axis=1),
        },
        index=query.index,
    ))
nuclei = nuclei.join(pd.concat(neighbour_counts))

nuclei["layer"] = np.select(
    [
        nuclei.x_um_qc.isna(),
        nuclei.n_decidual > nuclei.n_fibroblast,
        nuclei.n_fibroblast > nuclei.n_decidual,
    ],
    ["not placed", "inner", "outer"],
    default="tie",
)
has_layer = nuclei.layer.isin(["inner", "outer"])  # placed and not tied

print(nuclei.layer.value_counts().to_string())"""

CLUSTER_LAYER_MD = """\
## Layer of each cluster

A cluster's layer by position is the layer that holds most of its placed nuclei. For every `epithelial_inner` and
`epithelial_outer` cluster it must be the layer that the cluster's label in the decisions table names."""

CLUSTER_LAYER = """\
def count_per_cluster(is_counted: pd.Series) -> pd.Series:
    return nuclei[is_counted].groupby("cluster").size().reindex(decisions.index, fill_value=0)


clusters = pd.DataFrame({
    "label": decisions.label,
    "n_nuclei": size_in_object,
    "n_inner": count_per_cluster(nuclei.layer == "inner"),
    "n_outer": count_per_cluster(nuclei.layer == "outer"),
    "n_tie": count_per_cluster(nuclei.layer == "tie"),
})
clusters["pct_inner_layer"] = 100 * clusters.n_inner / (clusters.n_inner + clusters.n_outer)

# Only inner and outer clusters are checked against position; ciliated and excluded clusters keep their decision.
is_layer_cluster = clusters.label.isin(LAYER_LABELS)
has_no_votes = is_layer_cluster & (clusters.n_inner + clusters.n_outer == 0)
assert not has_no_votes.any(), f"clusters without untied placed nuclei: {list(clusters.index[has_no_votes])}"
is_tied = is_layer_cluster & (clusters.pct_inner_layer == 50)
assert not is_tied.any(), f"clusters tied between the layers: {list(clusters.index[is_tied])}"

# The layer by position must be the layer of the committed label; the error lists every cluster where it is not.
layer_by_position = np.where(clusters.pct_inner_layer > 50, "epithelial_inner", "epithelial_outer")
clusters["label_by_position"] = pd.Series(layer_by_position, index=clusters.index).where(is_layer_cluster)
disagrees = is_layer_cluster & (clusters.label_by_position != clusters.label)
assert not disagrees.any(), (
    "layer by position differs from the decisions table for these clusters:\\n"
    + clusters.loc[disagrees, ["label", "label_by_position", "n_inner", "n_outer", "pct_inner_layer"]].to_string()
)

# How clear-cut the split is: share of a cluster's placed nuclei in its own layer, and the same per nucleus.
layer_clusters = clusters[is_layer_cluster]
majority_pct = np.where(
    layer_clusters.label == "epithelial_inner",
    layer_clusters.pct_inner_layer,
    100 - layer_clusters.pct_inner_layer,
)
nucleus_label = nuclei.cluster.map(clusters.label)
is_in_layer_cluster = nucleus_label.isin(LAYER_LABELS)
own_layer = nucleus_label[is_in_layer_cluster & has_layer].str.removeprefix("epithelial_")
pct_nuclei_in_own_layer = 100 * (own_layer == nuclei.layer[is_in_layer_cluster & has_layer]).mean()
is_placed_in_layer_cluster = is_in_layer_cluster & (nuclei.layer != "not placed")
n_ties = (nuclei.layer[is_placed_in_layer_cluster] == "tie").sum()

print(f"layer by position matches the decisions table for all {is_layer_cluster.sum()} inner and outer clusters")
print(clusters.groupby("label").agg(clusters=("n_nuclei", "size"), nuclei=("n_nuclei", "sum")).to_string())
print(f"nuclei with equal decidual and fibroblast neighbours: {n_ties:,} of {is_placed_in_layer_cluster.sum():,} "
      f"placed in inner and outer clusters")
print(f"cluster majority: median {np.median(majority_pct):.1f}%, "
      f"range {majority_pct.min():.1f}-{majority_pct.max():.1f}%")
print(f"placed nuclei in the layer of their cluster's label: {pct_nuclei_in_own_layer:.1f}%")
clusters.sort_values(["label", "pct_inner_layer"])"""

CHECK_MD = """\
## Check against the two main expression groups at resolution 0.5

At resolution 0.5 the embedding splits into two large groups, plus ciliated and small clusters. The group whose placed
nuclei are more often in the inner layer is the inner group. A cluster agrees when most of its nuclei fall in the group
of its own layer. This check changes no label."""

CHECK = """\
check_group = adata.obs[CHECK_KEY].astype(str)
main_groups = check_group.value_counts().index[:2]

inner_share = (nuclei.layer == "inner").groupby(check_group).sum() / has_layer.groupby(check_group).sum()
inner_group = inner_share[main_groups].idxmax()
print("share of placed nuclei in the inner layer, two main groups:", inner_share[main_groups].round(3).to_dict())

clusters["pct_in_inner_group"] = 100 * (check_group == inner_group).groupby(nuclei.cluster).mean()
clusters["expression_layer"] = np.where(clusters.pct_in_inner_group > 50, "inner", "outer")
agrees = clusters.label.str.removeprefix("epithelial_") == clusters.expression_layer

print(f"{(agrees & is_layer_cluster).sum()} of {is_layer_cluster.sum()} inner and outer clusters agree")
clusters.loc[is_layer_cluster & ~agrees, ["n_nuclei", "pct_inner_layer", "pct_in_inner_group", "label"]]"""

ANNOTATE = """\
assert nuclei.index.equals(adata.obs_names)

adata.obs["cell_type"] = pd.Categorical(np.where(nucleus_label == "EXCLUDE", "EXCLUDE", "epithelial_cell"))
adata.obs["cell_subtype"] = pd.Categorical(nucleus_label)

# The evidence behind each nucleus's layer; NaN for nuclei without tissue coordinates.
adata.obs["layer_n_decidual"] = nuclei.n_decidual.to_numpy()
adata.obs["layer_n_fibroblast"] = nuclei.n_fibroblast.to_numpy()
adata.obs["layer_of_nucleus"] = pd.Categorical(nuclei.layer, categories=["inner", "outer", "tie", "not placed"])

pd.crosstab(adata.obs.cell_subtype, adata.obs.timepoint, margins=True)"""

FIGURES = """\
# Inner and outer take the two ends of a diverging palette; ciliated and excluded nuclei sit outside it.
SUBTYPE_COLORS = {
    "epithelial_outer": "#e66101",
    "epithelial_inner": "#7b3294",
    "epithelial_ciliated": "#1b9e77",
    "EXCLUDE": "#000000",
}
# Background tissue in pale tints of the two compartments the layers are defined by; smooth muscle marks the outer rim.
BACKGROUND_COLORS = {"decidual_cell": "#f7c6c6", "endometrial_fibroblast": "#c9dcf2", "smooth_muscle_cell": "#5a5a5a"}
OTHER_COLOR = "#ebebeb"

# UMAP of the lineage embedding, one colour per label.
umap = adata.obsm["X_umap_scVI"]
fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
for subtype, color in SUBTYPE_COLORS.items():
    is_subtype = (adata.obs.cell_subtype == subtype).to_numpy()
    ax.scatter(umap[is_subtype, 0], umap[is_subtype, 1], s=1, color=color, label=f"{subtype} ({is_subtype.sum():,})",
               rasterized=True)
ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")
ax.legend(markerscale=6, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
fig.savefig(out_dir / "umap_by_subtype.png", dpi=200)

# Tissue maps: the sections of one timepoint share a coordinate frame, so they are drawn together.
timepoints = sorted(adata.obs.timepoint.astype(str).unique())
fig, axes = plt.subplots(1, len(timepoints), figsize=(7 * len(timepoints), 7), layout="constrained")
for ax, timepoint in zip(axes, timepoints):
    sections = placed[placed.timepoint.astype(str) == timepoint]
    background = sections[~sections.index.isin(adata.obs_names)]
    background_colors = background.cell_type.map(BACKGROUND_COLORS).fillna(OTHER_COLOR)
    ax.scatter(background.x_um_qc, background.y_um_qc, s=1, c=background_colors, rasterized=True)

    epithelium = nuclei[(nuclei.timepoint == timepoint) & nuclei.x_um_qc.notna()]
    for subtype, color in SUBTYPE_COLORS.items():
        is_subtype = nucleus_label[epithelium.index] == subtype
        ax.scatter(epithelium.x_um_qc[is_subtype], epithelium.y_um_qc[is_subtype], s=4, color=color, rasterized=True)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(f"{timepoint}: {sections['sample'].nunique()} sections overlaid")
# One legend for both kinds of dots: the epithelial labels, then the tissue around them.
legend_colors = SUBTYPE_COLORS | BACKGROUND_COLORS | {"other or not in release": OTHER_COLOR}
handles = [Line2D([], [], marker="o", linestyle="", color=color, label=label) for label, color in legend_colors.items()]
axes[-1].legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
fig.savefig(out_dir / "tissue_maps_by_subtype.png", dpi=200)"""

DE_MD = """\
## Genes that distinguish the layers

Inner and outer nuclei are compared within each donor: two-sided Wilcoxon rank-sum test with tie correction and
Benjamini-Hochberg correction, on counts normalised to 10,000 per nucleus and log-transformed. A gene is listed as
higher in a layer when it is at least two-fold higher with adjusted P < 0.05 in both donors. Genes are ranked by the
smaller of the two donors' differences in the fraction of nuclei that detect them."""

DE = """\
TARGET_SUM = 1e4  # counts per nucleus after normalisation
MIN_LOG2FC = 1  # two-fold
MAX_PVAL_ADJ = 0.05

adata_expr = ad.AnnData(
    adata.layers["counts"],
    obs=adata.obs[["timepoint", "cell_subtype"]].copy(),
    var=pd.DataFrame(index=adata.var_names),
)
adata_expr = adata_expr[adata_expr.obs.cell_subtype != "EXCLUDE"].copy()
sc.pp.normalize_total(adata_expr, target_sum=TARGET_SUM)
sc.pp.log1p(adata_expr)

de_by_donor = {}
for donor in timepoints:
    is_donor_layered = (adata_expr.obs.timepoint == donor) & adata_expr.obs.cell_subtype.isin(
        ["epithelial_inner", "epithelial_outer"]
    )
    adata_donor = adata_expr[is_donor_layered].copy()
    adata_donor.obs["cell_subtype"] = adata_donor.obs.cell_subtype.astype("category")
    # A gene undetected in this donor has all values tied: scanpy divides 0 by 0 (the RuntimeWarning below) and
    # gives it score 0 and P = 1, so it is tested but never significant.
    sc.tl.rank_genes_groups(
        adata_donor,
        "cell_subtype",
        groups=["epithelial_inner"],
        reference="epithelial_outer",
        method="wilcoxon",
        tie_correct=True,
    )
    test = sc.get.rank_genes_groups_df(adata_donor, group="epithelial_inner").set_index("names")

    is_inner = (adata_donor.obs.cell_subtype == "epithelial_inner").to_numpy()
    is_detected = adata_donor.X > 0
    detection = pd.DataFrame(
        {
            "frac_inner": np.asarray(is_detected[is_inner].mean(axis=0)).ravel(),
            "frac_outer": np.asarray(is_detected[~is_inner].mean(axis=0)).ravel(),
        },
        index=adata_donor.var_names,
    )
    de_by_donor[donor] = (
        test[["logfoldchanges", "pvals_adj"]]
        .rename(columns={"logfoldchanges": "log2fc", "pvals_adj": "pval_adj"})
        .join(detection)
    )
    print(f"{donor}: {is_inner.sum():,} inner and {(~is_inner).sum():,} outer nuclei")

de = pd.concat({donor: table for donor, table in de_by_donor.items()}, axis=1)
log2fc = de.xs("log2fc", axis=1, level=1)
pval_adj = de.xs("pval_adj", axis=1, level=1)
detection_diff = de.xs("frac_inner", axis=1, level=1) - de.xs("frac_outer", axis=1, level=1)

is_significant = (pval_adj < MAX_PVAL_ADJ).all(axis=1)
higher_in_inner = is_significant & (log2fc >= MIN_LOG2FC).all(axis=1)
higher_in_outer = is_significant & (log2fc <= -MIN_LOG2FC).all(axis=1)
print(f"at least two-fold, adjusted P < {MAX_PVAL_ADJ} in both donors: {higher_in_inner.sum()} genes higher in inner, "
      f"{higher_in_outer.sum()} higher in outer nuclei")

ranking = pd.DataFrame({
    "min_detection_diff_inner": detection_diff.min(axis=1),
    "min_detection_diff_outer": (-detection_diff).min(axis=1),
    "higher_in_inner": higher_in_inner,
    "higher_in_outer": higher_in_outer,
})
top_inner = ranking.nlargest(10, "min_detection_diff_inner")
top_outer = ranking.nlargest(5, "min_detection_diff_outer")
pd.concat([top_inner, top_outer]).round(3)"""

DOTPLOT_MD = """\
Dot plot by cluster: the inner and outer genes named in the Methods, then epithelial and endothelial markers. Colour
is mean expression scaled per gene; dot size is the fraction of nuclei in which the gene is detected."""

DOTPLOT = """\
# Clusters in label order, inner clusters from most to least inner.
LABEL_ORDER = ["epithelial_inner", "epithelial_outer", "epithelial_ciliated"]
shown = clusters[clusters.label.isin(LABEL_ORDER)].copy()
shown["label_rank"] = shown.label.map(LABEL_ORDER.index)
shown = shown.sort_values(["label_rank", "pct_inner_layer"], ascending=[True, False])
cluster_names = {k: f"c{k} {label.removeprefix('epithelial_')}" for k, label in shown.label.items()}

adata_expr.obs["cluster"] = pd.Categorical(
    cluster[adata_expr.obs_names].map(cluster_names),
    categories=list(cluster_names.values()),
)
# The inner and outer genes named in the Methods, each higher in its layer in both donors; then epithelial and
# endothelial markers, which show that the PECAM1+ clusters 11, 22 and 44 are epithelium.
marker_genes = {
    "inner": ["DST", "TANC2", "ABCA1", "SVIL", "SOX4", "SKIL", "RDH10", "MAP2"],
    "outer": ["CTNNA2", "ST3GAL5"],
    "epithelial": ["EPCAM", "CDH1", "PAX8"],
    "endothelial": ["PECAM1", "FLT1", "VWF"],
}
assert ranking.higher_in_inner[marker_genes["inner"]].all()
assert ranking.higher_in_outer[marker_genes["outer"]].all()

dotplot = sc.pl.dotplot(adata_expr, marker_genes, groupby="cluster", standard_scale="var", return_fig=True)
dotplot.savefig(out_dir / "dotplot_layer_markers.png", dpi=200)"""

WRITE = """\
clusters.to_csv(out_dir / "human_epithelial_layers_per_cluster.csv")
de.to_csv(out_dir / "de_inner_vs_outer_by_donor.csv")
ranking.to_csv(out_dir / "de_inner_vs_outer_ranking.csv")

# Annotated objects are written once, never over an existing file.
out_object_path = repo / OUT_OBJECT_PATH
assert not out_object_path.exists(), f"{out_object_path} already exists"
adata.write_h5ad(out_object_path)
print("wrote", out_object_path)
print(sorted(path.name for path in out_dir.iterdir()))"""


def build() -> nbformat.NotebookNode:
    md = nbformat.v4.new_markdown_cell
    code = nbformat.v4.new_code_cell
    cells = [
        md(INTRO),
        # papermill injects REPO and any overrides directly after the cell tagged "parameters".
        code(PARAMETERS, metadata={"tags": ["parameters"]}),
        code(IMPORTS),
        md("## Inputs"), code(LOAD),
        md(NUCLEUS_LAYER_MD), code(NUCLEUS_LAYER),
        md(CLUSTER_LAYER_MD), code(CLUSTER_LAYER),
        md(CHECK_MD), code(CHECK),
        md("## Annotation"), code(ANNOTATE),
        md("## Figures"), code(FIGURES),
        md(DE_MD), code(DE),
        md(DOTPLOT_MD), code(DOTPLOT),
        md("## Outputs"), code(WRITE),
    ]

    # Ids from cell order rather than nbformat's random ones, so an unchanged rebuild gives an identical file.
    for position, cell in enumerate(cells):
        cell.id = f"cell-{position:02d}"

    notebook = nbformat.v4.new_notebook(cells=cells)
    notebook.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    return notebook


nbformat.write(build(), NOTEBOOK_NAME)
