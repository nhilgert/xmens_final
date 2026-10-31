"""Write 04c_mouse_epithelial_compartments.ipynb, which annotates the mouse epithelial lineage by compartment.

The notebook belongs in notebooks/public/ and runs under rule render_notebook (papermill -p REPO <repo> --cwd <repo>),
so its paths are parameters relative to the repository root. It recomputes the luminal and glandular labels from tissue
position, checks the proliferative clusters against cell-cycle scores, tests marker genes, and writes the annotated
lineage object, a per-cluster table, the differential-expression tables and four figures.
"""
import nbformat

NOTEBOOK_PATH = "04c_mouse_epithelial_compartments.ipynb"

INTRO = """\
# Mouse epithelium: luminal and glandular compartments

Reads the unannotated mouse epithelial lineage object (Leiden clusters at resolution 3.0) and writes it back with an
annotation per nucleus in `cell_type` and `cell_subtype`, taken cluster by cluster from the decisions table.

- `epithelial_luminal` and `epithelial_glandular` follow tissue position, recomputed here. A spatially placed nucleus
  votes **luminal** when lumen-side stroma outnumbers deeper fibroblasts among its 20 nearest placed nuclei in the same
  section, and **glandular** for the reverse; nuclei with equal counts cast no vote. Lumen-side stroma is the
  subluminal fibroblasts (eFib3) and the decidual cells (dec1, dec2, dec3, decAp), which fill the subluminal stroma
  after induction; deeper fibroblasts are eFib1, eFib2 and eFib2i. The rule has no tuned threshold. A cluster takes the
  compartment of most of its votes, and that compartment must be its committed label.
- `epithelial_proliferative` and `EXCLUDE` keep their decision.

Three descriptions follow that change no label: the share of nuclei in S or G2/M phase in each cluster, the genes that
distinguish luminal from glandular nuclei (Wilcoxon rank-sum test within uninduced and within induced uteri), and the
genes higher in proliferative nuclei.

Writes the annotated object, a per-cluster table, the differential-expression tables with their ranking, and four
figures."""

PARAMETERS = """\
# Paths are relative to the repository root; the report rule passes the root as REPO.
REPO = "."
EPITHELIAL_PATH = "data/interim/latent_mouse_epithelial.h5ad"  # unannotated lineage object: counts, UMAP, clusters
SPATIAL_PATH = "data/interim/qc_spatial_mouse.h5ad"  # QC-placed coordinates of every mouse nucleus
TISSUE_PATH = "data/deliverables/2026_10_28/Cevrim_XMens_All_Cells.h5ad"  # cell types around each nucleus
VOCAB_PATH = "resources/annotation_vocabulary.yaml"  # release timepoint names and their legacy aliases
DECISIONS_PATH = "resources/annotations/mouse_epithelial.csv"  # committed label of each cluster
OUT_DIR = "results/annotation/mouse_epithelial_compartments"
OUT_OBJECT_PATH = "results/annotation/mouse_epithelial_compartments/annotated_mouse_epithelial.h5ad"

CLUSTER_KEY = "leiden_scVI_3.0"  # resolution the annotation is made at
N_NEIGHBOURS = 20  # placed nuclei around an epithelial nucleus that decide its vote"""

SETUP = """\
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

# Repository-relative paths resolve against REPO; an absolute path, as in a test run, is used as given.
repo = Path(REPO)
out_dir = repo / OUT_DIR
out_object_path = repo / OUT_OBJECT_PATH

# The annotated object is written once and never over an existing file; checked here, before any work is done.
assert not out_object_path.exists(), f"{out_object_path} already exists"
out_dir.mkdir(parents=True, exist_ok=True)"""

LOAD = """\
def read_obs(path: Path, columns: list[str]) -> pd.DataFrame:
    \"\"\"obs columns of an h5ad file, without loading its matrices.\"\"\"
    adata_backed = ad.read_h5ad(path, backed="r")
    obs = adata_backed.obs[columns].copy()
    adata_backed.file.close()
    return obs


adata = ad.read_h5ad(repo / EPITHELIAL_PATH)
cluster = adata.obs[CLUSTER_KEY].astype(str).astype(int)

# One decision per cluster, made at this resolution, and the table's cluster sizes must be this object's.
decisions = pd.read_csv(repo / DECISIONS_PATH, index_col="cluster")
assert decisions.index.is_unique, "a cluster has more than one decision"
unmatched = sorted(set(cluster) ^ set(decisions.index))
assert not unmatched, f"clusters without a decision, or decisions without a cluster: {unmatched}"
assert (decisions.res == CLUSTER_KEY).all(), f"decisions made at a resolution other than {CLUSTER_KEY}"
size_in_object = cluster.value_counts().reindex(decisions.index)
assert (decisions.n_cells == size_in_object).all(), "decision table cluster sizes differ from the object"

spatial = read_obs(repo / SPATIAL_PATH, ["sample", "timepoint", "x_um_qc", "y_um_qc"])
tissue = read_obs(repo / TISSUE_PATH, ["timepoint", "cell_type", "cell_subtype"]).astype(str)
assert adata.obs_names.isin(spatial.index).all(), "lineage nuclei missing from the spatial table"

# The spatial table and the lineage object keep the timepoint names from before vocabulary v5. Its legacy aliases
# rename them to the release names, which must then match the release nucleus by nucleus.
vocabulary = yaml.safe_load((repo / VOCAB_PATH).read_text())
timepoint_aliases = vocabulary["legacy_aliases"]["timepoint"]
release_timepoints = vocabulary["categorical_order"]["timepoint"]["mouse"]

spatial["timepoint"] = spatial.timepoint.astype(str).replace(timepoint_aliases)
lineage_timepoint = adata.obs.timepoint.astype(str).replace(timepoint_aliases)
for source, timepoint in {"spatial table": spatial.timepoint, "lineage object": lineage_timepoint}.items():
    unknown = set(timepoint) - set(release_timepoints)
    assert not unknown, f"{source}: timepoints outside the vocabulary: {unknown}"
    shared = timepoint.index.intersection(tissue.index)
    n_differ = (timepoint[shared] != tissue.timepoint[shared]).sum()
    assert n_differ == 0, f"{source}: {n_differ} nuclei have another timepoint in the release"
    print(f"{source}: timepoint matches the release for all {len(shared):,} nuclei in both")

print(f"{adata.n_obs:,} epithelial lineage nuclei in {cluster.nunique()} clusters")
print(decisions.label.value_counts().to_string())"""

VOTE_MD = """\
## Vote of each placed nucleus

Neighbours are drawn from every QC-placed nucleus of the same section, whatever its cell type; each sample is one
section. Fibroblast and decidual neighbours are named by their release subtype. Nuclei of other cell types, and nuclei
outside the release (such as the excluded clusters), count as neither lumen-side nor deeper."""

VOTE = """\
# Stroma on the lumen side of an epithelial nucleus, and deeper in the endometrium. Decidual cells count as lumen side:
# after induction they fill the subluminal stroma that eFib3 occupies in uninduced uteri.
DECIDUAL = ["dec1", "dec2", "dec3", "decAp"]
LUMEN_SIDE = ["eFib3", *DECIDUAL]
DEEP = ["eFib1", "eFib2", "eFib2i"]
absent = set(LUMEN_SIDE + DEEP) - set(tissue.cell_subtype)
assert not absent, f"stromal subtypes of the rule missing from the release: {absent}"

# Fibroblast and decidual nuclei are named by subtype, the rest by cell type; nuclei outside the release stay NaN.
placed = spatial[spatial.x_um_qc.notna()].join(tissue[["cell_type", "cell_subtype"]])
is_stroma = placed.cell_type.isin(["endometrial_fibroblast", "decidual_cell"])
placed["neighbour_type"] = placed.cell_subtype.where(is_stroma, placed.cell_type)

nuclei = pd.DataFrame({"cluster": cluster, "timepoint": lineage_timepoint})
nuclei = nuclei.join(placed[["sample", "x_um_qc", "y_um_qc"]])

neighbour_counts = []
for sample, section in placed.groupby("sample", observed=True):
    query = nuclei[nuclei["sample"] == sample]
    tree = cKDTree(section[["x_um_qc", "y_um_qc"]].to_numpy())
    _, idx = tree.query(query[["x_um_qc", "y_um_qc"]].to_numpy(), k=N_NEIGHBOURS + 1)

    # The nearest hit is the nucleus itself; the next N_NEIGHBOURS are its neighbours.
    assert (section.index.to_numpy()[idx[:, 0]] == query.index.to_numpy()).all()
    neighbour_types = section.neighbour_type.to_numpy()[idx[:, 1:]]

    neighbour_counts.append(pd.DataFrame(
        {
            "n_lumen_side": np.isin(neighbour_types, LUMEN_SIDE).sum(axis=1),
            "n_deep": np.isin(neighbour_types, DEEP).sum(axis=1),
        },
        index=query.index,
    ))
nuclei = nuclei.join(pd.concat(neighbour_counts))

nuclei["vote"] = np.select(
    [
        nuclei.x_um_qc.isna(),
        nuclei.n_lumen_side > nuclei.n_deep,
        nuclei.n_deep > nuclei.n_lumen_side,
    ],
    ["not placed", "luminal", "glandular"],
    default="tie",
)
has_vote = nuclei.vote.isin(["luminal", "glandular"])  # placed and not tied

print(f"{len(placed):,} QC-placed nuclei in {placed['sample'].nunique()} sections")
print(nuclei.vote.value_counts().to_string())"""

CLUSTER_MD = """\
## Compartment of each cluster

A cluster's compartment by position is the vote of most of its placed nuclei. For every luminal and glandular cluster it
must be the committed label. Proliferative and excluded clusters keep their decision; their majority is shown for
information only."""

CLUSTER = """\
def count_per_cluster(is_counted: pd.Series) -> pd.Series:
    return nuclei[is_counted].groupby("cluster").size().reindex(decisions.index, fill_value=0)


clusters = pd.DataFrame({
    "label": decisions.label,
    "n_nuclei": size_in_object,
    "pct_uninduced": 100 * (nuclei.timepoint == "uninduced").groupby(nuclei.cluster).mean(),
    "n_luminal": count_per_cluster(nuclei.vote == "luminal"),
    "n_glandular": count_per_cluster(nuclei.vote == "glandular"),
    "n_tie": count_per_cluster(nuclei.vote == "tie"),
})
clusters["pct_luminal_votes"] = 100 * clusters.n_luminal / (clusters.n_luminal + clusters.n_glandular)
clusters["position_majority"] = np.select(
    [clusters.pct_luminal_votes > 50, clusters.pct_luminal_votes < 50],
    ["epithelial_luminal", "epithelial_glandular"],
    default="none",  # no votes, or as many luminal as glandular
)

# Luminal and glandular clusters follow position; proliferative and excluded clusters keep their decision.
by_position = clusters.label.isin(["epithelial_luminal", "epithelial_glandular"])
disagree = clusters.index[by_position & (clusters.position_majority != clusters.label)].tolist()
assert not disagree, f"tissue position disagrees with the committed label of clusters {disagree}"

# How clear-cut the split is: the share of a cluster's votes cast for its label, and the same over nuclei.
split = clusters[by_position]
majority_pct = np.where(split.label == "epithelial_luminal", split.pct_luminal_votes, 100 - split.pct_luminal_votes)
nucleus_label = nuclei.cluster.map(clusters.label)
is_in_split = nucleus_label.isin(["epithelial_luminal", "epithelial_glandular"])
own_compartment = nucleus_label[is_in_split & has_vote].str.removeprefix("epithelial_")
pct_votes_for_own = 100 * (own_compartment == nuclei.vote[is_in_split & has_vote]).mean()
is_placed_in_split = is_in_split & (nuclei.vote != "not placed")
n_ties = (nuclei.vote[is_placed_in_split] == "tie").sum()

print(f"tissue position matches the committed label in all {by_position.sum()} luminal and glandular clusters")
print(clusters.groupby("label").agg(clusters=("n_nuclei", "size"), nuclei=("n_nuclei", "sum")).to_string())
print(f"nuclei with as many lumen-side as deeper neighbours: {n_ties:,} of {is_placed_in_split.sum():,} placed")
print(f"cluster majority: median {np.median(majority_pct):.1f}%, "
      f"range {majority_pct.min():.1f}-{majority_pct.max():.1f}%")
print(f"placed nuclei voting for their cluster's compartment: {pct_votes_for_own:.1f}%")
clusters.sort_values(["label", "pct_luminal_votes"], ascending=[True, False]).round(1)"""

CYCLE_MD = """\
## Proliferative clusters: share of nuclei in S or G2/M

Every nucleus of the lineage is scored for S and G2/M with scanpy's `score_genes_cell_cycle`, using the S and G2/M
genes of Tirosh et al. (2016) in mouse capitalisation, on counts normalised to 10,000 per nucleus and log-transformed.
Scores are relative to genes of similar expression in this lineage, so the shares differ from scoring the whole atlas.
Each proliferative cluster must have a larger share of nuclei in S or G2/M than every luminal and glandular cluster.
This check changes no label."""

CYCLE = """\
TARGET_SUM = 1e4  # counts per nucleus after normalisation

# Tirosh et al. 2016, with the symbol updates of Seurat's cc.genes.updated.2019 (Cenpu, Pimreg, Jpt1).
S_GENES = [
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm7", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1",
    "Cenpu", "Hells", "Rfc2", "Polr1b", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2",
    "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1",
    "Chaf1b", "Mrpl36", "E2f8",
]
G2M_GENES = [
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo",
    "Cenpf", "Tacc3", "Pimreg", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b",
    "Gtse1", "Kif20b", "Hjurp", "Cdca3", "Jpt1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", "Dlgap5",
    "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2",
    "G2e3", "Gas2l3", "Cbx5", "Cenpa",
]

# Normalised expression of every lineage nucleus. The counts are copied because normalisation works in place, and the
# annotated object keeps its raw counts.
adata_expr = ad.AnnData(
    adata.layers["counts"].copy(),
    obs=pd.DataFrame({"cluster": cluster, "label": cluster.map(decisions.label), "timepoint": lineage_timepoint}),
    var=pd.DataFrame(index=adata.var_names),
)
sc.pp.normalize_total(adata_expr, target_sum=TARGET_SUM)
sc.pp.log1p(adata_expr)

missing_genes = sorted(set(S_GENES + G2M_GENES) - set(adata_expr.var_names))
assert not missing_genes, f"cell-cycle genes missing from the object: {missing_genes}"
sc.tl.score_genes_cell_cycle(adata_expr, s_genes=S_GENES, g2m_genes=G2M_GENES)

is_cycling = adata_expr.obs.phase.isin(["S", "G2M"])
is_detected = pd.DataFrame(
    adata_expr[:, ["Mki67", "Cenpf"]].X.toarray() > 0,
    index=adata_expr.obs_names,
    columns=["Mki67", "Cenpf"],
)
clusters["pct_S_or_G2M"] = 100 * is_cycling.groupby(cluster).mean()
clusters["pct_Mki67"] = 100 * is_detected.Mki67.groupby(cluster).mean()
clusters["pct_Cenpf"] = 100 * is_detected.Cenpf.groupby(cluster).mean()

is_proliferative = clusters.label == "epithelial_proliferative"
lowest_proliferative = clusters.pct_S_or_G2M[is_proliferative].min()
highest_by_position = clusters.pct_S_or_G2M[by_position].max()
assert lowest_proliferative > highest_by_position, "a luminal or glandular cluster cycles like a proliferative one"

print(f"nuclei in S or G2M: every proliferative cluster at least {lowest_proliferative:.1f}%, "
      f"every luminal and glandular cluster at most {highest_by_position:.1f}%")
cycle_columns = ["label", "n_nuclei", "pct_S_or_G2M", "pct_Mki67", "pct_Cenpf"]
clusters[cycle_columns].sort_values("pct_S_or_G2M", ascending=False).round(1)"""

ANNOTATE_MD = """\
## Annotation

`cell_type` is `epithelial_cell`, or `EXCLUDE` for the excluded clusters; `cell_subtype` is the committed label of the
nucleus's cluster. The neighbour counts and vote of each nucleus are kept beside them."""

ANNOTATE = """\
assert nuclei.index.equals(adata.obs_names)

adata.obs["cell_type"] = pd.Categorical(np.where(nucleus_label == "EXCLUDE", "EXCLUDE", "epithelial_cell"))
adata.obs["cell_subtype"] = pd.Categorical(nucleus_label)

# The evidence behind each nucleus's vote; NaN for nuclei without tissue coordinates.
adata.obs["position_n_lumen_side"] = nuclei.n_lumen_side.to_numpy()
adata.obs["position_n_deep"] = nuclei.n_deep.to_numpy()
adata.obs["position_vote"] = pd.Categorical(nuclei.vote, categories=["luminal", "glandular", "tie", "not placed"])

nuclei_per_label = pd.crosstab(adata.obs.cell_subtype, nuclei.timepoint, margins=True)
nuclei_per_label[[*release_timepoints, "All"]]"""

UMAP_MD = """\
## Figures

UMAP of the lineage embedding, one colour per label."""

UMAP = """\
# One hue per label (ColorBrewer Dark2); excluded nuclei in black.
LABEL_COLORS = {
    "epithelial_luminal": "#d95f02",
    "epithelial_glandular": "#1b9e77",
    "epithelial_proliferative": "#7570b3",
    "EXCLUDE": "#000000",
}

umap = adata.obsm["X_umap_scVI"]
fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
for label, color in LABEL_COLORS.items():
    is_label = (adata.obs.cell_subtype == label).to_numpy()
    ax.scatter(
        umap[is_label, 0],
        umap[is_label, 1],
        s=2,
        color=color,
        linewidths=0,
        label=f"{label} ({is_label.sum():,})",
        rasterized=True,
    )
ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")
ax.legend(markerscale=4, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
fig.savefig(out_dir / "umap_by_label.png", dpi=200)"""

MAPS_MD = """\
Tissue maps, one panel per timepoint. The sections of a timepoint share a coordinate frame, so they are drawn
together. Epithelial nuclei are coloured by label over the stroma the rule counts."""

MAPS = """\
# Stroma in pale tints: warm for the lumen side (eFib3, then decidual cells), blue for the deeper fibroblasts the rule
# compares them with. Smooth muscle in dark grey marks the myometrium at the outer edge.
BACKGROUND = {
    "eFib3, subluminal (lumen side)": (["eFib3"], "#eea3a3"),
    "decidual cells (lumen side)": (DECIDUAL, "#f3cfa6"),
    "eFib1, eFib2, eFib2i (deeper)": (DEEP, "#a9c6e6"),
    "smooth muscle": (["smooth_muscle_cell"], "#636363"),
}
OTHER_COLOR = "#ececec"
SCALE_BAR_UM = 1000  # 1 mm

background_color = pd.Series(OTHER_COLOR, index=placed.index)
for types, color in BACKGROUND.values():
    background_color[placed.neighbour_type.isin(types)] = color

fig, axes = plt.subplots(2, 3, figsize=(15, 10.5), layout="constrained")
for ax, timepoint in zip(axes.flat, release_timepoints):
    sections = placed[placed.timepoint == timepoint]
    background = sections[~sections.index.isin(adata.obs_names)]
    ax.scatter(
        background.x_um_qc,
        background.y_um_qc,
        s=2,
        c=background_color[background.index],
        linewidths=0,
        rasterized=True,
    )

    # Epithelial nuclei on top, larger than the stroma so that thin epithelial layers stay visible.
    epithelium = nuclei[(nuclei.timepoint == timepoint) & nuclei.x_um_qc.notna()]
    for label, color in LABEL_COLORS.items():
        is_label = nucleus_label[epithelium.index] == label
        ax.scatter(
            epithelium.x_um_qc[is_label],
            epithelium.y_um_qc[is_label],
            s=7,
            color=color,
            linewidths=0,
            rasterized=True,
        )

    # 1 mm scale bar at the lower left of every panel, because the panels are drawn at different scales.
    x0, y0 = sections.x_um_qc.min(), sections.y_um_qc.min()
    ax.plot([x0, x0 + SCALE_BAR_UM], [y0, y0], color="black", linewidth=2, solid_capstyle="butt")
    ax.text(x0 + SCALE_BAR_UM / 2, y0, "1 mm", ha="center", va="top")

    ax.set_aspect("equal")
    ax.set_axis_off()
    n_sections = sections["sample"].nunique()
    ax.set_title(f"{timepoint.replace('_', ' ')} ({n_sections} section{'s' if n_sections > 1 else ''})")

# One legend for both kinds of dots: the epithelial labels, then the tissue around them.
legend_colors = LABEL_COLORS | {name: color for name, (_, color) in BACKGROUND.items()}
legend_colors["other, or not in the release"] = OTHER_COLOR
handles = [Line2D([], [], marker="o", linestyle="", color=color, label=label) for label, color in legend_colors.items()]
fig.legend(handles=handles, frameon=False, loc="outside lower center", ncol=5)
fig.savefig(out_dir / "tissue_maps_by_label.png", dpi=200)"""

DE_MD = """\
## Genes that distinguish the compartments

Luminal and glandular nuclei are compared separately in uninduced and in induced uteri (all timepoints after induction
together): two-sided Wilcoxon rank-sum test with tie correction and Benjamini-Hochberg correction, on counts normalised
to 10,000 per nucleus and log-transformed. A gene is listed as higher in a compartment when it is at least two-fold
higher with adjusted P < 0.05 in both strata. Genes are ranked by the smaller of the two strata's differences in the
fraction of nuclei that detect them.

Proliferative nuclei are compared with the other epithelial nuclei of induced uteri, where most of them are; across all
uteri, the genes of induction would dominate the comparison. Excluded clusters take part in neither test."""

DE = """\
MIN_LOG2FC = 1  # two-fold
MAX_PVAL_ADJ = 0.05


def test_markers(adata_subset: ad.AnnData, key: str, group: str, reference: str) -> pd.DataFrame:
    \"\"\"Wilcoxon test of group against reference: log2 fold change, adjusted P and detection fractions per gene.\"\"\"
    adata_subset = adata_subset.copy()
    adata_subset.obs[key] = adata_subset.obs[key].astype("category")

    # A gene undetected in the subset has all values tied: scanpy divides 0 by 0 (the RuntimeWarning below) and gives
    # it score 0 and P = 1, so it is tested but never significant.
    sc.tl.rank_genes_groups(
        adata_subset,
        key,
        groups=[group],
        reference=reference,
        method="wilcoxon",
        tie_correct=True,
    )
    test = sc.get.rank_genes_groups_df(adata_subset, group=group).set_index("names").rename_axis("gene")

    in_group = (adata_subset.obs[key] == group).to_numpy()
    is_detected = adata_subset.X > 0
    detection = pd.DataFrame(
        {
            "frac_group": np.asarray(is_detected[in_group].mean(axis=0)).ravel(),
            "frac_reference": np.asarray(is_detected[~in_group].mean(axis=0)).ravel(),
        },
        index=adata_subset.var_names,
    )
    return (
        test[["logfoldchanges", "pvals_adj"]]
        .rename(columns={"logfoldchanges": "log2fc", "pvals_adj": "pval_adj"})
        .join(detection)
    )


# Luminal against glandular nuclei, separately in uninduced and in induced uteri.
adata_expr.obs["stratum"] = np.where(adata_expr.obs.timepoint == "uninduced", "uninduced", "induced")
is_compartment = adata_expr.obs.label.isin(["epithelial_luminal", "epithelial_glandular"])
de_by_stratum = {}
for stratum in ["uninduced", "induced"]:
    adata_stratum = adata_expr[is_compartment & (adata_expr.obs.stratum == stratum)]
    print(f"{stratum}:", adata_stratum.obs.label.value_counts().to_dict())
    de_by_stratum[stratum] = test_markers(adata_stratum, "label", "epithelial_luminal", "epithelial_glandular")

de_compartments = pd.concat(de_by_stratum, axis=1)
log2fc = de_compartments.xs("log2fc", axis=1, level=1)
pval_adj = de_compartments.xs("pval_adj", axis=1, level=1)
frac_luminal = de_compartments.xs("frac_group", axis=1, level=1)
frac_glandular = de_compartments.xs("frac_reference", axis=1, level=1)
detection_diff = frac_luminal - frac_glandular

is_significant = (pval_adj < MAX_PVAL_ADJ).all(axis=1)
higher_in_luminal = is_significant & (log2fc >= MIN_LOG2FC).all(axis=1)
higher_in_glandular = is_significant & (log2fc <= -MIN_LOG2FC).all(axis=1)
print(f"at least two-fold, adjusted P < {MAX_PVAL_ADJ} in both strata: {higher_in_luminal.sum()} genes higher in "
      f"luminal, {higher_in_glandular.sum()} higher in glandular nuclei")

ranking = pd.DataFrame({
    "min_detection_diff_luminal": detection_diff.min(axis=1),
    "min_detection_diff_glandular": (-detection_diff).min(axis=1),
    "higher_in_luminal": higher_in_luminal,
    "higher_in_glandular": higher_in_glandular,
})
top_luminal = ranking[ranking.higher_in_luminal].nlargest(10, "min_detection_diff_luminal")
top_glandular = ranking[ranking.higher_in_glandular].nlargest(10, "min_detection_diff_glandular")
pd.concat([top_luminal, top_glandular]).round(3)"""

DE_PROLIFERATIVE = """\
# Proliferative against the other epithelial nuclei of induced uteri.
is_proliferative_nucleus = adata_expr.obs.label == "epithelial_proliferative"
n_proliferative_induced = (is_proliferative_nucleus & (adata_expr.obs.stratum == "induced")).sum()
print(f"{n_proliferative_induced:,} of {is_proliferative_nucleus.sum():,} proliferative nuclei are from induced uteri")

is_induced_epithelial = (adata_expr.obs.stratum == "induced") & (adata_expr.obs.label != "EXCLUDE")
adata_induced = adata_expr[is_induced_epithelial].copy()
adata_induced.obs["group"] = np.where(
    adata_induced.obs.label == "epithelial_proliferative",
    "proliferative",
    "other",
)
print("induced:", adata_induced.obs.group.value_counts().to_dict())

de_proliferative = test_markers(adata_induced, "group", "proliferative", "other")
de_proliferative["detection_diff"] = de_proliferative.frac_group - de_proliferative.frac_reference
de_proliferative["higher_in_proliferative"] = (
    (de_proliferative.pval_adj < MAX_PVAL_ADJ) & (de_proliferative.log2fc >= MIN_LOG2FC)
)
print(f"at least two-fold, adjusted P < {MAX_PVAL_ADJ}: {de_proliferative.higher_in_proliferative.sum():,} genes "
      "higher in proliferative nuclei")
de_proliferative[de_proliferative.higher_in_proliferative].nlargest(10, "detection_diff").round(3)"""

DOTPLOT_MD = """\
## Dot plots

Genes, left to right: Wnt7a, then the six top-ranked named genes higher in luminal nuclei; Foxa2, then the six
top-ranked named genes higher in glandular nuclei (predicted genes, Gm, and RIKEN clones, Rik, are skipped); then five
G2/M genes. Wnt7a and Foxa2 pass the same test but rank lower, because their difference in detection between the
compartments is smaller (ranks printed below). Colour is mean expression scaled per gene; dot size is the fraction of
nuclei that detect the gene. Excluded clusters are not shown."""

DOTPLOT_GENES = """\
CANONICAL = {"luminal": "Wnt7a", "glandular": "Foxa2"}  # classical marker of each compartment, shown first
N_TOP_NAMED = 6  # test-ranked genes shown per compartment, after the classical marker
CELL_CYCLE_GENES = ["Mki67", "Top2a", "Cenpf", "Cdk1", "Birc5"]
LABEL_ORDER = ["epithelial_luminal", "epithelial_glandular", "epithelial_proliferative"]

for compartment, canonical in CANONICAL.items():
    is_higher = ranking[f"higher_in_{compartment}"]
    assert is_higher[canonical], f"{canonical} is not higher in {compartment} nuclei"
    rank = ranking.loc[is_higher, f"min_detection_diff_{compartment}"].rank(ascending=False)[canonical]
    print(f"{canonical}: rank {rank:.0f} of the {is_higher.sum()} genes higher in {compartment} nuclei")

# Predicted genes (Gm) and RIKEN clones (Rik) carry no name a reader can interpret.
is_named = ~ranking.index.str.contains(r"^Gm\\d+$|Rik$")
marker_genes = {}
for compartment, canonical in CANONICAL.items():
    is_shown = ranking[f"higher_in_{compartment}"] & is_named & (ranking.index != canonical)
    top = ranking[is_shown].nlargest(N_TOP_NAMED, f"min_detection_diff_{compartment}").index.tolist()
    marker_genes[compartment] = [canonical, *top]
marker_genes["cell cycle"] = CELL_CYCLE_GENES
print(marker_genes)

# Clusters in label order; within a label, from most to fewest luminal votes.
shown = clusters[clusters.label.isin(LABEL_ORDER)].copy()
shown["label_rank"] = shown.label.map(LABEL_ORDER.index)
shown = shown.sort_values(["label_rank", "pct_luminal_votes"], ascending=[True, False])
cluster_names = {k: f"c{k} {label.removeprefix('epithelial_')}" for k, label in shown.label.items()}

adata_shown = adata_expr[adata_expr.obs.label.isin(LABEL_ORDER)].copy()
adata_shown.obs["cluster_name"] = pd.Categorical(
    adata_shown.obs.cluster.map(cluster_names),
    categories=list(cluster_names.values()),
)
adata_shown.obs["label"] = pd.Categorical(adata_shown.obs.label, categories=LABEL_ORDER)"""

DOTPLOT_BY_CLUSTER = """\
dotplot = sc.pl.dotplot(adata_shown, marker_genes, groupby="cluster_name", standard_scale="var", return_fig=True)
dotplot.savefig(out_dir / "dotplot_by_cluster.png", dpi=200)"""

DOTPLOT_BY_LABEL = """\
dotplot = sc.pl.dotplot(adata_shown, marker_genes, groupby="label", standard_scale="var", return_fig=True)
dotplot.savefig(out_dir / "dotplot_by_label.png", dpi=200)"""

WRITE = """\
clusters.to_csv(out_dir / "mouse_epithelial_compartments_per_cluster.csv")
de_compartments.to_csv(out_dir / "de_luminal_vs_glandular_by_stratum.csv")
ranking.to_csv(out_dir / "de_luminal_vs_glandular_ranking.csv")
de_proliferative.to_csv(out_dir / "de_proliferative_vs_other_induced.csv")

adata.write_h5ad(out_object_path)
print("wrote", out_object_path)
print(sorted(path.name for path in out_dir.iterdir()))"""


def build() -> nbformat.NotebookNode:
    md = nbformat.v4.new_markdown_cell
    code = nbformat.v4.new_code_cell
    cells = [
        md(INTRO),
        code(PARAMETERS, metadata={"tags": ["parameters"]}),
        code(SETUP),
        md("## Inputs"), code(LOAD),
        md(VOTE_MD), code(VOTE),
        md(CLUSTER_MD), code(CLUSTER),
        md(CYCLE_MD), code(CYCLE),
        md(ANNOTATE_MD), code(ANNOTATE),
        md(UMAP_MD), code(UMAP),
        md(MAPS_MD), code(MAPS),
        md(DE_MD), code(DE), code(DE_PROLIFERATIVE),
        md(DOTPLOT_MD), code(DOTPLOT_GENES), code(DOTPLOT_BY_CLUSTER), code(DOTPLOT_BY_LABEL),
        md("## Outputs"), code(WRITE),
    ]
    notebook = nbformat.v4.new_notebook(cells=cells)
    notebook.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    return notebook


nbformat.write(build(), NOTEBOOK_PATH)
