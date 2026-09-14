"""
LIANA + COMMOT spatial cell-cell communication pipeline.

Run for a single timepoint + cell_type_key combination.

Usage (single):
    run_timepoint_analysis(adata_full, timepoint='3dpi', cell_type_key='cell_type_annotations')

Usage (loop):
    for tp in timepoints:
        for key in ['cell_type_annotations', 'cell_type_fine']:
            run_timepoint_analysis(adata_full, timepoint=tp, cell_type_key=key)

To choose bandwidth before running:
    choose_bandwidth(adata_spatial, spatial_key='X_spatial_coords')
"""

import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
# import commot as ct  # uncomment to enable Step 4 (COMMOT)
import squidpy as sq
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_helpers import plot_spatial_dark, save_dark_figure
import plotnine as p9
from pathlib import Path
from IPython.display import display
from statsmodels.stats.multitest import multipletests
from sklearn.neighbors import NearestNeighbors


# ============================================================
# CONFIG
# ============================================================

OUTPUT_BASE   = None     # legacy /n/eddy_lab default removed (2026-09-08 port); pass --output-base
SPATIAL_KEY   = 'X_spatial_coords'
BANDWIDTH     = 100      # µm — run choose_bandwidth() to verify for your data
RESOURCE      = 'consensus'  # 'consensus' for human, 'mouseconsensus' for mouse
TOP_N_LR      = None     # optional cap on discovered LR pairs after the specificity filter (--top-n-lr); None = no cap
TOP_N_DOTPLOT = 25       # interactions per dotplot
TOP_N_COMMOT  = 10       # pairs to run COMMOT on
COMMOT_K      = 30       # neighbors for communication_direction
EXPR_PROP     = 0.1
RANK_COLUMN   = 'magnitude_rank'  # L-R discovery column: magnitude_rank (Dimitrov et al. 2024: magnitude rank < 0.05) or specificity_rank
RANK_CUTOFF   = 0.05              # L-R discovery: RANK_COLUMN < cutoff (--rank-column / --rank-cutoff)
SKIP_SPATIAL  = False             # --skip-spatial: stop after rank_aggregate (non-spatial deliverable; Steps 2-4 not run)
MORAN_I_MIN   = 0.01        # Moran's I effect-size floor for spatially variable interactions (--moran-i-min)


# ============================================================
# BANDWIDTH HELPER
# ============================================================

def choose_bandwidth(adata, spatial_key='X_spatial_coords', start=10, end=300, interval_n=15,
                     chosen=None):
    """
    Plot neighbour count vs bandwidth to help choose an appropriate bandwidth.
    Wraps li.ut.query_bandwidth. Run this interactively on a spatial subset
    before committing to a BANDWIDTH value in CONFIG.

    For human decidualization data: ~100 µm works well (targeting ~50 neighbors)
    For mouse XMens data:          ~125 µm (see liana_pipeline_no_commot_mouse.py)
    """
    plot, df = li.ut.query_bandwidth(
        coordinates=adata.obsm[spatial_key],
        start=start,
        end=end,
        interval_n=interval_n
    )
    if chosen is not None:
        plot = (plot
                + p9.geom_vline(xintercept=chosen, linetype='dashed', color='red')
                + p9.annotate('text', x=chosen, y=df.neighbours.max(),
                              label=f'chosen={chosen}µm',
                              angle=90, va='bottom', ha='right', color='red'))
    print(df.to_string())
    return plot


# ============================================================
# HELPERS
# ============================================================

def ensure_log1p_layer(adata, counts_layer='counts', target_sum=1e4):
    """Recompute adata.layers['log1p'] = log1p(normalize_total(counts, target_sum)).

    Provenance note (2026-09-08 Snakemake port): the June 2026 LIANA inputs carried a
    pre-computed 'log1p' layer (median-target normalisation). The 2026_08_28 deposited
    objects carry only layers['counts'] (their X is 1e4-target log1p). Recomputing from
    counts here makes the normalisation explicit and identical across species/variants.
    Any pre-existing 'log1p' layer is overwritten.
    """
    if counts_layer not in adata.layers:
        raise KeyError(f"layers['{counts_layer}'] not found; available: {list(adata.layers)}")
    tmp = sc.AnnData(X=adata.layers[counts_layer].copy())
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    sc.pp.log1p(tmp)
    adata.layers['log1p'] = tmp.X
    print(f"  layers['log1p'] <- log1p(normalize_total(layers['{counts_layer}'], target_sum={target_sum}))")
    return adata


def add_pseudotime_bins(adata, pt_key, n_bins, pool_timepoints=False):
    """Add obs['pt_bin'] = equal-count quantile bins of a rank-normalised pseudotime column.

    Added 2026-09-11 for the fib/dec pseudotime lane (rule liana_pt, 17_liana_pseudotime.smk).
    The bins become the --cell-type-key, so LIANA scores every (source bin -> target bin)
    context: the diagonal is within-stage signalling and the off-diagonals are
    earlier->later (upper) vs later->earlier (lower) along the trajectory.

    pool_timepoints collapses obs['timepoint'] to a single 'all' level so the pseudotime axis
    is the only grouping; without it, bins are still global but the script loops per timepoint.
    """
    import pandas as pd
    if pt_key not in adata.obs:
        raise KeyError(f"--pt-key {pt_key!r} not in obs; available pseudotime-like columns: "
                       f"{[c for c in adata.obs.columns if 'pseudotime' in c or 'dpt' in c]}")
    pt = adata.obs[pt_key].astype(float)
    n_nan = int(pt.isna().sum())
    if n_nan:
        raise ValueError(f"{pt_key} has {n_nan} NaN of {len(pt)} cells; subset the object before binning")
    labels = [f"pt_q{i + 1}" for i in range(n_bins)]
    bins = pd.qcut(pt.rank(method="first"), n_bins, labels=labels)
    adata.obs["pt_bin"] = pd.Categorical(bins, categories=labels, ordered=True)
    print(f"  pt_bin <- {n_bins} equal-count quantile bins of obs[{pt_key!r}] "
          f"(range {pt.min():.3f}-{pt.max():.3f})")
    for lab in labels:
        m = adata.obs["pt_bin"] == lab
        rng = (pt[m].min(), pt[m].max())
        print(f"    {lab}: n={int(m.sum())}  {pt_key} {rng[0]:.3f}-{rng[1]:.3f}")
    if "cell_subtype" in adata.obs:
        print("  pt_bin x cell_subtype:")
        print(pd.crosstab(adata.obs["pt_bin"], adata.obs["cell_subtype"]).to_string())
    if "timepoint" in adata.obs:
        print("  pt_bin x timepoint:")
        print(pd.crosstab(adata.obs["pt_bin"], adata.obs["timepoint"]).to_string())
    if pool_timepoints:
        adata.obs["timepoint"] = pd.Categorical(["all"] * adata.n_obs)
        print("  timepoint collapsed to a single 'all' level (--pool-timepoints)")
    return adata


def apply_cell_type_map(adata, map_yaml):
    """Add obs['cell_type_merged'] and obs['cell_type_xspecies'] from a YAML mapping (2026-09-09 port of
    merge_immune_xspecies.ipynb). Every source label must be in base_map or passthrough."""
    import yaml
    m = yaml.safe_load(open(map_yaml))
    src_col = m['source_column']
    if src_col not in adata.obs:
        raise KeyError(f"cell-type map source column {src_col!r} not in obs; have {list(adata.obs.columns)}")
    src = adata.obs[src_col].astype(str)
    labels = set(src.unique())
    covered = set(m['base_map']) | set(m['passthrough'])
    unmapped = sorted(labels - covered)
    if unmapped:
        raise ValueError(f"cell-type map {map_yaml}: unmapped {src_col} labels: {unmapped}")
    merged = src.replace(m['base_map'])
    adata.obs['cell_type_merged'] = merged.astype('category')
    adata.obs['cell_type_xspecies'] = merged.replace(m.get('to_lcd', {})).astype('category')
    print(f"  cell_type_merged from {src_col} via {map_yaml}:")
    print(adata.obs['cell_type_merged'].value_counts().to_string())
    print(f"  cell_type_xspecies categories: {sorted(adata.obs['cell_type_xspecies'].unique())}")
    return adata


def setup_output_dirs(timepoint, cell_type_key, output_base=None):
    root = f'{output_base or OUTPUT_BASE}/{timepoint}/{cell_type_key}'
    dirs = {
        'root':        root,
        'dotplots':    f'{root}/01_dotplots',
        'lr_heatmaps': f'{root}/02_lr_heatmaps',
        'spatial':     f'{root}/03_spatial',
        'tables':      f'{root}/tables',
    }
    for d in dirs.values():
        Path(d).mkdir(parents=True, exist_ok=True)
    return dirs


def _filter_spatial_knn_top_pct_by_sample(
    adata,
    sample_key='sample',
    x_key='x_um_dbscan',
    y_key='y_um_dbscan',
    k=5,
    keep_pct=0.95,
):
    """Per-sample kNN outlier removal. Returns (inliers, all_with_coords)."""
    mask_xy = adata.obs[x_key].notna() & adata.obs[y_key].notna()
    adata_xy = adata[mask_xy].copy()
    adata_xy.obsm['X_spatial'] = adata_xy.obs[[x_key, y_key]].to_numpy(float)

    kept = []
    for s in adata_xy.obs[sample_key].unique():
        sub = adata_xy[adata_xy.obs[sample_key] == s].copy()
        X = sub.obsm['X_spatial']
        n = X.shape[0]
        if n <= k + 1:
            sub.obs['spatial_knn_keep'] = True
            kept.append(sub)
            continue
        nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
        dists, _ = nn.kneighbors(X)
        mean_dist = dists[:, 1:].mean(axis=1)
        q = np.quantile(mean_dist, keep_pct)
        sub.obs['spatial_knn_mean_dist'] = mean_dist
        sub.obs['spatial_knn_keep'] = mean_dist <= q
        sub.uns['spatial_knn_params'] = dict(k=k, keep_pct=float(keep_pct), threshold=float(q))
        kept.append(sub[mean_dist <= q].copy())

    adata_inliers = kept[0].concatenate(*kept[1:], join='outer', batch_key=None)
    adata_inliers.obsm['X_spatial_plot'] = adata_inliers.obsm['X_spatial']
    return adata_inliers, adata_xy


def filter_spatially(adata, sample_key='sample', x_key='x_um_dbscan', y_key='y_um_dbscan',
                     k=10, keep_pct=0.95):
    """
    Per-sample kNN outlier removal + dataset-specific coordinate exclusions.
    Matches the logic used in the bandwidth investigation notebooks.
    """
    adata_inliers, adata_xy = _filter_spatial_knn_top_pct_by_sample(
        adata, sample_key=sample_key, x_key=x_key, y_key=y_key,
        k=k, keep_pct=keep_pct,
    )
    n_total = len(adata)
    n_mapped = len(adata_xy)
    n_inliers = len(adata_inliers)

    adata_out = adata_inliers.copy()

    # Human data: no dataset-specific coordinate exclusions needed.
    # (Mouse version has exclusions for 2dpi and 1dpi_late outlier regions.)

    print(f"    filter_spatially: {n_total} total -> {n_mapped} mapped -> "
          f"{n_inliers} after outlier removal "
          f"(per-sample kNN, k={k}, keep_pct={keep_pct})")
    return adata_out


def safe_savefig(path):
    plt.savefig(path, bbox_inches='tight')
    plt.close()


def tic(tag, step):
    print(f"\n[{tag}] >>> {step}...")
    return time.time()


def toc(tag, step, t0):
    elapsed = time.time() - t0
    mins = int(elapsed // 60)
    secs = int(elapsed % 60)
    print(f"[{tag}] <<< {step} complete ({mins}m {secs}s)")
    return time.time()


# ============================================================
# MAIN PIPELINE
# ============================================================

def run_timepoint_analysis(adata_full, timepoint, cell_type_key='cell_type_annotations',
                           bandwidth=None, expr_prop=None, output_base=None, resource=None):
    """
    bandwidth:   if None, uses global BANDWIDTH config value.
    expr_prop:   if None, uses global EXPR_PROP config value.
    output_base: if None, uses global OUTPUT_BASE config value.
    resource:    if None, uses global RESOURCE config value.
    """

    bw = bandwidth if bandwidth is not None else BANDWIDTH
    ep = expr_prop if expr_prop is not None else EXPR_PROP
    res = resource if resource is not None else RESOURCE
    tag = f"{timepoint}/{cell_type_key}"
    pipeline_start = time.time()

    print(f"\n{'='*60}")
    print(f"  LIANA pipeline: {tag}")
    print(f"  Bandwidth: {bw}um")
    print(f"  expr_prop: {ep}")
    print(f"  resource:  {res}")
    print(f"  Started: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")

    dirs = setup_output_dirs(timepoint, cell_type_key, output_base=output_base)
    print(f"[{tag}] Output root: {dirs['root']}")

    # ----------------------------------------------------------
    # Subset to timepoint
    # ----------------------------------------------------------
    adata_tp = adata_full[adata_full.obs['timepoint'] == timepoint].copy()
    adata_tp.X = adata_tp.layers['log1p'].copy()
    cell_types = adata_tp.obs[cell_type_key].dropna().unique().tolist()

    n_mapped = adata_tp.obs['x_um_dbscan'].notna().sum()
    print(f"\n[{tag}] Dataset summary:")
    print(f"    Cells      : {len(adata_tp)}")
    print(f"    Genes      : {adata_tp.n_vars}")
    print(f"    Cell types : {len(cell_types)} — {sorted(cell_types)}")
    print(f"    Mapped     : {n_mapped}/{len(adata_tp)} ({100*n_mapped/len(adata_tp):.1f}%)")

    coverage = adata_tp.obs.groupby(cell_type_key).apply(
        lambda x: x['x_um_dbscan'].notna().mean()
    ).sort_values()
    print(f"\n[{tag}] Spatial coverage per cell type:")
    for ct_name, cov in coverage.items():
        n_ct = (adata_tp.obs[cell_type_key] == ct_name).sum()
        flag = ' <-- LOW' if cov < 0.5 else ''
        print(f"    {ct_name:<35} {cov*100:5.1f}%  (n={n_ct}){flag}")
    coverage.to_csv(f"{dirs['tables']}/spatial_coverage.csv")


    # ----------------------------------------------------------
    # STEP 1: LR DISCOVERY ON FULL (unfiltered) DATA
    # ----------------------------------------------------------
    t0 = tic(tag, "Step 1: rank_aggregate on full data")

    li.mt.rank_aggregate(
        adata_tp,
        groupby=cell_type_key,
        resource_name=res,
        expr_prop=ep,
        use_raw=False,
        verbose=True,
    )

    liana_res = adata_tp.uns['liana_res']
    liana_res.to_csv(f"{dirs['tables']}/liana_res_full.csv")
    print(f"[{tag}]   liana_res shape: {liana_res.shape}")
    print(f"[{tag}]   sources: {sorted(liana_res['source'].unique().tolist())}")
    print(f"[{tag}]   targets: {sorted(liana_res['target'].unique().tolist())}")

    # L-R discovery criterion (2026-09-09): significance = RANK_COLUMN < RANK_CUTOFF (default
    # magnitude_rank < 0.05, the rank_aggregate threshold used in Dimitrov et al. 2024); ordered by magnitude_rank;
    # optional TOP_N_LR cap. Previously: top 50 by magnitude_rank with no significance filter.
    sig_res = liana_res[liana_res[RANK_COLUMN] < RANK_CUTOFF].sort_values('magnitude_rank')
    sig_res.to_csv(f"{dirs['tables']}/liana_res_significant.csv")
    top_lr = sig_res[['ligand_complex', 'receptor_complex']].drop_duplicates()
    if TOP_N_LR is not None:
        top_lr = top_lr.head(TOP_N_LR)
    top_lr.to_csv(f"{dirs['tables']}/top_lr_pairs.csv")
    print(f"[{tag}]   LR pairs with {RANK_COLUMN} < {RANK_CUTOFF}: {len(top_lr)} unique pairs "
          f"({len(sig_res)} source->target rows of {len(liana_res)}; cap={TOP_N_LR}):")
    print(top_lr.to_string(index=False))

    t0 = toc(tag, "Step 1", t0)

    # Dotplot: each cell type as source
    t0 = tic(tag, "Step 1b: dotplots per source cell type")
    valid_sources = set(liana_res['source'].unique())
    valid_targets = set(liana_res['target'].unique())
    for source in cell_types:
        if source not in valid_sources:
            print(f"[{tag}]   Skipping dotplot for {source}: not in liana_res sources")
            continue
        target_labels = [c for c in cell_types if c != source and c in valid_targets]
        if not target_labels:
            print(f"[{tag}]   Skipping dotplot for {source}: no valid targets")
            continue
        source_clean = source.replace(' ', '_')
        print(f"[{tag}]   Dotplot: {source} -> {len(target_labels)} targets")
        try:
            p = li.pl.dotplot(
                adata=adata_tp,
                colour='lr_means',
                size='magnitude_rank',
                inverse_size=True,
                source_labels=[source],
                target_labels=target_labels,
                top_n=TOP_N_DOTPLOT,
                orderby='magnitude_rank',
                orderby_ascending=True,
                figure_size=(12, 12)
            )
            p = (p
                 + p9.ggtitle(f'{source} -> targets ({timepoint})')
                 + p9.theme(
                     axis_text_y=p9.element_text(size=8),
                     axis_text_x=p9.element_text(size=8, angle=45, ha='right'),
                     plot_title=p9.element_text(size=10),
                     plot_margin_left=0.3,
                 ))
            fig = p.draw()
            display(fig)
            fig.savefig(f"{dirs['dotplots']}/{source_clean}_dotplot.pdf",
                        bbox_inches='tight', dpi=150)
            plt.close(fig)
        except Exception as e:
            print(f"[{tag}]   Skipping dotplot for {source}: {e}")

    t0 = toc(tag, "Step 1b", t0)

    # Heatmaps: top LR pairs across cell type pairs
    t0 = tic(tag, "Step 1c: LR heatmaps")
    for _, row in top_lr.head(25).iterrows():
        lig, rec = row['ligand_complex'], row['receptor_complex']
        interaction_df = liana_res[
            (liana_res['ligand_complex'] == lig) &
            (liana_res['receptor_complex'] == rec)
        ].copy()
        if interaction_df.empty:
            print(f"[{tag}]   Skipping {lig}-{rec}: empty")
            continue
        print(f"[{tag}]   Heatmap: {lig} -> {rec} "
              f"({len(interaction_df)} source/target pairs)")
        try:
            fig, ax = plt.subplots(figsize=(8, 6))
            pivot_score = interaction_df.pivot_table(
                index='source', columns='target', values='lrscore'
            )
            pivot_annot = interaction_df.pivot_table(
                index='source', columns='target', values='specificity_rank'
            ).applymap(lambda x: '★' if x < 0.05 else '')
            sns.heatmap(pivot_score, annot=pivot_annot, fmt='',
                        square=True, cmap='viridis',
                        cbar_kws={'label': 'lrscore'}, ax=ax)
            ax.set_title(f'{lig} -> {rec} ({timepoint})')
            plt.tight_layout()
            safe_savefig(f"{dirs['lr_heatmaps']}/{lig}_{rec}_heatmap.pdf")
        except Exception as e:
            print(f"[{tag}]   Skipping heatmap {lig}-{rec}: {e}")
            plt.close()

    t0 = toc(tag, "Step 1c", t0)


    if SKIP_SPATIAL:
        print(f"[{tag}] --skip-spatial set: non-spatial deliverable complete after rank_aggregate; Steps 2-4 (spatial subset, Moran's I, bivariate) not run.")
        _summarize(tag, dirs, pipeline_start)
        return

    # ----------------------------------------------------------
    # STEP 2: SPATIAL SUBSET
    # ----------------------------------------------------------
    t0 = tic(tag, "Step 2: spatial subset + neighbors")

    adata_spatial = filter_spatially(adata_tp)
    adata_spatial.X = adata_spatial.layers['log1p'].copy()
    adata_spatial.obsm[SPATIAL_KEY] = np.column_stack([
        adata_spatial.obs['x_um_dbscan'].values,
        adata_spatial.obs['y_um_dbscan'].values
    ])
    adata_spatial.obsm['spatial'] = adata_spatial.obsm[SPATIAL_KEY].copy()

    spatial_cell_types = adata_spatial.obs[cell_type_key].dropna().unique().tolist()
    print(f"[{tag}]   Spatial cells: {len(adata_spatial)}")
    print(f"[{tag}]   Cell types in spatial subset ({len(spatial_cell_types)}): "
          f"{sorted(spatial_cell_types)}")
    dropped_ct = set(cell_types) - set(spatial_cell_types)
    if dropped_ct:
        print(f"[{tag}]   WARNING: cell types lost after spatial filter: {dropped_ct}")

    li.ut.spatial_neighbors(
        adata=adata_spatial,
        bandwidth=bw,
        spatial_key=SPATIAL_KEY
    )
    print(f"[{tag}]   Spatial neighbors computed (bandwidth={bw}um)")

    # Cell type spatial layout
    try:
        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.embedding(adata_spatial, basis=SPATIAL_KEY, color=cell_type_key,
                        s=10, ax=ax, show=False,
                        title=f'Cell types ({timepoint})')
        safe_savefig(f"{dirs['spatial']}/cell_types_spatial.pdf")
    except Exception as e:
        print(f"[{tag}]   Skipping cell type spatial plot: {e}")
        plt.close()

    # Cell type proximity heatmap
    print(f"[{tag}]   Computing cell type spatial proximity...")
    try:
        pair_proximity = li.ut.spatial_pair_proximity(
            adata_spatial,
            groupby=cell_type_key,
            kernel='gaussian',
            bandwidth=bw,
            spatial_key=SPATIAL_KEY,
            verbose=True
        )
        pair_proximity.to_csv(f"{dirs['tables']}/cell_type_proximity.csv")
        n_interacting = (pair_proximity['interacting'] == 1).sum()
        print(f"[{tag}]   Proximity computed: {n_interacting} interacting cell type pairs")
        fig, ax = plt.subplots(figsize=(7, 6))
        sns.heatmap(
            data=pair_proximity.pivot(index='source', columns='target', values='proximity'),
            annot=pair_proximity.pivot(index='source', columns='target', values='interacting')
                               .fillna(0).astype(int).replace({1: '.', 0: ''}),
            fmt='', cmap='cividis', square=True, ax=ax
        )
        ax.set_title(f'Cell type proximity ({timepoint})')
        plt.tight_layout()
        safe_savefig(f"{dirs['spatial']}/cell_type_proximity.pdf")
    except Exception as e:
        print(f"[{tag}]   Skipping proximity: {e}")
        plt.close()

    t0 = toc(tag, "Step 2", t0)


    # ----------------------------------------------------------
    # STEP 3: INFLOW + MORAN'S I
    # ----------------------------------------------------------
    t0 = tic(tag, "Step 3: inflow + Moran's I")

    resource = li.rs.select_resource(res)
    print(f"[{tag}]   Running inflow on {len(adata_spatial)} spatial cells...")
    lrdata = li.mt.inflow(
        adata_spatial,
        groupby=cell_type_key,
        resource=resource,
        use_raw=False
    )
    print(f"[{tag}]   lrdata shape: {lrdata.shape}")

    print(f"[{tag}]   Running Moran's I on {lrdata.n_vars} LR pairs...")
    sq.gr.spatial_autocorr(lrdata, mode='moran', use_raw=False, show_progress_bar=True)

    moran_df = lrdata.uns['moranI'].dropna(subset=['I'])
    print(f"[{tag}]   Moran's I computed: {len(moran_df)} pairs with valid I scores")
    print(f"[{tag}]   I score distribution: "
          f"min={moran_df['I'].min():.3f}, "
          f"median={moran_df['I'].median():.3f}, "
          f"max={moran_df['I'].max():.3f}")

    _, pvals_corr, _, _ = multipletests(moran_df['pval_norm'], method='fdr_bh')
    moran_df['pval_norm_fdr_bh'] = pvals_corr
    lrdata.uns['moranI'] = moran_df
    moran_df.to_csv(f"{dirs['tables']}/moranI_all.csv")

    svis = moran_df.index[
        (moran_df['pval_norm_fdr_bh'] <= 0.05) &
        (moran_df['I'] > MORAN_I_MIN)
    ]
    print(f"[{tag}]   Spatially variable interactions: {len(svis)} "
          f"(FDR<=0.05, I>{MORAN_I_MIN})")
    lrdata = lrdata[:, svis]

    # Intersection: rank_aggregate top hits that are also spatially variable
    top_lr_names = set(top_lr['ligand_complex'] + '^' + top_lr['receptor_complex'])
    spatial_hits = lrdata.uns['moranI'].sort_values('I', ascending=False).copy()
    spatial_hits['lr_pair'] = spatial_hits.index.map(
        lambda x: x.split('^')[1] + '^' + x.split('^')[2]
    )
    priority_hits = spatial_hits[spatial_hits['lr_pair'].isin(top_lr_names)]
    priority_hits.to_csv(f"{dirs['tables']}/priority_hits_intersection.csv")

    print(f"[{tag}]   Priority hits (rank_aggregate AND spatially variable): "
          f"{len(priority_hits)}")
    if len(priority_hits) > 0:
        print(priority_hits[['lr_pair', 'I', 'pval_norm_fdr_bh']].head(20).to_string())
    else:
        print(f"[{tag}]   WARNING: no intersection found — check rank_aggregate "
              f"and Moran's I thresholds")

    # Global specificity
    print(f"[{tag}]   Computing global specificity...")
    li.mt.compute_global_specificity(lrdata, groupby=cell_type_key, use_raw=False, verbose=True)
    gi = lrdata.uns['global_interactions']
    gi.to_csv(f"{dirs['tables']}/global_interactions.csv")
    print(f"[{tag}]   global_interactions shape: {gi.shape}")

    # Spatial embeddings: top inflow interactions
    print(f"[{tag}]   Plotting spatial embeddings for top 10 inflow interactions...")
    top_inflow_df = lrdata.uns['moranI'].sort_values('I', ascending=False).head(25)
    top_inflow = top_inflow_df.index.tolist()
    print(f"[{tag}]   Top inflow interactions (rank | Moran's I | interaction):")
    for rank, (idx, row) in enumerate(top_inflow_df.iterrows(), 1):
        print(f"    #{rank:02d}  I={row['I']:.3f}  {idx}")

    for rank, interaction in enumerate(top_inflow, 1):
        comp = interaction.split('^')
        cell_type_label, lig, rec = comp[0], comp[1], comp[2]
        moran_i = top_inflow_df.loc[interaction, 'I']
        receptor_genes = rec.split('_')
        genes_to_plot = [g for g in [lig] + receptor_genes if g in adata_spatial.var_names]
        safe_name = interaction.replace(' ', '_').replace('^', '_')
        missing = [g for g in [lig] + receptor_genes if g not in adata_spatial.var_names]
        if missing:
            print(f"[{tag}]     {lig}-{rec}: genes not in var_names: {missing}")

        try:
            inflow_vals = lrdata[:, interaction].X
            coords = lrdata.obsm[SPATIAL_KEY]
            fig, ax = plt.subplots(figsize=(7, 6))
            fig.patch.set_facecolor('black')
            plot_spatial_dark(
                ax, coords[:, 0], coords[:, 1], inflow_vals,
                title=f'#{rank} {lig}->{rec} | {cell_type_label}\n'
                      f"Moran's I={moran_i:.3f} | {timepoint}",
            )
            save_dark_figure(fig, f"{dirs['spatial']}/rank{rank:02d}_{safe_name}_inflow.pdf")
        except Exception as e:
            print(f"[{tag}]     Skipping inflow plot {interaction}: {e}")
            plt.close()

        if genes_to_plot:
            try:
                coords = adata_spatial.obsm[SPATIAL_KEY]
                fig, axes = plt.subplots(1, len(genes_to_plot),
                                         figsize=(7 * len(genes_to_plot), 6))
                fig.patch.set_facecolor('black')
                if len(genes_to_plot) == 1:
                    axes = [axes]
                for ax, gene in zip(axes, genes_to_plot):
                    plot_spatial_dark(
                        ax, coords[:, 0], coords[:, 1],
                        adata_spatial[:, gene].X,
                        title=f'#{rank} {gene} | {lig}->{rec} ({timepoint})',
                    )
                save_dark_figure(fig, f"{dirs['spatial']}/rank{rank:02d}_{safe_name}_expression.pdf")
            except Exception as e:
                print(f"[{tag}]     Skipping expression plot {interaction}: {e}")
                plt.close()

    t0 = toc(tag, "Step 3", t0)


    # ----------------------------------------------------------
    # STEP 4: COMMOT (disabled — uncomment to enable)
    # ----------------------------------------------------------
    # To re-enable: uncomment `import commot as ct` at top,
    # and restore this block from liana_pipeline-5.py

    _summarize(tag, dirs, pipeline_start)


def _summarize(tag, dirs, pipeline_start):
    total = time.time() - pipeline_start
    hours = int(total // 3600)
    mins  = int((total % 3600) // 60)
    secs  = int(total % 60)
    print(f"\n{'='*60}")
    print(f"  [{tag}] PIPELINE COMPLETE")
    print(f"  Total time : {hours}h {mins}m {secs}s")
    print(f"  Output     : {dirs['root']}")
    print(f"{'='*60}\n")


# ============================================================
# MAIN
# ============================================================

# Per-timepoint bandwidths from bandwidth investigation notebook.
# Human: ~50nn target for both timepoints.
HUMAN_BANDWIDTHS = {
    'menstrual':  112,
    'secretory':  101,
}

# HUMAN_H5AD legacy default removed (dead /n/eddy_lab path); --h5ad is required.


def main():
    import argparse
    global RANK_COLUMN, RANK_CUTOFF, TOP_N_LR, MORAN_I_MIN, SKIP_SPATIAL

    parser = argparse.ArgumentParser(description='LIANA spatial CCC pipeline (human, no COMMOT)')
    parser.add_argument('--timepoints', nargs='*', default=None,
                        help='Timepoints to run (default: all)')
    parser.add_argument('--h5ad', required=True,
                        help='Path to h5ad file')
    parser.add_argument('--cell-type-key', default='cell_type_annotations',
                        help='obs column for cell type labels')
    parser.add_argument('--expr-prop', type=float, default=EXPR_PROP,
                        help=f'Min fraction of cells expressing a gene (default: {EXPR_PROP})')
    parser.add_argument('--rank-column', default=RANK_COLUMN, choices=['magnitude_rank', 'specificity_rank'],
                        help=f'rank_aggregate column defining significance (default: {RANK_COLUMN})')
    parser.add_argument('--rank-cutoff', type=float, default=RANK_CUTOFF,
                        help=f'LR discovery: keep source->target rows with rank-column < cutoff (default: {RANK_CUTOFF})')
    parser.add_argument('--skip-spatial', action='store_true',
                        help='Stop after rank_aggregate + non-spatial plots (no spatial subset, Moran\'s I, or bivariate)')
    parser.add_argument('--top-n-lr', type=int, default=None,
                        help='Optional cap on discovered LR pairs after the specificity filter (default: no cap)')
    parser.add_argument('--moran-i-min', type=float, default=MORAN_I_MIN,
                        help=f"Moran's I floor for spatially variable interactions (default: {MORAN_I_MIN})")
    parser.add_argument('--cell-type-map', default=None,
                        help='YAML mapping (source_column, base_map, passthrough, to_lcd) that adds obs cell_type_merged / cell_type_xspecies before analysis')
    parser.add_argument('--output-base', required=True,
                        help='Output root; results land in <output-base>/<timepoint>/<cell-type-key>/')
    parser.add_argument('--counts-layer', default='counts',
                        help="layers key holding raw counts; log1p is recomputed from it (default: counts)")
    parser.add_argument('--target-sum', type=float, default=1e4,
                        help='normalize_total target_sum for the recomputed log1p layer (default: 1e4)')
    parser.add_argument('--pt-key', default=None,
                        help='obs column with a rank-normalised pseudotime; enables pt_bin as a cell-type key')
    parser.add_argument('--pt-bins', type=int, default=5,
                        help='number of equal-count pseudotime bins (default: 5)')
    parser.add_argument('--pool-timepoints', action='store_true',
                        help="collapse obs['timepoint'] to one 'all' level so pseudotime is the only axis")
    args = parser.parse_args()
    RANK_COLUMN, RANK_CUTOFF, TOP_N_LR, MORAN_I_MIN, SKIP_SPATIAL = (
        args.rank_column, args.rank_cutoff, args.top_n_lr, args.moran_i_min, args.skip_spatial)
    print(f"rank_column: {RANK_COLUMN}  rank_cutoff: {RANK_CUTOFF}  top_n_lr: {TOP_N_LR}  moran_i_min: {MORAN_I_MIN}  skip_spatial: {SKIP_SPATIAL}")

    timepoints = args.timepoints or list(HUMAN_BANDWIDTHS.keys())

    # Output root is owned by the caller (Snakemake encodes expr_prop in its own path/params).
    output_base = args.output_base

    print(f"Loading {args.h5ad} ...")
    print(f"expr_prop: {args.expr_prop}")
    print(f"resource:  {RESOURCE}")
    print(f"Output:    {output_base}")
    adata_full = sc.read_h5ad(args.h5ad)
    print(f"Loaded: {adata_full.n_obs} cells, {adata_full.n_vars} genes")
    ensure_log1p_layer(adata_full, counts_layer=args.counts_layer, target_sum=args.target_sum)
    if args.cell_type_map:
        apply_cell_type_map(adata_full, args.cell_type_map)
    if args.pt_key:
        add_pseudotime_bins(adata_full, args.pt_key, args.pt_bins, pool_timepoints=args.pool_timepoints)
    if args.cell_type_key not in adata_full.obs:
        raise KeyError(f"--cell-type-key {args.cell_type_key!r} not in obs (pass --cell-type-map to derive cell_type_merged/cell_type_xspecies)")

    for tp in timepoints:
        bw = HUMAN_BANDWIDTHS.get(tp, BANDWIDTH)
        print(f"\n>>> Starting {tp} (bandwidth={bw}µm, expr_prop={args.expr_prop})")
        run_timepoint_analysis(
            adata_full,
            timepoint=tp,
            cell_type_key=args.cell_type_key,
            bandwidth=bw,
            expr_prop=args.expr_prop,
            output_base=output_base,
        )

    print("\n=== ALL TIMEPOINTS DONE ===")


if __name__ == '__main__':
    main()
