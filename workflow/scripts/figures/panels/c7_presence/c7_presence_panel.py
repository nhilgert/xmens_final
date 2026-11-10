"""The C7-fibroblast figure (reviewer comment 1), panels only and annotated: Ulrich et al. 2024's Stro-7 markers and
Stro-7 similarity single out human Leiden cluster 3 (a-c); through the SAMap graph, mouse outer fibroblasts (eFib1)
draw their human neighbours from cluster 3 (d-f).

Called by rule panels_c7_presence (workflow/rules/41_fig_c7_presence.smk), on the tables of rule c7_tables.
A port of fig_c7_presence_v8.py (FIG_c7_presence v8, 2 Nov 2026), whose drawing code is kept line for line; changed
here: inputs and outputs from the command line, fonts from c7_common (Helvetica Neue, else Nimbus Sans), the mouse
population colors from --mouse-colors, the mouse sections turned by resources/figures/mouse_section_rotation.yaml
through workflow/scripts/figures/section_rotation.py with the windows of --mouse-crops, PDF + SVG + 300-dpi PNG per
variant, and one provenance record.
  a  dot plot of the eight Stro-7 markers over the six Leiden clusters (scanpy dotplot semantics, standard_scale='var')
  b  human UMAP: clusters | C7 | contrastive Stro-7 cosine
  c  human tissue (one section per donor): cluster 3 | contrastive Stro-7 cosine
  d  mean share of human SAMap neighbours in cluster 3 per mouse population
  e  mouse UMAP: populations | share of human SAMap neighbours in cluster 3 | the neighbours' Stro-7 cosine
  f  mouse tissue, one section per timepoint: the same two quantities
Marker genes: the 8 Stro-7 genes with the largest margin against all 39 subtypes in Ulrich Dataset S5 sheet 5B; all 8 are
in their published Stro-7 DE list (Dataset S6 sheet 6C). Contrastive cosine: rule c7_stro7_cosine.

usage: c7_presence_panel.py --human-table F --dotplot-table F --mouse-table F --background-table F --mouse-colors YAML
                            --mouse-crops YAML --section-rotation YAML --joint-name NAME
                            --panels-pdf F --panels-svg F --panels-png F --annotated-pdf F --annotated-svg F
                            --annotated-png F --provenance JSON
"""
import argparse, logging, os, re, sys, textwrap
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt, yaml
from matplotlib.lines import Line2D
import xml.etree.ElementTree as ET
import c7_common as common
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # workflow/scripts/figures: the mouse section rotation reader
import section_rotation  # noqa: E402

_P = argparse.ArgumentParser(description='The C7-fibroblast figure')
for _k in ('human-table', 'dotplot-table', 'mouse-table', 'background-table', 'mouse-colors', 'mouse-crops', 'section-rotation',
           'panels-pdf', 'panels-svg', 'panels-png', 'annotated-pdf', 'annotated-svg', 'annotated-png', 'provenance'):
    _P.add_argument(f'--{_k}', required=True, type=Path)
_P.add_argument('--joint-name', required=True, help='file name of the SAMap object, for the annotated footnote')
ARGS = _P.parse_args()
logging.basicConfig(level=logging.INFO, format='%(message)s')
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR); logging.getLogger('fontTools').setLevel(logging.WARNING)
LOG = logging.getLogger('panels_c7_presence')
OUTPATHS = {'panels': {'pdf': ARGS.panels_pdf, 'svg': ARGS.panels_svg, 'png': ARGS.panels_png},
            'annotated': {'pdf': ARGS.annotated_pdf, 'svg': ARGS.annotated_svg, 'png': ARGS.annotated_png}}
for _v in OUTPATHS.values():
    for _p in _v.values(): _p.parent.mkdir(parents=True, exist_ok=True)

FIG_STYLE = {
    'axes.labelsize': 10.5,
    'axes.linewidth': 0.6,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.titlesize': 10.5,
    'figure.dpi': 200.0,
    'font.size': 10.5,
    'legend.fontsize': 9.5,
    'legend.frameon': False,
    'lines.linewidth': 1.2,
    'patch.linewidth': 0.6,
    'ps.fonttype': 42,
    'savefig.bbox': 'tight',
    'savefig.dpi': 300.0,
    'xtick.labelsize': 9.0,
    'xtick.major.size': 3.0,
    'xtick.major.width': 0.6,
    'ytick.labelsize': 9.0,
    'ytick.major.size': 3.0,
    'ytick.major.width': 0.6,
}
mpl.rcParams.update(FIG_STYLE)
mpl.rcParams['axes.titlelocation'] = 'center'
FSC = 1.25                                               # v7.1: panel text (ticks, labels, legends, colorbars) scaled up from v7
FS_K, FS_TT, FS_L = 9 * FSC, 12.5, 17                    # labels | panel titles (was 9.5) | panel letters (was 15)
FS_H, FS_P, FS_F = 14, 9.5, 8.5                          # annotated version only: header titles (was 13), explanatory paragraphs (was 9), footnote (was 8)
RASTER_DPI = 600                                                                # rasterized point layers in the SVG and the PNG
PAL = ['#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#56B4E9']   # v7.4: Okabe-Ito (was scanpy default_20)        # scanpy default_20, first six
GREY = '#d3d3d3'
CM_C7 = mpl.colors.LinearSegmentedColormap.from_list('c7_grayfloor', ['#d3d3d3', '#21918c', '#440154'])   # v7.5: C7 on white, undetected nuclei light gray (as v7.3)
CM_FRAC = mpl.colors.LinearSegmentedColormap.from_list('plasma_from35', mpl.colormaps['plasma'](np.linspace(0.35, 1, 256)))  # v7.7: on black
ICE_GRAY = ['#bde7db', '#aed5d3', '#9fc4cc', '#90b4c6', '#7fa5c2', '#6f96bf', '#5f87bf', '#4f77c1', '#4167c6', '#4564ad', '#486098', '#495c85', '#4a5774', '#495266', '#484d59', '#46494d', '#454444', '#514746', '#5e4847', '#6d4848', '#7b4747', '#8a4546', '#9a4245', '#aa3c43', '#ba3741', '#c34d4d', '#cc615a', '#d57466', '#de8774', '#e69b82', '#efae90', '#f7c29f', '#ffd4ac']
CM_COS = mpl.colors.LinearSegmentedColormap.from_list('icefire_graycenter', ICE_GRAY)   # v7.4: seaborn 0.13.2 icefire at 0, .25, .75, 1 with a #444444 center, interpolated in CAM02-UCS
BLK, BG_DARK = 'black', '#262626'                                                         # v7.4: background of continuous panels; other mouse cells in e
C_OTHER = '#4a4a4a'                                                                       # v7.7: other clusters in the cluster-3 map (c), on black
def scalebar(ax, x0, y0, span):                                                           # v7.7: 1 mm bar just under a panel, in that panel's data units
    yb = y0 - 0.04 * span; ax.plot([x0, x0 + 1000], [yb, yb], color='0.2', lw=1.5, clip_on=False, scalex=False, scaley=False, solid_capstyle='butt')
    ax.text(x0 + 500, yb - 0.03 * span, '1 mm', ha='center', va='top', fontsize=8 * FSC, color='0.2')
REF_S, REF_IN, REF_DPI = 25.0, 6.0, 300                                                    # v7.8: dot size, panel size and dpi of cagri_mouse_spatial.py
def ref_um_per_in(x, y):                                                                  # that script's scale for these points (its padding rule)
    xr, yr = np.ptp(x), np.ptp(y); r_px = np.sqrt(REF_S / np.pi) / 72 * REF_DPI; w_px = REF_IN * REF_DPI
    px, py = max(.005 * xr, (r_px + 2) / w_px * xr), max(.005 * yr, (r_px + 2) / w_px * yr)
    return max(xr + 2 * px, yr + 2 * py) / REF_IN
matched_s = lambda x, y, ours_um_per_in: REF_S * (ref_um_per_in(x, y) / ours_um_per_in) ** 2   # same dot diameter in micrometers
DOT_S = {}
def blackbox(ax):                                                                         # v7.4: black panel without ticks or frame
    ax.set_facecolor(BLK); ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(False) for sp in ax.spines.values()]
CLN = lambda g: f'cluster {int(str(g).lstrip("c"))}'                                   # v7.2: display name of a cluster key (c3 or 3 -> 'cluster 3')

# ---------------- data: the 2026-11-02 final objects, all nuclei, no exclusions, as tables of rule c7_tables: human (47,533 nuclei:
#   leiden_0.5 of rule c7_leiden_labels, stro7_cosine of rule c7_stro7_cosine, C7, X_umap, QC'd positions), dot plot, mouse (43,197
#   cells: cell_subtype, X_umap, positions, and from all 47,533 human nuclei through the SAMap graph the share of neighbours per cluster
#   and the transferred Stro-7 cosine) and background (other positioned mouse cells).
HF = pd.read_csv(ARGS.human_table, index_col=0); HF.index = HF.index.astype(str)
assert len(HF) == 47533 and HF['leiden_0.5'].notna().all() and HF['stro7_cosine'].notna().all() and HF['C7_log'].notna().all()
GENES = ['C7', 'ADH1B', 'DPT', 'CCDC80', 'OGN', 'SERPINF1', 'DCN', 'MGP']
SPECIFIC4 = ['C7', 'ADH1B', 'DPT', 'CCDC80']
OURS = pd.read_csv(ARGS.dotplot_table)
OURS['group'] = 'c' + OURS.group.astype(int).astype(str)
GRP = sorted(OURS.group.unique(), key=lambda g: int(g[1:]))
FR = OURS.pivot_table(index='group', columns='human_gene', values='pct_expressing').reindex(index=GRP, columns=GENES) / 100
MN = OURS.pivot_table(index='group', columns='human_gene', values='mean_expr').reindex(index=GRP, columns=GENES)
MN = MN - MN.min(0); MN = (MN / MN.max(0)).fillna(0)                                   # standard_scale='var'
assert all(FR[g].idxmax() == 'c3' for g in SPECIFIC4), FR[SPECIFIC4].idxmax().to_dict()
DOT_MAX = np.ceil(FR.to_numpy().max() * 10) / 10; DOT_MIN = 0.0
size_of = lambda f: (np.clip((np.asarray(f) - DOT_MIN) / (DOT_MAX - DOT_MIN), 0, 1) ** 1.5) * 200.0
size_range = np.arange(DOT_MAX, DOT_MIN, -0.2)[::-1] if DOT_MAX - DOT_MIN > .6 else np.arange(DOT_MAX, DOT_MIN, -0.1)[::-1]

T = HF.reset_index().rename(columns={'umap_1': 'lin_umap_1', 'umap_2': 'lin_umap_2', 'stro7_cosine': 'cos_Stro-7'})
T['c3'] = T['leiden_0.5'] == 3; T['cos'] = T['cos_Stro-7']; T['timepoint'] = T['sample'].str.split('_').str[0]
U = T.dropna(subset=['lin_umap_1']).reset_index(drop=True); XY = T.dropna(subset=['x_um_qc'])
n_all, n_c3 = len(T), int(T.c3.sum()); n_c7 = int((U.C7_counts > 0).sum()); n_c7_c3 = int(((U.C7_counts > 0) & U.c3).sum())
pct = lambda g, k: 100 * FR.loc[k, g]
clim = float(np.quantile(np.abs(T.cos), .98)); cos_med = T.groupby('leiden_0.5').cos.median()
assert cos_med.idxmax() == 3
cos_oth = cos_med.drop(3)
pc3_tp = {tp: 100 * XY[XY.timepoint == tp].c3.mean() for tp in ('menstrual', 'secretory')}

# ---------------- mouse (v8: final object + SAMap transfer from all 47,533 human nuclei)
PALETTE = yaml.safe_load(ARGS.mouse_colors.read_text())
MX = pd.read_csv(ARGS.mouse_table, index_col=0); MX.index = MX.index.astype(str)
MX = MX.rename(columns={'umap_1': 'X_umap_1', 'umap_2': 'X_umap_2'})
MX['c3s'] = MX['share_c3']; MX['mcos'] = MX['stro7_cosine_transferred']; MX['pop'] = MX['cell_subtype'].astype(str)
MPOP = list(PALETTE['order']); assert set(MX['pop']) == set(MPOP), set(MX['pop'])
MCOL = dict(PALETTE['colors']); assert set(MCOL) == set(MPOP)
MTP = [('uninduced', 'uninduced'), ('1.5_dpi_early', '1.5 dpi early'), ('1.5_dpi_late', '1.5 dpi late'), ('2.5_dpi', '2.5 dpi'), ('3.5_dpi', '3.5 dpi'), ('4.5_dpi_bleeding', '4.5 dpi bleeding')]
MU = MX.dropna(subset=['X_umap_1', 'c3s']).copy()
# mouse sections: every placed nucleus of a timepoint turned as one body by the shared rotation table (section_rotation), then
#   the window of --mouse-crops applied on the raw QC coordinates (bounds inclusive), so it selects the same nuclei at any angle.
ANGLES = section_rotation.load_rotation_table(ARGS.section_rotation)
_crops = yaml.safe_load(ARGS.mouse_crops.read_text())['crop_um']
MSEC = {tp: _crops[tp] for tp, _ in MTP}
N_CROP = sum(v is not None for v in MSEC.values())
BGF = pd.read_csv(ARGS.background_table, index_col=0); BGF.index = BGF.index.astype(str)
SMM = pd.concat([MX[['sample', 'timepoint', 'x_um_qc', 'y_um_qc']], BGF[['sample', 'timepoint', 'x_um_qc', 'y_um_qc']]]).rename_axis('obs_name').reset_index()
assert SMM.obs_name.is_unique
_pos = SMM.dropna(subset=['x_um_qc', 'y_um_qc']).reset_index(drop=True)
assert set(_pos['timepoint']) <= set(MSEC), 'a positioned mouse timepoint has no window in --mouse-crops'
_pos['x_disp'], _pos['y_disp'] = section_rotation.rotate_by_timepoint(_pos.x_um_qc.to_numpy(float), _pos.y_um_qc.to_numpy(float), _pos.timepoint.to_numpy(), ANGLES)
_parts = []; MSEC_N = {}
for _tp, _lim in MSEC.items():
    _a = _pos[_pos.timepoint == _tp]; _n0 = len(_a)
    if _lim is not None:
        (_x0, _x1), (_y0, _y1) = _lim['x'], _lim['y']
        _a = _a[(_a.x_um_qc >= _x0) & (_a.x_um_qc <= _x1) & (_a.y_um_qc >= _y0) & (_a.y_um_qc <= _y1)]
    _parts.append(_a.set_index('obs_name')[['x_disp', 'y_disp']]); MSEC_N[_tp] = (_n0, len(_a))
DISP = pd.concat(_parts); assert DISP.index.is_unique
MXY = MX.join(DISP, how='inner').dropna(subset=['c3s'])
MBG = SMM[~SMM.obs_name.isin(MX.index)].join(DISP, on='obs_name', how='inner')                    # other mouse cells, tissue context
base_c3 = float(MX['source_frac_c3'].iloc[0])
mstat = MX.groupby('pop').agg(n=('c3s', 'size'), c3=('c3s', 'mean'), am=('argmax_cluster', lambda s: 100 * (s == 'c3').mean()), cosmed=('mcos', 'median')).reindex(MPOP)
msec = MX.groupby(['pop', 'sample']).c3s.agg(['mean', 'size']); msec = msec[msec['size'] >= 20].reset_index()
assert mstat.c3.idxmax() == 'eFib1' and mstat.cosmed.idxmax() == 'eFib1', (mstat.c3.idxmax(), mstat.cosmed.idxmax())
EBOX = {}                                                                                  # square limits per section (all shown cells, 4% margin)
for _tp in MXY.timepoint.unique():
    _xy = np.r_[MXY.loc[MXY.timepoint == _tp, ['x_disp', 'y_disp']].to_numpy(), MBG.loc[MBG.timepoint == _tp, ['x_disp', 'y_disp']].to_numpy()]
    _lo, _hi = _xy.min(0), _xy.max(0); _c, _h = (_lo + _hi) / 2, (_hi - _lo).max() / 2 * 1.04; EBOX[_tp] = (_c[0] - _h, _c[0] + _h, _c[1] - _h, _c[1] + _h)

# ---------------- commentary (annotated version only); every number from the data above
oth = lambda g: max(pct(g, k) for k in GRP if k != 'c3')
TXT = dict(
    a=('Ulrich et al.’s C7-fibroblast markers pick out one cluster, c3',
       'Genes: the eight with the largest Stro-7 specificity margin against all 39 Ulrich subtypes (Stro-7 centroid minus the highest other centroid; their Dataset S5, '
       'sheet 5B). All eight are also in their published Stro-7 marker list (Dataset S6, sheet 6C). '
       f'C7, ADH1B, DPT and CCDC80 are highest in c3 (C7 in {pct("C7", "c3"):.1f}% of c3 nuclei vs at most {oth("C7"):.1f}% in any other cluster); '
       'DCN, OGN, SERPINF1 and MGP are expressed across fibroblasts in both datasets. Dot size: % of nuclei with the gene detected; color: mean log-normalized expression, scaled 0–1 per gene.', 90),
    b=('c3 is a discrete region of the UMAP, and both C7 and the Stro-7 similarity peak there',
       f'{len(U):,} fibroblast/decidual nuclei from both donors on the UMAP of the final fibroblast/decidual object. Left: Leiden clusters (resolution 0.5). '
       f'Middle: C7, log-normalized; {n_c7:,} nuclei have at least one count, {n_c7_c3:,} of them in c3. '
       'Right: the contrastive Stro-7 cosine, an atlas-similarity score that uses none of the genes in a: the cosine between a nucleus’s expression (relative to the average of our nuclei) '
       'and Ulrich’s Stro-7 centroid (relative to the average of their 11 stromal centroids), over the 3,000 genes that vary most across those centroids after removing each subtype’s top 50 published markers. '
       f'Positive means more Stro-7-like than the average stromal subtype. Median {cos_med[3]:+.3f} in c3, the highest of the {len(GRP)} clusters ({cos_oth.min():+.3f} to {cos_oth.max():+.3f} in the others).', 150),
    c=('In tissue, c3 and high Stro-7 similarity sit toward the edge of each section',
       f'{len(XY):,} of {n_all:,} fibroblast/decidual nuclei with a QC’d position. Each timepoint is one section; its three capture lanes share one coordinate frame and are overlaid. '
       f'Left pair: c3 in red-orange over the other clusters in dark gray ({pc3_tp["menstrual"]:.0f}% of positioned nuclei in the menstrual section, {pc3_tp["secretory"]:.0f}% in the secretory). '
       'Right pair: the contrastive Stro-7 cosine on the same color scale as in b. The cosine is graded rather than all-or-none, highest toward the edge and lowest in the interior bands, '
       'so c3 marks the high end of a gradient. Each donor’s section fills its panels and has its own 1 mm scale bar. Which edge is the basalis side has not been checked against histology.', 190),
    foot=('Reading notes: detection rates come from single nuclei and are lower than in Ulrich et al.’s whole-cell data, so compare clusters rather than absolute values. '
          'The dot plot and UMAPs follow scanpy defaults (sc.pl.dotplot with standard_scale=\'var\'; sc.pl.umap point size 120,000/n, highest values drawn last; C7 from light gray (not detected) through teal to purple, the fraction of human neighbors in cluster 3 in plasma without its darkest 35%, on black). '
          f'Clusters: leiden_0.5 as used in the paper; every panel shows all {n_all:,} human fibroblast/decidual nuclei, with none excluded. ' + f'Cosine colors (icefire with a dark-gray center, on black) are centered at 0 and clipped at ±{clim:.3f} (98th percentile of |cosine|), with nuclei drawn in random order. Tissue dots (c, e) are sized so each covers the same tissue distance as in single-section plots at 25 pt² on a 6-inch square.\n'
          'Sources: the final objects Cevrim_Human_Only_Fib_Dec_Cells.h5ad, Cevrim_XMens_Only_Fib_Dec_Cells.h5ad and Cevrim_XMens_All_Cells.h5ad (2026-11-02), '
          'through the rules of workflow/rules/41_fig_c7_presence.smk (xmens_final): c7_leiden_labels, c7_stro7_cosine, c7_tables, panels_c7_presence; '
          'Ulrich et al. 2024 Datasets S5 (5B) and S6 (6C).'))

TXT['d'] = ('Through the SAMap graph, mouse outer fibroblasts (eFib1) draw their human neighbors from c3',
    f'Each of the {len(MX):,} mouse fibroblast/decidual cells receives the connectivity-weighted average of its human neighbors in the joint SAMap graph (the method used for the pseudotime transfer): '
    'the fraction of neighbor weight on c3 nuclei, and the neighbors’ contrastive Stro-7 cosine. '
    f'Left: mean fraction of human neighbors in c3 per mouse population, all cells pooled; for comparison, c3 is {100 * base_c3:.1f}% of the human nuclei. '
    f'eFib1 averages {100 * mstat.c3["eFib1"]:.0f}% ({mstat.c3["eFib1"] / base_c3:.1f}×), and c3 is the top cluster for {mstat.am["eFib1"]:.0f}% of eFib1 cells; no decidual population exceeds {100 * mstat.c3[["dec1", "dec2", "dec3", "decAp"]].max():.0f}%. '
    f'Their median transferred Stro-7 cosine, {mstat.cosmed["eFib1"]:+.3f}, is also the highest of the {len(MPOP)} mouse populations. '
    'Right: the mouse fibroblast/decidual UMAP colored by population, by the fraction of human neighbors in c3, and by the neighbors’ Stro-7 cosine.', 150)
TXT['e'] = ('In mouse tissue, the c3-like cells form the outer fibroblast layer',
    f'{len(MXY):,} of {len(MX):,} mouse fibroblast/decidual cells with a QC’d position in the sections shown: one section per timepoint (lanes overlaid; {N_CROP} of {len(MSEC)} cropped to a window), rotated for display and scaled to fill its square, with its own 1 mm scale bar; other mouse cells (epithelium, immune, vascular, smooth muscle) in dark gray. '
    'Top: fraction of human neighbors in c3. Bottom: the neighbors’ Stro-7 cosine, on the color scale of b and c. Transferred values are neighbor averages, so they vary more smoothly than per-nucleus values.', 190)
TXT['foot'] = TXT['foot'].replace('\nSources:', f' Mouse values (d, e) are transferred from all {n_all:,} human nuclei through the SAMap graph of the joint fit ({ARGS.joint_name}).\nSources:')
# v7.2: cluster names in all commentary
_cl = lambda t: re.sub(r'\bc(\d)\b', lambda m: CLN(m.group(1)), t)
TXT = {k: (tuple(_cl(x) if isinstance(x, str) else x for x in v) if isinstance(v, tuple) else _cl(v)) for k, v in TXT.items()}
TXT['a'] = ('Ulrich et al.’s C7-fibroblast markers pick out one cluster: cluster 3',) + TXT['a'][1:]
TXT['b'] = (TXT['b'][0][0].upper() + TXT['b'][0][1:],) + TXT['b'][1:]
TXT['e'] = ('In mouse tissue, cells with human neighbors in cluster 3 form the outer fibroblast layer',) + TXT['e'][1:]
# v7.9: cross-references for the new panel letters (bar chart d, mouse UMAPs e, mouse tissue f)
assert TXT['d'][1].count('Left: ') == 1 and TXT['d'][1].count('Right: the mouse') == 1 and TXT['foot'].count('(d, e)') == 1 and TXT['foot'].count('(c, e)') == 1
TXT['d'] = (TXT['d'][0], TXT['d'][1].replace('Left: ', 'd: ').replace('Right: the mouse', 'e: the mouse'), TXT['d'][2])
TXT['foot'] = TXT['foot'].replace('(d, e)', '(d–f)').replace('(c, e)', '(c, f)')
GIDS = ['a_dotplot', 'a_size_legend', 'a_color_legend', 'b_umap_leiden', 'b_umap_C7', 'b_colorbar_C7', 'b_umap_cosine', 'b_colorbar_cosine',
        'c_tissue_menstrual_c3', 'c_tissue_menstrual_cosine', 'c_tissue_secretory_c3', 'c_tissue_secretory_cosine', 'c_colorbar_cosine',
        'd_c3share_bars', 'd_umap_mouse_pop', 'd_umap_c3share', 'd_colorbar_c3share', 'd_umap_cosine', 'd_colorbar_cosine',
        'e_tissue_uninduced_c3share', 'e_tissue_1dpiE_c3share', 'e_tissue_1dpiL_c3share', 'e_tissue_2dpi_c3share', 'e_tissue_3dpi_c3share', 'e_tissue_4dpi_c3share', 'e_colorbar_c3share',
        'e_tissue_uninduced_cosine', 'e_tissue_1dpiE_cosine', 'e_tissue_1dpiL_cosine', 'e_tissue_2dpi_cosine', 'e_tissue_3dpi_cosine', 'e_tissue_4dpi_cosine', 'e_colorbar_cosine']

def build(ANN):
    rng = np.random.default_rng(0)
    B, OFF, Hh = (0.65, 1.55, 12.6) if ANN else (0.0, 0.0, 10.2)
    TOPA = 0.35 if ANN else 0.0                                                      # v7.1: extra room for the larger annotated header text
    hE, AE, AD = 2.45, (1.0 if ANN else 0.0), (1.5 if ANN else 0.0)                 # mouse block geometry (inches)
    yE2 = B + 0.45 + (0.6 if ANN else 0.0); yE1 = yE2 + hE + 0.3; eHead = yE1 + hE + 0.55 + AE; yD = eHead + 1.1; dHead = yD + 1.4 + 0.7 + AD
    MB = dHead + 0.45 - (0.5 + B); Hh = Hh + MB
    W = 17.0; fig = plt.figure(figsize=(W, Hh)); fx, fy = 1 / W, 1 / Hh
    def head(x, y, letter, key=None):
        fig.text(x * fx, y * fy, letter, fontsize=FS_L, fontweight='bold', va='top')
        if ANN and key:
            t, sub, wrap = TXT[key]
            fig.text((x + .35) * fx, y * fy, t, fontsize=FS_H, va='top', fontweight='medium')
            fig.text((x + .35) * fx, (y - .33) * fy, '\n'.join(textwrap.wrap(sub, int(wrap * 9 / FS_P))), fontsize=FS_P, va='top', color='0.25', linespacing=1.3)
    def frame(ax, title):
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('UMAP', fontsize=8.5 * FSC, labelpad=3)
        for sp in ax.spines.values(): sp.set_linewidth(.8)
        ax.set_title(title, fontsize=FS_TT, pad=4)
    def cbar_right(ax, sc, label, ticks=None):
        p = ax.get_position(); cax = fig.add_axes([p.x1 + .06 * fx, p.y0, .085 * fx, p.height])
        cb = fig.colorbar(sc, cax=cax, ticks=ticks); cb.ax.tick_params(labelsize=7.5 * FSC, length=2); cb.outline.set_linewidth(.5); cb.set_label(label, fontsize=8 * FSC)

    # a: dot plot
    yA = 6.35 + OFF + MB; head(0.6, (12.0 + TOPA if ANN else 9.15) + MB, 'a', 'a')
    ng, nk = len(GENES), len(GRP)
    axA = fig.add_axes([1.75 * fx, (yA + (6 - nk) * .35) * fy, (ng * .37 + .3) * fx, (nk * .35 + .35) * fy])        # top-aligned whatever the cluster count
    yy, xx = np.meshgrid(np.arange(nk) + .5, np.arange(ng) + .5, indexing='ij')
    axA.scatter(xx.ravel(), yy.ravel(), s=size_of(FR.to_numpy().ravel()), c=MN.to_numpy().ravel(), cmap='Reds', vmin=0, vmax=1, edgecolors='black', linewidths=.2, zorder=3)
    axA.set_xlim(-.3, ng + .3); axA.set_ylim(nk + .5, -.5); axA.grid(False)
    axA.set_xticks(np.arange(ng) + .5); axA.set_xticklabels(GENES, rotation=90, style='italic', fontsize=FS_K)
    axA.set_yticks(np.arange(nk) + .5); axA.set_yticklabels([CLN(g) for g in GRP], fontsize=FS_K); axA.tick_params(length=2)
    pA = axA.get_position()
    axS = fig.add_axes([pA.x1 + .3 * fx, pA.y1 - .75 * fy, 1.45 * fx, .3 * fy])
    axS.scatter(np.arange(len(size_range)) + .5, np.zeros(len(size_range)), s=size_of(size_range), color='gray', edgecolor='black', linewidth=.2, zorder=3)
    axS.set_xticks(np.arange(len(size_range)) + .5); axS.set_xticklabels([f'{int(round(v * 100))}' for v in size_range], fontsize=7.5 * FSC)
    axS.set_xlim(0, len(size_range)); axS.set_ylim(-.8, .8); axS.tick_params(axis='y', left=False, labelleft=False); axS.tick_params(axis='x', length=0, pad=3)
    for sp in axS.spines.values(): sp.set_visible(False)
    axS.set_title('Fraction of nuclei\nin group (%)', fontsize=8 * FSC, pad=6)
    axC = fig.add_axes([pA.x1 + .3 * fx, pA.y1 - 1.7 * fy, 1.45 * fx, .11 * fy])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap='Reds'), cax=axC, orientation='horizontal', ticks=[0, 1])
    cb.ax.set_xticklabels(['0', '1']); cb.ax.tick_params(labelsize=7.5 * FSC, length=2); cb.outline.set_linewidth(.5)
    axC.set_title('Mean expression in group\n(scaled 0–1 per gene)', fontsize=8 * FSC, pad=5)

    # b: three UMAPs, frames on an equal pitch
    head(6.95, (12.0 + TOPA if ANN else 9.15) + MB, 'b', 'b')
    wU, hU, yU, X0, GAP = 2.3, 1.4, 7.05 + OFF + MB, 7.4, 1.0
    xs = [X0 + k * (wU + GAP) for k in range(3)]
    ax = fig.add_axes([xs[0] * fx, yU * fy, wU * fx, hU * fy]); d = U.iloc[rng.permutation(len(U))]
    PS = 120000 / len(U)
    ax.scatter(d.lin_umap_1, d.lin_umap_2, s=PS, c=[PAL[k] for k in d['leiden_0.5']], linewidths=0, rasterized=True); frame(ax, 'Leiden clusters')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=PAL[k], ms=5 * FSC, label=CLN(k)) for k in sorted(U['leiden_0.5'].unique())],
              loc='center left', bbox_to_anchor=(1.0, .5), frameon=False, fontsize=8 * FSC, handletextpad=.2, labelspacing=.45, borderaxespad=.3, handlelength=1.0)
    ax = fig.add_axes([xs[1] * fx, yU * fy, wU * fx, hU * fy]); d = U.sort_values('C7_log', kind='stable')
    sc = ax.scatter(d.lin_umap_1, d.lin_umap_2, s=PS, c=d.C7_log, cmap=CM_C7, vmin=float(d.C7_log.min()), vmax=float(d.C7_log.max()), linewidths=0, rasterized=True)
    frame(ax, 'C7'); cbar_right(ax, sc, 'log(1 + CP10k)')
    ax = fig.add_axes([xs[2] * fx, yU * fy, wU * fx, hU * fy]); d = U.iloc[rng.permutation(len(U))]
    sc = ax.scatter(d.lin_umap_1, d.lin_umap_2, s=PS, c=d.cos, cmap=CM_COS, norm=mpl.colors.TwoSlopeNorm(0, -clim, clim), linewidths=0, rasterized=True)
    frame(ax, 'contrastive Stro-7 cosine'); ax.set_facecolor(BLK); cbar_right(ax, sc, 'cosine', ticks=[-.1, -.05, 0, .05, .1])

    # c: tissue by timepoint
    head(0.6, (6.5 if ANN else 5.0) + MB, 'c', 'c')
    TP = [('menstrual', 'menstrual donor'), ('secretory', 'secretory donor')]
    xc = [0.95, 4.65, 8.6, 12.3]; yC, hC, wC = 0.5 + B + MB, 4.1, 3.5
    for j, (tp, nm) in enumerate(TP):
        d = XY[XY.timepoint == tp]; o = d[~d.c3]; r = d[d.c3]
        pad = 0.03 * max(np.ptp(d.x_um_qc), np.ptp(d.y_um_qc)); lim = ((d.x_um_qc.min() - pad, d.x_um_qc.max() + pad), (d.y_um_qc.min() - pad, d.y_um_qc.max() + pad))
        sC = DOT_S[f'c_{tp}'] = matched_s(d.x_um_qc, d.y_um_qc, max(np.ptp(lim[0]) / wC, np.ptp(lim[1]) / hC))
        ax = fig.add_axes([xc[j] * fx, yC * fy, wC * fx, hC * fy])
        ax.scatter(o.x_um_qc, o.y_um_qc, s=sC, c=C_OTHER, linewidths=0, rasterized=True); ax.scatter(r.x_um_qc, r.y_um_qc, s=sC, c=PAL[3], linewidths=0, rasterized=True)
        ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1]); ax.set_aspect('equal'); blackbox(ax); ax.set_title(nm, fontsize=FS_TT, color='0.15', pad=4)
        scalebar(ax, lim[0][0], lim[1][0], max(np.ptp(lim[0]), np.ptp(lim[1])))                 # the two donors are drawn at different scales
        if j == 1:
            ax.legend(handles=[Line2D([], [], marker='o', ls='', color=PAL[3], ms=5 * FSC, label='cluster 3'), Line2D([], [], marker='o', ls='', color=C_OTHER, ms=5 * FSC, label='other clusters')],
                      loc='upper right', bbox_to_anchor=(1.0, -0.005), ncol=2, frameon=False, fontsize=8 * FSC, handletextpad=.2, columnspacing=1.0)   # v7.7: below the black box
        ax = fig.add_axes([xc[j + 2] * fx, yC * fy, wC * fx, hC * fy]); dd = d.iloc[rng.permutation(len(d))]
        sc = ax.scatter(dd.x_um_qc, dd.y_um_qc, s=sC, c=dd.cos, cmap=CM_COS, norm=mpl.colors.TwoSlopeNorm(0, -clim, clim), linewidths=0, rasterized=True)
        ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1]); ax.set_aspect('equal'); blackbox(ax); ax.set_title(nm, fontsize=FS_TT, color='0.15', pad=4); lastC = ax
    p = lastC.get_position(); cax = fig.add_axes([p.x1 + .25 * fx, p.y0 + p.height * .3, .1 * fx, p.height * .4])
    cb = fig.colorbar(sc, cax=cax, ticks=[-.1, -.05, 0, .05, .1]); cb.ax.tick_params(labelsize=7.5 * FSC, length=2); cb.outline.set_linewidth(.5); cb.set_label('contrastive Stro-7 cosine', fontsize=8 * FSC)

    # d: mouse, c3 share by population + three UMAPs on b's pitch
    head(0.6, dHead, 'd', 'd'); head(6.95, yD + (1.85 if ANN else 2.1), 'e')                                     # v7.9: mouse UMAPs are their own panel
    axB = fig.add_axes([1.75 * fx, (yD - 0.35) * fy, 3.4 * fx, 1.75 * fy]); yy_ = np.arange(len(MPOP))
    axB.barh(yy_, mstat.c3.values, color=[MCOL[p] for p in MPOP], height=.72, zorder=2)
    axB.set_ylim(len(MPOP) - .4, -.6); axB.set_yticks(yy_); axB.set_yticklabels(MPOP, fontsize=8 * FSC)
    axB.set_xlim(0, 1); axB.set_xticks([0, .25, .5, .75, 1]); axB.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'], fontsize=7.5 * FSC); axB.tick_params(length=2)
    axB.set_xlabel('mean fraction of human SAMap neighbors in cluster 3', fontsize=8.5 * FSC, labelpad=3)
    for sp in ('top', 'right'): axB.spines[sp].set_visible(False)
    PSm = 120000 / len(MU)
    ax = fig.add_axes([xs[0] * fx, yD * fy, wU * fx, hU * fy]); d = MU.iloc[rng.permutation(len(MU))]
    ax.scatter(d.X_umap_1, d.X_umap_2, s=PSm, c=d['pop'].map(MCOL), linewidths=0, rasterized=True); frame(ax, 'mouse populations')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=MCOL[p], ms=4.5 * FSC, label=p) for p in MPOP], loc='center left', bbox_to_anchor=(1.0, .5), frameon=False,
              fontsize=7.5 * FSC, handletextpad=.15, labelspacing=.2, borderaxespad=.3); axDleg = ax
    ax = fig.add_axes([xs[1] * fx, yD * fy, wU * fx, hU * fy]); d = MU.sort_values('c3s', kind='stable')
    sc = ax.scatter(d.X_umap_1, d.X_umap_2, s=PSm, c=d.c3s, cmap=CM_FRAC, vmin=0, vmax=1, linewidths=0, rasterized=True); frame(ax, 'human SAMap neighbors\nin cluster 3'); ax.set_facecolor(BLK); cbar_right(ax, sc, 'fraction', ticks=[0, .5, 1])
    ax = fig.add_axes([xs[2] * fx, yD * fy, wU * fx, hU * fy]); d = MU.iloc[rng.permutation(len(MU))]
    sc = ax.scatter(d.X_umap_1, d.X_umap_2, s=PSm, c=d.mcos, cmap=CM_COS, norm=mpl.colors.TwoSlopeNorm(0, -clim, clim), linewidths=0, rasterized=True)
    frame(ax, 'human SAMap neighbors’\nStro-7 cosine'); ax.set_facecolor(BLK); cbar_right(ax, sc, 'cosine', ticks=[-.1, -.05, 0, .05, .1])
    # e: mouse tissue, six sections, each scaled to fill its square (cropped and rotated per MSEC); top c3 share, bottom Stro-7 cosine
    head(0.6, eHead, 'f', 'e')                                                                  # v7.9: mouse tissue relabeled f (text key stays 'e')
    wE, pitch, x0E = 2.45, 2.55, 0.95
    for row, (yR, key) in enumerate(((yE1, 'c3s'), (yE2, 'mcos'))):
        for j, (tp, nm) in enumerate(MTP):
            g = MXY[MXY.timepoint == tp]; bg = MBG[MBG.timepoint == tp]; ex0, ex1, ey0, ey1 = EBOX[tp]
            sE = DOT_S[f'e_{tp}'] = matched_s(np.r_[g.x_disp, bg.x_disp], np.r_[g.y_disp, bg.y_disp], (ex1 - ex0) / wE)
            ax = fig.add_axes([(x0E + j * pitch) * fx, yR * fy, wE * fx, hE * fy])
            ax.scatter(bg.x_disp, bg.y_disp, s=sE, c=BG_DARK, linewidths=0, rasterized=True)
            if key == 'c3s':
                gg = g.sort_values('c3s', kind='stable'); sc = ax.scatter(gg.x_disp, gg.y_disp, s=sE, c=gg.c3s, cmap=CM_FRAC, vmin=0, vmax=1, linewidths=0, rasterized=True)
            else:
                gg = g.iloc[rng.permutation(len(g))]; sc = ax.scatter(gg.x_disp, gg.y_disp, s=sE, c=gg.mcos, cmap=CM_COS, norm=mpl.colors.TwoSlopeNorm(0, -clim, clim), linewidths=0, rasterized=True)
            ax.set_xlim(ex0, ex1); ax.set_ylim(ey0, ey1); ax.set_aspect('equal'); blackbox(ax)
            if row == 0: ax.set_title(nm, fontsize=FS_TT, color='0.15', pad=3)
            if row == 1: scalebar(ax, ex0, ey0, ex1 - ex0)                                        # each section has its own scale
            lastE = ax
        p = lastE.get_position(); cax = fig.add_axes([p.x1 + .2 * fx, p.y0 + p.height * .2, .09 * fx, p.height * .6])
        cb = fig.colorbar(sc, cax=cax, ticks=[0, .5, 1] if key == 'c3s' else [-.1, -.05, 0, .05, .1]); cb.ax.tick_params(labelsize=7.5 * FSC, length=2); cb.outline.set_linewidth(.5)
        cb.set_label('fraction of human SAMap\nneighbors in cluster 3' if key == 'c3s' else 'human SAMap neighbors’\nStro-7 cosine', fontsize=8 * FSC)
    if ANN:
        fig.text(0.6 * fx, -0.25 * fy, '\n'.join(sum((textwrap.wrap(par, int(230 * 8 / FS_F)) for par in TXT['foot'].split('\n')), [])), fontsize=FS_F, color='0.3', va='bottom', linespacing=1.3)

    # checks, ids, save
    assert len(fig.axes) == len(GIDS), len(fig.axes)
    assert [fig.axes[k].get_title() for k in (3, 4, 6, 8, 9, 10, 11)] == ['Leiden clusters', 'C7', 'contrastive Stro-7 cosine', 'menstrual donor', 'menstrual donor', 'secretory donor', 'secretory donor']
    for a_, g_ in zip(fig.axes, GIDS): a_.set_gid(g_)
    fig.canvas.draw(); r_ = fig.canvas.get_renderer()
    tx = [(t, t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    OVL = [(a.get_text()[:22], b.get_text()[:22]) for i, (a, ba) in enumerate(tx) for b, bb in tx[i + 1:] if ba.overlaps(bb) and a.get_text() != b.get_text()]
    CB = [a for a in fig.axes if a.get_label() == '<colorbar>']
    COL = [(i, j) for i, a in enumerate(CB) for j, b in enumerate(fig.axes) if b not in CB and a.get_window_extent(r_).overlaps(b.get_window_extent(r_))]
    um = [fig.axes[k].get_window_extent(r_) for k in (3, 4, 6)]; dpi = fig.dpi
    gaps_in = [round((um[k + 1].x0 - um[k].x1) / dpi, 3) for k in range(2)]
    # clearance between each UMAP's legend/colorbar block and the next frame
    acc = [fig.axes[3].get_legend().get_window_extent(r_), fig.axes[5].get_tightbbox(r_)]
    clear_in = [round((um[k + 1].x0 - acc[k].x1) / dpi, 3) for k in range(2)]
    dl = axDleg.get_legend().get_window_extent(r_); mclear = round((fig.axes[GIDS.index('d_umap_c3share')].get_window_extent(r_).x0 - dl.x1) / dpi, 3)
    BLACK = [a.get_window_extent(r_) for a in fig.axes if tuple(np.round(a.get_facecolor()[:3], 3)) == (0, 0, 0)]       # v7.7: text or legends must not sit on black panels
    ON_BLACK = [t.get_text()[:22] for t, bt in tx if any(bt.overlaps(b) for b in BLACK)] + [f'legend of axes {i}' for i, a in enumerate(fig.axes) if a.get_legend() and any(a.get_legend().get_window_extent(r_).overlaps(b) for b in BLACK)]
    tag = 'annotated' if ANN else 'panels'
    common.save_figure(fig, OUTPATHS[tag]); plt.close(fig)
    root = ET.parse(OUTPATHS[tag]['svg']).getroot(); ns = {'s': 'http://www.w3.org/2000/svg'}
    G = {g.get('id'): g for g in root.iter('{http://www.w3.org/2000/svg}g') if g.get('id') in GIDS}
    imgs = {k: len(G[k].findall('.//s:image', ns)) for k in GIDS if k.startswith(('b_umap', 'c_tissue', 'd_umap', 'e_tissue'))}
    return dict(version=tag, dot_s={k: round(float(v), 2) for k, v in DOT_S.items()}, on_black=ON_BLACK, old_cluster_names=sorted({t.get_text()[:40] for t, _ in tx if re.search(r'\bc\d\b', t.get_text())}), mouse_legend_to_next_frame_in=mclear, overlaps=OVL, colorbar_collisions=COL, umap_frame_gaps_in=gaps_in, legend_or_cbar_to_next_frame_in=clear_in, point_images_per_group=imgs,
                missing_groups=[k for k in GIDS if k not in G])


FONT = common.set_font()
REPORT = [build(False), build(True)]
for r in REPORT: LOG.info('%s', r)
inputs = {'human_table': ARGS.human_table, 'dotplot_table': ARGS.dotplot_table, 'mouse_table': ARGS.mouse_table, 'background_table': ARGS.background_table,
          'mouse_colors': ARGS.mouse_colors, 'mouse_crops': ARGS.mouse_crops, 'section_rotation': ARGS.section_rotation}
record = common.base_record(Path(__file__), inputs, {f'{v}_{k}': p for v, d in OUTPATHS.items() for k, p in d.items()})
record['font'] = FONT
record['mouse_rotation'] = section_rotation.table_record(ARGS.section_rotation)
record['mouse_sections'] = {tp: {'crop_um': MSEC[tp], 'rotation_deg_ccw': ANGLES[tp], 'nuclei_placed': MSEC_N[tp][0], 'nuclei_shown': MSEC_N[tp][1]} for tp in MSEC}
record['counts'] = {'human_nuclei': n_all, 'human_cluster3': n_c3, 'human_positioned': int(len(XY)), 'C7_detected': n_c7, 'C7_detected_in_cluster3': n_c7_c3,
                    'mouse_cells': int(len(MX)), 'mouse_on_umap': int(len(MU)), 'mouse_in_sections': int(len(MXY)), 'mouse_background_in_sections': int(len(MBG))}
record['values'] = {'cosine_median_by_cluster': {int(k): float(v) for k, v in cos_med.items()}, 'cosine_clim': clim,
                    'cluster3_share_of_positioned_by_donor_pct': pc3_tp, 'cluster3_fraction_of_human_nuclei': base_c3,
                    'mouse_by_population': mstat.round(6).reset_index().to_dict(orient='records')}
record['variants'] = REPORT
common.write_json(ARGS.provenance, record)
_bad = {r['version']: {k: r[k] for k in ('missing_groups', 'on_black', 'colorbar_collisions', 'old_cluster_names') if r[k]} for r in REPORT}
_bad = {k: v for k, v in _bad.items() if v}
if _bad:
    sys.exit(f'panels_c7_presence: layout checks failed: {_bad}')
LOG.info('wrote %s and %s', [str(p) for d in OUTPATHS.values() for p in d.values()], ARGS.provenance)
