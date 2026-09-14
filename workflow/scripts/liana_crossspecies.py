#!/usr/bin/env python
"""liana_crossspecies.py — cross-species LIANA conservation deliverable.

Snakemake port (2026-09-10) of the deliverable spine of
  decidualization/revisions/liana_pipeline_final/cross_species_lr_explorer_v8-2.ipynb (cells 5-31)
  + lr_analysis.py (orthology key, loading, conservation, Venn/bar figures).

Changes vs the June 2026 notebook
  * Orthology: Ensembl Compara table (config liana.orthology), homology_type == ortholog_one2one,
    no manual overrides. Was JAX 1:1 + 9 curated overrides via ortholog_utils.py (module no longer present).
  * Symmetric orthology test: a HUMAN subunit must also be in the one2one set for the row to be
    'has_ortholog'. The June code uppercased human symbols unconditionally, so every human row was
    flagged orthologous even when the gene has no mouse one2one partner.
  * Rank column / cutoff come from config (liana.rank_column, liana.rank_cutoff), default magnitude_rank < 0.05.
  * Cell-type vocabulary = 2026_08_28 objects (cell_type_merged from resources/liana_cell_type_map.yaml);
    compartment rollup + display names are read from the same YAML (compartment_map, display_names).
  * No human->mouse label harmonization needed (labels are shared by construction).
  * Figures: fig1 Venn, fig3 matched-vs-rewired bar, fig2 compartment heatmaps in pooled (immune pooled)
    and expanded (immune lineages separate) views: 5 conservation stats + per-species counts.
"""
from __future__ import annotations
import argparse, json, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.stats import hypergeom
from openpyxl.utils import get_column_letter


# ----------------------------------------------------------------------------- orthology
def load_one2one(tsv: Path, homology_type: str) -> tuple[dict, set, dict]:
    d = pd.read_csv(tsv, sep='\t')
    need = {'component_id', 'species', 'gene_symbol', 'homology_type'}
    assert need <= set(d.columns), f'orthology TSV missing {need - set(d.columns)}'
    o = d[d.homology_type == homology_type]
    n_nosym = int(o['gene_symbol'].isna().sum())
    o = o.dropna(subset=['gene_symbol'])          # Ensembl-ID-only genes cannot join symbol-keyed LIANA tables
    piv = o.pivot_table(index='component_id', columns='species', values='gene_symbol', aggfunc=list)
    pairs, bad = [], 0
    for _, row in piv.iterrows():
        hs, ms = row.get('human'), row.get('mouse')
        if not isinstance(hs, list) or not isinstance(ms, list) or len(hs) != 1 or len(ms) != 1:
            bad += 1; continue
        pairs.append((ms[0], hs[0]))
    # strict 1:1 at the SYMBOL level: a symbol present in >1 one2one component (patch / readthrough loci)
    # is ambiguous and is dropped on both sides rather than letting one component win silently.
    mcount = pd.Series([m for m, _ in pairs]).value_counts(); hcount = pd.Series([h for _, h in pairs]).value_counts()
    dup = [(m, h) for m, h in pairs if mcount[m] > 1 or hcount[h] > 1]
    m2h = {m: h for m, h in pairs if mcount[m] == 1 and hcount[h] == 1}
    human_set = set(m2h.values())
    meta = {'orthology_file': str(tsv), 'homology_type': homology_type,
            'source_release': sorted(d['source_release'].dropna().unique().tolist()) if 'source_release' in d else None,
            'n_components': len(pairs), 'n_pairs': len(m2h), 'n_genes_without_symbol': n_nosym, 'n_components_not_1_per_species': bad,
            'n_dropped_duplicate_symbol': len(dup), 'dropped_duplicate_symbols': sorted(set(dup))}
    return m2h, human_set, meta


def to_canonical(complex_str, species: str, m2h: dict, human_set: set):
    if pd.isna(complex_str):
        return None
    subs = str(complex_str).split('_')
    if species == 'mouse':
        out = [m2h.get(s) for s in subs]
    else:
        out = [s if s in human_set else None for s in subs]
    return None if any(x is None for x in out) else '_'.join(out)


def add_canonical_keys(df: pd.DataFrame, species: str, m2h: dict, human_set: set) -> pd.DataFrame:
    df['ligand_canonical'] = df['ligand_complex'].map(lambda x: to_canonical(x, species, m2h, human_set))
    df['receptor_canonical'] = df['receptor_complex'].map(lambda x: to_canonical(x, species, m2h, human_set))
    ok = df['ligand_canonical'].notna() & df['receptor_canonical'].notna()
    df['lr_canonical'] = np.where(ok, df['ligand_canonical'].astype(str) + '^' + df['receptor_canonical'].astype(str), None)
    return df


def unmappable_subunits(df: pd.DataFrame, species: str, m2h: dict, human_set: set) -> pd.Series:
    dropped = df[df['lr_canonical'].isna()]
    subs = []
    for col in ('ligand_complex', 'receptor_complex'):
        for v in dropped[col].dropna().astype(str):
            for s in v.split('_'):
                if (species == 'mouse' and s not in m2h) or (species == 'human' and s not in human_set):
                    subs.append(s)
    return pd.Series(subs, dtype=object).value_counts()


# ----------------------------------------------------------------------------- loading
def load_all_timepoints(root: Path, species: str, ctk: str, timepoints: list[str], log) -> pd.DataFrame:
    frames = []
    for tp in timepoints:
        p = root / species / tp / ctk / 'tables' / 'liana_res_full.csv'
        if not p.exists():
            sys.exit(f'missing LIANA table: {p}')
        df = pd.read_csv(p, index_col=0)
        df['timepoint'] = tp
        frames.append(df)
        log(f'  {species:6s} {tp:16s} {len(df):>7,} rows  {df.source.nunique():>2} sources  {df.target.nunique():>2} targets')
    return pd.concat(frames, ignore_index=True)


# ----------------------------------------------------------------------------- conservation
def sig_triples(df: pd.DataFrame, rank_col: str, cut: float) -> set:
    s = df[df[rank_col] < cut].dropna(subset=['lr_canonical'])
    return set(map(tuple, s[['lr_canonical', 'source', 'target']].itertuples(index=False, name=None)))


def hypergeometric(mouse_t, human_t, mouse_lr, human_lr, conserved) -> dict:
    universe = set(mouse_t['lr_canonical'].dropna()) & set(human_t['lr_canonical'].dropna())
    M = len(universe); n = len(mouse_lr & universe); N = len(human_lr & universe); k = len(conserved & universe)
    expected = n * N / M if M else float('nan')
    return {'M': M, 'n': n, 'N': N, 'k': k, 'expected': expected,
            'fold': (k / expected) if expected else float('nan'),
            'p_value': float(hypergeom.sf(k - 1, M, n, N)) if M else float('nan')}


def per_cell_type_pair_table(mouse_sig, human_sig) -> pd.DataFrame:
    def sets(df):
        return {k: set(g['lr_canonical'].dropna()) for k, g in df.groupby(['source', 'target'], sort=False)}
    ms, hs = sets(mouse_sig), sets(human_sig)
    rows = []
    for src in sorted(set(mouse_sig.source) & set(human_sig.source)):
        for tgt in sorted(set(mouse_sig.target) & set(human_sig.target)):
            m, h = ms.get((src, tgt), set()), hs.get((src, tgt), set())
            rows.append({'source': src, 'target': tgt, 'mouse_only': len(m - h), 'conserved': len(m & h),
                         'human_only': len(h - m)})
    df = pd.DataFrame(rows)
    df['mouse_total'] = df.mouse_only + df.conserved; df['human_total'] = df.human_only + df.conserved
    df['union'] = df.mouse_only + df.conserved + df.human_only
    df['min_total'] = df[['mouse_total', 'human_total']].min(axis=1)
    df['jaccard_pct'] = (100 * df.conserved / df['union'].clip(lower=1)).round(1)
    df['overlap_coef_pct'] = (100 * df.conserved / df.min_total.clip(lower=1)).round(1)
    return df


# ----------------------------------------------------------------------------- heatmaps
PCT_CMAP, COUNT_CMAP = 'Reds', 'magma_r'


def build_grids(m_trip, h_trip, matched, rollup: dict, order: list, drop: set, display: dict):
    fine = {x for tr in (m_trip, h_trip) for _, s, t in tr for x in (s, t)}
    missing = fine - set(rollup)
    if missing:
        raise KeyError(f'labels missing from compartment map: {sorted(missing)}')
    present = {rollup[l] for l in fine}
    groups = [g for g in order if g in present] + sorted(present - set(order))
    groups = [g for g in groups if g not in drop]
    ix = {g: i for i, g in enumerate(groups)}; ng = len(groups)
    ml, hl, cl = defaultdict(set), defaultdict(set), defaultdict(set)
    for lr, s, t in m_trip: ml[(rollup[s], rollup[t])].add(lr)
    for lr, s, t in h_trip: hl[(rollup[s], rollup[t])].add(lr)
    for lr, s, t in matched: cl[(rollup[s], rollup[t])].add(lr)
    Z = lambda: np.full((ng, ng), np.nan)
    Mc, Mj, Mo, Mm, Mh, Cm, Ch = np.zeros((ng, ng)), Z(), Z(), Z(), Z(), np.zeros((ng, ng)), np.zeros((ng, ng))
    for a in groups:
        for b in groups:
            i, j = ix[a], ix[b]
            c, mt, ht = len(cl[(a, b)]), len(ml[(a, b)]), len(hl[(a, b)])
            union = len(ml[(a, b)] | hl[(a, b)])
            Mc[i, j] = c; Cm[i, j] = mt; Ch[i, j] = ht
            if union: Mj[i, j] = 100 * c / union
            if min(mt, ht): Mo[i, j] = 100 * c / min(mt, ht)
            if mt: Mm[i, j] = 100 * c / mt
            if ht: Mh[i, j] = 100 * c / ht
    labels = [display.get(g, g) for g in groups]
    return labels, {'M_count': Mc, 'jaccard': Mj, 'overlap_coef': Mo, 'pct_mouse_in_human': Mm,
                    'pct_human_in_mouse': Mh, 'count_mouse': Cm, 'count_human': Ch}


def draw_heat(ax, data, labels, title, cmap, vmax, is_pct):
    m = len(labels)
    im = ax.imshow(data, cmap=cmap, vmin=0, vmax=vmax, aspect='equal')
    ax.set_xticks(range(m)); ax.set_yticks(range(m))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=13); ax.set_yticklabels(labels, fontsize=13)
    ax.set_xlabel('target', fontsize=14); ax.set_ylabel('source', fontsize=14); ax.set_title(title, fontsize=15, pad=10)
    hi = vmax if vmax else (np.nanmax(data) or 1)
    for i in range(m):
        for j in range(m):
            v = data[i, j]
            if np.isnan(v):
                ax.text(j, i, '·', ha='center', va='center', color='lightgray', fontsize=13); continue
            ax.text(j, i, f"{v:.0f}{'%' if is_pct else ''}", ha='center', va='center', fontsize=11,
                    fontweight='bold', color='white' if v > 0.55 * hi else 'black')
    return im


SPECS = [('M_count', 'Conserved L-R count', COUNT_CMAP, None, False, '# conserved pairs'),
         ('jaccard', 'Jaccard %  (symmetric, size-sensitive)', PCT_CMAP, 100, True, '|M∩H| / |M∪H|  (%)'),
         ('overlap_coef', 'Overlap coefficient %  (containment)', PCT_CMAP, 100, True, '|M∩H| / min(|M|,|H|)  (%)'),
         ('pct_mouse_in_human', '% mouse interactions detected in human', PCT_CMAP, 100, True, '|M∩H| / |M|  (%)'),
         ('pct_human_in_mouse', '% human interactions detected in mouse', PCT_CMAP, 100, True, '|M∩H| / |H|  (%)')]


def save_heatmaps(figdir: Path, views: dict, rank_col: str, cut: float):
    count_vmax = max(np.nanmax(v[1]['M_count']) for v in views.values()) or 1
    sp_vmax = max(max(np.nanmax(v[1]['count_mouse']), np.nanmax(v[1]['count_human'])) for v in views.values()) or 1
    for view, (labels, mats) in views.items():
        m = len(labels)
        for key, title, cmap, vmax, is_pct, cbl in SPECS:
            fig, ax = plt.subplots(figsize=(1.05 * m + 3.0, 0.95 * m + 2.2))
            im = draw_heat(ax, mats[key], labels, title, cmap, count_vmax if key == 'M_count' else vmax, is_pct)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label(cbl, fontsize=12)
            for ext in ('pdf', 'png'):
                fig.savefig(figdir / f'fig2_{view}_{key}.{ext}', dpi=200, bbox_inches='tight')
            plt.close(fig)
        fig, axes = plt.subplots(2, 3, figsize=(24, 15)); flat = axes.flatten()
        for k, (key, title, cmap, vmax, is_pct, cbl) in enumerate(SPECS):
            im = draw_heat(flat[k], mats[key], labels, title, cmap, count_vmax if key == 'M_count' else vmax, is_pct)
            cb = fig.colorbar(im, ax=flat[k], fraction=0.046, pad=0.04); cb.set_label(cbl, fontsize=11)
        flat[5].axis('off')
        fig.suptitle(f'Cross-species conservation by compartment ({view}; strict matched-context, {rank_col} < {cut})',
                     fontsize=17, y=1.0)
        fig.tight_layout(rect=[0, 0.02, 1, 0.96])
        fig.savefig(figdir / f'fig2_{view}_conservation_heatmaps_combined.pdf', bbox_inches='tight'); plt.close(fig)
        # per-species counts
        fig, axes = plt.subplots(1, 2, figsize=(2 * (1.05 * m + 3.0), 0.95 * m + 2.2))
        for ax, sp in zip(axes, ('mouse', 'human')):
            im = draw_heat(ax, mats[f'count_{sp}'], labels, f'{sp.capitalize()} L-R interactions by compartment',
                           COUNT_CMAP, sp_vmax, False)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label('# significant L-R pairs', fontsize=11)
            fig_s, ax_s = plt.subplots(figsize=(1.05 * m + 3.0, 0.95 * m + 2.2))
            im_s = draw_heat(ax_s, mats[f'count_{sp}'], labels, f'{sp.capitalize()} L-R interactions by compartment',
                             COUNT_CMAP, sp_vmax, False)
            cb_s = fig_s.colorbar(im_s, ax=ax_s, fraction=0.046, pad=0.04); cb_s.set_label('# significant L-R pairs', fontsize=12)
            for ext in ('pdf', 'png'):
                fig_s.savefig(figdir / f'fig2_count_{view}_{sp}.{ext}', dpi=200, bbox_inches='tight')
            plt.close(fig_s)
        fig.tight_layout()
        for ext in ('pdf', 'png'):
            fig.savefig(figdir / f'fig2_count_{view}_by_species.{ext}', dpi=200, bbox_inches='tight')
        plt.close(fig)


def plot_venn(n_m, n_c, n_h, rank_col, cut, hg, path):
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.add_patch(Circle((-0.6, 0), 1.0, alpha=0.55, facecolor='#4e79a7', edgecolor='black', linewidth=1.5))
    ax.add_patch(Circle((0.6, 0), 1.0, alpha=0.55, facecolor='#e15759', edgecolor='black', linewidth=1.5))
    ax.text(-1.15, 0, f'{n_m}\nmouse only', ha='center', va='center', fontsize=16, fontweight='bold')
    ax.text(1.15, 0, f'{n_h}\nhuman only', ha='center', va='center', fontsize=16, fontweight='bold')
    ax.text(0, 0, f'{n_c}\nconserved', ha='center', va='center', fontsize=18, fontweight='bold')
    ax.text(-1.0, 1.25, 'Mouse', ha='center', fontsize=18, fontweight='bold')
    ax.text(1.0, 1.25, 'Human', ha='center', fontsize=18, fontweight='bold')
    ax.set_xlim(-2, 2); ax.set_ylim(-1.5, 1.6); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f'Cross-species ligand-receptor conservation ({rank_col} < {cut})\n'
                 f"{hg['fold']:.2f}× enrichment, p = {hg['p_value']:.1e} (hypergeom, universe = {hg['M']:,})",
                 fontsize=22, pad=12)
    fig.tight_layout(); fig.savefig(path, dpi=300, bbox_inches='tight'); plt.close(fig)


def plot_matched_vs_rewired(matched, rewired, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    counts = [len(matched), len(rewired)]; total = sum(counts) or 1
    labels = [f'Matched-context\n({counts[0]}, {100*counts[0]/total:.0f}%)', f'Rewired\n({counts[1]}, {100*counts[1]/total:.0f}%)']
    bars = ax.bar(labels, counts, color=['#59a14f', '#f28e2b'], edgecolor='white', linewidth=1.5)
    for b, cnt in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1, str(cnt), ha='center', fontweight='bold')
    ax.set_ylabel('Conserved LR pairs'); ax.set_title(f'Matched-context vs rewired decomposition (total = {sum(counts)})')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout(); fig.savefig(path, dpi=300, bbox_inches='tight'); plt.close(fig)


# ----------------------------------------------------------------------------- deliverable tables
def build_deliverables(mouse, human, rank_col, cut, display: dict, log):
    spec_col = 'specificity_rank' if 'specificity_rank' in mouse.columns else None

    def augment(df, sp):
        s = df.copy()
        s['has_ortholog'] = s['lr_canonical'].notna()
        s['ligand'] = np.where(s['has_ortholog'], s['ligand_canonical'], s['ligand_complex'])
        s['receptor'] = np.where(s['has_ortholog'], s['receptor_canonical'], s['receptor_complex'])
        s['lr_key'] = np.where(s['has_ortholog'], s['lr_canonical'],
                               sp + '::' + s['ligand'].astype(str) + '^' + s['receptor'].astype(str))
        return s
    m_all, h_all = augment(mouse, 'mouse'), augment(human, 'human')
    m_h, h_h = m_all[m_all[rank_col] < cut].copy(), h_all[h_all[rank_col] < cut].copy()

    def raw(s, sp):
        cols = ['ligand', 'receptor', 'source', 'target', 'timepoint', 'magnitude_rank', 'has_ortholog'] + ([spec_col] if spec_col else [])
        r = s[cols].copy(); r.insert(0, 'species', sp)
        return r.rename(columns={'source': 'source_cell_type', 'target': 'target_cell_type'})
    raw_calls = pd.concat([raw(m_h, 'mouse'), raw(h_h, 'human')], ignore_index=True)

    KEYS = ['lr_key', 'ligand', 'receptor', 'source', 'target', 'has_ortholog']

    def agg(s, with_tp):
        g = s.groupby(KEYS, dropna=False)
        spec = {'min_magnitude_rank': ('magnitude_rank', 'min')}
        if with_tp: spec['timepoints'] = ('timepoint', lambda x: ', '.join(sorted(x.unique())))
        if spec_col: spec['min_specificity_rank'] = (spec_col, 'min')
        return g.agg(**spec).reset_index()
    m_ct, h_ct = agg(m_h, True), agg(h_h, True)
    m_fill, h_fill = agg(m_all, False), agg(h_all, False)

    hits = m_ct.merge(h_ct, on=KEYS, how='outer', suffixes=('_mouse', '_human'))
    hits['significant_species'] = np.select(
        [hits['min_magnitude_rank_mouse'].notna() & hits['min_magnitude_rank_human'].notna(),
         hits['min_magnitude_rank_mouse'].notna()], ['both', 'mouse'], default='human')
    idx = hits.set_index(KEYS).index
    for side, fill in (('mouse', m_fill), ('human', h_fill)):
        for base in ['min_magnitude_rank'] + (['min_specificity_rank'] if spec_col else []):
            col = f'{base}_{side}'
            lut = fill.set_index(KEYS)[base].reindex(idx).values
            mask = hits[col].isna().values
            hits.loc[mask, col] = lut[mask]
    _m = set(hits.loc[hits.significant_species.isin(['both', 'mouse']), 'lr_key'])
    _h = set(hits.loc[hits.significant_species.isin(['both', 'human']), 'lr_key'])
    conserved_keys = _m & _h
    hits['is_conserved_lr'] = hits['lr_key'].isin(conserved_keys)
    hits = hits.rename(columns={'source': 'source_cell_type', 'target': 'target_cell_type'})
    cols = ['ligand', 'receptor', 'source_cell_type', 'target_cell_type', 'significant_species', 'has_ortholog',
            'is_conserved_lr', 'min_magnitude_rank_mouse', 'min_magnitude_rank_human']
    if spec_col: cols += ['min_specificity_rank_mouse', 'min_specificity_rank_human']
    cols += ['timepoints_mouse', 'timepoints_human']
    hits = hits[cols].sort_values(['significant_species', 'has_ortholog', 'min_magnitude_rank_mouse',
                                   'min_magnitude_rank_human']).reset_index(drop=True)

    # pooled
    mp, hp = set(m_h.lr_key), set(h_h.lr_key)
    keyinfo = (pd.concat([m_h, h_h]).groupby('lr_key')
               .agg(ligand=('ligand', 'first'), receptor=('receptor', 'first'), has_ortholog=('has_ortholog', 'first')).reset_index())
    pooled = keyinfo.copy()
    pooled['conservation'] = np.select([pooled.lr_key.isin(mp & hp), pooled.lr_key.isin(mp - hp)],
                                       ['conserved', 'mouse_only'], default='human_only')
    for sp, sig, allc in (('mouse', m_h, m_all), ('human', h_h, h_all)):
        pooled[f'min_magnitude_rank_{sp}'] = pooled.lr_key.map(sig.groupby('lr_key')['magnitude_rank'].min())
        pooled[f'min_magnitude_rank_{sp}'] = pooled[f'min_magnitude_rank_{sp}'].fillna(
            pooled.lr_key.map(allc.groupby('lr_key')['magnitude_rank'].min()))
        pooled[f'{sp}_timepoints'] = pooled.lr_key.map(sig.groupby('lr_key')['timepoint'].apply(lambda x: ', '.join(sorted(x.unique())))).fillna('')
    pooled = pooled[['ligand', 'receptor', 'conservation', 'has_ortholog', 'min_magnitude_rank_mouse', 'min_magnitude_rank_human',
                     'mouse_timepoints', 'human_timepoints']].sort_values(
        ['conservation', 'has_ortholog', 'min_magnitude_rank_mouse', 'min_magnitude_rank_human']).reset_index(drop=True)

    def relabel(df, cols):
        df = df.copy()
        for c in cols:
            miss = set(df[c].dropna().unique()) - set(display)
            if miss: raise KeyError(f'{c}: labels missing from display_names: {sorted(miss)}')
            df[c] = df[c].map(display)
        return df
    raw_calls = relabel(raw_calls, ['source_cell_type', 'target_cell_type'])
    hits = relabel(hits, ['source_cell_type', 'target_cell_type'])

    # consistency checks (notebook cell 31)
    LR = ['ligand', 'receptor']
    S = lambda df, mask=None: set(map(tuple, (df if mask is None else df[mask])[LR].itertuples(index=False, name=None)))
    po = S(pooled, pooled.conservation == 'conserved')
    hi = S(hits, hits.significant_species.isin(['both', 'mouse'])) & S(hits, hits.significant_species.isin(['both', 'human']))
    checks = {'pooled_eq_hits': po == hi,
              'conserved_all_ortholog': bool(pooled.loc[pooled.conservation == 'conserved', 'has_ortholog'].all()),
              'no_nonortholog_conserved': not bool(hits.loc[~hits.has_ortholog, 'is_conserved_lr'].any()),
              'flag_eq_set': S(hits, hits.is_conserved_lr) == hi}
    for k, v in checks.items(): log(f"  [{'PASS' if v else 'FAIL'}] {k}")
    if not all(checks.values()):
        sys.exit('deliverable consistency check failed')
    return hits, pooled, raw_calls, spec_col


def autosize(ws, df, padding=2, max_width=60):
    for i, col in enumerate(df.columns, start=1):
        n = int(df[col].astype(str).map(len).max()) if len(df) else 0
        ws.column_dimensions[get_column_letter(i)].width = min(max(len(str(col)), n) + padding, max_width)


def format_ranks(ws, df, fmt='0.000E+00'):
    for c in [c for c in df.columns if 'magnitude_rank' in c or 'specificity_rank' in c]:
        for cell in ws[get_column_letter(df.columns.get_loc(c) + 1)][1:]:
            cell.number_format = fmt


def write_workbook(path, hits, pooled, rank_col, cut, ortho_meta, ctk):
    n = len(hits) + 1
    cl = lambda name: get_column_letter(hits.columns.get_loc(name) + 1)
    A, B, C, Sc = (f'lr_hits!${cl(c)}$2:${cl(c)}${n}' for c in ('ligand', 'receptor', 'is_conserved_lr', 'significant_species'))
    f_cons = f'=SUMPRODUCT(({C}=TRUE)/COUNTIFS({A},{A},{B},{B},{C},{C}))'
    f_tot = f'=SUMPRODUCT(1/COUNTIFS({A},{A},{B},{B}))'
    f_ctx = f'=COUNTIF({Sc}, "both")'
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        pd.DataFrame().to_excel(xw, sheet_name='README', index=False, header=False)
        hits.to_excel(xw, sheet_name='lr_hits', index=False); pooled.to_excel(xw, sheet_name='lr_pooled', index=False)
        autosize(xw.sheets['lr_hits'], hits); autosize(xw.sheets['lr_pooled'], pooled)
        format_ranks(xw.sheets['lr_hits'], hits); format_ranks(xw.sheets['lr_pooled'], pooled)
        ws = xw.sheets['README']
        lines = [
            ('CROSS-SPECIES LIGAND-RECEPTOR INTERACTIONS (LIANA)', None),
            (f'generated {pd.Timestamp.now():%Y-%m-%d %H:%M}  -  {rank_col} cutoff < {cut}  -  cell-type key: {ctk}', None),
            (f"orthology: {ortho_meta['source_release']} Compara, homology_type = {ortho_meta['homology_type']} "
             f"({ortho_meta['n_pairs']:,} human-mouse pairs), no manual overrides", None),
            ('', None), ('SHEETS', None),
            ('  lr_hits    L-R interactions significant in >=1 species; 1 row per (ligand, receptor,', None),
            ('             source, target), both species side by side. significant_species = both/mouse/human.', None),
            ('  lr_pooled  1 row per L-R pair, pooled across cell types, with the', None),
            ('             conserved / mouse_only / human_only verdict.', None),
            ('  (raw_calls.csv is a separate file: every significant call, per timepoint.)', None),
            ('', None), ('HEADLINE NUMBERS   (live - click the value to see the formula)', None),
            ('  conserved L-R pairs (both species)', f_cons), ('  conserved cell-type contexts', f_ctx),
            ('  total L-R pairs (either species)', f_tot),
            ('', None), ('KEY COLUMNS', None),
            ('  significant_species   both / mouse / human', None),
            ('  has_ortholog          TRUE  = every subunit has a 1:1 ortholog (cross-species comparable).', None),
            ('                        FALSE = no 1:1 ortholog in one species; one species only, never "conserved".', None),
            ('  is_conserved_lr       TRUE  = L-R pair significant in both species (pooled).', None),
            (f'  min_{rank_col}_*  strongest rank for that species. Populated even when NOT', None),
            ('                        significant (value > cutoff shows how close it came);', None),
            ('                        blank = pair was not testable in that species/context.', None),
            ('', None), ('NOTES', None),
            ('  - An L-R pair spans MULTIPLE rows in lr_hits (one per cell-type context).', None),
            ('    Do not count rows to count pairs - use the headline formulas above.', None),
            ('  - lr_pooled ranks are the best (min) across contexts, so they can be lower than', None),
            ('    the matched lr_hits rank for the same pair.', None),
            ('  - Cross-species comparison is restricted to Ensembl one-to-one orthologs.', None),
        ]
        for r, (label, formula) in enumerate(lines, start=1):
            ws.cell(row=r, column=1, value=label)
            if formula is not None: ws.cell(row=r, column=2, value=formula)
        ws.column_dimensions['A'].width = 82; ws.column_dimensions['B'].width = 22


def write_immune(path_xlsx, path_csv, hits, immune_labels: set, rank_col, cut):
    mask = hits.source_cell_type.isin(immune_labels) | hits.target_cell_type.isin(immune_labels)
    keep = [c for c in ['ligand', 'receptor', 'source_cell_type', 'target_cell_type', 'significant_species', 'has_ortholog',
                        'min_magnitude_rank_mouse', 'min_magnitude_rank_human', 'timepoints_mouse', 'timepoints_human'] if c in hits.columns]
    imm = hits.loc[mask, keep].copy()
    imm.to_csv(path_csv, index=False)
    with pd.ExcelWriter(path_xlsx, engine='openpyxl') as xw:
        pd.DataFrame().to_excel(xw, sheet_name='README', index=False, header=False)
        imm.to_excel(xw, sheet_name='lr_hits_immune', index=False)
        autosize(xw.sheets['lr_hits_immune'], imm); format_ranks(xw.sheets['lr_hits_immune'], imm)
        ws = xw.sheets['README']
        for r, txt in enumerate([
            'IMMUNE-INVOLVING LIGAND-RECEPTOR INTERACTIONS (LIANA)',
            'subset of lr_hits where source OR target is an immune cell type.', f'{rank_col} cutoff < {cut}.', '',
            'COLUMNS', '  ligand, receptor          interaction partners (canonical symbols where orthologous)',
            '  source_cell_type          sender (ligand-expressing)', '  target_cell_type          receiver (receptor-expressing)',
            '  significant_species       both / mouse / human',
            '  has_ortholog              TRUE = 1:1 ortholog in both species; FALSE = lineage-specific',
            '  min_magnitude_rank_*      strongest rank per species (lower = stronger). Populated even',
            '                            when not significant (shows how close it came); blank = untestable.',
            '  timepoints_*              timepoints where the pair was significant in that species', '',
            f'immune cell types: {", ".join(sorted(immune_labels))}'], start=1):
            ws.cell(row=r, column=1, value=txt)
        ws.column_dimensions['A'].width = 88
    return len(imm)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--liana-root', required=True, help='e.g. results/liana/no_commot')
    ap.add_argument('--cell-type-key', required=True)
    ap.add_argument('--mouse-timepoints', nargs='+', required=True)
    ap.add_argument('--human-timepoints', nargs='+', required=True)
    ap.add_argument('--orthology', required=True); ap.add_argument('--homology-type', default='ortholog_one2one')
    ap.add_argument('--cell-type-map', required=True, help='YAML with compartment_map, display_names, compartment_order, heatmap_drop')
    ap.add_argument('--rank-column', default='magnitude_rank'); ap.add_argument('--rank-cutoff', type=float, default=0.05)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()

    out = Path(a.outdir); TAB, FIG, DEL = out / 'tables', out / 'figures', out / 'deliverables'
    for d in (TAB, FIG, DEL): d.mkdir(parents=True, exist_ok=True)
    log = lambda *x: print(*x, flush=True)
    rc, cut = a.rank_column, a.rank_cutoff

    cfg = yaml.safe_load(open(a.cell_type_map))
    comp_map = cfg['compartment_map']; display = cfg['display_names']
    order = cfg.get('compartment_order', []); drop = set(cfg.get('heatmap_drop', []))
    immune_labels = {display[k] for k, v in comp_map.items() if v == 'immune' and k in display}

    m2h, human_set, ometa = load_one2one(Path(a.orthology), a.homology_type)
    log(f"orthology: {ometa['n_components']:,} one2one components -> {ometa['n_pairs']:,} symbol-unique pairs "
        f"({ometa['source_release']}); {ometa['n_dropped_duplicate_symbol']} dropped for duplicated symbols: "
        f"{ometa['dropped_duplicate_symbols']}; {ometa['n_components_not_1_per_species']} components not 1 gene/species")

    log('loading LIANA tables:')
    root = Path(a.liana_root)
    mouse = load_all_timepoints(root, 'mouse', a.cell_type_key, a.mouse_timepoints, log)
    human = load_all_timepoints(root, 'human', a.cell_type_key, a.human_timepoints, log)
    for df, sp in ((mouse, 'mouse'), (human, 'human')):
        add_canonical_keys(df, sp, m2h, human_set)
        n, t = int(df['lr_canonical'].notna().sum()), len(df)
        log(f'  {sp} rows with 1:1 orthologous L and R: {n:,} / {t:,} ({100*n/t:.1f}%)')
        um = unmappable_subunits(df, sp, m2h, human_set)
        um.rename('n_rows').rename_axis('subunit').to_csv(TAB / f'unmappable_subunits_{sp}.csv')
        log(f'  {sp} unmappable subunits (top 15): ' + ', '.join(f'{k}({v})' for k, v in um.head(15).items()))

    labels = set(mouse.source) | set(mouse.target) | set(human.source) | set(human.target)
    missing = labels - set(comp_map)
    if missing: sys.exit(f'labels missing from compartment_map: {sorted(missing)}')
    missing = labels - set(display)
    if missing: sys.exit(f'labels missing from display_names: {sorted(missing)}')

    mouse_t = mouse.dropna(subset=['lr_canonical']).copy(); human_t = human.dropna(subset=['lr_canonical']).copy()
    mouse_sig, human_sig = mouse_t[mouse_t[rc] < cut], human_t[human_t[rc] < cut]
    mouse_lr, human_lr = set(mouse_sig.lr_canonical), set(human_sig.lr_canonical)
    conserved, mouse_only, human_only = mouse_lr & human_lr, mouse_lr - human_lr, human_lr - mouse_lr
    hg = hypergeometric(mouse_t, human_t, mouse_lr, human_lr, conserved)
    log(f'significant orthologous LRs: mouse {len(mouse_lr)}, human {len(human_lr)}, conserved {len(conserved)}, '
        f'mouse-only {len(mouse_only)}, human-only {len(human_only)}; hypergeom fold {hg["fold"]:.2f}, p {hg["p_value"]:.2e}')

    m_trip, h_trip = sig_triples(mouse_t, rc, cut), sig_triples(human_t, rc, cut)
    matched = m_trip & h_trip
    matched_lrs = {t[0] for t in matched}; rewired_lrs = conserved - matched_lrs
    log(f'fine triples: mouse {len(m_trip):,} human {len(h_trip):,} matched {len(matched):,}; '
        f'matched-context LRs {len(matched_lrs)}, rewired {len(rewired_lrs)}')

    # tables
    per_cell_type_pair_table(mouse_sig, human_sig).to_csv(TAB / 'lr_conservation_by_pair.csv', index=False)
    long = pd.DataFrame([{'lr_canonical': lr, 'ligand': lr.partition('^')[0], 'receptor': lr.partition('^')[2],
                          'source_compartment': comp_map[s], 'target_compartment': comp_map[t],
                          'source_fine': s, 'target_fine': t} for lr, s, t in matched])
    long = long.sort_values(['source_compartment', 'target_compartment', 'lr_canonical']).reset_index(drop=True)
    long.to_csv(TAB / 'compartment_conserved_long.csv', index=False)
    expanded_map = {k: (k if v == 'immune' else v) for k, v in comp_map.items()}
    views = {'pooled': build_grids(m_trip, h_trip, matched, comp_map, order, drop, display),
             'expanded': build_grids(m_trip, h_trip, matched, expanded_map, order, drop, display)}
    for view, (lab, mats) in views.items():
        pd.DataFrame(mats['M_count'].astype(int), index=lab, columns=lab).to_csv(TAB / f'compartment_conserved_grid_{view}.csv')
    pd.DataFrame(sorted(conserved), columns=['lr_canonical']).assign(context=lambda d: np.where(d.lr_canonical.isin(matched_lrs), 'matched', 'rewired')).to_csv(TAB / 'conserved_lr_pairs.csv', index=False)

    # figures
    plot_venn(len(mouse_only), len(conserved), len(human_only), rc, cut, hg, FIG / 'fig1_conservation_venn.pdf')
    plot_matched_vs_rewired(matched_lrs, rewired_lrs, FIG / 'fig3_matched_vs_rewired_summary.pdf')
    save_heatmaps(FIG, views, rc, cut)

    # deliverables
    log('deliverable consistency checks:')
    hits, pooled, raw_calls, spec_col = build_deliverables(mouse, human, rc, cut, display, log)
    raw_calls.to_csv(DEL / 'raw_calls.csv', index=False); hits.to_csv(DEL / 'lr_hits.csv', index=False)
    pooled.to_csv(DEL / 'lr_pooled.csv', index=False)
    write_workbook(DEL / 'liana_table.xlsx', hits, pooled, rc, cut, ometa, a.cell_type_key)
    n_imm = write_immune(DEL / 'liana_immune.xlsx', DEL / 'lr_hits_immune.csv', hits, immune_labels, rc, cut)

    summary = {
        'rank_column': rc, 'rank_cutoff': cut, 'cell_type_key': a.cell_type_key, 'orthology': ometa,
        'timepoints': {'mouse': a.mouse_timepoints, 'human': a.human_timepoints},
        'rows': {'mouse_total': len(mouse), 'human_total': len(human), 'mouse_orthologous': len(mouse_t), 'human_orthologous': len(human_t)},
        'significant_orthologous_lr': {'mouse': len(mouse_lr), 'human': len(human_lr), 'conserved': len(conserved),
                                       'mouse_only': len(mouse_only), 'human_only': len(human_only),
                                       'pct_human_in_mouse': 100 * len(conserved) / len(human_lr) if human_lr else None,
                                       'pct_mouse_in_human': 100 * len(conserved) / len(mouse_lr) if mouse_lr else None,
                                       'jaccard_pct': 100 * len(conserved) / len(mouse_lr | human_lr) if (mouse_lr | human_lr) else None},
        'hypergeometric': hg,
        'context_decomposition': {'fine_triples_mouse': len(m_trip), 'fine_triples_human': len(h_trip), 'matched_triples': len(matched),
                                  'matched_context_lrs': len(matched_lrs), 'rewired_lrs': len(rewired_lrs)},
        'deliverables': {'lr_hits_rows': len(hits), 'conserved_contexts': int((hits.significant_species == 'both').sum()),
                         'unique_conserved_lr': int(hits.loc[hits.is_conserved_lr].drop_duplicates(['ligand', 'receptor']).shape[0]),
                         'non_ortholog_rows': int((~hits.has_ortholog).sum()), 'lr_pooled_rows': len(pooled),
                         'lr_pooled_conserved': int((pooled.conservation == 'conserved').sum()),
                         'raw_calls_rows': len(raw_calls), 'immune_rows': n_imm},
    }
    json.dump(summary, open(TAB / 'conservation_summary.json', 'w'), indent=1, default=float)
    log(json.dumps(summary['significant_orthologous_lr'], indent=1)); log(json.dumps(summary['deliverables'], indent=1))


if __name__ == '__main__':
    main()
