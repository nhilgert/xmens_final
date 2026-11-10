#!/usr/bin/env python
"""Supplementary Fig. S28: the two recurrence cutoffs behind the floor-6 SCENIC consensus regulons.

Port of regulon_recurrence_cutoffs_fig.py (Claude Science artifact 7bf91dd3, 2 Nov 2026). Panels:
  A  runs (of 100) in which each TF formed a regulon; bars at >= 80 runs (kept TFs) in species colour
  B  regulons kept vs TF threshold, per species and for TFs kept in both species (one-to-one orthologs)
  C  share of TF-target pairs, and of their occurrence weight, kept vs target threshold (log x)
  D  targets per regulon (median and IQR) vs target threshold (log-log); dotted = median single-run size

Inputs
  --report-mouse, --report-human  recurrence_report.tsv of the floor-6 consensus (per TF: n_present, n_targets_union,
                                  n_targets_kept, kept)
  --long-mouse, --long-human      consensus_regulons_long.tsv of the floor-2 consensus (s28_consensus_edges.py)
  --perrun                        perrun_sizes.json (s28_perrun_sizes.py)
  --orthologs                     Ensembl 116 component table; one-to-one pairs are homology_type ortholog_one2one
Outputs: <out>.svg (live text), <out>.png (300 dpi), <out>.provenance.json (sha256 of inputs and outputs)
Occurrence weight of a pair = runs containing the pair / runs in which its TF formed a regulon (the AUCell weight).
"""
import argparse
import hashlib
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import NullLocator  # noqa: E402

GATE, FLOOR, MIN_T = 80, 6, 10          # TF threshold (runs), target threshold (runs), minimum targets per regulon
N_ORTHO_PAIRS = 16275                    # one-to-one pairs with unique symbols, Ensembl 116 (as in Fig. 5E)
BASE, SMALL = 16, 14                     # label and legend/annotation sizes (pt); ticks use SMALL
COL = {"mouse": "#4d4d4d", "human": "#1b9e77"}
GREY, REF = "#c8c8c8", "0.55"
SPECIES = ("mouse", "human")
FONTS = ("Helvetica Neue", "Helvetica", "Arial", "Nimbus Sans", "DejaVu Sans")
F = np.arange(1, 101)
LOGT = [1, 2, 5, 10, 20, 50, 100]
TICKS_F = [1, 20, 40, 60, 80, 100]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def apply_figure_style(sizes=(BASE, SMALL, SMALL)):
    """Type ladder and axis style of the manuscript figures (local stand-in for the figure-style helper)."""
    base, legend, tick = sizes
    mpl.rcParams.update({
        "font.size": base, "axes.labelsize": base, "axes.titlesize": base, "legend.fontsize": legend,
        "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
        "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "xtick.major.size": 4, "ytick.major.size": 4,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })


def panel_letter(ax, letter, dx=-0.15, dy=1.07, fontsize=BASE + 4):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=fontsize, fontweight="bold", ha="left", va="bottom")


def ortholog_pairs(path):
    comp = pd.read_csv(path, sep="\t")
    one = comp[comp["homology_type"].eq("ortholog_one2one")]
    h = one[one["species"].eq("human")].groupby("component_id")["gene_symbol"].first()
    m = one[one["species"].eq("mouse")].groupby("component_id")["gene_symbol"].first()
    p = pd.DataFrame({"human": h, "mouse": m}).dropna().drop_duplicates()
    # unique gene symbols on both sides: SPACDR and TSC22D4 both map to mouse Tsc22d4, so both pairs are dropped
    p = p[~p.human.duplicated(keep=False) & ~p.mouse.duplicated(keep=False)]
    return p.reset_index(drop=True)


def load_species(sp, report, long_path, perrun):
    r = pd.read_csv(report, sep="\t").set_index("name")
    assert r.index.str.endswith("(+)").all()
    k = r[r.kept]
    assert len(k) == int(((r.n_present >= GATE) & (r.n_targets_kept >= MIN_T)).sum()) == int((r.n_present >= GATE).sum())
    d = pd.read_csv(long_path, sep="\t")
    assert set(d.regulon) == set(k.index)
    single = k.n_targets_union - d.groupby("regulon").size().reindex(k.index)
    assert (single >= 0).all()
    assert (d[d.n_seeds >= FLOOR].groupby("regulon").size().reindex(k.index) == k.n_targets_kept).all()
    ns = d.n_seeds.to_numpy()
    wv = (d.n_seeds / d.regulon.map(k.n_present)).to_numpy()
    w1 = float((single / k.n_present).sum())
    tot_n, tot_w = len(d) + int(single.sum()), wv.sum() + w1
    return dict(
        r=r, k=k,
        pairs=np.array([((ns >= f).sum() + (single.sum() if f == 1 else 0)) / tot_n for f in F]),
        weight=np.array([(wv[ns >= f].sum() + (w1 if f == 1 else 0)) / tot_w for f in F]),
        sizes=np.vstack([k.n_targets_union.to_numpy() if f == 1 else
                         d[d.n_seeds >= f].groupby("regulon").size().reindex(k.index).fillna(0).to_numpy()
                         for f in F]),
        perrun=float(np.median([perrun[sp]["per_tf_median_size"][n] for n in k.index])),
        gate=np.array([int(((r.n_present >= g) & (r.n_targets_kept >= MIN_T)).sum()) for g in F]),
    )


def passing(r, g):
    return {n[:-3] for n in r.index[(r.n_present >= g) & (r.n_targets_kept >= MIN_T)]}


def draw(D, PAIRED):
    apply_figure_style(sizes=(BASE, SMALL, SMALL))
    have = {f.name for f in font_manager.fontManager.ttflist}
    mpl.rcParams.update({"font.family": [f for f in FONTS if f in have] or ["DejaVu Sans"]})
    plt.close("all")
    fig = plt.figure(figsize=(11.0, 8.4), dpi=300)
    outer = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.30, left=0.085, right=0.975, top=0.925, bottom=0.085)
    ga = outer[0, 0].subgridspec(2, 1, hspace=0.22)
    axa1 = fig.add_subplot(ga[0])
    axa2 = fig.add_subplot(ga[1], sharex=axa1, sharey=axa1)
    axb, axc, axd = fig.add_subplot(outer[0, 1]), fig.add_subplot(outer[1, 0]), fig.add_subplot(outer[1, 1])

    # A: TF recurrence
    ymax = max(np.bincount(D[sp]["r"].n_present.to_numpy(), minlength=101)[1:].max() for sp in SPECIES) * 1.32
    for ax, sp in ((axa1, "mouse"), (axa2, "human")):
        cnt = np.bincount(D[sp]["r"].n_present.to_numpy(), minlength=101)[1:]
        ax.bar(F, cnt, width=1.0, color=np.where(F >= GATE, COL[sp], GREY), linewidth=0)
        ax.axvline(GATE - 0.5, color=REF, lw=0.9)
        ax.set_ylim(0, ymax)
        ax.set_xlim(0, 101.5)
        ax.text(0.015, 0.95, f"{sp.capitalize()}: {int(cnt.sum())} TFs with a regulon\nin at least one run",
                transform=ax.transAxes, ha="left", va="top", fontsize=SMALL, color=COL[sp])
        ax.text(GATE - 2, ymax * 0.40, f"{len(D[sp]['k'])} kept", ha="right", va="center", fontsize=SMALL,
                color=COL[sp])
        ax.set_ylabel("TFs")
    plt.setp(axa1.get_xticklabels(), visible=False)
    axa2.set_xticks(TICKS_F)
    axa2.set_xlabel("Runs (of 100) in which the TF formed a regulon")

    # B: regulons kept vs TF threshold
    for sp in SPECIES:
        axb.plot(F, D[sp]["gate"], color=COL[sp], lw=1.8, label=sp.capitalize())
    axb.plot(F, PAIRED, color="black", lw=1.4, ls=(0, (4, 2)), label="In both species")
    axb.axvline(GATE, color=REF, lw=0.9)
    for v, c_, dx, dy, ha, va in (
            (D["mouse"]["gate"][GATE - 1], COL["mouse"], 1.5, 6, "left", "bottom"),
            (D["human"]["gate"][GATE - 1], COL["human"], -1.5, 17, "right", "bottom"),
            (PAIRED[GATE - 1], "black", -1.5, -6, "right", "top")):
        axb.plot([GATE], [v], "o", ms=4.5, color=c_, zorder=4)
        axb.text(GATE + dx, v + dy, f"{v}", ha=ha, va=va, fontsize=SMALL, color=c_)
    axb.set_xlim(0, 101.5)
    axb.set_ylim(0, max(D["mouse"]["gate"].max(), D["human"]["gate"].max()) * 1.12)
    axb.set_xticks(TICKS_F)
    axb.set_xlabel("TF threshold (minimum runs, of 100)")
    axb.set_ylabel("Regulons kept")
    axb.legend(loc="lower left", frameon=False, fontsize=SMALL, handlelength=2.2)

    # C: pairs and occurrence weight kept vs target threshold
    for sp in SPECIES:
        axc.plot(F, 100 * D[sp]["pairs"], color=COL[sp], lw=1.8)
        axc.plot(F, 100 * D[sp]["weight"], color=COL[sp], lw=1.8, ls=(0, (4, 2)))
        for kind in ("pairs", "weight"):
            y = 100 * D[sp][kind][FLOOR - 1]
            fx_, dy, ha, va = (1.1, 2.5, "left", "bottom") if sp == "mouse" else (1 / 1.1, -2.0, "right", "top")
            axc.plot([FLOOR], [y], "o", ms=4.5, color=COL[sp], zorder=4)
            axc.text(FLOOR * fx_, y + dy, f"{y:.0f} %", ha=ha, va=va, fontsize=SMALL, color=COL[sp])
    axc.axvline(FLOOR, color=REF, lw=0.9)
    axc.set_xscale("log")
    axc.set_xticks(LOGT)
    axc.set_xticklabels([str(t) for t in LOGT])
    axc.xaxis.set_minor_locator(NullLocator())
    axc.set_xlim(0.92, 108)
    axc.set_ylim(0, 104)
    axc.set_xlabel("Target threshold (minimum runs)")
    axc.set_ylabel("Kept (%)")
    axc.text(0.98, 0.97, "100 % = all pairs\nseen in any run", transform=axc.transAxes, ha="right", va="top",
             fontsize=SMALL, color="0.35")
    axc.legend(handles=[Line2D([], [], color="0.3", lw=1.8, label="TF\u2013target pairs"),
                        Line2D([], [], color="0.3", lw=1.8, ls=(0, (4, 2)), label="Occurrence weight")],
               loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=SMALL,
               handlelength=2.2, columnspacing=1.6, borderaxespad=0.2)

    # D: regulon size vs target threshold
    for sp in SPECIES:
        s = D[sp]["sizes"]
        med = np.median(s, axis=1)
        q1, q3 = np.percentile(s, [25, 75], axis=1)
        ok = med > 0
        axd.fill_between(F[ok], np.clip(q1[ok], 1, None), q3[ok], color=COL[sp], alpha=0.16, lw=0)
        axd.plot(F[ok], med[ok], color=COL[sp], lw=1.8)
        axd.axhline(D[sp]["perrun"], color=COL[sp], lw=1.1, ls=":")
        m6 = med[FLOOR - 1]
        axd.plot([FLOOR], [m6], "o", ms=4.5, color=COL[sp], zorder=4)
        if sp == "mouse":
            axd.text(FLOOR * 1.12, m6 * 1.05, f"{m6:.0f}", ha="left", va="bottom", fontsize=SMALL, color=COL[sp])
        else:
            axd.text(FLOOR / 1.12, m6 / 1.05, f"{m6:.0f}", ha="right", va="top", fontsize=SMALL, color=COL[sp])
    axd.text(1.0, D["mouse"]["perrun"] * 1.08, f"one run: {D['mouse']['perrun']:.0f}", ha="left", va="bottom",
             fontsize=SMALL, color=COL["mouse"])
    axd.text(1.0, D["human"]["perrun"] / 1.08, f"one run: {D['human']['perrun']:.0f}", ha="left", va="top",
             fontsize=SMALL, color=COL["human"])
    axd.axvline(FLOOR, color=REF, lw=0.9)
    axd.set_xscale("log")
    axd.set_yscale("log")
    axd.set_xticks(LOGT)
    axd.set_xticklabels([str(t) for t in LOGT])
    axd.xaxis.set_minor_locator(NullLocator())
    axd.set_xlim(0.92, 108)
    axd.set_ylim(1, 6000)
    axd.set_yticks([1, 10, 100, 1000])
    axd.set_yticklabels(["1", "10", "100", "1k"])
    axd.yaxis.set_minor_locator(NullLocator())
    axd.set_xlabel("Target threshold (minimum runs)")
    axd.set_ylabel("Targets per regulon")

    for ax in fig.axes:
        ax.xaxis.label.set_fontweight("bold")
        ax.yaxis.label.set_fontweight("bold")
    for ax, lt in ((axa1, "A"), (axb, "B"), (axc, "C"), (axd, "D")):
        panel_letter(ax, lt, dx=-0.15, dy=1.07, fontsize=BASE + 4)
    return fig


def audit(fig):
    """Render-then-verify (figure-style 9.1): no overlapping text, no text outside the figure."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    offview = set()
    for a in fig.axes:
        for axis, lim in ((a.xaxis, a.get_xlim()), (a.yaxis, a.get_ylim())):
            lo, hi = sorted(lim)
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                loc = tk.get_loc()
                if loc is None or not (lo - 1e-9 <= loc <= hi + 1e-9):
                    offview.update({tk.label1, tk.label2})
    texts = [(t, t.get_window_extent(rend)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible() and t not in offview]
    spines = [(s, s.get_window_extent(rend)) for a in fig.axes for s in a.spines.values() if s.get_visible()]
    ticklabels = {a: set(a.get_xticklabels(which="both") + a.get_yticklabels(which="both")) for a in fig.axes}
    overlaps = [(a.get_text()[:30], b.get_text()[:30]) for i, (a, ba) in enumerate(texts)
                for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    overlaps += [(t.get_text()[:30], "spine") for t, bt in texts for s, bs in spines
                 if bt.overlaps(bs) and t not in ticklabels[s.axes]]
    fb = fig.bbox
    outside = [t.get_text()[:30] for t, b in texts if b.x0 < fb.x0 or b.x1 > fb.x1 or b.y0 < fb.y0 or b.y1 > fb.y1]
    return overlaps, outside


def main():
    ap = argparse.ArgumentParser(description="Fig. S28: SCENIC recurrence cutoffs")
    ap.add_argument("--report-mouse", required=True)
    ap.add_argument("--report-human", required=True)
    ap.add_argument("--long-mouse", required=True)
    ap.add_argument("--long-human", required=True)
    ap.add_argument("--perrun", required=True)
    ap.add_argument("--orthologs", required=True)
    ap.add_argument("--out", required=True, help="output prefix, without extension")
    a = ap.parse_args()

    with open(a.perrun) as fh:
        perrun = json.load(fh)
    O = ortholog_pairs(a.orthologs)
    assert len(O) == N_ORTHO_PAIRS, f"{len(O)} one-to-one pairs, expected {N_ORTHO_PAIRS}"
    D = {
        "mouse": load_species("mouse", a.report_mouse, a.long_mouse, perrun),
        "human": load_species("human", a.report_human, a.long_human, perrun),
    }
    PAIRED = np.array([int((O.human.isin(passing(D["human"]["r"], g)) &
                            O.mouse.isin(passing(D["mouse"]["r"], g))).sum()) for g in F])
    assert PAIRED[GATE - 1] == 111, PAIRED[GATE - 1]

    fig = draw(D, PAIRED)
    png, svg = a.out + ".png", a.out + ".svg"
    fig.savefig(png, metadata={"Software": None})
    with mpl.rc_context({"svg.hashsalt": "figS28"}):
        fig.savefig(svg, metadata={"Date": None})
    overlaps, outside = audit(fig)
    assert not overlaps and not outside, (overlaps, outside)

    inputs = {"report_mouse": a.report_mouse, "report_human": a.report_human, "long_mouse": a.long_mouse,
              "long_human": a.long_human, "perrun": a.perrun, "orthologs": a.orthologs}
    prov = {
        "figure": "Fig. S28 (supplementary): SCENIC recurrence cutoffs",
        "script_sha256": sha256(__file__),
        "inputs_sha256": {k: sha256(v) for k, v in inputs.items()},
        "outputs_sha256": {"svg": sha256(svg), "png": sha256(png)},
        "one_to_one_pairs": len(O),
        "kept_in_both_at_gate": int(PAIRED[GATE - 1]),
        "kept_at_gate": {sp: int(D[sp]["gate"][GATE - 1]) for sp in SPECIES},
        "matplotlib": matplotlib.__version__,
    }
    with open(a.out + ".provenance.json", "w") as fh:
        json.dump(prov, fh, indent=1)
    print("kept at gate", prov["kept_at_gate"], "in both", prov["kept_in_both_at_gate"], flush=True)
    print("wrote", png, svg, a.out + ".provenance.json", flush=True)


if __name__ == "__main__":
    main()
