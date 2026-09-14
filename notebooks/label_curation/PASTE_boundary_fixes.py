# ====================================================================================================
# §2b — marker selection (the cell starting `# A marker used as a THRESHOLD`)
# ====================================================================================================
# per-cluster detection fraction, ordered by the best decidual marker
# A marker used as a THRESHOLD needs two properties, and AUROC only tests one of them: it must rank
# cells well AND have low background in the negative class. FOXO1 has the best decidual AUROC (0.716)
# but is detected in ~65% of fibroblasts, so averaging it into a detection score shifts every cluster
# above any threshold calibrated on a low-background gene -- which is what emptied the fibroblast body
# in the first version of this notebook. Require both.
MAXBG = 0.35
cand = auroc[auroc.group.str.startswith("decidual")]
DECGENES = [g for g in cand.index if cand.auroc_dec_vs_fib[g] > 0.6 and cand.det_in_fib[g] < MAXBG]
print(f"decidual markers usable as thresholds (AUROC > 0.6 AND background in fibroblasts < {MAXBG}):")
for g in cand.index:
    keep = g in DECGENES
    why = "" if keep else ("AUROC too low" if cand.auroc_dec_vs_fib[g] <= 0.6 else f"background {cand.det_in_fib[g]:.2f} too high")
    print(f"   {'KEEP' if keep else 'drop'}  {g:8s} AUROC={cand.auroc_dec_vs_fib[g]:.3f} det_in_fib={cand.det_in_fib[g]:.3f}  {why}")
verify("at least one usable threshold marker", len(DECGENES) > 0, f"{DECGENES}")
BANDGENE = DECGENES[0]      # the band is cut on ONE gene so its thresholds keep their meaning
print("\nband will be cut on:", BANDGENE)

# two observations from the same table, worth having in front of you:
# The ECM axis (class B, section 4) must use fibroblast-panel genes that are actually fibroblast-enriched
# HERE. On this dataset most of the panel is not: AUROC-for-decidual > 0.5 means the gene is higher in
# DECIDUAL cells. Select rather than assume.
fibcand = auroc[auroc.group == "fibroblast"]
ECMGENES = [g for g in fibcand.index if fibcand.auroc_dec_vs_fib[g] <= 0.60]
inverted = [g for g in fibcand.index if g not in ECMGENES]
print(f"\nfibroblast-panel genes, direction on THIS dataset:")
for g in fibcand.index:
    print(f"   {'ECM axis' if g in ECMGENES else 'inverted'}  {g:8s} AUROC_dec={fibcand.auroc_dec_vs_fib[g]:.3f}  "
          f"det_in_dec={fibcand.det_in_dec[g]:.3f} det_in_fib={fibcand.det_in_fib[g]:.3f}")
if inverted: print(f"   -> {', '.join(inverted)} are detected MORE in decidual cells; excluded from the ECM axis.")
verify("at least one fibroblast-enriched ECM gene", len(ECMGENES) > 0, f"{ECMGENES}")
mp = auroc[auroc.group == "decidual_mousepanel"]
if len(mp): print(f"NOTE  mouse-panel orthologs do not transfer: AUROC {mp.auroc_dec_vs_fib.min():.3f}-{mp.auroc_dec_vs_fib.max():.3f} "
                  f"({', '.join(mp.index)}) -- the mouse decidual markers are not usable in human.")

cl = o[RES].astype(str)
frac = (E > 0).groupby(cl).mean(); mean = E.groupby(cl).mean()
ORD  = frac[BANDGENE].sort_values(ascending=False).index.tolist()

tab = pd.concat([cl.value_counts().rename("n"),
                 V[["label","purity"]].rename(columns={"label":"vote"}),
                 frac[[g for g in DECGENES]].add_prefix("det_").round(3),
                 mean[ECMGENES].mean(axis=1).rename("ecm_mean").round(3),
                 mean[PANEL["prolif"]].mean(axis=1).rename("mean_prolif").round(2),
                 mean[PANEL["contam"]].mean(axis=1).rename("mean_contam").round(2)], axis=1).reindex(ORD)
tab.index.name = "cluster"
print(f"\n=== per-cluster, ordered by {BANDGENE} detection ===")
print(tab.to_string())
tab.to_csv(OUT/"cluster_marker_table.csv")

fig, ax = plt.subplots(figsize=(11, 9))
sns.heatmap(frac.reindex(ORD)[DECGENES + PANEL["fibroblast"] + PANEL["prolif"]], cmap="viridis", vmin=0, vmax=1,
            cbar_kws=dict(label="fraction of cells with non-zero counts"), ax=ax)
for i, c_ in enumerate(ORD):
    ax.text(-0.4, i + 0.5, ("D" if V.label.get(c_) == "decidual_cell" else "F" if V.label.get(c_) == "endometrial_fibroblast" else "?"),
            ha="right", va="center", fontsize=8, color="0.3")
ax.set_title("human stromal clusters x panel markers (row label = committed vote)", fontsize=11)
ax.set_ylabel(f"{RES} cluster"); fig.tight_layout(); fig.savefig(OUT/"cluster_marker_heatmap.png", dpi=130, bbox_inches="tight")


# ====================================================================================================
# §4 — band + class B (the cell starting `# Class A = clusters`)
# ====================================================================================================
# Class A = clusters whose BANDGENE detection falls between the two bodies. Cut on one gene, so the
# thresholds mean what they say (LEFTY2 detection runs 0.03-0.70 across clusters).
det = frac[BANDGENE].reindex(tab.index)
LO, HI = 0.20, 0.35     # <- edit and re-run this cell to move the band
print(f"{BANDGENE} detection across clusters: {det.min():.3f} - {det.max():.3f}   band = [{LO}, {HI}]")
verify("thresholds are on the same scale as the score",
       (det < LO).any() and (det > HI).any(),
       f"{int((det<LO).sum())} clusters below LO, {int((det>HI).sum())} above HI")
band = det[(det >= LO) & (det <= HI)].index.tolist()
body_dec = det[det > HI].index.tolist(); body_fib = det[det < LO].index.tolist()
print(f"band (det {LO}-{HI}): {band}  -> {int(tab.n[band].sum()):,} cells")
print(f"decidual body: {body_dec}\nfibroblast body: {body_fib}")
print("\nvotes within the band:", V.label[band].value_counts().to_dict())

# --- class B: DOUBLE-POSITIVE clusters -- decidual marker present AND fibroblast/ECM markers retained.
# The audit found these are 82% of the low-purity votes (2, 5, 8, 11, 24) and a single-axis decidual
# threshold misses them entirely: their decidual detection sits in the dec body, so only the fibroblast
# axis distinguishes them. Defined by percentile so it is not tied to one dataset's absolute levels.
# ECM axis = mean expression of the fibroblast-enriched panel genes selected in section 2 (ECMGENES;
# the inverted ones are excluded there). Percentile-based, so not tied to absolute levels.
print("ECM axis genes:", ECMGENES)
ecmp = tab.ecm_mean.rank(pct=True)
decp = det.rank(pct=True)
DP_ECM_PCT, DP_DEC_PCT = 0.70, 0.55     # <- edit: "top 30% ECM AND upper 45% decidual"
dpos = tab.index[(ecmp >= DP_ECM_PCT) & (decp >= DP_DEC_PCT)].tolist()
print(f"\nclass B, double-positive (ECM pct >= {DP_ECM_PCT}, {BANDGENE} pct >= {DP_DEC_PCT}): {dpos}"
      f"  -> {int(tab.n[dpos].sum()):,} cells")
print("   ECM mean / dec det for these:",
      {k: (round(tab.ecm_mean[k],2), round(det[k],2)) for k in dpos})
# Being double-positive and being a low-purity vote are different things: report the intersection
# rather than tuning the cut until they coincide.
if V.purity.notna().any():
    _lp = set(V.index[V.purity < 0.8])
    print(f"   of these, the reference DISAGREED on: {sorted(_lp & set(dpos), key=lambda s: int(s))}"
          f" | agreed on: {sorted(set(dpos) - _lp, key=lambda s: int(s))}")
print("   votes:", V.label[dpos].value_counts().to_dict(),
      "| purity:", V.purity[dpos].round(2).to_dict())
print(f"A n B overlap: {sorted(set(band) & set(dpos))}")
if V.purity.notna().any():
    lowp = sorted(V.index[V.purity < 0.8], key=lambda s: int(s))
    n_low = int(tab.n.reindex(lowp).sum())
    print(f"\nlow-purity clusters {lowp} = {n_low:,} cells")
    for nm, S_ in (("A band", band), ("B double-pos", dpos)):
        ov = sorted(set(lowp) & set(S_), key=lambda s: int(s))
        print(f"   {nm:14s} explains {int(tab.n.reindex(ov).sum()):,} of them ({100*tab.n.reindex(ov).sum()/max(n_low,1):.0f}%)  {ov}")
    verify("the two classes together explain most low-purity votes",
           (f := tab.n.reindex(sorted(set(lowp) & (set(band) | set(dpos)))).sum() / max(n_low, 1)) > 0.8,
           f"{f:.2f} of low-purity cells")

o["band_cluster"] = cl.isin(band); o["doublepos_cluster"] = cl.isin(dpos)
grp = pd.Series(np.where(cl.isin(dpos), "double_positive",
                np.where(cl.isin(band), "band",
                np.where(cl.isin(body_dec), "dec_body", "fib_body"))), index=o.index)
o["group3"] = grp
print("\ngroups:", grp.value_counts().to_dict())
for g_ in ("fib_body", "dec_body", "band", "double_positive"):
    verify(f"group {g_} is non-empty", int((grp == g_).sum()) > 0, f"{int((grp == g_).sum())} cells")

rep = pd.crosstab(grp, o["sample"].astype(str))
print("\n=== replication: band cells per sample ===")
print(rep.to_string())
print("\nband as a fraction of each sample:")
print((rep.loc["band"] / rep.sum()).round(3).to_string())
verify("band replicates across samples", (rep.loc["band"] > 0).all(),
       f"present in {int((rep.loc['band']>0).sum())}/{rep.shape[1]} samples")
verify("band is not one sample's artifact", (rep.loc["band"] / rep.loc["band"].sum()).max() < 0.75,
       f"largest single-sample share {(rep.loc['band']/rep.loc['band'].sum()).max():.2f}")

axes_tab = o.groupby(grp).agg(n=("dec_score","size"), dec_score=("dec_score","median"),
                              isodepth=("isodepth_consensus","median"), pt=("dpt_pseudotime","median"),
                              frac_menstrual=("phase", lambda s: (s.astype(str)=="menstrual").mean())).round(3)
print("\n=== axis position by group ==="); print(axes_tab.reindex(["fib_body","band","double_positive","dec_body"]).to_string())
for col in ("isodepth_consensus", "dpt_pseudotime"):
    m = o.groupby(grp)[col].median()
    if not {"fib_body","dec_body"} <= set(m.index): continue
    for nm in ("band", "double_positive"):
        if nm not in m.index: continue
        mid = (m["fib_body"] < m[nm] < m["dec_body"]) or (m["fib_body"] > m[nm] > m["dec_body"])
        verify(f"{nm} is intermediate on {col}", mid,
               f"fib={m['fib_body']:.3f} {nm}={m[nm]:.3f} dec={m['dec_body']:.3f}")


# ====================================================================================================
# §6b — violins (replace the whole function cell)
# ====================================================================================================
def violins(by, cols=("dec_score","isodepth_consensus","dpt_pseudotime"), order=None, cmap="viridis", figsize=None):
    """Fill colour = that category's median of the column plotted, on that column's own scale
    (so a violin's colour is the colour the tissue takes at that category's typical position).
    Black bar = median; outline = distribution."""
    d = o[o[by].notna()].copy(); d[by] = d[by].astype(str)
    present = set(d[by])
    want = order or CALLORDER
    dropped = [s for s in want if s not in present]
    order = [s for s in want if s in present]
    if dropped: print(f"  (no cells in: {dropped} -- omitted)")
    assert order, f"no categories of {by} have cells"
    cols = [c for c in cols if c in d and pd.to_numeric(d[c], errors="coerce").notna().any()]
    fig, axes = plt.subplots(1, len(cols), figsize=figsize or (5.8*len(cols), 4.6))
    axes = np.atleast_1d(axes); meds = {}
    for ax, col in zip(axes, cols):
        d[col] = pd.to_numeric(d[col], errors="coerce")
        med = d.groupby(by, observed=True)[col].median().reindex(order)
        meds[col] = med
        nrm = Normalize(vmin=float(med.min()), vmax=float(med.max()))
        # a palette entry for EVERY category in `order` -- seaborn raises on a missing key, and an
        # all-NaN median for one group must not take the whole figure down
        pal = {k: (plt.get_cmap(cmap)(nrm(med[k])) if pd.notna(med.get(k)) else (0.82, 0.82, 0.82)) for k in order}
        sns.violinplot(data=d, x=by, y=col, order=order, hue=by, hue_order=order, palette=pal,
                       legend=False, cut=0, linewidth=.6, ax=ax)
        for k, cat in enumerate(order):
            if pd.notna(med.get(cat)): ax.plot([k-.35, k+.35], [med[cat]]*2, c="k", lw=1.8, zorder=5)
        ax.set_title(col, fontsize=11); ax.set_xlabel(""); ax.tick_params(axis="x", rotation=25, labelsize=8)
        fig.colorbar(ScalarMappable(norm=nrm, cmap=cmap), ax=ax, shrink=.85, label=f"median {col}")
    fig.suptitle(f"by {by}  (n={len(d):,})  —  fill = category median", y=1.03)
    fig.tight_layout(); return fig, pd.DataFrame(meds).round(4)

f1, M1 = violins("boundary_call"); f1.savefig(OUT/"violins_boundary_call.png", dpi=130, bbox_inches="tight")
print(M1.to_string())
f2, _ = violins("group3", order=["fib_body","band","double_positive","dec_body"]); f2.savefig(OUT/"violins_group3.png", dpi=130, bbox_inches="tight")


