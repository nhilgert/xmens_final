# Ulrich et al. 2024 stromal subtypes - definitions and high-information markers

Definitions are from Ulrich et al. 2024 (PNAS 121:e2404775121), main text. Markers and all
specificity numbers are computed from their own Dataset S5B centroid matrix (43,090 genes x 39
subtypes) and cross-checked against the Dataset S6C DE sheet.

## How 'high information' is defined here

Ranking by their published `avg_log2FC` answers 'is this gene up in this subtype', which is not the
same question as 'does this gene tell me which subtype I am looking at'. What is computed below is a
**specificity margin**: a gene's centroid expression in the subtype minus its highest expression in
any *other* stromal subtype, in log1p CP10K units. A gene with margin +1.5 is roughly a doubling
above the best competing subtype; a gene with margin +0.1 is shared.

Two margins are reported:

- `margin_stromal` - against the other 10 stromal subtypes. The relevant number when you already
  know a cell is stromal and want to know which subtype.
- `margin_all39` - against all 38 other atlas subtypes including epithelium, endothelium and immune.
  A gene that clears the first but fails the second is not a usable marker in whole tissue.

Genes are required to be expressed (centroid > 0.25) in the subtype they mark.

## The one-page summary

| subtype | Ulrich label | kind | markers margin > 0.5 | of those, specific vs all 39 | verdict |
|---|---|---|---|---|---|
| Stro-1 | Perivascular, MYH11+ FOXC1-high CTGF-low | cell type | 31 | 11 | USABLE WITH CARE |
| Stro-2 | Perivascular, MYH11+ FOXC1-low CTGF-high | cell type | 5 | 3 | WEAK |
| Stro-3 | Perivascular, STEAP4+ (pericyte-like) | cell type | 8 | 6 | USABLE |
| Stro-4 | Tissue-resident vascular progenitor | cell type | 38 | 0 | NOT USABLE AS A STROMAL MARKER SET |
| Stro-5 | RUNX3+ Fibroblast | cell type | 66 | 0 | NOT USABLE IN snRNA-seq |
| Stro-6 | Uterine smooth muscle cell (uSMC) | cell type | 12 | 10 | USABLE |
| Stro-7 | C7 Fibroblast | cell type | 51 | 36 | BEST-DEFINED FIBROBLAST IN THE ATLAS |
| Stro-8 | dS1 - early decidualised stroma | STAGE | 6 | 3 | USABLE ONLY FROM CENTROIDS, NOT FROM THE PUBLISHED DE LIST |
| Stro-9 | dS2 - mid decidualised stroma | STAGE | 0 | 0 | NOT MARKER-DEFINABLE |
| Stro-10 | dS3 - late / senescent decidualised stroma | STAGE | 0 | 0 | NOT MARKER-DEFINABLE |
| Stro-11 | Stromal progenitor | STATE | 157 | 14 | MOST MARKER-DEFINABLE, BUT IT IS A STATE |

Read that table before anything else. **Stro-9 and Stro-10 have no specific markers at all** - not
one gene in the atlas's own centroid matrix clears a 0.5 margin against the other stromal subtypes.
That is not a defect in the data; it follows from what they are. Ulrich define dS1/dS2/dS3 as three
*stages* of one decidualisation process and separate them by transcription-factor activity
(STAT3 falling, HOXA10/NR2F2 rising across the series), not by marker genes. Any method that tries
to identify dS2 or dS3 from a marker list is attempting something the reference cannot support.

## The trajectory, in the paper's own terms

```
Stro-11 (progenitor, proliferative phase)
    +-- dS1 -> dS2 -> dS3          (decidualisation stages, Stro-8/9/10)
            +-- Stro-7  C7 fibroblast
            +-- Stro-6  uSMC
                    +-- Stro-4  vascular progenitor -> blood endothelium
                    +-- Stro-1  PV --+-- Stro-2  PV terminal
                                     +-- Stro-3  PV STEAP4 terminal
```

## Stro-1 - Perivascular, MYH11+ FOXC1-high CTGF-low

**What it is.** Perivascular (PV) cell. One of three PV states (Stro-1/2/3). Sits upstream of Stro-2 and Stro-3 in the paper's trajectory, which places it on the PV branch downstream of uSMC (Stro-6).

**Cycle phase.** Peaks in early secretory phase.

**Transcription factors (paper).** PV-shared: EBF1, RARB, HES1, SATB2, NR3C1, FOXK1, BATF3. Stro-1-specific: BHLHE40, PRRX2, RFX2, FOXC1, FOXC2, ZBTB7A, TCF7L1.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**HSPA1A** (+0.95), **CH25H** (+0.92), **C11orf96** (+0.90), **MT1A** (+0.81), **NR4A1** (+0.80), **HSPA1B** (+0.79), **KLF2** (+0.78), **CEBPD** (+0.75), **BTG2** (+0.75), **CCN3** (+0.69)

**Of those, also specific against all 39 subtypes:** **CH25H** (+0.92), **C11orf96** (+0.90), **MT1A** (+0.81), **NR4A1** (+0.80), **CEBPD** (+0.75), **CCN3** (+0.69), **NFIL3** (+0.69), **ECRG4** (+0.65)

**With a strict 1:1 mouse ortholog:** **CH25H** (+0.92), **C11orf96** (+0.90), **NR4A1** (+0.80), **KLF2** (+0.78), **CEBPD** (+0.75), **BTG2** (+0.75), **CCN3** (+0.69), **NFIL3** (+0.69)

**Verdict - USABLE WITH CARE.** 31 genes clear the margin but the top of the list is heat-shock and immediate-early (HSPA1A/B, NR4A1, BTG2, KLF2, CEBPD). Stro-1 and Stro-2 are each other's nearest neighbour, so much of what separates them is activation/stress state rather than identity.

_31 genes clear margin > 0.5; 11 of those also clear it against all 39; best margin +0.95._

---

## Stro-2 - Perivascular, MYH11+ FOXC1-low CTGF-high

**What it is.** Perivascular cell, terminal state derived from Stro-1. Distinguished from Stro-1 by the FOXC1/CTGF axis in its name.

**Cycle phase.** Peaks in early secretory phase.

**Transcription factors (paper).** TBX2, HEYL, HEY2 (Notch effectors) on top of the shared PV set.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**GPX3** (+0.85), **PHLDA2** (+0.69), **NET1** (+0.60), **LBH** (+0.56), **C10orf10** (+0.54), **MAP3K7CL** (+0.48), **RERGL** (+0.45), **SERPINI1** (+0.43), **C2orf40** (+0.41), **PRKCDBP** (+0.41)

**Of those, also specific against all 39 subtypes:** **PHLDA2** (+0.69), **NET1** (+0.60), **LBH** (+0.56)

**With a strict 1:1 mouse ortholog:** **GPX3** (+0.85), **PHLDA2** (+0.69), **NET1** (+0.60), **LBH** (+0.56), **RERGL** (+0.45), **SERPINI1** (+0.43), **PLN** (+0.40), **TMEM38B** (+0.40)

**Verdict - WEAK.** Only 5 genes clear a 0.5 margin, best +0.85. Effectively a second state of the same PV population as Stro-1.

_5 genes clear margin > 0.5; 3 of those also clear it against all 39; best margin +0.85._

---

## Stro-3 - Perivascular, STEAP4+ (pericyte-like)

**What it is.** Perivascular cell, the second terminal state from Stro-1. Marker profile is pericyte: HIGD1B, GUCY1A2, basement-membrane collagens.

**Cycle phase.** Peaks in early secretory phase.

**Transcription factors (paper).** ZBTB41, ESRRA on top of the shared PV set.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**STEAP4** (+1.27), **GUCY1A2** (+0.94), **COL4A1** (+0.74), **HIGD1B** (+0.70), **PMEPA1** (+0.59), **FAM162B** (+0.57), **COL4A2** (+0.55), **CAMK2N1** (+0.54), **GGT5** (+0.48), **ABCC9** (+0.45)

**Of those, also specific against all 39 subtypes:** **STEAP4** (+1.27), **GUCY1A2** (+0.94), **HIGD1B** (+0.70), **PMEPA1** (+0.59), **FAM162B** (+0.57), **CAMK2N1** (+0.54)

**With a strict 1:1 mouse ortholog:** **STEAP4** (+1.27), **GUCY1A2** (+0.94), **COL4A1** (+0.74), **HIGD1B** (+0.70), **PMEPA1** (+0.59), **FAM162B** (+0.57), **COL4A2** (+0.55), **CAMK2N1** (+0.54)

**Verdict - USABLE.** 8 genes, 6 of them specific against all 39 subtypes. STEAP4 (+1.27) and the pericyte pair HIGD1B/GUCY1A2 are clean. Note: its DE sheet rows are split across ~93 auto-numbered labels, so pool them before use.

_8 genes clear margin > 0.5; 6 of those also clear it against all 39; best margin +1.27._

---

## Stro-4 - Tissue-resident vascular progenitor

**What it is.** Explicitly a TRANSITIONAL subtype in the paper - a bridge from stroma toward blood endothelium, described as a stromal-to-endothelial transition. It carries stromal and endothelial lineage TFs simultaneously.

**Cycle phase.** Peaks in early secretory phase.

**Transcription factors (paper).** MECOM, MEOX1 (stromal) together with FLI1, SOX18 (endothelial) - the dual identity is the point.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**ACKR1** (+1.79), **CLDN5** (+1.56), **ENPP2** (+1.33), **IFI27** (+1.26), **CCL14** (+1.15), **AQP1** (+1.13), **VWF** (+1.09), **LIFR** (+1.08), **EGFL7** (+1.03), **RNASE1** (+0.96)

**Of those, also specific against all 39 subtypes:** _none_

**With a strict 1:1 mouse ortholog:** **ACKR1** (+1.79), **CLDN5** (+1.56), **ENPP2** (+1.33), **AQP1** (+1.13), **VWF** (+1.09), **LIFR** (+1.08), **EGFL7** (+1.03), **RNASE1** (+0.96)

**Verdict - NOT USABLE AS A STROMAL MARKER SET.** 38 genes clear the stromal margin but ZERO survive against all 39 subtypes - every one is shared with endothelium (ACKR1, CLDN5, VWF, EGFL7, AQP1). Exactly what the paper's stromal-to-endothelial framing predicts.

_38 genes clear margin > 0.5; 0 of those also clear it against all 39; best margin +1.79._

---

## Stro-5 - RUNX3+ Fibroblast

**What it is.** Labelled a fibroblast by Ulrich, but the paper itself notes that the Imm-2 fibrocyte population most resembles Stro-5. Its centroid is dominated by NK/T and mast-cell genes.

**Cycle phase.** Not phase-restricted in the paper's cycle analysis.

**Transcription factors (paper).** RUNX3 - an NK/T-lineage factor.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**GNLY** (+2.67), **TPSB2** (+2.44), **TPSAB1** (+2.06), **CCL5** (+1.97), **NKG7** (+1.84), **GZMA** (+1.66), **CCL4** (+1.52), **CD69** (+1.52), **PTPRC** (+1.38), **HCST** (+1.35)

**Of those, also specific against all 39 subtypes:** _none_

**With a strict 1:1 mouse ortholog:** **CCL5** (+1.97), **NKG7** (+1.84), **GZMA** (+1.66), **CD69** (+1.52), **PTPRC** (+1.38), **HCST** (+1.35), **LTB** (+1.33), **CD52** (+1.30)

**Verdict - NOT USABLE IN snRNA-seq.** 66 genes clear the margin, but they are NK/T and mast genes - GNLY, TPSB2, TPSAB1, CCL5, NKG7, GZMA, plus PTPRC (pan-leukocyte). Zero survive against all 39. Your data has ~1% immune-panel-positive nuclei, so this subtype cannot map, and that is a property of the reference, not a limitation of your assay.

_66 genes clear margin > 0.5; 0 of those also clear it against all 39; best margin +2.67._

---

## Stro-6 - Uterine smooth muscle cell (uSMC)

**What it is.** Smooth muscle. In the trajectory it sits between the dS states and the PV/Stro-4 branches, i.e. upstream of the perivascular compartment.

**Cycle phase.** Not phase-restricted.

**Transcription factors (paper).** Smooth-muscle programme; contractile genes dominate the centroid.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**DES** (+0.98), **ACTG2** (+0.94), **MYLK** (+0.82), **CNN1** (+0.81), **LEFTY2** (+0.75), **ATP1B1** (+0.62), **PDLIM7** (+0.60), **RAMP1** (+0.56), **PAGE4** (+0.55), **SVIL** (+0.55)

**Of those, also specific against all 39 subtypes:** **DES** (+0.98), **ACTG2** (+0.94), **MYLK** (+0.82), **CNN1** (+0.81), **LEFTY2** (+0.75), **PDLIM7** (+0.60), **RAMP1** (+0.56), **PAGE4** (+0.55)

**With a strict 1:1 mouse ortholog:** **DES** (+0.98), **ACTG2** (+0.94), **MYLK** (+0.82), **CNN1** (+0.81), **ATP1B1** (+0.62), **PDLIM7** (+0.60), **RAMP1** (+0.56), **SVIL** (+0.55)

**Verdict - USABLE.** 12 genes, 10 specific against all 39. Clean contractile signature (DES, ACTG2, MYLK, CNN1).

_12 genes clear margin > 0.5; 10 of those also clear it against all 39; best margin +0.98._

---

## Stro-7 - C7 Fibroblast

**What it is.** The structural, ECM-organising endometrial/myometrial fibroblast. Branches off the dS trajectory alongside uSMC. This is the subtype your reviewer knows.

**Cycle phase.** Not phase-restricted.

**Transcription factors (paper).** ECM and complement programme.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**C7** (+1.87), **CCDC80** (+1.82), **ADH1B** (+1.54), **OGN** (+1.31), **DCN** (+1.27), **DPT** (+1.22), **SERPINF1** (+1.18), **MGP** (+1.16), **APOD** (+1.15), **SRPX** (+1.06)

**Of those, also specific against all 39 subtypes:** **C7** (+1.87), **CCDC80** (+1.82), **ADH1B** (+1.54), **OGN** (+1.31), **DCN** (+1.27), **DPT** (+1.22), **SERPINF1** (+1.18), **MGP** (+1.16)

**With a strict 1:1 mouse ortholog:** **C7** (+1.87), **CCDC80** (+1.82), **OGN** (+1.31), **DCN** (+1.27), **DPT** (+1.22), **SERPINF1** (+1.18), **MGP** (+1.16), **APOD** (+1.15)

**Verdict - BEST-DEFINED FIBROBLAST IN THE ATLAS.** 51 genes clear the margin and 36 hold against all 39 subtypes - by far the highest tissue-level specificity of any stromal subtype. This is why the C7 mapping is the one that works everywhere.

_51 genes clear margin > 0.5; 36 of those also clear it against all 39; best margin +1.87._

---

## Stro-8 - dS1 - early decidualised stroma

**What it is.** STAGE, not a cell type: the earliest of three decidualisation states. Ulrich order them dS1 -> dS2 -> dS3 as a temporal progression.

**Cycle phase.** Peaks in early secretory phase.

**Transcription factors (paper).** STAT3 highest here and in dS2; FOXO1 high across all three; HOXA10 and NR2F2 LOWEST here. Stage-specific: PBX3, NFKB2.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**IGFBP1** (+0.77), **MMP10** (+0.67), **PLIN2** (+0.66), **CXCL8** (+0.61), **WTAP** (+0.58), **MMP3** (+0.52), **S100A6** (+0.48), **TWISTNB** (+0.45), **TGM2** (+0.44), **TMEM45A** (+0.44)

**Of those, also specific against all 39 subtypes:** **IGFBP1** (+0.77), **MMP10** (+0.67), **MMP3** (+0.52)

**With a strict 1:1 mouse ortholog:** **IGFBP1** (+0.77), **PLIN2** (+0.66), **WTAP** (+0.58), **MMP3** (+0.52), **S100A6** (+0.48), **TGM2** (+0.44), **RND3** (+0.44), **APOE** (+0.41)

**Verdict - USABLE ONLY FROM CENTROIDS, NOT FROM THE PUBLISHED DE LIST.** Six genes clear a 0.5 margin, in order: IGFBP1 (+0.77), MMP10 (+0.67), PLIN2 (+0.66), CXCL8 (+0.61), WTAP (+0.58), MMP3 (+0.52). Three of those six appear NOWHERE in Dataset S6C - IGFBP1 and MMP10 at ranks 1 and 2, and MMP3 at rank 6. Those same three are also exactly the three of the six that clear the all-39 margin (IGFBP1 +0.77, MMP10 +0.67, MMP3 +0.52) - WTAP comes closest at +0.48, just under the 0.5 bar used throughout this document, while the two S6C-listed genes that outrank MMP3 on stromal margin fail the whole-tissue test outright (PLIN2 -0.18, CXCL8 -1.90) and are unusable outside a purely stromal object. PRL, the other canonical decidualisation marker, IS in the centroid matrix and does peak in Stro-8 (highest of all 39 subtypes) but only at 0.089 log1p CP10K - barely detected anywhere in this atlas, which is expected for cycling rather than pregnant endometrium. So PRL is directionally right but too sparse to use as a marker here. Build a dS1 panel from S5B centroids using IGFBP1/MMP10/MMP3; the published DE list omits all three.

_6 genes clear margin > 0.5; 3 of those also clear it against all 39; best margin +0.77._

---

## Stro-9 - dS2 - mid decidualised stroma

**What it is.** STAGE. The middle decidualisation state.

**Cycle phase.** Peaks mid/late secretory, declines in proliferative phase.

**Transcription factors (paper).** STAT3 high; HOXA10 and NR2F2 intermediate. Stage-specific: MEIS1, MAFB, XBP1.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**SCGB1D4** (+0.39), **IER2** (+0.32), **RASD1** (+0.27), **CTD-3252C9.4** (+0.19), **FOSB** (+0.17), **IER3** (+0.15), **HNRNPU** (+0.15), **DDX5** (+0.14), **STAT3** (+0.13), **TPT1** (+0.11)

**Of those, also specific against all 39 subtypes:** _none_

**With a strict 1:1 mouse ortholog:** **IER2** (+0.32), **RASD1** (+0.27), **FOSB** (+0.17), **IER3** (+0.15), **HNRNPU** (+0.15), **DDX5** (+0.14), **STAT3** (+0.13), **RHOU** (+0.10)

**Verdict - NOT MARKER-DEFINABLE.** ZERO genes clear a 0.5 margin; best is SCGB1D4 at +0.39, followed by immediate-early genes (IER2, FOSB, IER3, RASD1) and housekeeping transcripts. dS2 cannot be identified by markers - 'unresolved' is the correct result, not a failure of your data.

_0 genes clear margin > 0.5; 0 of those also clear it against all 39; best margin +0.39._

---

## Stro-10 - dS3 - late / senescent decidualised stroma

**What it is.** STAGE. The latest decidualisation state, described in the paper as senescent.

**Cycle phase.** Accumulates through mid/late secretory and peaks in the proliferative phase.

**Transcription factors (paper).** STAT3 LOW; HOXA10 and NR2F2 HIGHEST here. Stage-specific: CERS4, MXD4, ING4.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**COL3A1** (+0.34), **HOXA11** (+0.30), **HOXA10** (+0.27), **GABPB1-AS1** (+0.23), **APCDD1** (+0.23), **MARCKSL1** (+0.22), **WNT2** (+0.18), **CHCHD10** (+0.18), **PIK3R1** (+0.18), **PGR** (+0.17)

**Of those, also specific against all 39 subtypes:** _none_

**With a strict 1:1 mouse ortholog:** **COL3A1** (+0.34), **HOXA11** (+0.30), **HOXA10** (+0.27), **APCDD1** (+0.23), **MARCKSL1** (+0.22), **WNT2** (+0.18), **CHCHD10** (+0.18), **PIK3R1** (+0.18)

**Verdict - NOT MARKER-DEFINABLE.** ZERO genes clear a 0.5 margin; best is COL3A1 at +0.34. The list is canonical endometrial-stromal identity (HOXA11, HOXA10, PGR, APCDD1, WNT2) at graded rather than specific levels - consistent with the paper distinguishing dS3 by HOXA10/NR2F2 being HIGHEST, not by any gene being exclusive.

_0 genes clear margin > 0.5; 0 of those also clear it against all 39; best margin +0.34._

---

## Stro-11 - Stromal progenitor

**What it is.** Described as a likely multipotent, cross-compartment progenitor and placed at the ROOT of the stromal trajectory: Stro-11 -> dS1/dS2/dS3 -> {C7 fibroblast, uSMC} -> {vascular progenitor, PV}.

**Cycle phase.** Enriched in endometrium; peaks in proliferative phase.

**Transcription factors (paper).** FOXM1, E2F2/E2F7/E2F8, MYBL1/MYBL2 - all proliferation regulators.

**High-information markers** (specificity margin vs the other 10 stromal subtypes):

**TOP2A** (+2.10), **CENPF** (+1.99), **HMGB2** (+1.90), **UBE2C** (+1.82), **MKI67** (+1.79), **STMN1** (+1.76), **KIAA0101** (+1.75), **TUBA1B** (+1.72), **TYMS** (+1.71), **HIST1H4C** (+1.61)

**Of those, also specific against all 39 subtypes:** **TOP2A** (+2.10), **UBE2C** (+1.82), **KIAA0101** (+1.75), **PTTG1** (+1.59), **TUBB** (+1.49), **CDKN3** (+1.34), **MXD3** (+1.05), **CKAP2** (+1.04)

**With a strict 1:1 mouse ortholog:** **TOP2A** (+2.10), **CENPF** (+1.99), **HMGB2** (+1.90), **UBE2C** (+1.82), **MKI67** (+1.79), **STMN1** (+1.76), **TUBA1B** (+1.72), **TYMS** (+1.71)

**Verdict - MOST MARKER-DEFINABLE, BUT IT IS A STATE.** 157 genes clear the margin, 36 above 1.0 - the most of any stromal subtype. Every one is cell cycle (TOP2A, MKI67, CENPF, UBE2C, TYMS). A match to Stro-11 means shared proliferation, not shared identity.

_157 genes clear margin > 0.5; 14 of those also clear it against all 39; best margin +2.10._

---

## Five things that follow from this

1. **Stro-7 is the only stromal subtype with a large, tissue-specific signature.** 36 of its 51
   margin-clearing genes hold up against all 39 subtypes. This is the structural reason the C7
   mapping reproduces everywhere and the others do not.
2. **Stro-5 is not a fibroblast reference.** Its centroid is NK/T and mast (GNLY, TPSB2, TPSAB1,
   GZMA, NKG7) plus PTPRC, and the paper itself notes the Imm-2 fibrocytes resemble it. Failing to
   recover it in nuclei is the correct outcome.
3. **Stro-4 has no stromal-specific markers.** All 38 of its margin-clearing genes are shared with
   endothelium, which is precisely the transitional identity the paper describes - which is also
   why it is a sound negative control for a fibroblast object.
4. **The dS series is a gradient, not three types.** Zero specific markers for dS2 and dS3; the
   separating axis is HOXA10/NR2F2 rising and STAT3 falling. Report dS assignments as positions on
   a decidualisation axis, not as discrete subtype calls.
5. **Do not build a dS1 panel from Dataset S6C.** Of the six genes clearing a 0.5 margin for
   Stro-8, three are absent from the published DE sheet - IGFBP1 and MMP10 at ranks 1 and 2, and
   MMP3 at rank 6 - and those three are exactly the three of the six that also clear the all-39
   margin. The two S6C-listed genes that outrank MMP3 on stromal margin (PLIN2, CXCL8) both fail
   the whole-tissue test outright.
   PRL is present in the centroid matrix and peaks in Stro-8 as expected, but at 0.089 log1p
   CP10K it is barely detected anywhere in this atlas and is not usable as a marker from it.
   Use S5B and take IGFBP1/MMP10/MMP3.
