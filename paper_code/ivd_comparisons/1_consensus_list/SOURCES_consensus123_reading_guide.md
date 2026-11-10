# Where the 123-gene consensus comes from: documents and lines to read

**Status, 2026-10-22: correction applied.** The line numbers below describe the original pipeline, which is unchanged. The corrected list (Stadtmauer table genes called with the authors' criteria) was rebuilt by `rebuild_consensus123_stadtmauer_thresholds.py`; see `README_consensus123_stadtmauer_fix.md`.

Repository root (cluster): `/n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro_final_defense_20260812/`
Pipeline that uses the list: `/n/lab_storage/eddy_lab/users/nhilgert/cagri_ivd_end_to_end_20260815/`
Line numbers were read from the files on the cluster; page numbers are journal pages (PDF page in parentheses).

## Stadtmauer and Wagner 2022 (`papers/stadtmauer/`)

**`stadtmauer.pdf`**
- p. 158 (PDF 4): Methods, "Cell culture and in vitro decidualization" and "Bulk RNA sequencing": immortalized ESF (ATCC CRL-4003); base medium, cAMP + MPA, PGE2 + MPA, MPA alone, PGE2 alone; 6 days.
- p. 159 (PDF 5): the differential expression thresholds (edgeR; FDR < 10^-6, fold change > 1.5, > 3 TPM). This is the criterion to hold the table values against.
- p. 162 (PDF 8): k-means clusters. Cluster 2 = genes associated with MPA alone and PGE2 alone, no distinguishing functional pattern (KLF4's cluster). Clusters 9 and 5 = genes up with the decidualizing treatments; Supplementary Table S2 cited here.
- p. 163 (PDF 9): Supplementary Table S3 cited; ADAMTS5, CLU and SOX4 described as elevated only in cAMP + MPA-treated cells.

**`tables2.docx`** (Supplementary Table S2: selected genes related to decidual cell identity, first table, or decidual phenotype, second table). Columns: gene, cluster, TPM in base / MPA+PGE2 / MPA+cAMP / MPA / PGE2, then FDR base vs MPA+PGE2, base vs MPA+cAMP, MPA+PGE2 vs MPA+cAMP. Row numbers count the header as row 1.
- First table, row 4, KLF4: cluster 2; TPM 2.76 / 3.75 / 2.84 / 5.84 / 4.4; FDR 7.4e-3 / 0.40. Fold change 1.36 (PGE2 arm) and 1.03 (cAMP arm): neither passes p. 159.
- First table, row 13, TSHZ1: cluster 8; TPM 30.2 / 34.2 / 45.1 / 12.9 / 15.7; FDR 0.23 / 4.4e-3. Fold change 1.13 and 1.49: neither passes.
- Second table (no header row of its own): the caption's "phenotype" group, the effector genes BMP2, SST, IGFBP1, LEFTY2 and PRL. The first table is the "identity" group (17 genes, mostly transcription factors). The original parser skipped this headerless table until the BMP2 fix (methods.md section 2.4).

**`tables3.docx`** (Supplementary Table S3: selected markers of the senescent decidual cell state described by Lucas et al. 2020; same columns as S2)
- Row 3, ADAMTS5: TPM 146 / 127 / 226; FDR 0.49 / 7.8e-7. Significantly up in the cAMP arm (1.55-fold), unchanged in the PGE2 arm.
- Row 5, CLU: TPM 162 / 142 / 454; FDR 0.75 / 3.7e-13. Significantly up in the cAMP arm (2.8-fold), unchanged in the PGE2 arm.
- Row 9, SOX4: TPM 117 / 103 / 184; FDR 0.72 / 3.4e-3. Not significant in either arm.

**`tables4.docx`** (Supplementary Table S4): DREMI scores (single-cell co-expression) between the senescence genes and SOX4 at day 2 and day 6, cited on p. 166 (PDF 12). There is no comparison with untreated cells, so it gives no up/down direction; the parser does not read it.

Supplementary Table S1 is not in the folder; the paper cites it (p. 157, PDF 3) for references on PGE2 production at the fetal-maternal interface, not for expression data.

**`stadtmauer_manual_overlay.csv`** (hand-entered directions, 38 rows)
- Line 12: `KLF4,up,boundary FC~1 curation call`, the judgment call. The other 37 rows are directions transcribed from the Results text.

## The parser and the vote (repository root)

**`derive_signatures.py`**
- Line 98: `STAD_FLOOR = 3.0`, base TPM raised to 3 before comparing.
- Line 99: `STAD_ARMS = [3, 4]`, only the MPA+PGE2 and MPA+cAMP columns vote.
- Line 167: `def parse_stadtmauer`.
- Line 181: `fb = max(base, STAD_FLOOR)`.
- Line 186: `"up" if treat > fb else "down"`. Every tabulated gene gets a direction from the raw TPM difference; FDR and fold change are never checked. This is the root of both problems below.
- Line 200: `if g in overlay`, the hand-entered file overrides the table (how KLF4's computed call was replaced).
- Lines 88-93: `RYT_SHEETS`, the four Rytkönen sheets used (3 d and 8 d vs untreated).
- Line 233: `AGH_TP`, Aghajanova 6 h and 48 h from Table 4, 14 d from Table 5.
- Line 311: Aghajanova rows kept only for women without endometriosis (`nonendo`).

**`signatures/stadtmauer2022_hits.csv`** (the parser's output for Stadtmauer)
- Lines 6-7 ADAMTS5, 9-10 CLU, 18-19 SOX4: one `down` and one `up` row each (the non-significant PGE2-arm dip counted as `down`), so all three are marked conflicting and excluded.
- Line 33: TSHZ1 `up` from the table.
- Line 44: KLF4 `up` from the override.

**`build_consensus.py`**
- Line 77: `PAPER_SPEC`, the three voting studies (Lucas = annotation only).
- Line 86: `MIN_MATURE = 2`.
- Line 111: within a study, `up` only if all rows are up, `down` only if all are down, otherwise conflicting (`biphasic`) and excluded.
- Lines 153-154: the vote, at least two clean calls, all in the same direction.

**`methods.md`** (the repository's own notes)
- Lines 31-112 (section 2): how the consensus was built. Section 2.3 (line 66) covers the Stadtmauer overlay; lines 91-95 are the KLF4 judgment; line 396 lists it among the judgment calls.
- Line 81 says 23 of the 38 hand-entered genes are in the 123; the current files give 22 (the note predates BMP2 moving to the parsed table).

## Aghajanova et al. 2011 (`papers/aghajanova/aghajanova.pdf`)
- p. 802 (PDF 2): treatments (vehicle, E2, E2 + P4), harvest times (0 h, 6 h, 48 h, 14 d), duplicate cultures, and the normalization giving the "pure P4 response" (E2 + P4 relative to E2).
- p. 803 (PDF 3): statistics (two-way ANOVA, Benjamini-Hochberg, > 1.5-fold).
- p. 805 (PDF 5): Table 3, endometriosis vs control at time 0 (not used).
- p. 807 (PDF 7): Table 4, E2 + P4 vs E2 at 6 h and 48 h.
- pp. 808-810 (PDF 8-10): Table 5, E2 + P4 vs E2 at 14 days.
- p. 814 (PDF 14): comparison with an 8-Br-cAMP decidualization list (Supplemental Table S7).
- `signatures/aghajanova2011_invitro_decidualization_nonendo.csv`: line 44 OSR2 and line 51 NPR3, both listed only at 48 h.

## Rytkönen et al. 2019 (`papers/rytkonen/`)
- `rytkonen.pdf`, PDF p. 2: Methods; T HESC (ATCC CRL-4003); 3 days (n = 3) and 8 days (n = 2) of 8-Br-cAMP + MPA; edgeR (FDR < 0.05, |fold change| > 2, TPM > 2).
- `Supplemental_Table_S1_AtoH.xlsx`: sheets S1A-S1D used (untreated to 3 d up/down, untreated to 8 d up/down); S1E (intersection), S1F-S1H (3 d to 8 d) and S1I (transcription factors) not used.
- `signatures/rytkonen2019_complete_signature.csv`: line 567 KLF4 up (3 d); line 1137 TSHZ1 up (3 d).

## The frozen list the X-Mens pipeline reads (`cagri_ivd_end_to_end_20260815/`)
- `resources/consensus123.csv`: line 59 KLF4 and line 112 TSHZ1 (both Rytkönen + Stadtmauer); line 79 NPR3; line 83 OSR2 (all three studies).
- `config.yaml` lines 23-25: the consensus file the pipeline reads, pinned by checksum.

## What the correction changed (applied 2026-10-22)
- KLF4 and TSHZ1 lose their Stadtmauer vote and drop out (Rytkönen only).
- ADAMTS5 and CLU gain a clean Stadtmauer `up` (cAMP arm significant, PGE2 arm unchanged) and enter, with Rytkönen `up`.
- The list stays at 123 genes. X-Mens: ADAMTS5 is called and up (log2 fold change +1.86, adjusted p 0.012); CLU, KLF4 and TSHZ1 are not called. Concordance is now 20 of 50 (40.0%; 45.6% expected by chance) instead of 19 of 49 (38.8%; 45.2%).
