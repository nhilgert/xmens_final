# Staged inputs of Fig. 5E (rule panels_programs_5e)

Copied byte for byte on 2026-10-30 from
/n/lab_storage/eddy_lab/users/nhilgert/xmens_figures/scenic_quadrants/ (the 2026-10-26 quadrant figure folder):

- regulon_quadrants_input.csv
  sha256 0b93fd193dc8a17e3fad72a47c196f9c97e7b80a45992dafa4b0d859a28e0d40 (pinned in config/panels_programs.yaml)
- scenic_quadrant_rho.provenance.json
  sha256 fbd1ac9628bda1c8abca0c11be2a530e03d88d4c69652adc7c370ec32dbf2130

regulon_quadrants_input.csv has one row per transcription factor with a floor-6 consensus SCENIC regulon in
both species (Ensembl one-to-one pair; 111 rows): human_tf, mouse_tf, rho_human, rho_mouse, label. rho is the
Spearman correlation of regulon activity (pySCENIC aucell, weighted) with pt_harmonized_rank over the
fibroblast/decidual nuclei. It was written by scenic_quadrant_rho.py in the same folder (sha256
19996b6d824fc0d7ea394b4080bb5df63a7b80b2823cefe6e76deae60ad0f0b3); the paths and sha256 of that script's
inputs are in scenic_quadrant_rho.provenance.json. The panel rule does not use the label column: the named
regulons come from config/panels_programs.yaml.

## 2026-10-30: _fibdec_latent table (pinned instead of the table above once installed)

- regulon_quadrants_input_fibdec_latent.csv
  sha256 6c3c99e07445bb3142af26f8900ea88e67eeb4f4a261de8dbce42a50ad4ae317
- scenic_quadrant_rho_fibdec_latent.provenance.json
  sha256 c232109d12a3aef3a1dcba9624adfa7774e1783d346edc7cfd2e76b608d63fbb

Same 111 rows and columns. rho_human is recomputed with the human pseudotime of 30 Oct 2026: pt_rank of the side table
built on the human stromal scVI latent refit on the 47,533 fibroblast/decidual nuclei only (as for mouse), joined on
obs_names (run4_a40/fibdec/pt_human_pca_released/cells.parquet, sha256 956a07cd...f94f). rho_mouse is unchanged
(byte-identical per-regulon table). Written by scenic_quadrant_rho_fibdec_latent.py, which differs from
scenic_quadrant_rho.py only in where the human pseudotime is read (diff in the provenance folder). Quadrant classes
at |rho| >= 0.2 are unchanged for all 111 regulons (18 / 8 / 14 / 0, 71 below).
