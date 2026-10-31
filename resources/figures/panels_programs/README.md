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
