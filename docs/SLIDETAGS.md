# Slide-tags mapping code

Spatial coordinates (step 02, workflow/rules/02_slidetags.smk) come from a fork of the Slide-tags pipeline, run
separately from Snakemake:

- repository: https://github.com/nhilgert/slide-tags.git (tag v0.1-fasrc)
- commit used: 0c1630520ffcfd7ea5d324c4db52660733f3cb40

Until 31 Oct 2026 the fork was linked as a git submodule at external/slidetags. The checkout stays on disk there
(ignored by git); its FASRC_SETUP.md describes how it was run.
