# Bringing the remaining Methods code into the pipeline

paper_code/ and the code listed here are cited in the Methods but are not steps of the pipeline yet. Each item is
ported the same way: its code goes into workflow/scripts/, a rule with declared inputs and outputs goes into
workflow/rules/, it runs into a scratch folder, and its outputs are compared with the existing ones (identical, or
every difference explained). Existing outputs are never overwritten. The item then leaves paper_code/.

Order (number of Methods statements in brackets; from the Methods-vs-code check of 31 Oct):
1. With the v10 release: the 30 Oct human pseudotime refit and its MI [10]; the cross-species MI heatmap
   selection [11]; the 22 Oct harmonized pseudotime, replaced by the 30 Oct refit [22].
2. Epithelial notebooks: mouse compartments [35], human layers [31] (paper_code/epithelial_notebooks).
3. IVD and Wang comparisons [35] (code outside the repository; IVD DE code to be added).
4. Ulrich/HECA mapping [4] (paper_code/ulrich_mapping; also analysis/hammoud/WRITEUP_2026_09_20 code on disk).
5. Stadtmauer consensus bundle [20], X-Mens fit with uninduced fibroblasts [16], hand-staged inputs [9]
   (resources/methods_numbers/artifact_inputs).
6. SCENIC and SAMap scripts that live outside the repository (24_scenic, 19_homology); SCENIC stages outside X [2].
7. Small items: proteome release check [1], zone-module seed-1 run [1], plot_pt_spatial_sections.py [1],
   library defaults [2].
Upstream of the pipeline and kept as documented provenance: Cell Ranger, CellBender, Slide-tags [12] (steps 01-06).
