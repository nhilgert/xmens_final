# De novo DRAFT-v2 notebook regeneration report

Date: 2026-08-12  
Repository: `/n/lab_storage/eddy_lab/users/nhilgert/xmens_final`  
Vocabulary: `resources/annotation_vocabulary.yaml` at SHA-256 `d02b68120a62d1764b6f0473976fcd72a85f5c30a081b3a857db9215e6a38847`  
Certified source run: `xmens_annotation_v2_20260811` (`MIGRATION_STATUS.json` status `PASS`)

## Outcome

PASS. All 11 parallel `DRAFT_v2_*_denovo.ipynb` notebooks were generated mechanically from the committed decision CSVs, executed end-to-end in isolated scratch, and independently audited as `CLEAN_PASS`.

- 11/11 generated.
- 11/11 executed with no captured error output.
- 66/66 code cells executed with exact execution counts 1–6 per notebook.
- 33/33 expected inline PNG plots rendered and passed PNG structural validation.
- 11/11 embedded `ANNOTATION` maps exactly equal their committed CSV maps.
- No annotation label or CSV decision was invented, changed, or inferred from marker plots.

## Per-notebook results

| Species | Lineage | New public DRAFT | Generated | Execution | Map = CSV | PNGs | Public SHA-256 prefix | Executed SHA-256 prefix |
|---|---|---|---:|---|---:|---:|---|---|
| `mouse` | `epithelial` | `DRAFT_v2_04c_annotate_mouse_epithelial_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `26d5487ea2ce…` | `3ea483a328a1…` |
| `mouse` | `perimetrium` | `DRAFT_v2_04c_annotate_mouse_perimetrium_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `fdb1b812a141…` | `21b46de0ffbf…` |
| `mouse` | `stromal` | `DRAFT_v2_04c_annotate_mouse_stromal_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `56296f1dce8d…` | `de9e416aded1…` |
| `mouse` | `immune` | `DRAFT_v2_04c_annotate_mouse_immune_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `c68225094f9d…` | `c396237942ae…` |
| `mouse` | `vascular` | `DRAFT_v2_04c_annotate_mouse_vascular_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `fcecfb44fa0c…` | `b12fde2e5e93…` |
| `mouse` | `smc` | `DRAFT_v2_04c_annotate_mouse_smc_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `9ad474ec4873…` | `e98bb4e7abe4…` |
| `human` | `epithelial` | `DRAFT_v2_04c_annotate_human_epithelial_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `88b69a802dd4…` | `b5e01f918142…` |
| `human` | `stromal` | `DRAFT_v2_04c_annotate_human_stromal_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `5aa81671fd92…` | `429b3a044377…` |
| `human` | `immune` | `DRAFT_v2_04c_annotate_human_immune_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `299a2594b126…` | `480cd9fc26b2…` |
| `human` | `perivascular` | `DRAFT_v2_04c_annotate_human_perivascular_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `abd3ea19285d…` | `11a946ec6186…` |
| `human` | `smc` | `DRAFT_v2_04c_annotate_human_smc_denovo.ipynb` | Yes | `CLEAN_PASS` | Yes | 3 | `2b295ca5552e…` | `a46063e37079…` |

Full hashes and per-cell/per-image evidence are in `data/interim/_draft_v2_exec_20260812/final_audit/DENOVO_V2_EXECUTION_AUDIT.json` (SHA-256 `97b11a8c70302b8e93b4e7ef5991455cd41687878cd6ee928c333c45ddda9d2f`).

## Generator changes

`tools/compile_denovo.py` was mechanically updated and independently reviewed at SHA-256 `82c7709ecfead2dd17ffe10ea4eb346c6e7cbcec8f048d4a2b259f2cabfb5a11`.

Changes were limited to plumbing, provenance, validation, and presentation scaffolding:

- emit only the parallel `DRAFT_v2_04c_annotate_{species}_{lineage}_denovo.ipynb` names and refuse overwrite;
- bind the exact new vocabulary, lineage specification, producer, v2 PASS/READY records, 11 committed CSVs, 11 frozen lineage latents, and 11 certified v2 reannotated objects by SHA-256;
- reject blank/duplicate/noncanonical CSV clusters, invalid labels, resolution drift, and incomplete marker panels;
- preserve numeric cluster order and the committed CSV row order;
- compute the existing mechanical per-cell marker-group mean scores for display only, print cluster score tables, and render the raw dotplot;
- apply the committed CSV map in memory with the shared producer, compare the retained cells/genes/types/subtypes exactly to the certified v2 object, render verification plots, echo the unchanged CSV, and print the positive freeze PASS line;
- stream large-file hashes, validate generated Python cells, and publish each new notebook with exclusive no-overwrite semantics.

The marker scores and plots never feed back into `ANNOTATION`; the map is read directly from each committed CSV.

Active stale paths fixed: none. Every active path in files touched for this run already uses `/n/lab_storage/eddy_lab/`. The literal `/n/eddy_lab/` occurs only as a forbidden-token sentinel in the preparation and audit code, not as an executable path.

## Locations and provenance

- New unexecuted public drafts: `notebooks/public/DRAFT_v2_04c_annotate_*_denovo.ipynb`
- New generation manifest: `notebooks/public/DRAFT_v2_04c_MANIFEST.json`
  - SHA-256: `061f38730cda723df8a3f76f117d2fbceb355702f584ea0989335fcfe3231b88`
- Executed review copies: `data/interim/_draft_v2_exec_20260812/executed/*.executed.ipynb`
- Per-notebook audits: `data/interim/_draft_v2_exec_20260812/inspect/*.json`
- Aggregate JSON/TSV: `data/interim/_draft_v2_exec_20260812/final_audit/`
- Generation log: `data/interim/_draft_v2_exec_20260812/generation.log`
  - SHA-256: `e5111007f89313c4655c1d44af213a2444dbbea1a9937e59c1ac022fe74ece55`

The original 11 `DRAFT_04c_*` notebooks and original `DRAFT_04c_MANIFEST.json` were retained byte-for-byte. The final auditor rechecked the original manifest SHA-256 `4cebfbe28a01e6029e48e561101a6036ed09628d5273f3ec2a865bb90b416fc2` and every original notebook hash.

## Execution

- Generation job `38498944`: `COMPLETED`, exit `0:0`, 4:55, CPU-only `eddy` node.
- Notebook array `38500442` (indices 0–10, maximum concurrency 3): every task `COMPLETED`, exit `0:0`, CPU-only `eddy` nodes.
- Aggregate/protected-live job `38501284`: `COMPLETED`, exit `0:0`, 9:01, CPU-only `eddy` node.
- Execution controller: `rsc_uv` `jupyter-nbconvert`.
- Notebook kernel: the requested `envs/samap/.venv/bin/python`, selected through the exact audited `python3` kernelspec after loading `gcc/14.2.0-fasrc01`.
- Matplotlib backend: `module://matplotlib_inline.backend_inline`.
- `allow_errors=True` was used only to preserve complete notebook diagnostics; the auditor required zero error outputs and failed closed otherwise.
- HOME, XDG, matplotlib, numba, Jupyter, Python bytecode, TMP, Torch, CUDA, Triton, Hugging Face, Keras, Dask, joblib, and Weights & Biases locations were all redirected under the execution scratch root.

Two non-content plumbing issues were encountered and resolved without changing annotations:

1. The requested SAMap environment contains the scientific/kernel stack but not `nbformat`/`nbconvert`; `rsc_uv` therefore controlled execution while the kernelspec launched the requested SAMap Python.
2. The first generation submission was rejected before entering the queue because its GPU-node exclusion list contained retired node names. The list was refreshed from Slurm; no job or scientific artifact resulted from that rejected submission.

## Integrity and safety

- Aggregate notebook audit: `CLEAN_PASS` (11 notebooks, 66 executed code cells, 33 PNGs, 0 errors).
- Protected-live content rehash: `PASS`, `n_protected=26`, `changes=[]`; report SHA-256 `0b6b24e3bb71ec53a109799736432c21a4c3e7e7b9ad9c8e57c7a476462d4b78`.
- Original public DRAFT set: unchanged.
- Certified v2 inputs and decision CSVs: unchanged and matched all reviewed hashes.
- Home-directory write audit: `PASS`. A completed recursive scan found zero files created or modified under `/n/home12/nhilgert` since 2026-08-12 00:00, and the Snakemake runtime-cache `tmp*` count was zero.
- Live `reannotated_*`, `refined_*`, `final_*`, deposits, and SAMap objects were not written by notebook generation or execution.
- Nothing was promoted. No promotion helper was invoked.
- No SAMap job or downstream workflow was run.

## STOP items

Content STOP items: none. The committed CSVs, generator maps, vocabulary, frozen latents, and certified v2 reannotated objects agreed for all 11 identities. No discrepancy was resolved by biological judgment.

These remain DRAFT notebooks for human review; this run only regenerated and verified their mechanical presentation of already committed annotation decisions.
