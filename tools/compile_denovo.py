#!/usr/bin/env python3
"""Compile a read-only DRAFT presentation notebook from a frozen 04c CSV.

The generated notebook embeds the committed cluster map verbatim, loads the
frozen latent object, shows raw and verification dotplots, applies the shared
canonical producer in memory, and finally proves that the CSV has not changed.
It contains no object/CSV write calls and makes no annotation judgments.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import sys

import nbformat as nbf
from nbformat.v4 import new_code_cell as code, new_markdown_cell as md, new_notebook
import yaml


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


ap = argparse.ArgumentParser()
ap.add_argument("--species", choices=["human", "mouse"], required=True)
ap.add_argument("--lineage", required=True)
ap.add_argument("--res", required=True)
ap.add_argument("--csv", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--level", choices=["fine"], default="fine")
ap.add_argument("--repo", default="/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
ap.add_argument("--spec", default="resources/lineage_spec.yaml")
ap.add_argument("--vocab", default="resources/annotation_vocabulary.yaml")
args = ap.parse_args()

repo = Path(args.repo)
csv_path = Path(args.csv)
spec_path = Path(args.spec)
vocab_path = Path(args.vocab)
rows = list(csv.DictReader(csv_path.open()))
if not rows or not {"cluster", "label"} <= set(rows[0]):
    raise SystemExit(f"{csv_path}: expected nonempty cluster,label CSV")
annotation = {row["cluster"].strip(): row["label"].strip() for row in rows}
spec = yaml.safe_load(spec_path.read_text())
lineage = spec["species"][args.species]["lineages"][args.lineage]
markers = lineage.get("markers", {})
csv_hash = sha256(csv_path)
vocab_hash = sha256(vocab_path)
stamp = datetime.date.today().isoformat()

annotation_literal = json.dumps(annotation, indent=4, sort_keys=True)
markers_literal = json.dumps(markers, indent=4)
cells = [
    md(f"""# DRAFT — 04c {args.species} {args.lineage} de novo annotation

**Presentation scaffold only.** Generated {stamp} from the committed decision
record `{csv_path}` and locked vocabulary `{vocab_path}`. The map below is
frozen from the CSV; this notebook does not infer, edit, or write annotations.

- resolution: `{args.res}`
- decision CSV SHA-256: `{csv_hash}`
- vocabulary SHA-256: `{vocab_hash}`
"""),
    md("## 1. Frozen configuration and annotation map"),
    code(f'''SPECIES = {args.species!r}
LINEAGE = {args.lineage!r}
RES = {args.res!r}
REPO = {str(repo)!r}
CSV_PATH = {str(csv_path)!r}
VOCAB_PATH = {str(vocab_path)!r}
CSV_SHA256 = {csv_hash!r}
VOCAB_SHA256 = {vocab_hash!r}
ANNOTATION = {annotation_literal}
MARKERS = {markers_literal}'''),
    md("## 2. Load the frozen lineage latent and verify complete cluster coverage"),
    code('''import hashlib, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import pandas as pd, scanpy as sc

REPO = Path(REPO)
LATENT = REPO / "data/interim" / f"latent_{SPECIES}_{LINEAGE}.h5ad"
adata = sc.read_h5ad(LATENT)
assert RES in adata.obs, f"{RES} absent from {LATENT}"
clusters = set(adata.obs[RES].astype(str).unique())
assert clusters == set(ANNOTATION), (
    f"frozen map mismatch: missing={sorted(clusters-set(ANNOTATION))}, "
    f"extra={sorted(set(ANNOTATION)-clusters)}")
print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes; {len(clusters)} clusters")'''),
    md("## 3. Raw marker dotplot by frozen cluster"),
    code('''plot_sets = {name: [gene for gene in genes if gene in adata.var_names]
             for name, genes in MARKERS.items()}
plot_sets = {name: genes for name, genes in plot_sets.items() if genes}
if plot_sets:
    sc.pl.dotplot(adata, plot_sets, groupby=RES, standard_scale="var",
                  figsize=(10, max(6, 0.22 * len(clusters))), show=True)
else:
    print("No configured marker genes are present.")'''),
    md("## 4. Apply the committed map in memory using the shared producer"),
    code('''import sys
sys.path.insert(0, str(REPO / "workflow/scripts"))
from apply_annotation import (apply_map, build_vocabulary_rollup,
                              validate_canonical_annotations)
import yaml

vocab = yaml.safe_load(Path(VOCAB_PATH).read_text())
rollup = build_vocabulary_rollup(vocab)
summary = apply_map(adata, RES, ANNOTATION, rollup=rollup, output_mode="fine")
kept = adata[~adata.obs.excluded].copy()
validate_canonical_annotations(kept, vocab, context=f"{SPECIES}/{LINEAGE} draft")
kept.obs["display_label"] = kept.obs.cell_subtype.astype(str).where(
    kept.obs.cell_subtype.astype(str).ne(""), kept.obs.cell_type.astype(str))
print(summary.to_string())
print("\\ncanonical cell types:")
print(kept.obs.cell_type.value_counts().to_string())
print("\\ncanonical subtypes (nonempty):")
print(kept.obs.loc[kept.obs.cell_subtype.ne(""), "cell_subtype"].value_counts().to_string())'''),
    md("## 5. Verification dotplot ordered by the frozen annotation"),
    code('''if plot_sets:
    sc.pl.dotplot(kept, plot_sets, groupby="display_label", standard_scale="var",
                  figsize=(10, max(6, 0.35 * kept.obs.display_label.nunique())), show=True)
sc.pl.embedding(kept, "X_umap_scVI", color=[RES, "display_label"],
                legend_loc="right margin", show=True)'''),
    md("## 6. Freeze check — assert the committed inputs are unchanged"),
    code('''def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
assert _sha(CSV_PATH) == CSV_SHA256, "decision CSV changed; regenerate this draft"
assert _sha(VOCAB_PATH) == VOCAB_SHA256, "vocabulary changed; regenerate this draft"
now = pd.read_csv(CSV_PATH, dtype=str)
now_map = dict(zip(now.cluster.str.strip(), now.label.str.strip()))
assert now_map == ANNOTATION, "embedded map differs from committed CSV"
print("PASS: CSV and vocabulary hashes match; embedded map is exact.")'''),
]

notebook = new_notebook(cells=cells)
notebook.metadata["kernelspec"] = {
    "display_name": "Python 3", "language": "python", "name": "python3"}
notebook.metadata["language_info"] = {"name": "python"}
notebook.metadata["xmens"] = {
    "status": "DRAFT", "schema_version": "xmens.annotation.v1",
    "decision_csv_sha256": csv_hash, "vocabulary_sha256": vocab_hash,
}
nbf.validate(notebook)
out = Path(args.out)
if "DRAFT" not in out.name or out.parent.name != "public":
    raise SystemExit("refusing: output must be notebooks/public/DRAFT_*.ipynb")
if out.exists():
    raise FileExistsError(f"refusing to overwrite draft notebook: {out}")
out.parent.mkdir(parents=True, exist_ok=True)
nbf.write(notebook, out)
print(f"compiled DRAFT -> {out} ({len(cells)} cells)")
