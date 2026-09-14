#!/usr/bin/env python3
"""Compile a read-only DRAFT presentation notebook from a frozen 04c CSV.

The generated notebook embeds the committed cluster map verbatim, loads the
frozen latent object, shows raw and verification dotplots, applies the shared
canonical producer in memory, and finally proves that the CSV has not changed.
It contains no object/CSV write calls and makes no annotation judgments.
"""
from __future__ import annotations

import argparse
import ast
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import stat
import sys

import nbformat as nbf
from nbformat.v4 import new_code_cell as code, new_markdown_cell as md, new_notebook
import yaml


EXPECTED_GLOBAL_SHA256 = {
    "vocabulary": "73efff69b3cfd3327d06b17ec05fd197da3ffefd5a363ceb40e8c4d31fef665f",
    "lineage_spec": "23be2688f07ba05e9f0dc63d6852fe844daa0a651a1e16660a5c55e88af82d53",
    "producer": "e57cd6e6b66ece138e770048060e7eece2eb6020395bc21b73aab4e78812b480",
    "migration_status": "def88e5240f653f98c8d78adf7f8dddcdb05bfda30d08bd0fada7df489969208",
    "migration_ready": "c26de83abdc9496cd1301470918ec39ecca1cf389ef0ae1c6504da1800d1c431",
}
EXPECTED_CSV_SHA256 = {
    ("human", "epithelial"): "2eb061d30d26b2455ddd7b8d894fc3f8bb30956f4aff987b60f68253f2e03f6e",
    ("human", "immune"): "8c24e61758fbc590e1bd8393fad26d028dbd6239edd024b065c57162474f2433",
    ("human", "perivascular"): "28fcac896c6b51fcc4bf2c3d89a366d730702e4cae27d2243b7da373497fd494",
    ("human", "smc"): "0886ff7dc5e94d3e31219ffbcca02b7bfe5f5e35ad44d49641e46102a6f79c5c",
    ("human", "stromal"): "aae1ea276571ee37bb56c0b704d7dcbc2747babb5b615c5671ee0eab015f2aef",
    ("mouse", "epithelial"): "0b47052c19eb2065d3266cb3b7d4aadc4298acd69c4dd171589a4e8f131d5be3",
    ("mouse", "immune"): "5ae30a059a3590c111bf4aebf8b9979250a78cf6938820ebda951db0e4c0d8f3",
    ("mouse", "perimetrium"): "dadb4d283ef511dce341e34c56fe16806b5434feea5dcdbd94c5c490f57d066c",
    ("mouse", "smc"): "25dec472cc51544300b3a8593ce493d56cabe2634f161b8e9638ebea9c47b7c6",
    ("mouse", "stromal"): "84fe2ed614de86e00403067f343980d8fd6966647d2c0313b760b9aeed7f11d2",
    ("mouse", "vascular"): "61b69436c526f4c91dfacd77112a3de17dc7ea3620d182c073e4edee38cb6adc",
}
EXPECTED_LATENT_SHA256 = {
    ("human", "epithelial"): "207d1f8542a918edb7755a416c9d273ef494e3877973f2a605edd9f36f505f75",
    ("human", "immune"): "50c0801938ed18e61a26f9d3f6a241106f7e9a5f2fc194984ac7785d0397b2a5",
    ("human", "perivascular"): "bb205ddfdc2f701d3f7babccbc8a9633c24925ad9014b93ff5e23638622671ab",
    ("human", "smc"): "cc9b1a30cce15b57632e6c32334ec6a888d26d7856d9ee2f68ab0f6bc216398f",
    ("human", "stromal"): "3e6a0e90699be4f86072d64ef896c139cd7cb3a44f226c46b70825beaa9acadd",
    ("mouse", "epithelial"): "75243c93d329466424aa338c466ac8b4a4df51ff4421795b2769bb7573526429",
    ("mouse", "immune"): "cac0499bd87bd314ef828cf14d2e776704389b5760018c0f5cb6594a6afe7363",
    ("mouse", "perimetrium"): "535ad3833eeb21c3b015554a90ba67c425ead35befb2feac36b37cfd637a66da",
    ("mouse", "smc"): "72512112c6eea4e477d4e29073942ad2f993efce3b3c03fc61bf0771556086c3",
    ("mouse", "stromal"): "79e8f3aa8e3b17f57548629c108900296cd9a6ced382a6ce03bad7ca02ff900e",
    ("mouse", "vascular"): "2ab3d7de2c48d5396ad3d9b1641b9baf5a37422ba3ebd6814a02c1d527c0bbc4",
}
EXPECTED_V2_SHA256 = {
    ("human", "epithelial"): "c3a873ca8ecff7341adec891834e0a7c2243341bba7413de5e3eea50e6a08182",
    ("human", "immune"): "810f2a391cac2d3065c06f7690841bc7b7ce9a437de75838e294c0c279f59575",
    ("human", "perivascular"): "af73161ee5c9cfc9f2cf43c38a41d366165a58a5b21a1facbcd3690c153d6dcc",
    ("human", "smc"): "d62f00723bd2a6d233aee1d4ba2d3129179ef49eb51fa5fa6e5406e723ab0d3c",
    ("human", "stromal"): "0d33b6176ce9573375bd1f681a2914a5971714d40d702233781f999a455e9c40",
    ("mouse", "epithelial"): "fe649c3a4520862bce9d240a8923d30a7b371e1a0ddf265a58dfa0ea162c4030",
    ("mouse", "immune"): "08778b953559738f0afc23ffbe3ec3e4921fae39f7ad23e88507ec594da8cdfa",
    ("mouse", "perimetrium"): "fc29226f7c0bc6d5a4ff20313b13cdc4d1051ced631ffb8db2c5183714cf7d12",
    ("mouse", "smc"): "f8052a23a401bfc75205ada85cb74fd60ba30caf4e4a317775cae8af7b335448",
    ("mouse", "stromal"): "496162972869d8babf72acf04b696b18bb07d5d7e5841f3678061102ff031d88",
    ("mouse", "vascular"): "9966b9f57f0d1424189b8803ed34859444d3ed2eba2f00ae200e0888b6ca3b79",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reject_symlink_components(path: Path) -> None:
    """Reject a missing or symlinked component in an absolute path."""
    if not path.is_absolute():
        raise SystemExit(f"refusing non-absolute path: {path}")
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        try:
            mode = os.lstat(current).st_mode
        except FileNotFoundError as exc:
            raise SystemExit(f"refusing missing path component: {current}") from exc
        if stat.S_ISLNK(mode):
            raise SystemExit(f"refusing symlinked path component: {current}")


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
ap.add_argument("--v2-root", default="data/interim/migrated/xmens_annotation_v2_20260811")
args = ap.parse_args()
identity = (args.species, args.lineage)
if identity not in EXPECTED_CSV_SHA256:
    raise SystemExit(f"unsupported reviewed lineage identity: {identity}")

repo_arg = Path(args.repo)
repo = repo_arg.resolve()
csv_path = Path(args.csv).resolve()
spec_path = Path(args.spec).resolve()
vocab_path = Path(args.vocab).resolve()
expected_repo = Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
if repo != expected_repo:
    raise SystemExit(f"refusing: repo must resolve exactly to {expected_repo}")
reject_symlink_components(expected_repo)
expected_csv = repo / f"resources/annotations/{args.species}_{args.lineage}.csv"
if csv_path != expected_csv:
    raise SystemExit(f"refusing: CSV must be exactly {expected_csv}")
if spec_path != repo / "resources/lineage_spec.yaml":
    raise SystemExit("refusing: spec must be the repo lineage specification")
if vocab_path != repo / "resources/annotation_vocabulary.yaml":
    raise SystemExit("refusing: vocabulary must be the repo annotation vocabulary")
v2_root = Path(args.v2_root)
if not v2_root.is_absolute():
    v2_root = repo / v2_root
v2_root = v2_root.resolve()
expected_v2_root = repo / "data/interim/migrated/xmens_annotation_v2_20260811"
if v2_root != expected_v2_root:
    raise SystemExit(f"refusing: v2 root must be exactly {expected_v2_root}")
v2_path = v2_root / f"reannotated_migrated/{args.species}_{args.lineage}.h5ad"
producer_path = repo / "workflow/scripts/apply_annotation.py"
latent_path = repo / f"data/interim/latent_{args.species}_{args.lineage}.h5ad"
status_path = v2_root / "manifests/MIGRATION_STATUS.json"
ready_path = v2_root / "MIGRATION_REVIEW_READY"
rows = list(csv.DictReader(csv_path.open()))
if not rows or not {"cluster", "label"} <= set(rows[0]):
    raise SystemExit(f"{csv_path}: expected nonempty cluster,label CSV")
annotation = {}
for row in rows:
    cluster, label = row["cluster"].strip(), row["label"].strip()
    if not cluster or not label:
        raise SystemExit(f"{csv_path}: blank cluster or label")
    if cluster in annotation:
        raise SystemExit(f"{csv_path}: duplicate cluster {cluster!r}")
    if not cluster.isdigit() or str(int(cluster)) != cluster:
        raise SystemExit(f"{csv_path}: noncanonical numeric cluster {cluster!r}")
    annotation[cluster] = label
spec = yaml.safe_load(spec_path.read_text())
lineage = spec["species"][args.species]["lineages"][args.lineage]
expected_res = lineage.get("res", spec["defaults"]["res"])
if args.res != expected_res:
    raise SystemExit(f"resolution mismatch: requested={args.res}, spec={expected_res}")
if "res" in rows[0]:
    csv_res_values = [row.get("res", "").strip() for row in rows]
    if any(not value for value in csv_res_values):
        raise SystemExit(f"{csv_path}: blank CSV resolution value")
    csv_res = set(csv_res_values)
    if csv_res != {expected_res}:
        raise SystemExit(f"{csv_path}: CSV resolution differs from spec: {sorted(csv_res)}")
markers = lineage.get("markers", {})
if not isinstance(markers, dict) or not markers:
    raise SystemExit(f"{identity}: marker panel must be a nonempty mapping")
for marker_group, genes in markers.items():
    if (not isinstance(marker_group, str) or not marker_group.strip()
            or not isinstance(genes, list) or not genes
            or any(not isinstance(gene, str) or not gene.strip() for gene in genes)):
        raise SystemExit(f"{identity}: invalid marker panel entry {marker_group!r}")
csv_hash = sha256(csv_path)
vocab_hash = sha256(vocab_path)
spec_hash = sha256(spec_path)
producer_hash = sha256(producer_path)
latent_hash = sha256(latent_path)
v2_hash = sha256(v2_path)
status_hash = sha256(status_path)
ready_hash = sha256(ready_path)
generator_hash = sha256(Path(__file__).resolve())
observed_global = {
    "vocabulary": vocab_hash, "lineage_spec": spec_hash,
    "producer": producer_hash, "migration_status": status_hash,
    "migration_ready": ready_hash,
}
if observed_global != EXPECTED_GLOBAL_SHA256:
    raise SystemExit(f"reviewed global input drift: {observed_global}")
if csv_hash != EXPECTED_CSV_SHA256[identity]:
    raise SystemExit(f"reviewed decision CSV drift: {csv_path}")
if latent_hash != EXPECTED_LATENT_SHA256[identity]:
    raise SystemExit(f"reviewed latent drift: {latent_path}")
if v2_hash != EXPECTED_V2_SHA256[identity]:
    raise SystemExit(f"reviewed certified v2 object drift: {v2_path}")
vocabulary = yaml.safe_load(vocab_path.read_text())
valid_labels = set(vocabulary["allowed"])
valid_labels.update(
    subtype
    for subtypes in vocabulary["allowed"].values()
    for subtype in (subtypes or [])
)
invalid_labels = sorted(
    {label for label in annotation.values() if label != "EXCLUDE"} - valid_labels)
if invalid_labels:
    raise SystemExit(f"{csv_path}: labels absent from locked vocabulary: {invalid_labels}")
if ready_path.read_text() != "PASS\n":
    raise SystemExit(f"v2 readiness marker is not exact PASS: {ready_path}")
v2_status = json.loads(status_path.read_text())
if (v2_status.get("status") != "PASS"
        or v2_status.get("run_id") != "xmens_annotation_v2_20260811"
        or v2_status.get("schema_version") != vocabulary["schema_version"]):
    raise SystemExit("v2 migration status/schema is not the expected PASS contract")
certified = []
for report in v2_status.get("annotation_checks", []):
    for record in report.get("objects", []):
        record_path = Path(record["path"])
        if not record_path.is_absolute():
            record_path = repo / record_path
        if record_path.resolve() == v2_path:
            certified.append(record)
if len(certified) != 1:
    raise SystemExit(f"expected one v2 certification record for {v2_path}; found {len(certified)}")
certified = certified[0]
if (certified.get("status") != "PASS" or certified.get("checker_exit_code") != 0
        or certified.get("sha256") != v2_hash
        or certified.get("size") != v2_path.stat().st_size):
    raise SystemExit(f"certified v2 record does not match current object: {v2_path}")
stamp = datetime.date.today().isoformat()

annotation_literal = json.dumps(annotation, indent=4)
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
SPEC_PATH = {str(spec_path)!r}
LATENT_PATH = {str(latent_path)!r}
V2_PATH = {str(v2_path)!r}
PRODUCER_PATH = {str(producer_path)!r}
STATUS_PATH = {str(status_path)!r}
READY_PATH = {str(ready_path)!r}
CSV_SHA256 = {csv_hash!r}
VOCAB_SHA256 = {vocab_hash!r}
SPEC_SHA256 = {spec_hash!r}
LATENT_SHA256 = {latent_hash!r}
V2_SHA256 = {v2_hash!r}
PRODUCER_SHA256 = {producer_hash!r}
STATUS_SHA256 = {status_hash!r}
READY_SHA256 = {ready_hash!r}
GENERATOR_SHA256 = {generator_hash!r}
ANNOTATION = {annotation_literal}
MARKERS = {markers_literal}'''),
    md("## 2. Load the frozen lineage latent and verify complete cluster coverage"),
    code('''import hashlib, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import pandas as pd, scanpy as sc

REPO = Path(REPO)
def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
assert _sha(LATENT_PATH) == LATENT_SHA256, (
    "latent changed; regenerate this draft")
adata = sc.read_h5ad(LATENT_PATH)
assert RES in adata.obs, f"{RES} absent from {LATENT_PATH}"
clusters = set(adata.obs[RES].astype(str).unique())
assert clusters == set(ANNOTATION), (
    f"frozen map mismatch: missing={sorted(clusters-set(ANNOTATION))}, "
    f"extra={sorted(set(ANNOTATION)-clusters)}")
def _cluster_order(value):
    try:
        return (0, float(value), str(value))
    except ValueError:
        return (1, 0.0, str(value))
cluster_order = sorted(clusters, key=_cluster_order)
adata.obs[RES] = pd.Categorical(
    adata.obs[RES].astype(str), categories=cluster_order, ordered=True)
assert list(adata.obs[RES].cat.categories) == cluster_order
print("numeric cluster order:", cluster_order)
print(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes; {len(clusters)} clusters")'''),
    md("## 3. Mechanical marker-group scores and RAW marker dotplot by frozen cluster"),
    code('''import numpy as np
plot_sets = {name: [gene for gene in genes if gene in adata.var_names]
             for name, genes in MARKERS.items()}
plot_sets = {name: genes for name, genes in plot_sets.items() if genes}
missing_groups = sorted(set(MARKERS) - set(plot_sets))
assert not missing_groups, f"configured marker groups absent from latent: {missing_groups}"
for name, genes in MARKERS.items():
    missing_genes = [gene for gene in genes if gene not in adata.var_names]
    if missing_genes:
        print(f"marker genes absent for {name}: {missing_genes}")
# Existing 04c convention: each marker-group score is the per-cell mean
# of its present marker genes. Scores are printed for review and never
# used to alter the CSV-derived ANNOTATION map.
group_scores = pd.DataFrame(index=adata.obs_names)
for name, genes in plot_sets.items():
    group_scores[f"score_{name}"] = np.asarray(
        adata[:, genes].X.mean(axis=1)).ravel()
cluster_scores = group_scores.join(
    adata.obs[[RES]].astype(str)).groupby(RES, sort=False).mean()
cluster_scores = cluster_scores.reindex(cluster_order)
print("marker-group mean scores by numeric cluster order:")
print(cluster_scores.round(3).to_string())
sc.pl.dotplot(adata, plot_sets, groupby=RES, standard_scale="var",
              figsize=(10, max(6, 0.22 * len(clusters))), show=True)'''),
    md("## 4. Apply the committed map in memory using the shared producer"),
    code('''import sys
assert _sha(PRODUCER_PATH) == PRODUCER_SHA256, (
    "annotation producer changed; regenerate this draft")
sys.path.insert(0, str(REPO / "workflow/scripts"))
from apply_annotation import (apply_map, build_vocabulary_rollup,
                              validate_canonical_annotations)
import yaml

vocab = yaml.safe_load(Path(VOCAB_PATH).read_text())
rollup = build_vocabulary_rollup(vocab)
summary = apply_map(adata, RES, ANNOTATION, rollup=rollup, output_mode="fine")
kept = adata[~adata.obs.excluded].copy()
validate_canonical_annotations(kept, vocab, context=f"{SPECIES}/{LINEAGE} draft")
kept.obs["display_label"] = kept.obs.cell_subtype.astype(str)
v2 = sc.read_h5ad(V2_PATH, backed="r")
try:
    assert kept.shape == v2.shape, "in-memory retained shape differs from certified v2"
    assert kept.obs_names.equals(v2.obs_names), "in-memory retained IDs/order differ from certified v2"
    assert kept.var_names.equals(v2.var_names), "in-memory genes/order differ from certified v2"
    assert kept.obs.cell_type.astype(str).tolist() == v2.obs.cell_type.astype(str).tolist(), (
        "in-memory cell_type differs from certified v2")
    assert kept.obs.cell_subtype.astype(str).tolist() == v2.obs.cell_subtype.astype(str).tolist(), (
        "in-memory cell_subtype differs from certified v2")
    assert not kept.obs.cell_subtype.astype(str).str.strip().eq("").any(), "empty canonical subtype"
finally:
    if getattr(v2, "file", None) is not None:
        v2.file.close()
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
    code('''assert _sha(CSV_PATH) == CSV_SHA256, "decision CSV changed; regenerate this draft"
assert _sha(VOCAB_PATH) == VOCAB_SHA256, "vocabulary changed; regenerate this draft"
assert _sha(SPEC_PATH) == SPEC_SHA256, "lineage spec changed; regenerate this draft"
assert _sha(LATENT_PATH) == LATENT_SHA256, "latent changed; regenerate this draft"
assert _sha(PRODUCER_PATH) == PRODUCER_SHA256, "annotation producer changed; regenerate this draft"
assert _sha(V2_PATH) == V2_SHA256, "certified v2 lineage changed; regenerate this draft"
assert _sha(STATUS_PATH) == STATUS_SHA256, "v2 status changed; regenerate this draft"
assert _sha(READY_PATH) == READY_SHA256, "v2 readiness changed; regenerate this draft"
now = pd.read_csv(CSV_PATH, dtype=str)
assert not now[["cluster", "label"]].isna().any().any(), "blank CSV value"
assert not now.cluster.str.strip().duplicated().any(), "duplicate CSV cluster"
now_map = dict(zip(now.cluster.str.strip(), now.label.str.strip()))
assert now_map == ANNOTATION, "embedded map differs from committed CSV"
print(now.to_string(index=False))
print("PASS: CSV and vocabulary hashes match; embedded map is exact.")'''),
]

notebook = new_notebook(cells=cells)
notebook.metadata["kernelspec"] = {
    "display_name": "Python 3", "language": "python", "name": "python3"}
notebook.metadata["language_info"] = {"name": "python"}
notebook.metadata["xmens"] = {
    "status": "DRAFT", "schema_version": vocabulary["schema_version"],
    "decision_csv_sha256": csv_hash, "vocabulary_sha256": vocab_hash,
    "lineage_spec_sha256": spec_hash, "latent_sha256": latent_hash,
    "producer_sha256": producer_hash, "v2_reannotated_sha256": v2_hash,
    "migration_status_sha256": status_hash, "migration_ready_sha256": ready_hash,
    "generator_sha256": generator_hash,
}
for cell_index, cell in enumerate(notebook.cells):
    if cell.cell_type == "code":
        ast.parse(cell.source, filename=f"generated-cell-{cell_index}", mode="exec")
nbf.validate(notebook)
out = Path(args.out)
if not out.is_absolute():
    out = repo / out
public_dir = repo / "notebooks/public"
reject_symlink_components(public_dir)
if public_dir.resolve(strict=True) != public_dir:
    raise SystemExit(f"refusing noncanonical public directory: {public_dir}")
expected_out = public_dir / (
    f"DRAFT_v2_04c_annotate_{args.species}_{args.lineage}_denovo.ipynb")
if out != expected_out:
    raise SystemExit(f"refusing: output must be exactly {expected_out}")
payload = nbf.writes(notebook, version=4).encode("utf-8")
directory_fd = os.open(
    public_dir, os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0))
raw_fd = None
owned_identity = None
try:
    raw_fd = os.open(
        out.name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o644,
        dir_fd=directory_fd,
    )
    owned_stat = os.fstat(raw_fd)
    if not stat.S_ISREG(owned_stat.st_mode) or owned_stat.st_nlink != 1:
        raise RuntimeError(f"unexpected newly created output identity: {out}")
    owned_identity = (owned_stat.st_dev, owned_stat.st_ino)
    with os.fdopen(raw_fd, "wb") as handle:
        raw_fd = None
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.fsync(directory_fd)
except BaseException:
    if raw_fd is not None:
        os.close(raw_fd)
    # Fail closed: never unlink a public pathname in an error handler. O_EXCL
    # proved the name was absent when this process created it; if a write is
    # interrupted, retain that exact partial artifact for explicit inspection.
    if owned_identity is not None:
        os.fsync(directory_fd)
    raise
finally:
    os.close(directory_fd)
print(f"compiled DRAFT -> {out} ({len(cells)} cells)")
