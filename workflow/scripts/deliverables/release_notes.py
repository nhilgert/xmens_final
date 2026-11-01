"""Release notes of the release of record, written from computed values only (no typed counts).

04g release_notes. Reads the release's MANIFEST.json and objects, the objects of the previous release (config
releases; none for the first) and the post-cut check reports (check_release_object); writes one Markdown file. Every
number in it is counted here from those files:
  objects       nuclei and sha256 from the manifest; composition by cell_type from the all-cells objects;
  changes       per species, against the previous release: nuclei added and removed, label changes of shared nuclei
                (largest transitions), spatial QC changes of shared nuclei split into the round-2 lineages and the
                nuclei outside them (obs origin_lineage), and the obs columns that differ outside the round-2 lineages;
  embedding     reused (X_scVI and X_umap_scVI identical to the previous release for every nucleus) or re-fitted, with
                the scVI fit record; embedding outliers now, before and in both;
  fib/dec       nuclei, the obs columns that changed, and the pseudotime record (uns pt_harmonized) with its rank
                correlation to the previous release;
  checks        the status and failed gates of every post-cut check report.
"""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

snk = snakemake  # noqa: F821 (injected by Snakemake)

SPATIAL_OBS = ("spatial_qc_reason", "x_um_qc", "y_um_qc", "tissue_section")
EMBEDDING_OBSM = ("X_scVI", "X_umap_scVI")
TOP_TRANSITIONS = 5                     # label transitions listed per species
OUTSIDE = "outside the round-2 lineages"


def count(n) -> str:
    return f"{int(n):,}"


def same_values(a: pd.Series, b: pd.Series) -> np.ndarray:
    """Element-wise equality, NaN matching NaN; strings compared as strings."""
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        x, y = a.to_numpy(float), b.to_numpy(float)
        return (x == y) | (np.isnan(x) & np.isnan(y))
    return (a.astype(object).where(a.notna(), "<NA>").astype(str).to_numpy()
            == b.astype(object).where(b.notna(), "<NA>").astype(str).to_numpy())


def transitions(before: pd.Series, after: pd.Series, top: int = None) -> str:
    """'a -> b n; ...' for the values that differ, largest first."""
    changed = ~same_values(before, after)
    pairs = pd.Series(list(zip(before[changed].astype(str), after[changed].astype(str)))).value_counts()
    if top:
        pairs = pairs.head(top)
    return "; ".join(f"{a} -> {b} {count(n)}" for (a, b), n in pairs.items()) or "none"


def by_label(values: pd.Series) -> str:
    return "; ".join(f"{label} {count(n)}" for label, n in values.astype(str).value_counts().items()) or "none"


def species_changes(new, previous, round2_lineages: list, lines: list):
    """What changed in one species' all-cells object, against the previous release."""
    shared = new.obs_names.intersection(previous.obs_names)
    added = new.obs_names.difference(previous.obs_names)
    removed = previous.obs_names.difference(new.obs_names)
    after, before = new.obs.loc[shared], previous.obs.loc[shared]
    group = after["origin_lineage"].astype(str).where(after["origin_lineage"].astype(str).isin(round2_lineages),
                                                      OUTSIDE)
    relabeled = ~same_values(before["cell_subtype"], after["cell_subtype"])
    lines += [
        f"- Nuclei: {count(new.n_obs)} ({count(len(shared))} shared with the previous release, {count(len(added))} "
        f"added, {count(len(removed))} removed).",
        f"- Added by label: {by_label(new.obs.loc[added, 'cell_subtype'])}.",
        f"- Removed by label: {by_label(previous.obs.loc[removed, 'cell_subtype'])}.",
        f"- Label changes of shared nuclei: {count(relabeled.sum())}; largest: "
        f"{transitions(before['cell_subtype'], after['cell_subtype'], TOP_TRANSITIONS)}.",
    ]

    # Spatial QC changes, split by lineage group, so outside-lineage counts are never mixed with the lineages.
    for name, members in ((OUTSIDE, group == OUTSIDE), ("in the round-2 lineages", group != OUTSIDE)):
        reason_before = before.loc[members.to_numpy(), "spatial_qc_reason"]
        reason_after = after.loc[members.to_numpy(), "spatial_qc_reason"]
        n_changed = int((~same_values(reason_before, reason_after)).sum())
        lines.append(f"- Spatial QC of shared nuclei {name}: {count(n_changed)} change "
                     f"({transitions(reason_before, reason_after)}).")

    outside = group == OUTSIDE
    skipped = set(SPATIAL_OBS) | {"embedding_outlier"}
    differing = {c: int((~same_values(before.loc[outside.to_numpy(), c], after.loc[outside.to_numpy(), c])).sum())
                 for c in after.columns if c in before.columns and c not in skipped}
    differing = {c: n for c, n in differing.items() if n}
    lines.append(f"- Obs columns that differ for shared nuclei {OUTSIDE} (spatial QC and outlier flags aside): "
                 f"{differing or 'none'}.")

    # Embedding: reused when identical for every nucleus, re-fitted otherwise.
    same_nuclei = len(added) == 0 and len(removed) == 0
    rows = previous.obs_names.get_indexer(new.obs_names)
    reused = same_nuclei and all(
        np.array_equal(np.asarray(new.obsm[key]), np.asarray(previous.obsm[key])[rows]) for key in EMBEDDING_OBSM)
    fit = dict(new.uns.get("scvi_fit", {}))
    fit_text = ", ".join(f"{k} {v}" for k, v in fit.items()) or "no scvi_fit record"
    flags_new = new.obs["embedding_outlier"].astype(bool)
    flags_previous = previous.obs["embedding_outlier"].astype(bool)
    in_both = int((flags_new.loc[shared] & flags_previous.loc[shared]).sum())
    embedding = "reused from the previous release (identical for every nucleus)" if reused else "re-fitted"
    lines += [
        f"- Whole-tissue embedding: {embedding}; scVI fit record: {fit_text}.",
        f"- Embedding outliers: {count(flags_new.sum())} (previous release {count(flags_previous.sum())}; "
        f"{count(in_both)} flagged in both).",
    ]


def fibdec_changes(new, previous, lines: list):
    shared = new.obs_names.intersection(previous.obs_names)
    after, before = new.obs.loc[shared], previous.obs.loc[shared]
    differing = {c: int((~same_values(before[c], after[c])).sum())
                 for c in after.columns if c in before.columns}
    differing = {c: n for c, n in differing.items() if n}
    lines.append(f"- Nuclei: {count(new.n_obs)} ({count(len(shared))} shared with the previous release); obs columns "
                 f"that changed for shared nuclei: {differing or 'none'}.")
    record = dict(new.uns.get("pt_harmonized", {}))
    if record:
        rank_now = after["pt_harmonized_rank"].astype(float)
        rank_before = before["pt_harmonized_rank"].astype(float)
        if record.get("table"):
            source = f"run table {record['table']} (sha256 {str(record.get('table_sha256', ''))[:16]}...)"
        else:
            source = "the release's own run (no run table recorded)"
        lines.append(f"- Pseudotime (obs pt_harmonized_rank, uns pt_harmonized): {source}, root {record.get('root')}; "
                     f"Spearman with the previous release {rank_now.corr(rank_before, method='spearman'):.4f}.")


def main():
    manifest = json.loads(Path(snk.input.manifest).read_text())
    release, version = str(snk.params.release), str(snk.params.version)
    previous_release = snk.params.previous_release or None
    round2 = dict(snk.params.round2_lineages)
    objects = {Path(path).stem: path for path in snk.input.objects}
    previous_objects = {Path(path).stem: path for path in snk.input.previous}

    lines = [
        f"# Release {release} ({version}): release notes",
        "",
        f"`data/deliverables/{release}`, MANIFEST.json created {manifest.get('created')}, vocabulary "
        f"`{manifest.get('vocabulary_schema')}`. Previous release: {previous_release or 'none'}. Every number below is "
        "counted by workflow/scripts/deliverables/release_notes.py from the files it names; none is typed.",
        "",
        f"MANIFEST note: {manifest.get('note')}",
        "",
        "| object | nuclei | sha256 |",
        "|---|---|---|",
    ]
    for name, entry in manifest["objects"].items():
        lines.append(f"| {entry['file']} | {count(entry['n_obs'])} | `{entry['sha256'][:16]}...` |")

    all_cells = {name: ad.read_h5ad(path, backed="r") for name, path in objects.items() if "_All_Cells" in name}
    composition = pd.DataFrame({name: a.obs["cell_type"].astype(str).value_counts() for name, a in all_cells.items()})
    composition = composition.fillna(0).astype(int)
    lines += ["", "## Composition (cell_type)", "", "| cell_type | " + " | ".join(composition.columns) + " |",
              "|---|" + "---|" * len(composition.columns)]
    for cell_type, row in composition.iterrows():
        lines.append(f"| {cell_type} | " + " | ".join(count(n) for n in row) + " |")
    lines.append("| total | " + " | ".join(count(n) for n in composition.sum()) + " |")

    if previous_release:
        lines += ["", f"## What changed from {previous_release}"]
        for name, new in all_cells.items():
            species = str(snk.params.species_of_tag[name.split("_")[1]])
            previous = ad.read_h5ad(previous_objects[name], backed="r")
            lines += ["", f"### {name} ({species}; round-2 lineages {round2.get(species, [])})", ""]
            species_changes(new, previous, list(round2.get(species, [])), lines)
        for name, path in objects.items():
            if "_Only_Fib_Dec_Cells" not in name:
                continue
            new = ad.read_h5ad(path, backed="r")
            previous = ad.read_h5ad(previous_objects[name], backed="r")
            lines += ["", f"### {name}", ""]
            fibdec_changes(new, previous, lines)

    lines += ["", "## Post-cut checks (tools/check_annotations.py)", "", "| object | status | gates not passed |",
              "|---|---|---|"]
    for path in snk.input.checks:
        report = json.loads(Path(path).read_text())
        failed = [g["check"] for g in report.get("count_gates", []) if g["status"] == "SUSPECT"]
        failed += list(report.get("schema_failures", []))
        lines.append(f"| {Path(path).stem} | {report['status']} | {', '.join(failed) or 'none'} |")

    Path(snk.output.notes).parent.mkdir(parents=True, exist_ok=True)
    Path(snk.output.notes).write_text("\n".join(lines) + "\n")
    Path(snk.log[0]).write_text(f"wrote {snk.output.notes}: {len(lines)} lines\n")


main()
