#!/usr/bin/env python3
"""Apply a frozen cluster->label choices CSV to an object.

THE single implementation of the apply step. Three ways in, one code path:

  * snakemake  `script:` directive   (the pipeline -- rules 04/04c)
  * CLI                             `python apply_annotation.py --level atlas ...`
  * import                          `from apply_annotation import apply_map`
                                    (the notebooks -- so they cannot drift)

Levels:
  atlas : latent_{sp}.h5ad + {sp}_lineage.csv
          -> annotated_{sp}.h5ad, lineages/{sp}_{slug}.h5ad, the leiden ladder
  fine  : latent_{sp}_{lin}.h5ad + {sp}_{lin}.csv
          -> reannotated/{sp}_{lin}.h5ad
"""
import sys
from pathlib import Path
import pandas as pd, numpy as np

__all__ = ["load_choices", "apply_map", "write_atlas", "write_fine", "build_rollup"]


# ----------------------------------------------------------------- choices ---
def load_choices(csv) -> dict:
    """cluster -> label, from the committed two-column record."""
    df = pd.read_csv(csv, dtype=str)
    need = {"cluster", "label"}
    if not need <= set(df.columns):
        raise ValueError(f"{csv}: needs columns {need}, has {list(df.columns)}")
    return dict(zip(df["cluster"].str.strip(), df["label"].str.strip()))


def build_rollup(fine_cfg: dict, coarse: str) -> dict:
    """fine label -> coarse label, from a spec `fine:` block."""
    if not fine_cfg:
        return {}
    if fine_cfg.get("mode") == "markers":
        return {lab: coarse for lab in fine_cfg.get("from_groups", {})}
    # Any mode declaring an explicit `labels` list rolls up the same way: the
    # rollup is a function of the LABELS, not of where they were sourced from.
    # (Gating on mode == "csv" silently returned {} when mouse stromal moved to
    # mode: h5ad -- which wrote fine labels into coarse_refined AND emptied
    # fine_refined, with the pipeline still exiting 0.)
    if fine_cfg.get("labels"):
        out = {}
        for lab in fine_cfg.get("labels", []):
            if "rollup_all" in fine_cfg:
                out[lab] = fine_cfg["rollup_all"]; continue
            for pref, tgt in fine_cfg.get("rollup_prefix", {}).items():
                if str(lab).startswith(pref):
                    out[lab] = tgt; break
            else:
                out[lab] = coarse
        return out
    return {}


# ------------------------------------------------------------------- apply ---
def apply_map(adata, res: str, choices: dict, rollup: dict | None = None,
              label_col: str = "label", strict: bool = True):
    """Write label / coarse_refined / fine_refined / excluded onto adata.obs.

    Returns a per-cluster summary frame. Pure: no I/O, no plotting -- this is
    what the pipeline and the notebooks share.
    """
    rollup = rollup or {}
    if res not in adata.obs:
        raise KeyError(f"{res} not in obs; have "
                       f"{[c for c in adata.obs if c.startswith('leiden')]}")
    clusters = set(adata.obs[res].astype(str).unique())
    missing = clusters - set(choices)
    extra = set(choices) - clusters
    if missing and strict:
        raise ValueError(f"clusters with no row in the choices CSV: "
                         f"{sorted(missing, key=lambda c: (len(c), c))}")
    if extra:
        print(f"  note: CSV has rows for {len(extra)} absent cluster(s) (ignored)")

    key = adata.obs[res].astype(str)
    adata.obs[label_col] = key.map(choices).astype(str)
    adata.obs["excluded"] = adata.obs[label_col].eq("EXCLUDE").values
    adata.obs["coarse_refined"] = np.where(
        adata.obs.excluded, "EXCLUDE",
        adata.obs[label_col].map(lambda l: rollup.get(l, l)))
    adata.obs["fine_refined"] = np.where(
        adata.obs.excluded, "",
        np.where(adata.obs[label_col].isin(rollup), adata.obs[label_col], ""))

    summary = pd.DataFrame({
        "n": key.value_counts(),
        "label": pd.Series(choices).reindex(sorted(clusters, key=lambda c: (len(c), c))),
    }).dropna(subset=["label"])
    summary["excluded"] = summary["label"].eq("EXCLUDE")
    summary["coarse"] = summary["label"].map(lambda l: rollup.get(l, l))
    summary.index.name = "cluster"
    n_ex = int(adata.obs.excluded.sum())
    print(f"  applied {len(choices)} choices at {res}: "
          f"{adata.n_obs - n_ex:,} kept, {n_ex:,} EXCLUDE")
    return summary


# ------------------------------------------------------------------- write ---
def write_atlas(adata, out_h5ad, lineage_dir, ladder_dir, slug: dict, species: str,
                label_col: str = "coarse_refined"):
    """Annotated object + the lineage subsets 04b consumes + the leiden ladder."""
    out_h5ad = Path(out_h5ad); out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    lineage_dir, ladder_dir = Path(lineage_dir), Path(ladder_dir)
    lineage_dir.mkdir(parents=True, exist_ok=True)
    ladder_dir.mkdir(parents=True, exist_ok=True)

    labels = set(adata.obs.loc[~adata.obs.excluded, label_col])
    unknown = labels - set(slug)
    if unknown:
        raise ValueError(f"labels in the CSV with no slug in the spec: {sorted(unknown)}. "
                         f"Add them to atlas.<species>.labels or fix the CSV.")

    ann = adata[~adata.obs.excluded].copy()
    for c in [c for c in ann.obs if c.startswith(("score_", "ref_"))] + ["label"]:
        if c in ann.obs: del ann.obs[c]
    ann.obs[label_col] = pd.Categorical(
        ann.obs[label_col], categories=[l for l in slug if l in labels])
    ann.write_h5ad(out_h5ad)
    print(f"  wrote {out_h5ad}  ({ann.n_obs:,} cells)")

    rows = []
    for lin in ann.obs[label_col].cat.categories:
        sub = ann[ann.obs[label_col] == lin].copy()
        if sub.n_obs == 0: continue
        p = lineage_dir / f"{species}_{slug[lin]}.h5ad"
        sub.write_h5ad(p); rows.append({"lineage": lin, "n": sub.n_obs, "file": p.name})
    print(pd.DataFrame(rows).set_index("lineage").to_string() if rows else "  (no subsets)")

    base = pd.DataFrame({"cell": adata.obs_names})
    for c in ("sample", "cell_barcode", "timepoint"):
        if c in adata.obs: base[c] = adata.obs[c].astype(str).values
    res_cols = sorted([c for c in adata.obs if c.startswith("leiden_scVI_")],
                      key=lambda c: float(c.rsplit("_", 1)[1]))
    for col in res_cols:
        base.assign(cluster=adata.obs[col].astype(str).values).to_csv(
            ladder_dir / f"leiden_res{col.rsplit('_', 1)[1]}.csv", index=False)
    wide = base.copy()
    for col in res_cols: wide[col] = adata.obs[col].astype(str).values
    wide[label_col] = adata.obs[label_col].astype(str).values
    wide.to_csv(ladder_dir / "leiden_all_resolutions.csv", index=False)
    print(f"  wrote {len(res_cols)} ladder files -> {ladder_dir}")
    return ann


def write_fine(adata, out_h5ad):
    """Cleaned per-lineage object for 04d_merge (EXCLUDE dropped)."""
    out_h5ad = Path(out_h5ad); out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    out = adata[~adata.obs.excluded].copy()
    for c in ["ref_type", "ref_lineage", "ref_fine", "label", "asn_lbl", "cluster_lbl"]:
        if c in out.obs: del out.obs[c]
    out.write_h5ad(out_h5ad)
    print(f"  wrote {out_h5ad}  ({out.n_obs:,} cells)")
    return out


# ------------------------------------------------------------------ entries --
def _run(level, species, lineage, res, in_h5ad, csv, spec_path, repo, out_h5ad):
    import scanpy as sc, yaml
    spec = yaml.safe_load(open(spec_path))
    repo = Path(repo)
    print(f"apply_annotation [{level}] {species}"
          f"{'/' + lineage if lineage else ''} | {Path(in_h5ad).name} | {res}")

    adata = sc.read_h5ad(in_h5ad)
    choices = load_choices(csv)

    if level == "atlas":
        cfg = spec["atlas"][species]
        apply_map(adata, res, choices, rollup={})
        write_atlas(adata, out_h5ad,
                    repo / "data/interim/lineages",
                    repo / "results/tables/clusters" / species,
                    slug=cfg["labels"], species=species)
    else:
        lin = spec["species"][species]["lineages"][lineage]
        rollup = build_rollup(lin.get("fine", {}), lin["coarse"])
        apply_map(adata, res, choices, rollup=rollup)
        write_fine(adata, out_h5ad)


def _from_snakemake(sm):
    p = sm.params
    _run(level=p.level, species=p.species, lineage=getattr(p, "lineage", None),
         res=p.res, in_h5ad=sm.input.h5ad, csv=sm.input.csv,
         spec_path=getattr(sm.input, "spec", None)
                   or "resources/lineage_spec.yaml",
         repo=p.repo, out_h5ad=sm.output.h5ad)


def _from_cli(argv):
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", choices=["atlas", "fine"], required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--lineage")
    ap.add_argument("--res", required=True)
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--spec", default="resources/lineage_spec.yaml")
    ap.add_argument("--repo", default=".")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    _run(a.level, a.species, a.lineage, a.res, a.in_h5ad, a.csv, a.spec, a.repo, a.out)


if "snakemake" in globals():          # script: directive
    _from_snakemake(globals()["snakemake"])   # noqa: F821
elif __name__ == "__main__":
    _from_cli(sys.argv[1:])
