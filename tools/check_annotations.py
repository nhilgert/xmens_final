#!/usr/bin/env python3
"""Pre-flight: every committed choices CSV must be consistent with the spec.

Three failure modes, all of which have silently shipped at least once:

  1. RESOLUTION MISMATCH -- the CSV's cluster IDs were decided at a different
     resolution than the spec declares. apply_annotation.py only raises when a
     cluster is MISSING from the CSV, so a finer-resolution CSV (a superset of
     IDs) passes and mislabels everything.
  2. ORPHANED LABEL -- a label in the CSV that build_rollup doesn't know. It
     becomes its own cell_type AND leaves fine_refined empty.
  3. Both exit 0 and produce plausible-looking output.

Usage:  python3 tools/check_annotations.py [--species mouse]
Exit 1 on any failure, so it can gate a run.
"""
import argparse, pathlib, sys
import pandas as pd, yaml

sys.path.insert(0, "workflow/scripts")
from apply_annotation import build_rollup

COARSE_TARGETS = {
    "epithelial cells", "immune cells", "vascular cells", "perivascular cells",
    "smooth muscle cells", "perimetrium cells", "stromal",
    "decidual cells", "endometrial fibroblasts",
}
# NB do NOT add subtype names here. This set exempts labels that are legitimately
# cell-type level (incl. cross-lineage throwbacks). Adding a subtype -- e.g.
# "ciliated epithelium" -- silently disables the orphan check for that label,
# which is the failure this script exists to catch.

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="*", default=["mouse", "human"])
    ap.add_argument("--spec", default="resources/lineage_spec.yaml")
    a = ap.parse_args()
    spec = yaml.safe_load(open(a.spec))
    fails = []

    for sp in a.species:
        print(f"=== {sp} ===")
        atlas_res = spec["atlas"][sp]["res"]
        pairs = [("(atlas)", atlas_res, f"resources/annotations/{sp}_lineage.csv",
                  f"data/interim/latent_{sp}.h5ad", None)]
        for lin, cfg in spec["species"][sp]["lineages"].items():
            pairs.append((lin, cfg["res"], f"resources/annotations/{sp}_{lin}.csv",
                          f"data/interim/latent_{sp}_{lin}.h5ad", cfg))

        for lin, res, csv, h5, cfg in pairs:
            cp, hp = pathlib.Path(csv), pathlib.Path(h5)
            if not cp.exists():
                print(f"  {lin:14} SKIP (no CSV)"); continue
            d = pd.read_csv(cp, dtype=str)
            if not {"cluster", "label"} <= set(d.columns):
                fails.append(f"{sp}_{lin}: CSV lacks cluster/label"); continue
            rows = set(d["cluster"].str.strip())
            labs = set(d["label"].str.strip()) - {"EXCLUDE"}

            # (1) resolution: exact cluster-set equality, not just counts
            if hp.exists():
                import anndata as ad
                obs = ad.read_h5ad(hp, backed="r").obs
                if res not in obs:
                    fails.append(f"{sp}_{lin}: spec res {res} not in {hp.name}")
                    print(f"  {lin:14} FAIL res {res} absent from object")
                else:
                    have = set(obs[res].astype(str).unique())
                    if have != rows:
                        fails.append(f"{sp}_{lin}: cluster sets differ at {res}")
                        print(f"  {lin:14} FAIL {res}: object={len(have)} csv={len(rows)}")
                        print(f"                 missing from csv: {sorted(have-rows)[:6]}")
                        print(f"                 extra in csv:     {sorted(rows-have)[:6]}")
                    else:
                        print(f"  {lin:14} ok  {res} ({len(rows)} clusters)")
            else:
                print(f"  {lin:14} ??  {res} ({len(rows)} clusters) -- no object to verify")

            # (2) orphaned labels
            if cfg is not None:
                r = build_rollup(cfg.get("fine", {}), cfg["coarse"])
                orphan = sorted(labs - set(r) - COARSE_TARGETS - {cfg["coarse"]})
                if orphan:
                    fails.append(f"{sp}_{lin}: orphaned labels {orphan}")
                    print(f"                 FAIL orphaned labels: {orphan}")
                    print(f"                 -> own cell_type + empty fine_refined")
        print()

    if fails:
        print("FAILURES:"); [print(f"  - {f}") for f in fails]
        return 1
    print("ALL OK")
    return 0

if __name__ == "__main__":
    sys.exit(main())
