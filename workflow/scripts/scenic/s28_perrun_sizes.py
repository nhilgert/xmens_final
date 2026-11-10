#!/usr/bin/env python
"""Regulon size within single runs, both species (Fig. S28 panel D reference lines).

For each per-seed cisTarget output (regulons_seed<S>.pkl) of a species: the number of regulons and the number of
targets each has. Writes one JSON with the per-TF median regulon size across single runs. Read-only.
Port of perrun_regulon_sizes.py (Claude Science artifact d27de97f, 2 Nov 2026). Run in the pyscenic environment.
"""
import argparse
import glob
import json
import pickle

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lane", required=True, help="scenic_seeds folder holding <species>/<regulon-dir>")
    ap.add_argument("--regulon-dir", default="regulons_sketch15k")
    ap.add_argument("--n-expected", type=int, default=100)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()
    out = {}
    for species in ("mouse", "human"):
        pkls = sorted(glob.glob(f"{a.lane}/{species}/{a.regulon_dir}/*.pkl"))
        assert len(pkls) == a.n_expected, f"{species}: {len(pkls)} pickles, expected {a.n_expected}"
        n_reg, sizes, per_tf = [], [], {}
        for p in pkls:
            with open(p, "rb") as fh:
                regs = pickle.load(fh)
            seen, k = set(), 0
            for r in regs:
                if r.name in seen:
                    continue
                seen.add(r.name)
                k += 1
                s = len(r.gene2weight)
                sizes.append(s)
                per_tf.setdefault(r.name, []).append(s)
            n_reg.append(k)
        sizes = np.array(sizes)
        med_per_tf = {k: float(np.median(v)) for k, v in per_tf.items()}
        out[species] = {
            "n_pickles": len(pkls),
            "regulons_per_run": {"min": int(min(n_reg)), "median": float(np.median(n_reg)), "max": int(max(n_reg))},
            "target_count_per_regulon_all_runs": {q: float(np.quantile(sizes, q / 100)) for q in (10, 25, 50, 75, 90)},
            "target_count_mean": float(sizes.mean()),
            "target_count_max": int(sizes.max()),
            "n_distinct_tfs_any_run": len(per_tf),
            "median_of_per_tf_median_size": float(np.median(list(med_per_tf.values()))),
            "per_tf_median_size": med_per_tf,
        }
        print(species, json.dumps({k: v for k, v in out[species].items() if k != "per_tf_median_size"}, indent=1),
              flush=True)
    with open(a.out_json, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote", a.out_json, flush=True)


if __name__ == "__main__":
    main()
