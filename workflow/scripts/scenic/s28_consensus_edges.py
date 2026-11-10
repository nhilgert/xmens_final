#!/usr/bin/env python
"""Long edge table of a SCENIC consensus (Fig. S28 input): one row per TF-target pair.

Reads consensus_regulons.pkl written by build_consensus_aucell.py: a list of ctxcore Regulon objects, one per kept TF.
gene2weight is the AUCell weight of the pair; gene2occurrence is the number of runs (of the seed ensemble) that
contain it (n_seeds). Writes regulon, tf, target, n_seeds, weight; weight rounded to 4 decimals, as the table of record.
Run in the pyscenic environment (unpickling needs ctxcore). Port of the one-off flattening of 10 Sep 2026.
"""
import argparse
import pickle

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.pkl, "rb") as fh:
        regs = pickle.load(fh)
    rows, seen = [], set()
    for r in regs:
        if r.name in seen:
            raise ValueError(f"duplicate regulon {r.name}")
        seen.add(r.name)
        tf = r.name[:-3] if r.name.endswith("(+)") else r.name
        for gene, weight in r.gene2weight.items():
            rows.append((r.name, tf, gene, int(r.gene2occurrence[gene]), round(float(weight), 4)))
    df = pd.DataFrame(rows, columns=["regulon", "tf", "target", "n_seeds", "weight"])
    df = df.sort_values(["regulon", "target"]).reset_index(drop=True)
    df.to_csv(a.out, sep="\t", index=False)
    print(f"{len(regs)} regulons, {len(df)} pairs -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
