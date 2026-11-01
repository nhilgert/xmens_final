#!/usr/bin/env python
"""Write pickled pySCENIC regulons to YAML for the pySCENIC command line (rule scenic_aucell_cli, 24_scenic.smk).

`pyscenic aucell` reads regulons from YAML (or GMT), not from the pickle that build_consensus_aucell.py writes. This
mirrors the 26 Oct 2026 step that produced consensus_regulons.yaml in the *_pyscenic_cli folders:
pyscenic.utils.save_to_yaml,
then a reload check that every regulon's name, targets and weights round-trip exactly. Writes atomically.
"""
import argparse, os, pickle
from pyscenic.utils import save_to_yaml, load_from_yaml


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pkl", required=True, help="consensus_regulons.pkl from build_consensus_aucell.py")
    ap.add_argument("--yaml", required=True, help="output YAML for `pyscenic aucell`")
    a = ap.parse_args()
    with open(a.pkl, "rb") as fh:
        regulons = pickle.load(fh)
    tmp = f"{a.yaml}.partial.{os.getpid()}"
    save_to_yaml(regulons, tmp)
    back = load_from_yaml(tmp)
    if len(back) != len(regulons):
        raise SystemExit(f"reload: {len(back)} regulons, wrote {len(regulons)}")
    for r, b in zip(regulons, back):
        if r.name != b.name or dict(r.gene2weight) != dict(b.gene2weight):
            raise SystemExit(f"reload differs for {r.name}")
    os.replace(tmp, a.yaml)
    print(f"wrote {len(regulons)} regulons -> {a.yaml} (reload identical)")


if __name__ == "__main__":
    main()
