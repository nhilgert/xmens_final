#!/usr/bin/env python3
"""
build_consensus.py - deterministically build the cross-paper decidualization consensus
("best-hits") from the four per-paper in-vitro signature tables.

WHAT THIS IS
------------
The consensus gene list is NOT hand-curated: membership falls out of a fixed rule applied
to the per-paper differential-expression signatures. This script encodes that rule so the
list is reproducible and auditable - re-run it on the same inputs and you always get the
same genes. Ship it with the paper alongside the four input tables.

INPUTS  (per-paper signature tables; the supp-derived intermediates)
------
Each is a CSV with at least a `gene` column and a `direction` column ('up'/'down'); a gene
may appear on several rows (timepoints / arms / contrasts). Direction is assumed already
derived from each paper's supplementary statistics (logFC sign at that paper's thresholds) -
see generate_panel.py / verify_panel.py for the supp -> direction derivation and its specs.

  Aghajanova 2011   E2P4-vs-E2, 6h/48h/14d          (mature vote)
  Rytkonen 2019     8-Br-cAMP+MPA, ESF->3d/8d       (mature vote)
  Stadtmauer 2022   PGE2+MPA & cAMP+MPA, 6d          (mature vote)
  Lucas 2020        S2-vs-S3 onset (D2, early)       (annotation only, NOT a vote)

THE RULE  (all editorial choices are explicit and live here, not in typed rows)
--------
  1. Symbol harmonization: uppercase + a small HGNC alias map (ALIASES).
  2. Per-paper collapse: a gene's rows within one paper collapse to a single direction -
     'up' if all its rows are up, 'down' if all down, else 'biphasic'.
  3. Mature vote: consider only MATURE papers (Aghajanova, Rytkonen, Stadtmauer). A gene is
     a best-hit iff >= MIN_MATURE of them give a CLEAN (up/down) call and those calls AGREE.
     Rationale: Lucas' D2 onset is the earliest point in the pool and reads the opposite sign
     from the mature program for late-induced genes, so it annotates rather than votes.
  4. Tier: 'core_KofK' if all mature papers call it, else 'supported_nofK'. This reflects
     COVERAGE (a gene can be robust yet miss the smallest panel), not confidence.
  5. Lucas onset is attached as an annotation: its collapsed direction and whether it agrees
     with the mature call (the disagreements are the late-induced/early-repressed genes).

Because generation is deterministic, verification is `--verify`: rebuild and diff against a
reference best-hits CSV; any mismatch is a surfaced, rule-traceable discrepancy.

USAGE
-----
  python build_consensus.py --supp-dir DIR --out best_hits.csv
        Build the best-hits list from the four tables in DIR.

  python build_consensus.py --supp-dir DIR --out best_hits.csv --consensus-out consensus.csv
        Also write the full per-gene consensus table (every gene, with the vote tally).

  python build_consensus.py --supp-dir DIR --verify decidualization_best_hits.csv
        Rebuild and assert it reproduces the reference list exactly (== verification).

  python build_consensus.py --selftest
        Synthetic per-paper tables -> assert the pooling rule yields the expected calls.
"""

import argparse
import os
import sys
import tempfile
from collections import defaultdict

import pandas as pd

# --------------------------------------------------------------------------------------
# CONFIG  (the editorial layer - explicit and version-controlled)
# --------------------------------------------------------------------------------------

# HGNC alias / previous-symbol map for cross-paper matching (incl. Excel-corruption fixes).
ALIASES = {
    "KIAA1199": "CEMIP", "GILZ": "TSC22D3", "SNAIL": "SNAI1", "KI67": "MKI67",
    "P65": "RELA", "NIK": "MAP3K14", "TAK1": "MAP3K7", "NRF2": "NFE2L2",
    "MARC1": "MTARC1", "MARC2": "MTARC2",
}

# Per-paper signature tables: (label, filename, role). role in {"mature", "annotation"}.
PAPER_SPEC = [
    ("Aghajanova", "aghajanova2011_invitro_decidualization_nonendo.csv", "mature"),
    ("Rytkonen",   "rytkonen2019_complete_signature.csv",                "mature"),
    ("Stadtmauer", "stadtmauer2022_hits.csv",                            "mature"),
    ("Lucas",      "lucas2020_S2vsS3_decidualization_onset.csv",         "annotation"),
]

GENE_COL = "gene"
DIR_COL = "direction"
MIN_MATURE = 2           # a best-hit needs >= this many agreeing clean mature calls
ANNOTATION_PAPER = "Lucas"

MATURE = [lab for lab, _, role in PAPER_SPEC if role == "mature"]


# --------------------------------------------------------------------------------------
# core
# --------------------------------------------------------------------------------------

def norm(sym):
    s = str(sym).upper().strip()
    return ALIASES.get(s, s)


def collapse_paper(path):
    """One paper's table -> {gene: 'up'|'down'|'biphasic'} (collapsed over all its rows)."""
    df = pd.read_csv(path)
    seen = defaultdict(set)
    for g, d in zip(df[GENE_COL], df[DIR_COL]):
        d = str(d).strip().lower()
        if d in ("up", "down"):
            seen[norm(g)].add(d)
    out = {}
    for g, ds in seen.items():
        out[g] = "up" if ds == {"up"} else "down" if ds == {"down"} else "biphasic"
    return out


def build(supp_dir):
    """Return (best_hits_df, consensus_df). Deterministic function of the input tables."""
    per_paper = {}
    for lab, fn, _role in PAPER_SPEC:
        path = os.path.join(supp_dir, fn)
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing input table for {lab}: {path}")
        per_paper[lab] = collapse_paper(path)

    # gene -> {paper: collapsed_dir}
    gd = defaultdict(dict)
    for lab in per_paper:
        for g, d in per_paper[lab].items():
            gd[g][lab] = d

    best_rows, cons_rows = [], []
    for g in sorted(gd):
        calls = gd[g]
        mature_calls = {p: calls[p] for p in MATURE if p in calls}
        clean = {p: d for p, d in mature_calls.items() if d in ("up", "down")}
        n_bip = sum(1 for d in mature_calls.values() if d == "biphasic")
        n_up = sum(1 for d in clean.values() if d == "up")
        n_dn = sum(1 for d in clean.values() if d == "down")

        onset = calls.get(ANNOTATION_PAPER, "absent")

        # full consensus row (every gene)
        if len(clean) == 0:
            cdir = "biphasic" if n_bip else "n.s."
        elif n_up and n_dn:
            cdir = "CONFLICT"
        else:
            cdir = "up" if n_up else "down"
        cons_rows.append(dict(gene=g, n_mature_clean=len(clean), n_up=n_up, n_down=n_dn,
                              n_mature_biphasic=n_bip, consensus_direction=cdir,
                              mature_papers="; ".join(sorted(clean.keys())),
                              lucas_onset_dir=onset))

        # best-hit rule: >= MIN_MATURE clean mature calls, all agreeing
        if len(clean) >= MIN_MATURE and n_up * n_dn == 0 and len(clean) > 0:
            cdir = "up" if n_up else "down"
            agrees = (onset == cdir) if onset in ("up", "down") else ""
            tier = (f"core_{len(clean)}of{len(MATURE)}" if len(clean) == len(MATURE)
                    else f"supported_{len(clean)}of{len(MATURE)}")
            best_rows.append(dict(gene=g, consensus_direction=cdir, tier=tier,
                                  n_mature=len(clean), mature_papers="; ".join(sorted(clean.keys())),
                                  mature_biphasic=n_bip, lucas_onset_dir=onset,
                                  onset_agrees_mature=agrees))

    best = pd.DataFrame(best_rows).sort_values(
        ["tier", "consensus_direction", "gene"],
        key=lambda s: s.map({"core_3of3": 0}).fillna(1) if s.name == "tier" else s
    ).reset_index(drop=True)
    cons = pd.DataFrame(cons_rows)
    return best, cons


def verify(best, ref_csv):
    ref = pd.read_csv(ref_csv)
    got = dict(zip(best.gene, best.consensus_direction))
    exp = dict(zip(ref.gene, ref.consensus_direction))
    gs, es = set(got), set(exp)
    only_new = sorted(gs - es)
    only_ref = sorted(es - gs)
    dir_mm = sorted(g for g in gs & es if got[g] != exp[g])
    ok = not (only_new or only_ref or dir_mm)
    print(f"\n=== VERIFY against {os.path.basename(ref_csv)} ===")
    print(f"rebuilt: {len(got)}   reference: {len(exp)}   identical set: {gs == es}")
    if only_new: print(f"  only in rebuilt ({len(only_new)}): {only_new}")
    if only_ref: print(f"  only in reference ({len(only_ref)}): {only_ref}")
    if dir_mm:   print(f"  direction mismatches ({len(dir_mm)}): {dir_mm}")
    print("VERIFY", "PASSED" if ok else "FAILED")
    return ok


# --------------------------------------------------------------------------------------
# self-test: the pooling rule
# --------------------------------------------------------------------------------------

def selftest():
    d = tempfile.mkdtemp()

    def w(fn, rows):  # rows: list of (gene, direction)
        pd.DataFrame(rows, columns=["gene", "direction"]).to_csv(os.path.join(d, fn), index=False)

    # HIT_CORE : up in all 3 mature                       -> core_3of3
    # HIT_2    : up in 2 mature                            -> supported_2of3
    # CONFLICT : up in one mature, down in another         -> excluded (clean votes disagree)
    # BIP_1CLN : up in ONE mature, biphasic in another     -> excluded (only 1 CLEAN vote)
    # BIP_2CLN : up in TWO mature, biphasic in a third     -> supported_2of3 (biphasic recorded,
    #                                                         doesn't count as a clean vote)
    # LUCAS_ONLY: only Lucas                               -> excluded (annotation, not a vote)
    # FLIP     : up in 2 mature, down at Lucas onset       -> included, onset_agrees_mature False
    w("aghajanova2011_invitro_decidualization_nonendo.csv",
      [("HIT_CORE", "up"), ("HIT_2", "up"), ("CONFLICT", "up"),
       ("BIP_1CLN", "up"), ("BIP_2CLN", "up"), ("FLIP", "up")])
    w("rytkonen2019_complete_signature.csv",
      [("HIT_CORE", "up"), ("HIT_2", "up"), ("CONFLICT", "down"),
       ("BIP_1CLN", "up"), ("BIP_1CLN", "down"),      # -> biphasic within Rytkonen
       ("BIP_2CLN", "up"), ("BIP_2CLN", "down"),      # -> biphasic within Rytkonen
       ("FLIP", "up")])
    w("stadtmauer2022_hits.csv",
      [("HIT_CORE", "up"), ("BIP_2CLN", "up")])       # 2nd clean vote for BIP_2CLN
    w("lucas2020_S2vsS3_decidualization_onset.csv",
      [("LUCAS_ONLY", "up"), ("FLIP", "down")])

    best, _ = build(d)
    got = {r.gene: (r.consensus_direction, r.tier, r.onset_agrees_mature, r.mature_biphasic)
           for r in best.itertuples()}
    expect = {
        "HIT_CORE": ("up", "core_3of3", "", 0),
        "HIT_2":    ("up", "supported_2of3", "", 0),
        "BIP_2CLN": ("up", "supported_2of3", "", 1),   # biphasic Rytkonen recorded, 2 clean votes qualify it
        "FLIP":     ("up", "supported_2of3", False, 0),
    }
    excluded = {"CONFLICT", "LUCAS_ONLY", "BIP_1CLN"}
    ok = True
    print(f"{'gene':11s} {'expected':32s} got")
    for g, e in expect.items():
        good = got.get(g) == e
        ok &= good
        print(f"{g:11s} {str(e):32s} {got.get(g)}  {'ok' if good else 'XX <-- mismatch'}")
    for g in excluded:
        good = g not in got
        ok &= good
        print(f"{g:11s} {'(excluded)':32s} {'absent' if good else 'PRESENT'}  {'ok' if good else 'XX <-- should be excluded'}")
    print("\nSELFTEST", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--supp-dir", help="directory holding the four per-paper signature CSVs")
    ap.add_argument("--out", help="write the best-hits list here")
    ap.add_argument("--consensus-out", help="also write the full per-gene consensus table here")
    ap.add_argument("--verify", metavar="REF_CSV", help="rebuild and diff against this best-hits CSV")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(selftest())
    if not args.supp_dir:
        ap.error("--supp-dir is required (or use --selftest)")

    best, cons = build(args.supp_dir)
    n_up = int((best.consensus_direction == "up").sum())
    n_dn = int((best.consensus_direction == "down").sum())
    n_core = int(best.tier.str.startswith("core").sum())
    print(f"best-hits: {len(best)}  ({n_up} up / {n_dn} down; {n_core} core, {len(best)-n_core} supported)")

    if args.out:
        best.to_csv(args.out, index=False)
        print(f"wrote {args.out}")
    if args.consensus_out:
        cons.to_csv(args.consensus_out, index=False)
        print(f"wrote {args.consensus_out}  ({len(cons)} genes)")

    ok = True
    if args.verify:
        ok = verify(best, args.verify)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
