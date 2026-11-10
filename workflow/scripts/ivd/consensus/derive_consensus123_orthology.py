"""Derive a frozen human->mouse orthology table for the 123 consensus IVD genes.

Pipeline copy (2026-11-02) of cagri_ivd_end_to_end_20260815/derive_consensus123_orthology.py (SHA-256 d4b22175...).
The only changes: the upstream inputs, the consensus list, the frozen human mapping and the output path are
command-line arguments instead of paths relative to the script and a hard-coded lab_storage folder, and the
"exactly 123 genes" guard is kept but can be relaxed (--allow-any-length) for sensitivity lists. The resolution
logic, column names and output format are unchanged, so the table it writes for the released list is
byte-identical to resources/ivd/consensus/released_consensus123_mouse_orthology.csv (checked by build_consensus123.py).

Two independent sources, both required to agree on the mouse copy number:

  1. Ensembl 116 human/mouse orthology components (ensembl_human_mouse_components.tsv, homology_type per edge)
  2. OMA All.May2026 Euarchontoglires ancestral HOGs (oma_human_mouse_ancestral_hogs.tsv, copy_number_class)

Legacy gene symbols in the consensus list are resolved to current symbols with the project's own alias audit
(77,453 rows, Entrez-backed) rather than by hand.

Usage:
    python derive_consensus123_orthology.py --consensus CONSENSUS.csv --heca-mapping consensus123_heca_mapping.csv \
        --components ensembl_human_mouse_components.tsv --oma oma_human_mouse_ancestral_hogs.tsv \
        --annotation ensembl_116_homo_sapiens_gene_annotation.tsv.gz --alias aghajanova_alias_to_current_audit.csv \
        --out consensus123_mouse_orthology.csv [--tested-support-rows support_rows.csv] [--check]
"""

import argparse
import collections
import csv
import gzip
import hashlib
import sys
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--consensus", required=True, type=Path, help="consensus list (columns gene, consensus_direction, ...)")
    ap.add_argument("--heca-mapping", required=True, type=Path, help="frozen consensus -> human Ensembl ID mapping")
    ap.add_argument("--components", required=True, type=Path, help="Ensembl 116 human/mouse orthology components (TSV)")
    ap.add_argument("--oma", required=True, type=Path, help="OMA All.May2026 Euarchontoglires HOGs (TSV)")
    ap.add_argument("--annotation", required=True, type=Path, help="Ensembl 116 human gene annotation (TSV, gzipped)")
    ap.add_argument("--alias", required=True, type=Path, help="alias -> current symbol audit (CSV)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--tested-support-rows", type=Path, default=None,
                    help="optional consensus_ortholog_xmens_broad_support_rows.csv of a release, for the consistency gate")
    ap.add_argument("--check", action="store_true", help="verify only, no write")
    ap.add_argument("--allow-any-length", action="store_true", help="accept lists that are not exactly 123 genes")
    args = ap.parse_args(argv)
    check_only = args.check

    consensus = read_csv(args.consensus)
    heca_map = {r["human_gene"]: r for r in read_csv(args.heca_mapping)}
    if len(consensus) != 123 and not args.allow_any_length:
        raise SystemExit(f"expected 123 consensus genes, got {len(consensus)}")

    # ---- alias -> current symbol (Entrez-backed, from the project's own audit) ----
    alias = {}
    for row in read_csv(args.alias):
        alias.setdefault(row["alias"].upper(), row["current_gene"])

    # ---- Ensembl 116 human annotation: current symbol -> gene id ----
    symbol_to_id = collections.defaultdict(list)
    for row in read_tsv(args.annotation, gz=True):
        symbol_to_id[row["gene_symbol"]].append(row["ensembl_gene_id"])

    # ---- Ensembl 116 orthology components ----
    components = read_tsv(args.components)
    by_component = collections.defaultdict(list)
    component_of = {}
    for row in components:
        by_component[row["component_id"]].append(row)
        component_of[row["ensembl_gene_id"]] = row

    # ---- OMA ancestral HOGs ----
    # Group on `oma_hog_id`, the Euarchontoglires-level HOG, NOT `root_hog_id`.
    # The root spans deeper ancestral groups, so pooling on it merges lineage-specific
    # expansions that do not descend from the same Euarchontoglires ancestor: human
    # IFIT1's root pulls in mouse Ifit1/Ifit1bl1/Ifit1bl2 even though OMA scores it
    # mouse_copy_number=0. `mouse_copy_number` is the authority and is asserted below.
    hogs = read_tsv(args.oma)
    hog_of = collections.defaultdict(list)
    by_hog = collections.defaultdict(list)
    for row in hogs:
        hog_of[row["ensembl_gene_id"]].append(row)
        by_hog[row["oma_hog_id"]].append(row)

    rows = []
    disagreements = []
    for entry in consensus:
        gene = entry["gene"]
        mapped = heca_map.get(gene, {})
        current = mapped.get("current_gene_symbol") or alias.get(gene.upper()) or gene
        gene_id = mapped.get("human_ensembl_gene_id") or ""
        if not gene_id:
            candidates = symbol_to_id.get(current, [])
            gene_id = candidates[0] if len(candidates) == 1 else ""

        ens_mouse, ens_type, ens_mouse_ids = [], "", {}
        component = component_of.get(gene_id)
        if component:
            ens_type = component["homology_type"]
            members = [x for x in by_component[component["component_id"]] if x["species"] == "mouse"]
            ens_mouse = sorted(x["gene_symbol"] for x in members)
            ens_mouse_ids = {x["gene_symbol"]: x["ensembl_gene_id"] for x in members}

        oma_mouse, oma_class, oma_n = [], "", ""
        for hog in hog_of.get(gene_id, []):
            oma_class = hog["copy_number_class"]
            oma_n = hog["mouse_copy_number"]
            oma_mouse = sorted(
                x["gene_symbol"] for x in by_hog[hog["oma_hog_id"]] if x["species"] == "mouse"
            )
            if bool(oma_mouse) != (int(oma_n) > 0):
                raise SystemExit(
                    f"{gene}: OMA hog {hog['oma_hog_id']} lists mouse members {oma_mouse} but "
                    f"mouse_copy_number={oma_n}"
                )
            break

        # Ensembl 116 is primary; OMA is an advisory cross-check.
        #
        # OMA All.May2026 infers its HOGs from Ensembl *102* proteins and uses 116 only as
        # "an identifier/display bridge" (its own attrition note), so its mouse coverage
        # lags. It scores ABI3BP, ADRA2C, EFHD1, TUBA1C and TUBB4B as `human_only` even
        # though Ensembl 116 gives each a one-to-one mouse ortholog and mouse plainly has
        # the gene. Two of those (ABI3BP, ADRA2C) are tested against Abi3bp/Adra2c in the
        # release, so an OMA veto would also contradict the shipped mapping.
        #
        # OMA is therefore consulted where Ensembl has no component (it recovers CFH) and
        # recorded per row, but it cannot overturn a positive Ensembl call.
        mouse_set = ens_mouse or oma_mouse
        if ens_mouse and oma_n and int(oma_n) == 0:
            oma_agrees = "False_oma_lags_ensembl_102"
        elif not ens_mouse and not oma_mouse:
            oma_agrees = "True_neither_has_mouse_member"
        elif ens_mouse and oma_mouse:
            oma_agrees = "True" if set(ens_mouse) & set(oma_mouse) else "False_disjoint"
        elif oma_mouse:
            oma_agrees = "oma_only_no_ensembl_component"
        else:
            oma_agrees = "ensembl_only_no_oma_row"
        # Only a STRICT one-to-one Ensembl 116 ortholog counts as "the" mouse ortholog.
        # Anything else has no unambiguous comparator: human PRL sits in a component with
        # all 26 mouse Prl-family genes, and choosing `Prl` from it by symbol identity
        # names the pituitary gene while the mouse decidual program lives in Prl8a2 /
        # Prl3c1. Resolving a one-to-many family by symbol is the exact error this
        # resource exists to prevent, so one2many, many2many and OMA-only calls report no
        # ortholog and keep their members in the audit columns instead.
        one_to_one = ens_type == "ortholog_one2one" and len(ens_mouse) == 1
        mouse = ens_mouse[0] if one_to_one else ""

        has_ortholog = one_to_one
        has_any_homolog = bool(mouse_set)
        # Only a genuine content conflict is fatal: both sources name mouse members and
        # they share none. Coverage gaps in either direction are recorded, not fatal.
        if ens_mouse and oma_mouse and not (set(ens_mouse) & set(oma_mouse)):
            disagreements.append((gene, f"disjoint mouse sets ensembl={ens_mouse}", f"oma={oma_mouse}"))
        if not has_any_homolog and oma_class and oma_class != "human_only":
            disagreements.append((gene, "no mouse member found but OMA class is", oma_class))

        rows.append(
            {
                "human_gene": gene,
                "current_human_symbol": current,
                "human_ensembl_gene_id": gene_id,
                "mouse_ortholog": mouse,
                "mouse_ensembl_gene_id": ens_mouse_ids.get(mouse, "") if one_to_one else "",
                "has_mouse_ortholog": str(has_ortholog),
                "strict_one_to_one": str(one_to_one),
                "has_any_mouse_homolog": str(has_any_homolog),
                "ensembl116_homology_type": ens_type,
                "ensembl116_mouse_members": ";".join(ens_mouse),
                "oma_copy_number_class": oma_class,
                "oma_mouse_copy_number": oma_n,
                "oma_mouse_members": ";".join(oma_mouse),
                "oma_agrees_with_ensembl": oma_agrees,
            }
        )

    if disagreements:
        for d in disagreements:
            print(f"DISAGREEMENT {d}", file=sys.stderr)
        raise SystemExit("Ensembl 116 and OMA disagree on ortholog existence; resolve before writing")

    # ---- consistency gate: every mouse gene a release tested must sit inside the resolved
    # orthology set for its human gene (optional; the upstream script found the release by glob). ----
    tested = args.tested_support_rows
    if tested is not None:
        resolved = {r["human_gene"]: r for r in rows}
        checked = skipped = 0
        for row in read_csv(tested):
            r = resolved[row["human_gene"]]
            members = set(filter(None, r["ensembl116_mouse_members"].split(";"))) | set(
                filter(None, r["oma_mouse_members"].split(";"))
            )
            if not members:
                skipped += 1
                continue
            if row["mouse_gene"] not in members:
                raise SystemExit(
                    f"{row['human_gene']}: tested mouse gene {row['mouse_gene']!r} is not in the "
                    f"resolved orthology set {sorted(members)}"
                )
            checked += 1
        print(f"consistency gate: {checked} tested mouse genes confirmed orthologous"
              f" ({skipped} unresolved) against {tested}")

    n_orth = sum(1 for r in rows if r["strict_one_to_one"] == "True")
    n_any = sum(1 for r in rows if r["has_any_mouse_homolog"] == "True")
    print(f"resolved {len(rows)} genes: {n_orth} strict one-to-one orthologs, {len(rows) - n_orth} not")
    print(f"  ({n_any} have some mouse homolog; the extra {n_any - n_orth} are one2many/many2many"
          f" with no unambiguous comparator)")
    for r in rows:
        if r["strict_one_to_one"] != "True":
            members = r["ensembl116_mouse_members"] or r["oma_mouse_members"] or "-"
            n = len([m for m in members.split(";") if m and m != "-"])
            print(f"    {r['human_gene']:<10} {r['ensembl116_homology_type'] or '(no ensembl component)':<24}"
                  f" mouse_n={n:<3} {members if n <= 4 else members.split(';')[0] + f' ...+{n - 1} more'}")

    for path in (args.components, args.oma, args.annotation, args.alias):
        print(f"  upstream {path.name}: sha256 {sha256(path)}")

    if check_only:
        return
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {args.out}")
    print(f"  sha256 {sha256(args.out)}")


if __name__ == "__main__":
    main(sys.argv[1:])
