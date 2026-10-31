"""Round-2 lineage subsets: the lineages that gain nuclei when atlas clusters excluded in July are
assigned to their lineages (rule split_lineages_round2_{species}, 04b2_reintegrate_round2.smk).

Reads the whole-tissue object data/interim/latent_{species}.h5ad (the July fit, pinned by sha256), the
round-2 atlas decision table, the July atlas decision table and the July lineage latents. Writes
{outdir}/lineages/{species}_{lineage}.h5ad for each round-2 lineage, and a report with one row per lineage.

The subsets are built the way stage 04 builds them (apply_annotation.apply_map, then the subset steps of
apply_annotation.write_atlas), so each one is what stage 04 would have written had the round-2 table been
used in July. One difference: the atlas Leiden ladder is dropped and only the atlas cluster is kept, as
obs['atlas_cluster']. In July the lineage fit overwrote that ladder with its own; round 2 clusters at one
resolution only, so atlas clusterings would otherwise sit under the lineage column names.
The run stops unless each subset holds exactly its July nuclei plus the nuclei of the atlas clusters newly
assigned to it.
"""
import hashlib
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from apply_annotation import apply_map, load_choices  # noqa: E402

sm = snakemake  # noqa: F821
species = sm.params.species
lineages = list(sm.params.lineages)
res = sm.params.res          # atlas clustering column, lineage_spec atlas.<species>.res
slug = sm.params.slug        # lineage label -> file slug, lineage_spec atlas.<species>.labels


def sha256(path, block_size=1 << 24):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def july_nuclei(path):
    """obs_names of a July lineage latent, without reading its matrices."""
    latent = ad.read_h5ad(path, backed="r")
    names = latent.obs_names.copy()
    latent.file.close()
    return names


with open(sm.log[0], "w") as log:
    def say(message):
        print(f"[split_lineages_round2:{species}] {message}", file=log, flush=True)

    found = sha256(sm.input.atlas)
    if found != sm.params.atlas_sha256:
        raise ValueError(f"{sm.input.atlas}: sha256 {found}, expected the July fit {sm.params.atlas_sha256}")
    say(f"atlas object is the July fit (sha256 {found[:12]})")

    # The only change round 2 may make: an atlas cluster that was EXCLUDE in July now goes to a round-2 lineage.
    july = load_choices(sm.input.july_decisions)
    round2 = load_choices(sm.input.decisions)
    if set(july) != set(round2):
        raise ValueError(f"the two atlas tables list different clusters: {sorted(set(july) ^ set(round2))}")
    changed = {c: (july[c], round2[c]) for c in round2 if july[c] != round2[c]}
    added = {c: new for c, (old, new) in changed.items() if old == "EXCLUDE" and new in lineages}
    if set(added) != set(changed):
        unexpected = {c: changed[c] for c in changed if c not in added}
        raise ValueError(f"round-2 table changes more than EXCLUDE -> round-2 lineage: {unexpected}")
    say(f"atlas clusters newly assigned: {added}")

    adata = sc.read_h5ad(sm.input.atlas)
    say(f"loaded {adata.n_obs:,} nuclei x {adata.n_vars:,} genes")
    apply_map(adata, res, round2, rollup={}, output_mode="atlas")
    atlas_cluster = adata.obs[res].astype(str)

    # Subset steps of apply_annotation.write_atlas, for the round-2 lineages only.
    kept = adata[~adata.obs.excluded].copy()
    for column in [c for c in kept.obs if c.startswith(("score_", "ref_"))] + ["label"]:
        if column in kept.obs:
            del kept.obs[column]
    kept.obs["atlas_cluster"] = kept.obs[res].astype(str)
    atlas_ladder = [c for c in kept.obs if c.startswith("leiden_scVI_")]
    kept.obs = kept.obs.drop(columns=atlas_ladder)
    for key in [k for k in kept.uns if k.startswith("leiden_scVI_")]:
        del kept.uns[key]
    say(f"dropped the atlas ladder ({len(atlas_ladder)} columns); atlas cluster kept as obs['atlas_cluster']")
    present = set(kept.obs["lineage"])
    kept.obs["lineage"] = pd.Categorical(kept.obs["lineage"], categories=[l for l in slug if l in present])

    rows = []
    for lineage, july_latent, out in zip(lineages, sm.input.july_latents, sm.output.subsets):
        name = f"{species}_{slug[lineage]}"
        if not (Path(july_latent).name == f"latent_{name}.h5ad" and Path(out).name == f"{name}.h5ad"):
            raise ValueError(f"inputs and outputs out of order for {lineage}: {july_latent}, {out}")
        subset = kept[kept.obs["lineage"] == lineage].copy()

        from_july = july_nuclei(july_latent)
        new_clusters = [c for c, l in added.items() if l == lineage]
        from_added = adata.obs_names[atlas_cluster.isin(new_clusters).to_numpy()]
        if len(from_july.intersection(from_added)):
            raise ValueError(f"{name}: newly assigned nuclei are already in the July lineage")
        expected = from_july.union(from_added)
        if subset.n_obs != len(expected) or not subset.obs_names.isin(expected).all():
            raise ValueError(f"{name}: subset has {subset.n_obs:,} nuclei, expected the {len(from_july):,} "
                             f"July nuclei plus {len(from_added):,} from atlas clusters {new_clusters}")

        subset.write_h5ad(out)
        rows.append(dict(species=species, lineage=lineage, july_nuclei=len(from_july),
                         added_atlas_clusters=";".join(new_clusters), added_nuclei=len(from_added),
                         round2_nuclei=subset.n_obs, file=out))
        say(f"{name}: {len(from_july):,} July + {len(from_added):,} added = {subset.n_obs:,} -> {out}")

    pd.DataFrame(rows).to_csv(sm.output.report, sep="\t", index=False)
