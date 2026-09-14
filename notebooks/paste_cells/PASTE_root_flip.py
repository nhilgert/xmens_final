# ==========================================================================
# PASTE THIS over the root-cell-selection cell (the one with dists.argmax()).
#
# MOUSE NEEDS "near". With "far" the root landed on a 3_dpi cell at the LATE end
# of the diffusion axis and pseudotime ran backwards — headless run 6a7ac749 gave
# rho(pt, day) = -0.4933 and rho(pt, subtype_rank) = -0.7693, with
# decidual_mature scoring the LOWEST pseudotime.
#
# |rho| = 0.77 means the trajectory is strong and coherent; only its sign is wrong.
# Flip the ROOT, never negate dpt_pseudotime — downstream code assumes [0, 1].
# ==========================================================================
ROOT_DIRECTION = "near"       # "near" = argmin (mouse) | "far" = argmax (original)

dc23 = adata.obsm["X_diffmap"][:, [2, 3]]
cluster_key = CONFIG["cluster_key"]
anchor = CONFIG["root_candidate_label"]

labels = adata.obs[cluster_key].astype(str)
if anchor in set(labels.unique()):
    root_mask = (labels == anchor).values
    root_source = f"{cluster_key}={anchor}"
else:
    earliest_day = adata.obs["day"].astype(float).min()
    root_mask = (adata.obs["day"].astype(float) == earliest_day).values
    root_source = f"earliest day ({earliest_day}) — {anchor!r} absent after filtering"
    print(f"WARNING: {anchor!r} not found; using earliest-timepoint cells")

print(f"Root source: {root_source}")
print(f"Candidate pool: {int(root_mask.sum()):,} cells")
assert root_mask.sum() > 0, "empty root candidate pool"

target_centroid = dc23[~root_mask].mean(axis=0)
candidate_indices = np.where(root_mask)[0]
dists = np.linalg.norm(dc23[root_mask] - target_centroid, axis=1)

assert ROOT_DIRECTION in ("far", "near"), ROOT_DIRECTION
pick = dists.argmax() if ROOT_DIRECTION == "far" else dists.argmin()
root_idx = int(candidate_indices[pick])

print(f"root_direction = {ROOT_DIRECTION!r}   "
      f"dist to non-candidate centroid = {dists[pick]:.4f}   "
      f"pool range [{dists.min():.4f}, {dists.max():.4f}]")
print(f"\nRoot cell: obs index {root_idx}")
print(f"  obs_name  = {adata.obs_names[root_idx]}")
print(f"  timepoint = {adata.obs['timepoint'].iloc[root_idx]}")
print(f"  day       = {adata.obs['day'].iloc[root_idx]}")
print(f"  subtype   = {labels.iloc[root_idx]}")
print(f"  DC2, DC3  = {dc23[root_idx]}")

adata.obs["is_root"] = "other"
adata.obs.loc[adata.obs.index[root_idx], "is_root"] = "root"
adata.obs["is_root"] = adata.obs["is_root"].astype("category")
