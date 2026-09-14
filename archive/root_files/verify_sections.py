
import scanpy as sc, json
EXPECTED = json.load(open("expected_sections.json"))
a = sc.read_h5ad("data/interim/qc_spatial_mouse.h5ad", backed="r"); o = a.obs
o[["sample","timepoint","x_um_dbscan","y_um_dbscan","x_um_qc","y_um_qc","spatial_qc_reason","tissue_section","coarse_refined","fine_refined"]].to_csv("results/tables/spatial_qc_obs_mouse.csv")
kept = o.spatial_qc_reason.astype(str) == "kept"
assert o.tissue_section.notna().sum() == kept.sum(), "section labels must exist exactly on kept cells"
assert o.tissue_section[~kept].isna().all()
vc = o.tissue_section.value_counts()
print("sections:", vc[vc > 0].to_dict())
for k, v in EXPECTED.items():
    assert int(vc.get(k, 0)) == v, f"{k}: got {vc.get(k,0)} expected {v}"
print("uninduced per-puck counts match offline expectation")
print("n sections:", int((vc > 0).sum()), "(expect 11)")
print("uns sections_json:", a.uns["spatial_qc"]["sections_json"])
