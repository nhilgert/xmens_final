
import pathlib, shutil, yaml
X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
for rel in ("data/interim/qc_spatial_mouse.h5ad", "results/tables/spatial_qc_mouse.csv",
            "results/tables/spatial_qc_obs_mouse.csv", "logs/04f_spatial_qc_mouse.log", "config/config.yaml"):
    p = X/rel; b = p.with_name(p.name + ".pre_boxes2.bak")
    assert not b.exists(), f"backup exists: {b}"
    shutil.copy2(p, b); print("  backed up", rel)
cp = X/"config/config.yaml"; y = cp.read_text()
anchor = "      - {timepoint: 3_dpi,       y_lt: 1500, x_gt: 6500}\n"
assert y.count(anchor) == 1
new = (anchor +
       "      # added 2026-09-02 (user-specified on the 5-MAD post-QC maps; see manual_box_queries_v2.png)\n"
       "      - {timepoint: uninduced,   x_gt: 5400, y_gt: 3000}   # 8 off-tissue stragglers top-right, both pucks (paired positions in 8a/8b)\n"
       "      - {timepoint: 3_dpi,       x_gt: 9000, y_gt: 8000}   # 8-cell off-tissue clump ~1.3 mm NE of mck10a\n")
y = y.replace(anchor, new); cp.write_text(y)
cfg = yaml.safe_load(y); ex = cfg["spatial_qc"]["exclusions"]["mouse"]
print("  exclusions now:", len(ex)); print("  last two:", ex[-2:])
assert cfg["spatial_qc"]["n_mads"] == 5.0 and cfg["spatial_qc"]["mad_scale"] == 1.0
