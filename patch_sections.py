
import pathlib, shutil, yaml, ast, json
X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final")
for rel in ("data/interim/qc_spatial_mouse.h5ad", "results/tables/spatial_qc_mouse.csv", "results/tables/spatial_qc_obs_mouse.csv",
            "logs/04f_spatial_qc_mouse.log", "workflow/scripts/spatial_qc.py", "workflow/rules/04f_spatial_qc.smk", "config/config.yaml"):
    p = X/rel; b = p.with_name(p.name + ".pre_sections.bak")
    assert not b.exists(), f"backup exists: {b}"; shutil.copy2(p, b); print("  backed up", rel)

def rep1(s, old, new):
    assert s.count(old) == 1, f"anchor count {s.count(old)} for: {old[:70]!r}"; return s.replace(old, new)

sp = X/"workflow/scripts/spatial_qc.py"; s = sp.read_text()
s = rep1(s, "  boxed      hand-drawn exclusion regions the statistics miss\n",
            "  boxed      hand-drawn exclusion regions the statistics miss\n"
            "  sections   label each kept cell's physical cross-section: '{sample}__A' by default;\n"
            "             timepoints listed in config sections get split by a straight line\n"
            "             (y > slope*x + intercept -> 'above' label, else 'below')\n")
s = rep1(s, "        QC'd + oriented coords added alongside)\n",
            "        QC'd + oriented coords + tissue_section added alongside)\n")
s = rep1(s, 'ROT      = P.get("rotations", {}).get(SPECIES, {}) or {}\n',
            'ROT      = P.get("rotations", {}).get(SPECIES, {}) or {}\n'
            'SECT     = P.get("sections", {}).get(SPECIES, []) or []   # line splits for multi-section pucks\n')
s = rep1(s, '    f"{ {k: v for k, v in ROT.items()} }")\n',
            '    f"{ {k: v for k, v in ROT.items()} }; {len(SECT)} section split(s)")\n')
STAGE = open("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final/section_stage.txt").read()
s = rep1(s, "# ---- rotate (per puck; 0 = identity) ---------------------------------------\n",
            STAGE + "# ---- rotate (per puck; 0 = identity) ---------------------------------------\n")
s = rep1(s, 'adata.obs["x_um_rot"], adata.obs["y_um_rot"] = x_rot, y_rot\n',
            'adata.obs["x_um_rot"], adata.obs["y_um_rot"] = x_rot, y_rot\n'
            'adata.obs["tissue_section"] = pd.Categorical(section.astype(object).where(section.notna(), None))\n')
s = rep1(s, '    "rotations_json": json.dumps({k: float(v) for k, v in ROT.items()}),\n',
            '    "rotations_json": json.dumps({k: float(v) for k, v in ROT.items()}),\n'
            '    "sections_json": json.dumps([dict(r) for r in SECT], default=str),\n')
s = rep1(s, 'say("  obs: spatial_qc_reason, x_um_qc/y_um_qc, x_um_rot/y_um_rot")\n',
            'say("  obs: spatial_qc_reason, x_um_qc/y_um_qc, x_um_rot/y_um_rot, tissue_section")\n')
ast.parse(s); sp.write_text(s); print("  spatial_qc.py patched")

rp = X/"workflow/rules/04f_spatial_qc.smk"; t = rp.read_text()
t = rep1(t, "# hand-drawn exclusion boxes). Annotations are present on the input object but\n",
            "# hand-drawn exclusion boxes) and labels each kept cell's physical cross-section\n"
            "# (`tissue_section`, split by a config line for multi-section pucks). Annotations\n"
            "# are present on the input object but\n")
rp.write_text(t); print("  04f rule header patched")

cp = X/"config/config.yaml"; y = cp.read_text()
anchor = "  # Hand-drawn exclusion regions (um, on x_um_dbscan/y_um_dbscan) the distance\n"
y = rep1(y, anchor,
'  # Physical cross-sections per puck. Every kept cell gets tissue_section = "{sample}__A"\n'
"  # unless its timepoint is listed here, in which case a straight line in the coordinate\n"
"  # frame splits it: y > slope*x + intercept -> `above` label, else `below`. Uninduced is\n"
"  # the only mouse timepoint with two cross-sections on one puck (both pucks share a frame;\n"
"  # line chosen 2026-09-02, agrees with k-means(2) on all 6,596 cells, margin >= 127 um).\n"
"  sections:\n"
"    mouse:\n"
"      - {timepoint: uninduced, line: {slope: 1.0, intercept: -500}, above: A, below: B}\n"
+ anchor)
cp.write_text(y)
cfg = yaml.safe_load(y); print("  config sections:", cfg["spatial_qc"]["sections"])
