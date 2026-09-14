
import pathlib, yaml, re
X = pathlib.Path("/n/lab_storage/eddy_lab/users/nhilgert/xmens_final"); cp = X/"config/config.yaml"
assert (X/"config/config.yaml.pre_boxes3.bak").exists()
y = cp.read_text()
# the mouse exclusion block = every '      - {timepoint: ...}' line (and my comment lines) between '    mouse:' under 'exclusions:' and the next dedent
m = re.search(r"(  exclusions:\n    mouse:\n)((?:      .*\n)+)", y)
assert m, "exclusions.mouse block not found"
old = m.group(2)
old_entries = [yaml.safe_load(l.strip()[2:]) for l in old.splitlines() if l.strip().startswith("- {")]
print("  old entries:", len(old_entries))
others = [e for e in old_entries if e["timepoint"] != "uninduced"]
assert len(others) == 5, others   # 1dpi_early(1), 2dpi(1), 3dpi(3)
def fmt(e):
    tp = e.pop("timepoint"); rest = ", ".join(f"{k}: {v}" for k, v in e.items())
    return "      - {timepoint: " + f"{tp:<12s}" + ", " + rest + "}"
new = (
"      # --- uninduced: union of all exclusions in minimal form (rewritten 2026-09-02). Each box verified to cover every\n"
"      #     previously-boxed cell and to add only the user-named stragglers. Both pucks (mck8a/8b) share one frame.\n"
"      - {timepoint: uninduced,   y_gt: 5000}                 # above both sections\n"
"      - {timepoint: uninduced,   x_lt: 500}                  # left margin\n"
"      - {timepoint: uninduced,   x_lt: 3000, y_lt: 1950}     # below-left; absorbs old (x<3000,y<1500) + (2500<x<3000,1500<y<2000) + clump [1300,1600]x[1700,1900]\n"
"      - {timepoint: uninduced,   x_lt: 3400, y_lt: 700}      # low strip for the isolated pair at [3000,3400]x[500,600]; right-section rim (y>=1036) untouched\n"
"      - {timepoint: uninduced,   x_gt: 5400, y_gt: 3000}     # top-right stragglers (paired positions in 8a/8b)\n"
"      # --- induced timepoints\n"
+ "\n".join(fmt(dict(e)) for e in others) + "\n")
y = y.replace(old, new); cp.write_text(y)
cfg = yaml.safe_load(y); ex = cfg["spatial_qc"]["exclusions"]["mouse"]
from collections import Counter; print("  new entries:", len(ex), dict(Counter(e["timepoint"] for e in ex)))
for e in ex: print("    ", e)
assert cfg["spatial_qc"]["n_mads"] == 5.0
