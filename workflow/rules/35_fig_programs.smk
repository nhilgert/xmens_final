# workflow/rules/35_fig_programs.smk
# Panels of the shared transcriptional programs (Fig. 5), added 2026-10-30:
#
#   panels_programs_5e   Fig. 5E: correlation of SCENIC regulon activity with pseudotime, human (x) against
#                        X-Mens (y), one point per TF with a regulon in both species, drawn in the style of
#                        the manuscript panel (quadrant labels, dotted lines at |rho| = 0.2, named regulons)
#
# Input: the regulon table of record (regulon_quadrants_input_fibdec_latent.csv, written by
# scenic_quadrant_rho_fibdec_latent.py from the floor-6 consensus regulons, the mouse pt_harmonized_rank and the
# human pseudotime of record), staged in resources/figures/panels_programs/ with its upstream provenance; its sha256
# is pinned in config/panels_programs.yaml and checked by the plot script. The human pseudotime its provenance names
# must be config pseudotime_of_record.human (sha256), checked when the workflow is read.
# Outputs: results/figures/panels/programs/panel_5E_one.{svg,png,provenance.json}; SVG and PNG are
# byte-reproducible (fixed SVG ids, no dates). Target: panels_programs_all (not in rule all).

import json

include: "91_panels_versions.smk"

configfile: "config/panels_programs.yaml"

PPG = config["panels_programs"]

# Decision G: Fig. 5E's human rho must be computed on the human pseudotime of record.
with open(PPG["regulon_table"]["upstream_provenance"]) as _ppg_handle:
    _PPG_HUMAN_PT_SHA256 = json.load(_ppg_handle)["human_pt_side_table"]["sha256"]
if _PPG_HUMAN_PT_SHA256 != config["pseudotime_of_record"]["human"]["sha256"]:
    raise ValueError(
        f"{PPG['regulon_table']['upstream_provenance']}: human pseudotime sha256 {_PPG_HUMAN_PT_SHA256} is not "
        f"config pseudotime_of_record.human ({config['pseudotime_of_record']['human']['sha256']})"
    )
PPG_OUT = "results/figures/panels/programs"
PPG_SCRIPTS = "workflow/scripts/figures/panels/programs"
PPG_EXTS = [".svg", ".png", ".provenance.json"]

# Matplotlib's cache must not land in HOME, which is full.
PPG_ENV = (
    "export XDG_CACHE_HOME={cache} MPLCONFIGDIR={cache}/mpl PYTHONDONTWRITEBYTECODE=1; "
    "mkdir -p {cache}/mpl"
).format(cache=PANEL_CACHE)


# Every panel of this rule file.
rule panels_programs_all:
    input:
        multiext(f"{PPG_OUT}/panel_5E_one", *PPG_EXTS),


# Fig. 5E: regulon-pseudotime correlation in human and X-Mens, with the manuscript's regulons named.
rule panels_programs_5e:
    input:
        table = PPG["regulon_table"]["path"],
        table_provenance = PPG["regulon_table"]["upstream_provenance"],
        script = f"{PPG_SCRIPTS}/plot_regulon_quadrants.py",
    output:
        multiext(f"{PPG_OUT}/panel_5E_one", *PPG_EXTS),
    params:
        py = PANEL_PY,
        env = PPG_ENV,
        sha256 = PPG["regulon_table"]["sha256"],
        source = PPG["regulon_table"]["source_of_record"],
        threshold = PPG["fig5e"]["threshold"],
        named = ",".join(PPG["fig5e"]["named"]),
        stem = f"{PPG_OUT}/panel_5E_one",
    log:
        "logs/panels_programs/5e.log",
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 10,
    shell:
        r"""
        {params.env}
        {params.py} -u {input.script} \
            --table {input.table} \
            --table-sha256 {params.sha256} \
            --table-provenance {input.table_provenance} \
            --source-of-record {params.source} \
            --threshold {params.threshold} \
            --named {params.named} \
            --out {params.stem} \
            > {log} 2>&1
        """
