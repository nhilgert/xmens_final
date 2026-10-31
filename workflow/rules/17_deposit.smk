# workflow/rules/17_deposit.smk
# Deposit copies of the release objects: each release object with obs columns, uns entries and obsm embeddings
# pruned and renamed as listed in config/deposit_objects.yaml (user decisions of 2026-10-30), for two sets:
#   with_new_cell_types      release 2026_10_30 (annotation vocabulary v9)
#   without_new_cell_types   release 2026_10_29b (vocabulary v8)
# The mouse objects carry the zone-module AUCell scores of the run of record (Fig. S24) for six modules, on the
# 3.5 dpi nuclei only. Objects listed under pseudotime_of_record (the human fibroblast/decidual objects, refit on
# 30 Oct) take their pseudotime, its diffusion map and its uns record from that run instead of the release. Release
# objects are read, never modified; counts, embeddings, graphs and var are copied as they are. Target:
# deposit_objects_all (not in rule all).

configfile: "config/deposit_objects.yaml"

DEP = config["deposit_objects"]
DEP_OUT = DEP["out_dir"]
DEP_SETS = list(DEP["sets"])
DEP_OBJECTS = DEP["objects"]


def deposit_release_object(wildcards):
    # ancient(): a release object whose mtime moved must not trigger a rebuild of its deposit copy on its own.
    return ancient(f"data/deliverables/{DEP['sets'][wildcards.set]['release']}/{wildcards.obj}.h5ad")


def deposit_zone_aucell(wildcards):
    # Only the mouse objects carry zone-module scores.
    return DEP["sets"][wildcards.set]["zone_aucell"] if wildcards.obj.startswith("Cevrim_XMens") else []


def deposit_pseudotime_of_record(wildcards):
    # An object whose pseudotime is newer than its release depends on that run's table, run record and script.
    entry = DEP.get("pseudotime_of_record", {}).get(wildcards.obj)
    return [entry["table"], entry["run_record"], entry["script"]] if entry else []


rule deposit_object:
    input:
        h5ad=deposit_release_object,
        zone_aucell=deposit_zone_aucell,
        pseudotime=deposit_pseudotime_of_record,
        settings="config/deposit_objects.yaml",
        script="workflow/scripts/deposit/prune_release_object.py",
    output:
        # protected(): deposit copies are never overwritten; a new build needs a new out_dir.
        h5ad=protected(DEP_OUT + "/{set}/{obj}.h5ad"),
        provenance=protected(DEP_OUT + "/{set}/{obj}.provenance.json"),
    params:
        py=DEP["python"],
        release=lambda wildcards: DEP["sets"][wildcards.set]["release"],
    wildcard_constraints:
        set="|".join(DEP_SETS),
        obj="|".join(DEP_OBJECTS),
    log:
        "logs/deposit_objects/{set}/{obj}.log",
    threads: 1
    resources:
        mem_mb=30000,
        runtime=90,
    shell:
        r"""
        {params.py} -u {input.script} \
            --input {input.h5ad} \
            --release {params.release} \
            --settings {input.settings} \
            --zone-aucell '{input.zone_aucell}' \
            --output {output.h5ad} \
            --provenance {output.provenance} \
            > {log} 2>&1
        """


# Both sets, every object, with its provenance record.
rule deposit_objects_all:
    input:
        expand(DEP_OUT + "/{set}/{obj}.{ext}", set=DEP_SETS, obj=DEP_OBJECTS, ext=["h5ad", "provenance.json"]),
