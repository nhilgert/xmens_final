# workflow/rules/17_deposit.smk
# Deposit copies of the release objects: each release object with obs columns, uns entries and obsm embeddings
# pruned and renamed as listed in config/deposit_objects.yaml (user decisions of 2026-10-30), for two sets:
#   with_new_cell_types      the release of record (config deliverables.version)
#   without_new_cell_types   the release before the new cell types (sets.without_new_cell_types.version)
# into the deposit folder of the release of record (out_dir, builds). The mouse objects carry the zone-module AUCell
# scores of the run of record (Fig. S24) for six modules, on the 3.5 dpi nuclei only. Objects listed under
# pseudotime_of_record (the human fibroblast/decidual objects, refit on 30 Oct) carry that run's pseudotime, its
# diffusion map and its uns record: copied from a release that already has them, substituted for an older release's.
# Release objects are read, never modified; counts, embeddings, graphs and var are copied as they are.
# Target: deposit_objects_all (in rule all).

configfile: "config/deposit_objects.yaml"

DEP = config["deposit_objects"]
DEP_OUT = DEP["out_dir"].format(release=RELEASE, build=DEP.get("builds", {}).get(RELEASE_VERSION, ""))
DEP_SETS = list(DEP["sets"])
DEP_OBJECTS = DEP["objects"]


def deposit_version(set_name: str) -> str:
    """Version label of a deposit set's release; `record` is the release of record."""
    version = str(DEP["sets"][set_name]["version"])
    return RELEASE_VERSION if version == "record" else version


def deposit_release(set_name: str) -> str:
    return RELEASES[deposit_version(set_name)]


def deposit_zone_table(set_name: str) -> str:
    """The zone-module AUCell table of the run of record made on the set's release."""
    version = deposit_version(set_name)
    recorded = DEP.get("zone_aucell_of_version", {}).get(version)
    return recorded or DEP["zone_aucell"].format(run=DEP["zone_aucell_run"], release=RELEASES[version])


RELEASE_READERS["config/deposit_objects.yaml sets.with_new_cell_types"] = deposit_release("with_new_cell_types")

# The deposit substitutes the same pseudotime run that 04g2 writes into the release (config pseudotime_of_record):
# a deposit entry naming a run table of record must agree with it on every field the two records share.
_PT_RECORD_OF_TABLE = {record["path"]: record for record in config.get("pseudotime_of_record", {}).values()}
for _obj, _entry in DEP.get("pseudotime_of_record", {}).items():
    _record = _PT_RECORD_OF_TABLE.get(_entry["table"])
    if _record is None:
        continue
    _pairs = [("table_sha256", "sha256"), ("run_record", "run_record"), ("script", "script"),
              ("latent_copy", "latent_copy"), ("latent_sha256", "latent_sha256"), ("note", "note")]
    _differ = [deposit_key for deposit_key, record_key in _pairs
               if str(_entry.get(deposit_key)) != str(_record.get(record_key))]
    if _differ:
        raise ValueError(f"config/deposit_objects.yaml pseudotime_of_record.{_obj}: {_differ} differ from config "
                         f"pseudotime_of_record for {_entry['table']}")


def deposit_files(set_name: str) -> list:
    return [f"{DEP_OUT}/{set_name}/{obj}.{ext}" for obj in DEP_OBJECTS for ext in ("h5ad", "provenance.json")]


# A set whose deposit copies all exist is complete: its files are protected and never rebuilt, so deposit_object is
# defined only for the sets still missing a file, and a later edit of the settings or the script cannot schedule one.
DEP_OPEN_SETS = [name for name in DEP_SETS if not all(Path(p).exists() for p in deposit_files(name))]


def deposit_release_object(wildcards):
    # ancient(): a release object whose mtime moved must not trigger a rebuild of its deposit copy on its own.
    return ancient(f"data/deliverables/{deposit_release(wildcards.set)}/{wildcards.obj}.h5ad")


def deposit_zone_aucell(wildcards):
    # Only the mouse objects carry zone-module scores.
    return deposit_zone_table(wildcards.set) if wildcards.obj.startswith("Cevrim_XMens") else []


def deposit_pseudotime_of_record(wildcards):
    # An object whose pseudotime is newer than its release depends on that run's table, run record and script.
    entry = DEP.get("pseudotime_of_record", {}).get(wildcards.obj)
    return [entry["table"], entry["run_record"], entry["script"]] if entry else []


if DEP_OPEN_SETS:

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
            release=lambda wildcards: deposit_release(wildcards.set),
        wildcard_constraints:
            set="|".join(DEP_OPEN_SETS),
            obj="|".join(DEP_OBJECTS),
        log:
            f"logs/deposit_objects/{Path(DEP_OUT).name}/{{set}}/{{obj}}.log",
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
        [path for set_name in DEP_SETS for path in deposit_files(set_name)],
