#!/usr/bin/env python3
"""Ulrich et al. 2024 pseudobulk DE (decidual stroma vs C7-Fibroblast and the other configured contrasts), human arm only.

Pipeline driver (2026-11-02) for the human part of ulrich_ivd_pipeline.py (in this folder, imported and not edited). The
record of 28 Oct 2026 (release_20261028_v7_dS_vs_C7fib_luteal_follicular) was produced by that script, which also
refits X-Mens with the 19 Oct design and writes the decidual-zone outputs; neither is needed for the IVD comparisons
(X-Mens is fitted with the uninduced fibroblasts by fit_mouse_de.py), so this driver runs only:

  prepare_ulrich_base -> per contrast: contrast_pseudobulks / donor_pseudobulks -> arm_aware_eligible -> fit_deseq2
  -> assemble_result, with the acceptance checks of the config (every pinned count must reproduce, unless --record),
  -> the consensus support rows per contrast (upstream build_heca_support), ulrich_contrasts.csv, ivd_support_summary.csv
  (human rows), validation_checks.csv, run_metadata.json and release_manifest.json.

The functions are the pipeline's own, so the tables are those of the record; the file names and columns are unchanged.
The output folder is written as <out>.inprogress and renamed when complete; an existing output is never overwritten.

Usage: python ulrich_human_de.py --config CONFIG.yaml --out OUT_DIR [--n-cpus 8] [--record]
"""
from __future__ import annotations

import argparse
import atexit
import importlib.util
import json
import platform
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True, help="output folder (created; must not exist)")
    ap.add_argument("--n-cpus", type=int, default=8)
    ap.add_argument("--record", action="store_true", help="write observed acceptance values instead of checking them")
    args = ap.parse_args()

    UP = load_module("ulrich_ivd_pipeline", HERE / "ulrich_ivd_pipeline.py")
    fail, sha256_file = UP.fail, UP.sha256_file
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text())
    if config.get("schema_version") != "1.0" or config.get("analysis_id") != "cagri_ivd_ulrich_deseq2_20261019":
        fail("Unexpected config schema/analysis_id")
    if platform.python_version() != str(config["runtime"]["python"]):
        fail(f"Expected Python {config['runtime']['python']}, observed {platform.python_version()}")
    import pydeseq2
    if pydeseq2.__version__ != str(config["runtime"]["pydeseq2"]):
        fail(f"Expected PyDESeq2 {config['runtime']['pydeseq2']}, observed {pydeseq2.__version__}")
    RP = UP.load_upstream(config)
    record = bool(args.record)
    acc = config.get("acceptance") or {}

    execution_inputs = {"config": config_path, "pipeline": (HERE / "ulrich_ivd_pipeline.py").resolve(),
                        "driver": Path(__file__).resolve(), "upstream_pipeline": Path(config["upstream"]["repo"]) / "run_pipeline.py"}
    exec_hash = {k: sha256_file(p) for k, p in execution_inputs.items()}
    res = {k: RP.resolve_resource(config_path.parent, spec) for k, spec in config["resources"].items()}
    ulrich_obj = Path(config["objects"]["ulrich_cells"]["path"])
    if sha256_file(ulrich_obj) != config["objects"]["ulrich_cells"]["sha256"]:
        fail(f"Object hash mismatch: {ulrich_obj}")

    output = args.out.resolve()
    staging = output.with_name(output.name + ".inprogress")
    if output.exists() or staging.exists():
        fail(f"Refusing existing output/staging path: {output}")
    staging.mkdir(parents=True)
    atexit.register(lambda: staging.exists() and shutil.rmtree(staging, ignore_errors=True))
    inter = staging / "intermediate"; inter.mkdir()
    checks: list[dict[str, Any]] = []
    calls = config["support_calls"]
    observed_acc: dict[str, Any] = {"human": {}, "status_partitions": {}}

    # ---- human (verbatim logic of ulrich_ivd_pipeline.main)
    human = config["human"]
    features_h, base_counts, base_meta = UP.prepare_ulrich_base(RP, ulrich_obj, config["objects"]["ulrich_cells"], res["ulrich_feature_manifest"], human, inter)
    human_results: dict[str, pd.DataFrame] = {}
    human_diag: dict[str, Any] = {}
    for con in human["contrasts"]:
        cid = con["id"]
        between = con.get("kind") == "between_donor"
        if between:
            meta, counts = UP.donor_pseudobulks(base_counts, base_meta, con, int(human["cell_floor"]))
            num, den, block = con["numerator"], con["denominator"], "block"
            if meta.groupby("block")["state"].nunique().lt(2).all():
                fail(f"{cid}: factor and covariate are confounded")
        else:
            meta, counts = UP.contrast_pseudobulks(base_counts, base_meta, con, int(con.get("cell_floor", human["cell_floor"])))
            num, den, block = "decidual", "reference", "donor"
        is_num = meta["state"].eq(num).to_numpy()
        eligible = UP.arm_aware_eligible(counts, is_num, float(human["eligibility_cpm_threshold"]))
        fitted, diag = UP.fit_deseq2(counts.loc[:, eligible.to_numpy()], meta, block=block, factor="state",
                                     numerator=num, denominator=den, n_cpus=args.n_cpus)
        result = UP.assemble_result(features_h, counts, eligible, fitted, meta, "state", num,
                                    source_id=f"ulrich_{cid}", positive_meaning=f"{num}_minus_{den}", calls=calls)
        obs_c = {"pseudobulks": int(len(meta)), "donors": int(meta.donor.nunique()), "eligible_genes": int(eligible.sum()),
                 "count_sum": int(counts.to_numpy().sum()), "strict_calls": int(result["significant_fdr05_lfc1"].sum()),
                 "residual_df": diag["residual_df"]}
        for k, v in obs_c.items():
            UP.check_equal(checks, f"human_{cid}_{k}", v, (acc.get("human", {}).get(cid) or {}).get(k), record)
        observed_acc["human"][cid] = obs_c
        human_diag[cid] = {**diag, **obs_c, "label": con["label"], "kind": "between_donor" if between else "within_donor",
                           "numerator": con["labels"] if between else con["decidual"], "denominator": con["denominator"] if between else con["reference"],
                           "numerator_value": num, "denominator_value": den, "unit": "donors" if between else "donor pairs"}
        out_meta = meta.copy(); RP.write_tsv(out_meta, inter / f"ulrich_{cid}_metadata.tsv")
        out_c = counts.copy(); out_c.insert(0, "pseudobulk_id", meta["pseudobulk_id"]); RP.write_tsv(out_c, inter / f"ulrich_{cid}_counts.tsv")
        RP.write_csv(result, staging / f"ulrich_{cid}_all_genes.csv")
        result.to_parquet(staging / f"ulrich_{cid}_all_genes.parquet", index=False)
        human_results[cid] = result
        print(f"{cid}: {obs_c}", flush=True)

    # ---- consensus support rows per contrast (upstream builder, the released 123-gene list and its human mapping)
    primary = human["primary"]
    support_h = {}
    for cid, result in human_results.items():
        s = RP.build_heca_support(res["consensus123"], res["consensus123_human_mapping"], result, calls)
        s["source_id"] = f"ulrich_{cid}"
        support_h[cid] = s
        RP.write_csv(s, staging / ("consensus123_ulrich_support_rows.csv" if cid == primary else f"consensus123_ulrich_{cid}_support_rows.csv"))
    cats = ["supported", "reversed", "tested_not_called", "ineligible_or_unavailable"]
    for cid, s in support_h.items():
        name = f"consensus123_vs_ulrich_{cid}"
        observed = RP.normalized_partition(s["support_status"], cats)
        if set(s["support_status"].astype(str)) - set(cats):
            fail(f"{name}: unexpected status values")
        UP.check_equal(checks, f"partition_{name}", observed, (acc.get("status_partitions") or {}).get(name), record)
        observed_acc["status_partitions"][name] = observed
    summary = pd.DataFrame([RP.comparison_summary_row(f"consensus123_vs_ulrich_{cid}", f"123 human consensus IVD genes vs Ulrich {human_diag[cid]['label']}",
                                                      s["support_status"], "ineligible_or_unavailable") for cid, s in support_h.items()])
    summary.insert(1, "is_primary_human", summary["comparison_id"].eq(f"consensus123_vs_ulrich_{primary}"))
    RP.write_csv(summary, staging / "ivd_support_summary.csv")
    contrasts = pd.DataFrame([{"contrast_id": cid, "primary": cid == primary, **{k: (";".join(v) if isinstance(v, list) else v) for k, v in d.items()
                                                                                 if k not in ("design_columns",)}} for cid, d in human_diag.items()])
    RP.write_csv(contrasts, staging / "ulrich_contrasts.csv")

    RP.write_csv(pd.DataFrame(checks), staging / "validation_checks.csv")
    if {k: sha256_file(p) for k, p in execution_inputs.items()} != exec_hash:
        fail("Config or code changed during the run")
    run_meta = {"analysis_id": config["analysis_id"], "created_utc": pd.Timestamp.utcnow().isoformat(), "mode": "record" if record else "verify",
                "arm": "human_only", "python_executable": sys.executable, "python_version": platform.python_version(),
                "packages": RP.package_versions(), "n_cpus": args.n_cpus, "execution_sha256": exec_hash,
                "primary_human_contrast": primary, "human_fits": human_diag,
                "support_rule": "eligible and estimable and padj<0.05 and abs(unshrunk log2FC)>1; supported requires same direction"}
    (staging / "run_metadata.json").write_text(json.dumps(run_meta, indent=2, sort_keys=True, default=str) + "\n")
    if record:
        (staging / "acceptance_observed.yaml").write_text(yaml.safe_dump(observed_acc, sort_keys=False))
    if sha256_file(ulrich_obj) != config["objects"]["ulrich_cells"]["sha256"]:
        fail(f"Source object changed during the run: {ulrich_obj}")
    UP.write_manifest(staging, config, {"ulrich_cells": ulrich_obj}, exec_hash)
    staging.rename(output)
    print(json.dumps({"status": "complete_immutable", "mode": run_meta["mode"], "output": str(output),
                      "summary": summary[["comparison_id", "supported_n", "reversed_n", "denominator_n"]].to_dict("records")}, indent=1))


if __name__ == "__main__":
    main()
