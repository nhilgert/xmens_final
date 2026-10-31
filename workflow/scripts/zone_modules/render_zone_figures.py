#!/usr/bin/env python
"""Render the configured decidual-zone figures for one score run (rule zone_figures, workflow/rules/16_zone_modules.smk).

Runs fig_zone_modules.py once per entry of --jobs (config zone_modules.figures: [{name, args}]) on the run's plot table
(--plot-table, the only data input) into --out, then writes MANIFEST.json there: every file with its md5 and size, the
invocations, and the inputs with md5. Written 2026-10-24. --out must not exist or must be empty (Snakemake removes a
stale directory output before re-running the rule). --section-rotation (added 2026-10-30) is passed on to every job:
the spatial panels take their display rotation from that table (resources/figures/mouse_section_rotation.yaml).
"""
import argparse, datetime, hashlib, json, os, shlex, subprocess, sys


def md5(path):
    h = hashlib.md5(); f = open(path, "rb")
    for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--script", required=True)      # fig_zone_modules.py
    ap.add_argument("--plot-table", required=True)  # zone_plot_table.csv (+ its .provenance.json next to it)
    ap.add_argument("--jobs", required=True)        # JSON list of {"name": str, "args": str}
    ap.add_argument("--out", required=True)
    ap.add_argument("--section-rotation", default="")   # mouse section rotation table, passed on to every job
    a = ap.parse_args()
    if os.path.exists(a.out) and os.listdir(a.out):
        sys.exit(f"refusing to write into non-empty {a.out}")
    side = a.plot_table[:-len(".csv")] + ".provenance.json"
    assert a.plot_table.endswith(".csv") and os.path.exists(side), f"no sidecar {side}"
    jobs = json.loads(a.jobs)
    names = [j["name"] for j in jobs]
    assert len(set(names)) == len(names), f"duplicate figure job names: {names}"
    os.makedirs(a.out, exist_ok=True)
    shared = ["--plot-table", a.plot_table, "--out", a.out]
    if a.section_rotation:
        shared += ["--section-rotation", a.section_rotation]
    runs = []
    for j in jobs:
        cmd = [a.python, a.script] + shared + shlex.split(j["args"])
        print("RUN", j["name"], ":", " ".join(shlex.quote(x) for x in cmd), flush=True)
        r = subprocess.run(cmd, capture_output=True, text=True)
        sys.stdout.write(r.stdout); sys.stderr.write(r.stderr)
        if r.returncode != 0:
            sys.exit(f"figure job {j['name']} failed with exit code {r.returncode}")
        wrote = [x for line in r.stdout.splitlines() if line.startswith("wrote ") for x in line.split()[1:]]
        runs.append({"name": j["name"], "args": j["args"], "files": [os.path.basename(x) for x in wrote]})
    files = sorted(f for f in os.listdir(a.out) if f != "MANIFEST.json")
    man = {"written": datetime.datetime.now().isoformat(timespec="seconds"), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
           "script": {"path": os.path.abspath(a.script), "md5": md5(a.script)},
           "driver": {"path": os.path.abspath(__file__), "md5": md5(os.path.abspath(__file__))},
           "inputs": {"plot_table": {"path": os.path.abspath(a.plot_table), "md5": md5(a.plot_table)},
                      "plot_table_provenance": {"path": os.path.abspath(side), "md5": md5(side)},
                      **({"section_rotation": {"path": os.path.abspath(a.section_rotation), "md5": md5(a.section_rotation)}}
                         if a.section_rotation else {})},
           "jobs": runs,
           "files": {f: {"md5": md5(os.path.join(a.out, f)), "bytes": os.path.getsize(os.path.join(a.out, f))} for f in files}}
    json.dump(man, open(os.path.join(a.out, "MANIFEST.json"), "w"), indent=1)
    print(f"{len(files)} files in {a.out}", flush=True)


if __name__ == "__main__":
    main()
