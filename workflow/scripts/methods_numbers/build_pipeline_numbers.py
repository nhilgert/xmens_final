#!/usr/bin/env python
"""Build pipeline_numbers.csv: every pipeline-derived number in the computational methods, recomputed from the
current outputs, with the file and step that produce it and the words around it in the draft.

Usage:
  python build_pipeline_numbers.py --facts DIR --facts2 DIR --small-reads JSON --artifact-inputs DIR --proteome-check TXT \\
                                   --draft TXT --out CSV [--report MD] [--extra-inputs DIR]

--facts / --facts2 / --small-reads are written by collect_facts.py, collect_facts2.py and small_reads.py (cluster jobs);
--artifact-inputs holds inputs kept only as project artifacts (MANIFEST.json lists ids and sha256), including the corrected
in vitro consensus (2026-10-22 Stadtmauer correction); --proteome-check is the output of check_proteome_release.sh; --draft is the methods text
(extracted from the docx/PDF). --extra-inputs holds further artifact-held inputs staged outside the repository (with
their own MANIFEST.json): the SAMap pseudotime permutation-null statistics and the X-Mens fit with the uninduced
fibroblasts; without it, the rows that read them are no_source. The CSV is written fresh (refuses to overwrite).

Columns: key, value, unit, description, source_file, source_step, doc_phrase, updated, then
value_raw (unrounded), draft_value (as written in the draft), status, note.
status: match | MISMATCH (recomputed value differs from the draft) | no_source (not recomputed; note says why) |
phrase_found_Nx / value_not_in_phrase (the draft text moved; update doc_phrase).
"""
import os, sys, argparse, json, time
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import registry_core as rc
import registry_rows as rr
import registry_methods as rm

ap = argparse.ArgumentParser()
ap.add_argument("--facts", required=True); ap.add_argument("--facts2", required=True); ap.add_argument("--small-reads", required=True)
ap.add_argument("--artifact-inputs", required=True, help="directory with MANIFEST.json and the artifact-held inputs")
ap.add_argument("--proteome-check", required=True, help="stdout of check_proteome_release.sh")
ap.add_argument("--draft", required=True); ap.add_argument("--out", required=True); ap.add_argument("--report", default=None)
ap.add_argument("--extra-inputs", default=None, help="directory with MANIFEST.json and artifact-held inputs staged outside the repository")
a = ap.parse_args()
assert not os.path.exists(a.out), f"refusing to overwrite {a.out}"

draft = rc.normalize_draft(open(a.draft, encoding="utf-8").read())
FIXZ = os.path.join(a.artifact_inputs, "consensus123_stadtmauer_fix_bundle.zip")
FIXL = "artifact consensus123_stadtmauer_fix_bundle.zip (version %s), job_outputs" % json.load(open(os.path.join(a.artifact_inputs, "MANIFEST.json")))["consensus123_stadtmauer_fix_bundle.zip"]["latest_version_id"]
ctx = rm.Ctx2(a.facts, a.facts2, a.small_reads, FIXZ, FIXL, a.artifact_inputs, a.proteome_check, extra_inputs=a.extra_inputs)
rows = rr.rows(); ov = rm.overrides()
for r in rows:
    if r.key in ov: r.fn = ov[r.key]
res = rc.evaluate(rows, ctx, draft)

T = pd.DataFrame(res)
# a note comparing with superseded runs only matters where the draft disagrees
T.loc[(T.status == "match") & T.note.str.startswith("legacy runs"), "note"] = ""
T["value_raw"] = T["value_raw"].fillna("") if "value_raw" in T else ""      # set by registry_core.raw_value
T["value"] = T.display_value.fillna("")
cols = ["key", "value", "unit", "description", "source_file", "source_step", "doc_phrase", "updated", "value_raw", "draft_value", "status", "note"]
T = T[cols]
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
T.to_csv(a.out, index=False)
print(f"{len(T)} rows -> {a.out} | " + ", ".join(f"{k} {v}" for k, v in T.status.value_counts().items()), flush=True)

if a.report:
    assert not os.path.exists(a.report), f"refusing to overwrite {a.report}"
    L = [f"# pipeline_numbers: flagged rows ({time.strftime('%Y-%m-%d %H:%M')})", "",
         f"{len(T)} numbers; " + ", ".join(f"{v} {k}" for k, v in T.status.value_counts().items()) + ".", ""]
    for st in ("MISMATCH", "no_source"):
        S = T[T.status == st]
        if S.empty: continue
        L += [f"## {st} ({len(S)})", "", "| key | draft | recomputed | phrase | note |", "|---|---|---|---|---|"]
        for _, r in S.iterrows():
            L.append(f"| {r.key} | {r.draft_value} | {r.value} | {r.doc_phrase[:70]} | {str(r.note)[:160]} |".replace("\n", " "))
        L.append("")
    open(a.report, "w").write("\n".join(L))
