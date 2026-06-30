"""STUB — replace with real logic for this stage.

For now it just creates the declared output file(s) so the full DAG runs end to
end. Implement incrementally: read `snakemake.input`, do the work, write
`snakemake.output`. Params live in `snakemake.params` / `snakemake.config`.
"""
from pathlib import Path

snk = snakemake  # noqa: F821  (injected by Snakemake)
for f in snk.output:
    p = Path(str(f)); p.parent.mkdir(parents=True, exist_ok=True); p.touch()
print("STUB:", snk.rule, "->", list(snk.output))
