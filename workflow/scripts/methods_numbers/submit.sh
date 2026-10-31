#!/bin/bash
# Rebuild pipeline_numbers.csv from the current pipeline outputs.
#
# Usage (login node):  bash workflow/scripts/methods_numbers/submit.sh RUN_DIR DRAFT_TXT [EXTRA_INPUTS]
#   RUN_DIR       new directory for this run (refuses to reuse one), e.g. results/methods_numbers/2026_11_02
#   DRAFT_TXT     the methods text as plain text (e.g. extracted from the docx); copied into RUN_DIR
#   EXTRA_INPUTS  optional directory of artifact-held inputs staged outside the repository, with a MANIFEST.json
#                 (xmens_with_uninduced_20261027.zip); without it those rows are no_source
# Submits three read-only collection jobs and a build job that starts when all three succeed.
# Outputs: RUN_DIR/pipeline_numbers.csv and RUN_DIR/pipeline_numbers_flags.md
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
RES="$REPO/resources/methods_numbers"
RUN=$(realpath -m "$1"); DRAFT=$2; EXTRA=${3:-}
case "$RUN" in /n/home12/*) echo "refusing to write under home"; exit 3;; esac
test ! -e "$RUN" || { echo "$RUN exists; choose a new run directory"; exit 2; }
test -f "$DRAFT" || { echo "draft text $DRAFT not found"; exit 2; }
EXTRA_ARG=""
if [ -n "$EXTRA" ]; then
  test -f "$EXTRA/MANIFEST.json" || { echo "extra inputs $EXTRA have no MANIFEST.json"; exit 2; }
  EXTRA_ARG="--extra-inputs $(realpath "$EXTRA")"
fi
mkdir -p "$RUN"; cp "$DRAFT" "$RUN/draft.txt"
PY=/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/python
EXCL=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
COMMON="--partition=eddy --account=eddy_lab --exclude=$EXCL --chdir=$RUN --cpus-per-task=2"
# Lab cache locations for matplotlib and numba; HDF5 file locking off for reads from lab_storage and netscratch.
CACHE=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch
ENVSET="export PYTHONDONTWRITEBYTECODE=1 HDF5_USE_FILE_LOCKING=FALSE TMPDIR=$RUN/tmp MPLCONFIGDIR=$CACHE/mpl XDG_CACHE_HOME=$RUN/.cache XDG_CONFIG_HOME=$RUN/.config NUMBA_CACHE_DIR=$CACHE/numba; mkdir -p $RUN/tmp $RUN/.cache $RUN/.config"
j1=$(sbatch --parsable $COMMON --mem=64G --time=0-00:45 --job-name=numFacts1 --output="$RUN/slurm-facts1-%j.out" \
     --wrap="$ENVSET; $PY -u $HERE/collect_facts.py $RUN/facts")
j2=$(sbatch --parsable $COMMON --mem=80G --time=0-01:00 --job-name=numFacts2 --output="$RUN/slurm-facts2-%j.out" \
     --wrap="$ENVSET; $PY -u $HERE/collect_facts2.py $RUN/facts2")
j3=$(sbatch --parsable $COMMON --mem=8G --time=0-00:30 --job-name=numReads --output="$RUN/slurm-reads-%j.out" \
     --wrap="$ENVSET; $PY -u $HERE/small_reads.py $RUN/small_reads.json")
j4=$(sbatch --parsable $COMMON --mem=8G --time=0-00:20 --job-name=numBuild --dependency=afterok:$j1:$j2:$j3 --output="$RUN/slurm-build-%j.out" \
     --wrap="$ENVSET; $PY -u $HERE/build_pipeline_numbers.py --facts $RUN/facts --facts2 $RUN/facts2 --small-reads $RUN/small_reads.json --artifact-inputs $RES/artifact_inputs --proteome-check $RES/proteome_release_check.txt --draft $RUN/draft.txt --out $RUN/pipeline_numbers.csv --report $RUN/pipeline_numbers_flags.md $EXTRA_ARG")
echo "collection jobs $j1 $j2 $j3; build job $j4 (runs after all three succeed)"
echo "outputs: $RUN/pipeline_numbers.csv, $RUN/pipeline_numbers_flags.md"
