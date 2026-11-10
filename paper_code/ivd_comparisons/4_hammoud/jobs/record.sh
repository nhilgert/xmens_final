#!/bin/bash
#SBATCH --account=eddy_lab
#SBATCH --partition=eddy
#SBATCH --exclude=holygpu7c0522,holygpu7c0530,holygpu7c0920,holygpu8a26504,holygpu8a26505,holygpu8a26506
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1:00:00
set -euo pipefail
export XDG_CACHE_HOME=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch
export NUMBA_CACHE_DIR=$XDG_CACHE_HOME/numba MPLCONFIGDIR=$XDG_CACHE_HOME/mpl TMPDIR=/n/lab_storage/eddy_lab/users/nhilgert/tmp HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $XDG_CACHE_HOME $NUMBA_CACHE_DIR $MPLCONFIGDIR $TMPDIR
mkdir -p resources && mv ulrich_feature_manifest.tsv consensus123_ulrich_mapping.csv resources/
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111/bin/python ulrich_ivd_pipeline.py --config config_ulrich_dS_vs_C7fib_luteal_follicular.yaml --record --n-cpus 8
cp /n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular_record/acceptance_observed.yaml /n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular_record/ivd_support_summary.csv /n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular_record/ulrich_contrasts.csv /n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular_record/run_metadata.json /n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular_record/validation_checks.csv .
mkdir -p donors && cp /n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular_record/intermediate/ulrich_*_metadata.tsv donors/
tar czf release_tables.tgz -C /n/netscratch/eddy_lab/Everyone/nhilgert/cagri_ivd_ulrich_deseq2_runs/release_20261028_v7_dS_vs_C7fib_luteal_follicular_record --exclude="*.h5ad" --exclude="*.parquet" --exclude=intermediate .
