#!/bin/bash
#SBATCH -J fibnb
#SBATCH -p eddy
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 0-02:00
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibnb_%j.out
#SBATCH -e /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibnb_%j.err

source /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/activate
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=1800 \
    hammoud_fib_centroid_mapping.ipynb
echo "NBCONVERT EXIT: $?"
