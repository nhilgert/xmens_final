#!/bin/bash
#SBATCH -J fibmap
#SBATCH -p eddy
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 0-04:00
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibmap_%j.out
#SBATCH -e /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibmap_%j.err

source /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/activate
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping
python -u run_mapping.py
