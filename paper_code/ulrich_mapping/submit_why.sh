#!/bin/bash
#SBATCH -J fibwhy
#SBATCH -p eddy
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 0-04:00
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibwhy_%j.out
#SBATCH -e /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibwhy_%j.err

source /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/activate
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping
python -u why_no_match.py
