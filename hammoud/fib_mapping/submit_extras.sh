#!/bin/bash
#SBATCH -J fibextra
#SBATCH -p eddy
#SBATCH --account=eddy_lab
#SBATCH -c 8
#SBATCH --mem=96G
#SBATCH -t 0-02:00
#SBATCH -o /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibextra_%j.out
#SBATCH -e /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping/fibextra_%j.err

source /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv/bin/activate
cd /n/lab_storage/eddy_lab/users/nhilgert/xmens_final/hammoud/fib_mapping
python -u run_mapping_mouse_extras.py
echo "EXTRAS EXIT: $?"
