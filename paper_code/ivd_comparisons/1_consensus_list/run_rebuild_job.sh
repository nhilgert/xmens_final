#SBATCH --account=eddy_lab
#SBATCH --partition=eddy
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
set -euo pipefail
export HDF5_USE_FILE_LOCKING=FALSE XDG_CACHE_HOME=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch MPLCONFIGDIR=/n/lab_storage/eddy_lab/users/nhilgert/.cache_scratch/mpl
mkdir -p out
echo "49bf53d4537243d7ca5184fd94f20032bb2c60697fdcdeb1db182251f1bcdee9  rebuild_consensus123_stadtmauer_thresholds.py" | sha256sum -c -
/n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111/bin/python rebuild_consensus123_stadtmauer_thresholds.py --ivd-repo /n/lab_storage/eddy_lab/users/nhilgert/human_ivd_repro_final_defense_20260812 --pipeline-repo /n/lab_storage/eddy_lab/users/nhilgert/cagri_ivd_end_to_end_20260815 --xmens xmens_broad_all_genes.csv --released-support consensus_ortholog_xmens_broad_support_rows.csv --out out 2>&1 | tee out/run_log.txt
cd out && /n/lab_storage/eddy_lab/users/nhilgert/.uv/envs/rsc_uv_py31111/bin/python -c "import shutil; shutil.make_archive('../stadtmauer_fix_outputs', 'zip', '.')"