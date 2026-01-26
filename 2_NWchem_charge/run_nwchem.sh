#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=nwchem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x.out
#SBATCH --error=%x.err

set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: sbatch run_nwchem.sh path/to/job.nw"
  exit 1
fi

NWFILE="$1"
if [ ! -f "${NWFILE}" ]; then
  echo "ERROR: NWChem input not found: ${NWFILE}"
  exit 1
fi

# ===== 環境（乾淨版）=====
module --force purge

# conda 初始化時暫時關掉 -u，避免 /etc/bashrc 類問題
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nwchem_env
set -u

WORKDIR="$(dirname "${NWFILE}")"
BASENAME="$(basename "${NWFILE}" .nw)"
cd "${WORKDIR}"

echo "Running NWChem on ${BASENAME}.nw"
echo "Start time: $(date)"
nwchem "${BASENAME}.nw" > "${BASENAME}.out"
echo "End time: $(date)"
