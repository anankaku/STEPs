#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --job-name=xtb_one
#SBATCH --output=logs/xtb_%j.out
#SBATCH --error=logs/xtb_%j.err

set -eo pipefail
cd "$SLURM_SUBMIT_DIR" || exit 1
mkdir -p logs

# conda（batch 安全）
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "[ERROR] conda.sh not found" >&2
  exit 2
fi
conda activate xtb

: "${XYZ:?XYZ not set}"
: "${OUTDIR:?OUTDIR not set}"

if [ ! -s "$XYZ" ]; then
  echo "[ERROR] xyz missing/empty: $XYZ" >&2
  exit 3
fi

name=$(basename "$XYZ" .xyz)
jobdir="$OUTDIR/$name"
mkdir -p "$jobdir"

echo "[INFO] xyz    = $XYZ"
echo "[INFO] jobdir = $jobdir"
echo "[INFO] xtb    = $(which xtb)"
xtb --version

xtb "$XYZ" --gfn 2 --opt > "$jobdir/xtb.log" 2>&1
echo "[INFO] DONE"
