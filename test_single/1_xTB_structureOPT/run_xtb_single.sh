#!/bin/bash
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=xtb_opt
#SBATCH --output=xtb_%j.out
#SBATCH --error=xtb_%j.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR" || exit 1
echo "[INFO] submit_dir = $SLURM_SUBMIT_DIR"
echo "[INFO] pwd        = $(pwd)"

# --- conda (non-interactive safe) ---
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "[ERROR] conda.sh not found under miniconda3/anaconda3"
  exit 2
fi

conda activate xtb

echo "[INFO] xtb = $(which xtb)"
xtb --version

# --- input ---
INPUT="${1:-}"

if [ -z "$INPUT" ]; then
  echo "[ERROR] Usage: sbatch run_xtb_single.sh /full/path/to/file.xyz"
  exit 3
fi

if [ ! -f "$INPUT" ]; then
  echo "[ERROR] File not found: $INPUT"
  exit 4
fi

# --- output ---
OUT="$SLURM_SUBMIT_DIR/runs"
mkdir -p "$OUT"

name=$(basename "$INPUT" .xyz)
outdir="$OUT/$name"
mkdir -p "$outdir"

echo "[INFO] INPUT  = $INPUT"
echo "[INFO] OUTDIR = $outdir"

echo "[INFO] running $name"

(
  cd "$outdir" || exit 1
  xtb "$INPUT" --gfn 2 --opt > xtb.log 2>&1
)

echo "[INFO] DONE"