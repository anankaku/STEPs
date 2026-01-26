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

# --- input / output paths ---
SRC="$SLURM_SUBMIT_DIR/../0_openBable_SMILES/out/xyz"
OUT="$SLURM_SUBMIT_DIR/runs"

echo "[INFO] SRC = $SRC"
echo "[INFO] OUT = $OUT"

shopt -s nullglob
xyz_files=( "$SRC"/*.xyz )
echo "[INFO] xyz count = ${#xyz_files[@]}"

if [ ${#xyz_files[@]} -eq 0 ]; then
  echo "[ERROR] No xyz found in $SRC"
  ls -la "$SRC" || true
  exit 3
fi

mkdir -p "$OUT"

for f in "${xyz_files[@]}"; do
  name=$(basename "$f" .xyz)
  outdir="$OUT/$name"
  mkdir -p "$outdir"

  echo "[INFO] running $name (input: $f)"
  (
    cd "$outdir" || exit 1
    xtb "$f" --gfn 2 --opt > xtb.log 2>&1
  )
done

echo "[INFO] DONE"
