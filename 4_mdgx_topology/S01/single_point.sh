#!/bin/bash
#SBATCH --job-name=g09_array
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --array=0-9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=normal

cd $SLURM_SUBMIT_DIR

module purge
module load gaussian/g09-C01

export GAUSS_SCRDIR=/scratch/$USER/$SLURM_JOB_ID
mkdir -p $GAUSS_SCRDIR

file=$(ls *.com | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")

if [ -z "$file" ]; then
  echo "No .com file for task $SLURM_ARRAY_TASK_ID"
  exit 1
fi

echo "Running $file"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"

g09 $file

echo "End: $(date)"

rm -rf $GAUSS_SCRDIR
