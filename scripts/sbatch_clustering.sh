#!/bin/bash
#SBATCH -J clustering
#SBATCH -p workq
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/clustering-%j.out

set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_clustering.sh <unintegrated|harmony|cca> <input.rds>" >&2
  exit 1
fi
method="${1,,}"
case "$method" in
  unintegrated|harmony|cca) ;;
  *) echo "Unknown clustering method: $1" >&2; exit 1 ;;
esac
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-clustering-${method}-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/clustering/main.R "$method" "$2"
