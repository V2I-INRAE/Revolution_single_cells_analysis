#!/bin/bash
#SBATCH -J find-markers
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/find-markers-%j.out

set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_find_markers.sh <clustering/analysis.rds> <0.1,0.2,...>" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-find-markers-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/find_markers/main.R "$@"
