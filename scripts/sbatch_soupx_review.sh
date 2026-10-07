#!/usr/bin/env bash
#SBATCH --job-name=soupx-review
#SBATCH --partition=workq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/soupx-review-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
module purge
module load statistics/R/4.6.1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
Rscript scripts/soupx/review_soupx.R "$@"
