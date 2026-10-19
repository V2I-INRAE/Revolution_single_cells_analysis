#!/bin/bash
#SBATCH -J clean-annotations
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/clean-annotations-%j.out

# From the project root:
#   mkdir -p logs && sbatch scripts/sbatch_clean_annotations.sh
# Uses the input/output paths defined in R/clusters/clean_annotations.R.
# Resources match the existing full-parent annotation job; peak RAM is unmeasured.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-clean-annotations-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
bash scripts/run_clean_annotations.sh "$@"
