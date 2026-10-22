#!/bin/bash
#SBATCH -J refit-markers
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/refit-markers-%j.out

# Marker discovery + dotplot/heatmap for the refitted cleaned object:
# sbatch scripts/sbatch_refit_markers.sh data/job-44890936/refit_clean_annotated.rds 0.4
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_refit_markers.sh <data/job-XXX/refit_clean_annotated.rds> <resolution>" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-refit-markers-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/clusters/refit_markers.R "$1" "$2"
