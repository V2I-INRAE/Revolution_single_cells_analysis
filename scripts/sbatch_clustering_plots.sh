#!/bin/bash
#SBATCH -J clustering-plots
#SBATCH -p workq
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/clustering-plots-%j.out

set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_clustering_plots.sh <new_output_dir> <analysis.rds> [<analysis.rds> ...]" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-clustering-plots-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/clustering/plot_main.R "$@"
