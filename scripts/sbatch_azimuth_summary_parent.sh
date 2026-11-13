#!/bin/bash
#SBATCH -J azimuth-summary-parent
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/azimuth-summary-parent-%j.log

# Run after prediction succeeds:
# sbatch --dependency=afterok:<prediction_job> scripts/sbatch_azimuth_summary_parent.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

Rscript R/annotations/azimuth_summary_parent.R
