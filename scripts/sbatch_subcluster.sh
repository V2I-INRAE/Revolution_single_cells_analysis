#!/bin/bash
#SBATCH -J subcluster
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/subcluster-%j.out

# Example: sbatch scripts/sbatch_subcluster.sh \
#   data/job-44493014/annotated_clustered.rds "T cells & NK" T_cell_nk
set -euo pipefail
if [[ $# -ne 3 ]]; then
  echo 'Usage: sbatch scripts/sbatch_subcluster.sh <annotated.rds> "<subset label>" <folder>' >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-subcluster-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript R/clusters/subcluster.R "$1" "$2" "$3"
