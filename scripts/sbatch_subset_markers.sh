#!/bin/bash
#SBATCH -J subset-markers
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o logs/subset-markers-%j.out

# One Slurm job per subset, resolution passed as argument:
# sbatch scripts/sbatch_subset_markers.sh \
#   data/job-44542861/subsets/B_cells/job-44614229/subclustered.rds 0.7
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_subset_markers.sh <data/job-XXX/subsets/<Folder>/<subrun>/subclustered.rds> <resolution>" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-subset-markers-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/clusters/subset_markers.R "$1" "$2"
