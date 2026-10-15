#!/bin/bash
#SBATCH -J group-distribution
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/group-distribution-%j.out

# Parent: sbatch scripts/sbatch_group_distribution.sh data/job-44542861/annotated_clustered.rds 0.7
# Subset: sbatch scripts/sbatch_group_distribution.sh data/job-44542861/subsets/B_cells/job-44614229/subclustered.rds <resolution>
# Resolution must already exist as harmony_snn_res.<resolution>; no reclustering.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_group_distribution.sh <input.rds> <resolution>" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-group-distribution-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript R/plots/run_group_distribution.R "$1" "$2"
