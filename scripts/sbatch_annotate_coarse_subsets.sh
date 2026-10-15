#!/bin/bash
#SBATCH -J annotate-subsets
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/annotate-subsets-%j.out

set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_annotate_coarse_subsets.sh <clustered.rds> <annotated_clustered.rds>" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-annotate-subsets-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/clusters/annotate.R "$1" "$2"
