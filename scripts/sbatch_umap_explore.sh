#!/bin/bash
#SBATCH -J umap-explore
#SBATCH -p workq
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/umap-explore-%j.out

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/sbatch_umap_explore.sh <harmony_checkpoint.rds>" >&2
  exit 1
fi

cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1
export UMAP_LOG_FILE="logs/$(date +%Y%m%d-%H%M%S)-umap-explore-${SLURM_JOB_ID:-$$}.log"
set -o noclobber
exec > "$UMAP_LOG_FILE" 2>&1 || exit 1
set +o noclobber

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/integration/umap_explore.R "$@"
