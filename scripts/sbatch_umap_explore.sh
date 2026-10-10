#!/bin/bash
#SBATCH -J umap-explore
#SBATCH -p workq
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH -t 24:00:00
#SBATCH -o logs/umap-explore-%j.out

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: sbatch [--array=...] scripts/sbatch_umap_explore.sh <harmony_or_cca_checkpoint.rds> [test_name]" >&2
  echo "Array ranges: Harmony 0-13%2; CCA 0-14%2 (includes min.dist 0.5)." >&2
  exit 1
fi

cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1
job_tag="${SLURM_JOB_ID:-$$}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  job_tag="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
fi
export UMAP_LOG_FILE="logs/$(date +%Y%m%d-%H%M%S)-umap-explore-${job_tag}.log"
set -o noclobber
exec > "$UMAP_LOG_FILE" 2>&1 || exit 1
set +o noclobber

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/integration/umap_explore.R "$@"
