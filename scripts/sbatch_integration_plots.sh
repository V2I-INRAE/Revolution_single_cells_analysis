#!/bin/bash
#SBATCH -J integration-plots
#SBATCH -p workq
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/integration-plots-%j.out

set -euo pipefail
if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: sbatch scripts/sbatch_integration_plots.sh <scvi.rds> <harmony.rds> <cca.rds> <new_output_dir> [--embedding=umap|tsne|both]" >&2
  exit 1
fi
if [[ $# -eq 5 ]]; then
  case "$5" in
    --embedding=umap|--embedding=tsne|--embedding=both) ;;
    *) echo "Unknown embedding option: $5" >&2; exit 1 ;;
  esac
fi

cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-integration-plots-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

Rscript R/integration/plot_main.R "$@"
