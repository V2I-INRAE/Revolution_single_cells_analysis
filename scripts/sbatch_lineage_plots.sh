#!/bin/bash
#SBATCH -J lineage-plots
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -o logs/lineage-plots-%j.out

set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch scripts/sbatch_lineage_plots.sh <clustering/lognorm.rds> <resolution>" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
export LINEAGE_PLOTS_LOG="logs/$(date +%Y%m%d-%H%M%S)-lineage-plots-${SLURM_JOB_ID}.log"
set -o noclobber
exec > "$LINEAGE_PLOTS_LOG" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/find_markers/plot_features_main.R "$1" "$2" all
