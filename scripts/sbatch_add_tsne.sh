#!/bin/bash
#SBATCH -J add-tsne
#SBATCH -p workq
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/add-tsne-%j.out

set -euo pipefail
if [[ $# -ne 3 ]]; then
  echo "Usage: sbatch scripts/sbatch_add_tsne.sh <integration|clustering> <lognorm.rds|analysis.rds> <new_output_dir>" >&2
  exit 1
fi
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-add-tsne-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript R/integration/add_tsne.R "$@"
