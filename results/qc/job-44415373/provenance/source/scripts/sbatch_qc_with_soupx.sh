#!/bin/bash
#SBATCH -J qc-with-soupx
#SBATCH -p workq
#SBATCH -c 6
#SBATCH --mem=128G
#SBATCH -t 05:00:00
#SBATCH -o logs/qc-with-soupx-%j.out

# Sequential QC on the saved SoupX-corrected samples, Scrublet only.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-qc-with-soupx-${SLURM_JOB_ID:-$$}.log" 2>&1 || exit 1
set +o noclobber

module purge
module load statistics/R/4.6.1

source .venv-scrublet/bin/activate || exit 1
export RETICULATE_PYTHON="$VIRTUAL_ENV/bin/python"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"

Rscript R/qc/main_soupx.R
