#!/usr/bin/env bash
#SBATCH --job-name=cellbender-preparation
#SBATCH --partition=workq
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cellbender-preparation-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
.venv-scvi/bin/python scripts/cellbender/prepare_cellbender.py "$@"
