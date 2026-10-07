#!/usr/bin/env bash
#SBATCH --job-name=cellbender-cpu
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cellbender-cpu-%A_%a.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export NUMEXPR_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
export PATH="$PWD/.venv-scvi/bin:$PATH"
# Use this for a single pilot, or after QC with a reviewed TSV and --array=0-21%8.
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  [[ $# -eq 1 ]] || { echo "Array usage: sbatch --array=...%8 scripts/sbatch_cellbender.sh <reviewed.tsv>" >&2; exit 1; }
  read -r sample < <(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$1")
  [[ -n "$sample" ]] || exit 1
  .venv-scvi/bin/python scripts/cellbender/run_cellbender.py "$sample"
else
  .venv-scvi/bin/python scripts/cellbender/run_cellbender.py "$@"
fi
