#!/usr/bin/env bash
#SBATCH --job-name=soupx
#SBATCH --partition=workq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/soupx-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
[[ $# -eq 1 ]] || { echo "Usage: sbatch scripts/sbatch_soupx.sh <sample>; with --array: <sample-list>" >&2; exit 1; }
if [[ -n ${SLURM_ARRAY_TASK_ID:-} ]]; then
  [[ -f $1 ]] || { echo "Sample list not found: $1" >&2; exit 1; }
  sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$1")
  SOUPX_SAMPLE_MANIFEST=$(realpath "$1")
  export SOUPX_SAMPLE_MANIFEST
else
  sample=$1
fi
[[ -n $sample && -d data/raw_data/$sample ]] || { echo "Sample not found: $sample" >&2; exit 1; }
echo "SoupX sample: $sample; array task: ${SLURM_ARRAY_TASK_ID:-none}"
module purge
module load statistics/R/4.6.1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
Rscript scripts/soupx/run_soupx.R "$sample"
