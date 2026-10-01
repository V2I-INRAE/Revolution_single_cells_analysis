#!/bin/bash
#SBATCH -J norm-feat-pipeline
#SBATCH -p workq
#SBATCH -c 6
#SBATCH --mem=256G
#SBATCH -t 05:00:00
#SBATCH -o logs/norm-feat-pipeline-%j.out

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/sbatch_norm_feat.sh <clean_input.rds>" >&2
  exit 1
fi

# Use the submission directory: Slurm runs a spool copy of this script.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

# Keep timestamped logs without overwriting another job's output.
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-norm-feat-pipeline-${SLURM_JOB_ID:-$$}.log" 2>&1 || exit 1
set +o noclobber

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/norm_feat/main.R "$@"
