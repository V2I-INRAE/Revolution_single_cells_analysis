#!/bin/bash
#SBATCH -J integration-pipeline
#SBATCH -p workq
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH -o logs/integration-pipeline-%j.out

if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/sbatch_integration.sh <harmony|cca|scvi>" >&2
  exit 1
fi
method="${1,,}"
case "$method" in
  harmony|cca|scvi) ;;
  *) echo "Unknown integration method: $1" >&2; exit 1 ;;
esac

# Use the submission directory: Slurm runs a spool copy of this script.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

# Keep timestamped logs, separating methods/jobs without truncating earlier logs.
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-integration-${method}-${SLURM_JOB_ID:-$$}.log" 2>&1 || exit 1
set +o noclobber

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/integration/main.R "$@"
