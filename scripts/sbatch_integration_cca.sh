#!/bin/bash
#SBATCH -J integration-cca
#SBATCH -p workq
#SBATCH -c 8
#SBATCH --mem=512G
#SBATCH -t 72:00:00
#SBATCH -o logs/integration-cca-%j.out

set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch scripts/sbatch_integration_cca.sh <lognorm_input.rds>" >&2
  exit 1
fi

# Use the submission directory: Slurm runs a spool copy of this script.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"

# Keep timestamped logs, separating methods/jobs without truncating earlier logs.
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-integration-cca-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/integration/main.R cca "$1"
