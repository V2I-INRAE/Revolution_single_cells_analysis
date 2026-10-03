#!/bin/bash
#SBATCH -J integration-plots
#SBATCH -p workq
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o logs/integration-plots-%j.out

set -e
if [[ $# -ne 3 || ! "$1" =~ ^[0-9]+$ || ! "$2" =~ ^[0-9]+$ || ! "$3" =~ ^[0-9]+$ ]]; then
  echo "Usage: sbatch scripts/sbatch_integration_plots.sh <scvi_job> <harmony_job> <cca_job>" >&2
  exit 1
fi
scvi_job="$1"
harmony_job="$2"
cca_job="$3"

cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
exec > "logs/$(date +%Y%m%d-%H%M%S)-integration-plots-${SLURM_JOB_ID:-$$}.log" 2>&1

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/integration/plot_main.R \
  "data/integration/scvi/job-${scvi_job}/lognorm.rds" \
  "data/integration/harmony/job-${harmony_job}/lognorm.rds" \
  "data/integration/cca/job-${cca_job}/lognorm.rds" \
  "results/integration/comparison/scvi-job-${scvi_job}_harmony-job-${harmony_job}_cca-job-${cca_job}"
