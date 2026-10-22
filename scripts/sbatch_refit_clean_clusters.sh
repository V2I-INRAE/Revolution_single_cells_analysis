#!/bin/bash
#SBATCH -J refit-clean
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=256G
#SBATCH -t 72:00:00
#SBATCH -o logs/refit-clean-%j.out

# Self-contained refit of the cleaned global object:
#   sbatch scripts/sbatch_refit_clean_clusters.sh
# Input (read-only): data/job-44542861/clean_annotated_object.rds
# Output: data/job-44542861/refit_clean_annotated.rds + .INFO, figures under
# results/job-44542861/clusters_plots/refit_clean_annotated/.
# Refuses to run if the output RDS, INFO, tmp file or figure directory exists.
# Resources: oracle-recommended 256 GB envelope; parent peak is unmeasured.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-refit-clean-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript R/clusters/refit_clean_clusters.R
