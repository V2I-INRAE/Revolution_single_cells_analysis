#!/bin/bash
#SBATCH -J revo-pipeline
#SBATCH -p workq
#SBATCH -c 8
#SBATCH --mem=256G
#SBATCH -t 72:00:00
#SBATCH -o logs/pipeline-%j.out

# Usage: sbatch scripts/sbatch_pipeline.sh
# Runs the whole pipeline (QC -> normalization -> integration -> clustering)
# exactly as configured in R/config.R. No arguments: edit config.R to change
# input source, samples, integration method or any threshold.
# Slurm-level messages (before the redirect) go to logs/pipeline-<jobid>.out.

# main.R uses repo-relative paths (data/, R/, results/), so run from the
# analysis root. $0 cannot be trusted here: sbatch executes the script
# from an internal spool copy, so use the submission directory instead.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

# Everything below goes to the timestamped log (R startup, renv activation
# messages, pipeline output and warnings included).
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-pipeline-${SLURM_JOB_ID:-$$}.log" 2>&1 || exit 1
set +o noclobber

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

# Doublet detection is scDblFinder (Bioconductor, R-native): the QC stage no
# longer needs a Python environment or reticulate.
# Thread policy: everything is pinned to one thread. scDblFinder's xgboost
# step respects OMP, so one thread keeps the run deterministic and
# reproducible.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

Rscript main.R
