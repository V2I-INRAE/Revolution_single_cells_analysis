#!/bin/bash
#SBATCH -J azimuth-parent
#SBATCH -p workq
#SBATCH -c 4
#SBATCH --mem=256G
#SBATCH -t 12:00:00
#SBATCH -o logs/azimuth-parent-%j.log

# Run the existing prediction script unchanged on all cleaned parent cells.
# Usage: sbatch scripts/sbatch_azimuth_parent.sh
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export TF_NUM_INTRAOP_THREADS="${SLURM_CPUS_PER_TASK:-4}" TF_NUM_INTEROP_THREADS=1

Rscript R/annotations/azimuth.R --input=data/job-44542861/clean_annotated_object.rds
