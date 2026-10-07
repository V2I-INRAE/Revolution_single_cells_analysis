#!/bin/bash
#SBATCH -J qc-soupx-verify
#SBATCH -p workq
#SBATCH -c 2
#SBATCH --mem=96G
#SBATCH -t 01:00:00
#SBATCH -o logs/qc-soupx-verify-%j.log
set -e
module purge
module load statistics/R/4.6.1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export MKL_NUM_THREADS="$SLURM_CPUS_PER_TASK"
Rscript .amp/in/qc-mad3-run/verify.R
