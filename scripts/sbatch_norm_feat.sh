#!/bin/bash
#SBATCH -J norm-feat-pipeline
#SBATCH -p workq
#SBATCH -c 6
#SBATCH --mem=128G
#SBATCH -t 02:30:00
#SBATCH -o logs/norm-feat-pipeline-%j.out

# Use the submission directory: Slurm runs a spool copy of this script.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

# Match the QC pipeline's timestamped log naming convention.
exec > "logs/$(date +%Y%m%d-%H%M)-norm-feat-pipeline.log" 2>&1

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/norm_feat/main.R
