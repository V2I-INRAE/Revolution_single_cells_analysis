#!/bin/bash
#SBATCH -J integration-pipeline
#SBATCH -p workq
#SBATCH -c 4
#SBATCH --mem=128G
#SBATCH -t 06:00:00
#SBATCH -o logs/integration-pipeline-%j.out

# Use the submission directory: Slurm runs a spool copy of this script.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

# Match the existing pipelines' timestamped log naming convention.
exec > "logs/$(date +%Y%m%d-%H%M)-integration-pipeline.log" 2>&1

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript R/integration/main.R
