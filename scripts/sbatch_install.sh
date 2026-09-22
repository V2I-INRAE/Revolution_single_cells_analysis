#!/bin/bash
#SBATCH -J renv-install
#SBATCH -p workq
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o logs/renv_install_%j.log

# Install the single-cell R stack (Seurat v5 + skill packages) into the
# project renv library, per Genotoul doc: load R module, run Rscript.

module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

Rscript scripts/install_packages.R
