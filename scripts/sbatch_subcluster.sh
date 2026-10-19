#!/bin/bash
#SBATCH -J subcluster
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o logs/subcluster-%j.out

# Full run (all overrides are optional): sbatch scripts/sbatch_subcluster.sh \
#   data/job-44542861/clean_annotated_object.rds DC DC \
#   --hvg=1500 --dims=10 --k-param=50 --resolutions=0.1,0.2,0.3,0.4,0.5,0.6,0.7
# --hvg=N selects N HVGs; normalized RNA is reused, not normalized again.
# --dims=N uses components 1:N for integration, both UMAPs and clustering; computed PCs are unchanged.
# --k-param=N sets the clustering neighbour count, not UMAP n.neighbors.
# --resolutions=... retains each listed partition; no optimal resolution is selected automatically.
# Overrides affect only this run; without them, use R/config.R defaults.
# Submit all approved cleaned compartments with: bash scripts/submit_clean_subclusters.sh
# PCA only: sbatch scripts/sbatch_subcluster.sh \
#   data/job-44542861/clean_annotated_object.rds DC DC --only-pca
# --only-pca saves the subset RDS and PCA diagnostics with configured defaults;
# it cannot be combined with --hvg, --dims, --k-param or --resolutions.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}"
set -o noclobber
exec > "logs/$(date +%Y%m%d-%H%M%S)-subcluster-${SLURM_JOB_ID:-$$}.log" 2>&1
set +o noclobber
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript R/clusters/subcluster.R "$@"
