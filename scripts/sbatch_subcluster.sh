#!/bin/bash
#SBATCH -J subcluster
#SBATCH -p workq
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o logs/subcluster-%j.out

# Full run (both overrides are optional): sbatch scripts/sbatch_subcluster.sh \
#   data/job-44542861/annotated_clustered.rds "T cells & NK" T_cell_nk --hvg=1000 --dims=20
# --hvg=N selects N HVGs; normalized RNA is reused, not normalized again.
# --dims=N uses components 1:N for integration, both UMAPs and clustering; computed PCs are unchanged.
# Without overrides, use R/config.R defaults for this job.
# PCA only: sbatch scripts/sbatch_subcluster.sh \
#   data/job-44542861/annotated_clustered.rds "T cells & NK" T_cell_nk --only-pca
# --only-pca saves the subset RDS and PCA diagnostics with configured defaults;
# it cannot be combined with --hvg, --dims or any other optional argument.
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
