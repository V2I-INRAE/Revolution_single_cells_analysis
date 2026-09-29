#!/bin/bash
#SBATCH -J doublet-monitor
#SBATCH -p workq
#SBATCH -c 6
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -o logs/doublet-monitor-%j.out

# Run the standalone all-sample doublet investigation from the analysis root.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

exec > "logs/$(date +%Y%m%d-%H%M)-doublet-monitor.log" 2>&1

module purge
module load statistics/R/4.6.1

source .venv-scrublet/bin/activate || exit 1
export RETICULATE_PYTHON="$VIRTUAL_ENV/bin/python"

Rscript R/doublets/monitor.R "$@"
