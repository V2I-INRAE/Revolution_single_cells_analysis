#!/usr/bin/env bash
# Start the local Shiny app on a Genotoul compute node.
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd -- "$project_dir"

if ! type module >/dev/null 2>&1; then
  source /etc/profile.d/modules.sh
fi
module purge
module load compilers/gcc/12.2.0
module load statistics/R/4.6.1

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export RENV_CONFIG_NAMESPACES_CHECK=FALSE

printf 'Compute node: %s\nApp: 127.0.0.1:3838 (SSH tunnel required from your computer)\n' "$(hostname -s)"
exec Rscript app/run.R
