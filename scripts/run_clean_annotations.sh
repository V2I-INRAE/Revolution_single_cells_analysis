#!/bin/bash
# Local use: bash scripts/run_clean_annotations.sh
# Requires Rscript and the project's R packages in your local environment.
set -euo pipefail
cd "$(dirname "$0")/.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec Rscript R/clusters/clean_annotations.R "$@"
