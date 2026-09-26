#!/bin/bash
#SBATCH -J qc-pipeline
#SBATCH -p workq
#SBATCH -c 6
#SBATCH --mem=64G
#SBATCH -t 02:30:00
#SBATCH -o logs/qc-pipeline-%j.out

# Run the QC pipeline (R/qc/main.R) as a batch job. All output goes to
# logs/<YYYYMMDD>-<HHMM>-qc-pipeline.log, same timestamp format as the
# previous manual runs. Slurm-level messages (before the redirect) go to
# logs/qc-pipeline-<jobid>.out.

# main.R uses repo-relative paths (data/, R/, results/), so run from the
# analysis root. $0 cannot be trusted here: sbatch executes the script
# from an internal spool copy, so use the submission directory instead.
cd "${SLURM_SUBMIT_DIR:-/work/project/revo-pig-sc/analysis}" || exit 1

# everything below goes to the timestamped log (R startup, renv
# activation messages, pipeline output and warnings included). If this
# redirect fails, the error lands in logs/qc-pipeline-<jobid>.out.
exec > "logs/$(date +%Y%m%d-%H%M)-qc-pipeline.log" 2>&1

module purge
module load statistics/R/4.6.1

Rscript R/qc/main.R
