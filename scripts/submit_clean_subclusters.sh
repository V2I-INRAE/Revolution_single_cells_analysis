#!/bin/bash
# Submit the approved 11 cleaned compartments at k=100 plus DC at k=50.
# Run only when submission is authorized: bash scripts/submit_clean_subclusters.sh
# This is a submission script, not an sbatch job itself.
set -euo pipefail
cd "$(dirname "$0")/.."

input="data/job-44542861/clean_annotated_object.rds"
resolutions="0.1,0.2,0.3,0.4,0.5,0.6,0.7"
input_info="${input%.rds}.INFO"
if [[ ! -f "$input" || ! -f "$input_info" ]]; then
  printf '%s\n' "Cleaned input or its INFO is missing: $input" >&2
  exit 1
fi
if ! grep -qx 'Status: completed' "$input_info"; then
  printf '%s\n' "Cleaned input is not recorded as completed: $input_info" >&2
  exit 1
fi

mkdir -p logs
record="logs/$(date +%Y%m%d-%H%M%S)-clean-subclusters-submission-$$.INFO"
set -o noclobber
exec 3> "$record"
set +o noclobber
trap 'result=$?; printf "\nSubmission exit status: %s\nFinished: %s\n" "$result" "$(date --iso-8601=seconds)" >&3' EXIT
printf '%s\n' \
  '# Cleaned-compartment subclustering submission' \
  "Started: $(date --iso-8601=seconds)" \
  "Input: $input" \
  "Input record: $input_info; cleanup producer 44827588; pipeline ancestor 44542861." \
  'Method: reuse normalized RNA; refit subset features, scaling/regression, PCA, Harmony, neighbours, Leiden and diagnostic UMAPs.' \
  'Seed: 1234; k=100 baseline for all compartments; DC k=50 sensitivity comparison.' \
  "Resolutions: $resolutions; final reporting resolution remains unselected." \
  'Resources per job: workq; 1 CPU; 32 GB RAM; 2 hours.' \
  'Execution and saved-output validation: pending; see each job .INFO and original logs.' \
  'Job | Compartment | Folder | HVGs | Downstream dimensions | k_param' >&3

submit_subset() {
  local label="$1" folder="$2" hvg="$3" dims="$4" k="$5" job
  job=$(sbatch --parsable --job-name="sub-${folder}-k${k}" \
    scripts/sbatch_subcluster.sh "$input" "$label" "$folder" \
    --hvg="$hvg" --dims="$dims" --k-param="$k" --resolutions="$resolutions")
  printf '%s | %s | %s | %s | 1:%s | %s\n' \
    "$job" "$label" "$folder" "$hvg" "$dims" "$k" >&3
  printf 'Submitted %s: %s (k=%s)\n' "$job" "$label" "$k"
}

# Feature counts/dimensions come from the original job .INFO records,
# not the conflicting narrative tables or current shared defaults.
submit_subset "B cells"                  B_cells                    1500 15 100
submit_subset "CD4_T_cells"              CD4_T_cells                2000 15 100
submit_subset "CD8_T_cells"              CD8_T_cells                2000 15 100
submit_subset "DC"                       DC                         1500 10 100
submit_subset "Endothelial"              Endothelial                2000 20 100
submit_subset "Epithelial"               Epithelial                 2000 20 100
submit_subset "Granulocytes"             Granulocytes                500 10 100
submit_subset "Macrophage_Monocytes"     Macrophage_Monocytes       2000 20 100
submit_subset "NK_cells"                 NK_cells                   2000 15 100
submit_subset "Stromal"                  Stromal                    2000 20 100
submit_subset "Unconventional_T_and_ILCs" Unconventional_T_and_ILCs  2000 15 100
submit_subset "DC"                       DC                         1500 10  50

printf 'Submitted 12 jobs. Submission record: %s\n' "$record"
