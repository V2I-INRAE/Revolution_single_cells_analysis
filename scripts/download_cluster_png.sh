#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 user@remote-host" >&2
  exit 2
fi

remote_host=$1
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
relative_path="results/job-44542861/clusters_plots/subset/"
remote_repo_root="/work/project/revo-pig-sc/analysis"
remote_path="${remote_repo_root}/${relative_path}"
local_path="${repo_root}/${relative_path}"

mkdir -p "$local_path"
rsync -avm --include='*/' --include='*clustree*' --include='*umap*' --exclude='*' \
  "${remote_host}:${remote_path}" "$local_path"
