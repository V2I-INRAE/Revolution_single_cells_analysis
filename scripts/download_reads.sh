#!/bin/bash
# Download all REVO sample folders from Genotoul to your local computer.

SRC="genobioinfo:/work/project/revo-pig-sc/Revolution_single_cells_analysis/raw_data"
DEST="/Users/fbraza/Documents/PROJECTS/REVOLUTION/RESULTS/SINGLE_CELLS/reads"
mkdir -p "$DEST"

FOLDERS=(
  REVO26-C REVO27-C REVO29-C REVO30-C REVO31-C
  REVO26-N4 REVO26-N10 REVO27-N4 REVO27-N10
  REVO29-N4 REVO30-N4 REVO31-N4 REVO31-N8
  REVO26-P4 REVO26-P10 REVO27-P4 REVO27-P10
  REVO29-P4 REV30-P4 REV31-P4 REVO30-P10 REVO31-P10
)

# --- Option A (default): sequential, one folder after another ---
# --- USE THIS IF YOU HAVE NOT SETUP THE SSH KEY WITH THE CLUSTER ---
for f in "${FOLDERS[@]}"; do
  echo "=== $(date '+%H:%M:%S') downloading $f ==="
  rsync -avh \
    --partial \
    --partial-dir=.rsync-partial \
    --progress \
    "$SRC/$f/" "$DEST/$f/"
done

# --- Option B (tested 2026-09-22): 6 parallel streams, aggregate ~1.5 MB/s ---
# --- Use this if you ssh key with cluster because rsync opens several ssh connections ---
# printf '%s\n' "${FOLDERS[@]}" | xargs -P 6 -I{} \
#   "$RSYNC" -avh \
#     --partial \
#     --partial-dir=.rsync-partial \
#     --progress \
#     "$SRC/{}/" "$DEST/{}/"

echo "=== $(date '+%H:%M:%S') finished — rerun this script to fill any gaps ==="
