#!/bin/bash
# Download all REVO sample folders from Genotoul to the local Mac, in priority order.
# Prevent your Mac from sleeping before running: caffeinate -i bash download_reads.sh
# Re-run anytime — rsync resumes interrupted folders (--partial) and skips done files.

SRC="genobioinfo:/work/project/revo-pig-sc/Revolution_single_cells_analysis/raw_data"
DEST="/Users/fbraza/Documents/PROJECTS/REVOLUTION/RESULTS/SINGLE_CELLS/reads"
mkdir -p "$DEST"

# Prefer the Homebrew rsync 3.x — the macOS system rsync (openrsync) is limited
# (no --info=progress2, and flaky with other rsync-3 options).
if [ -x /opt/homebrew/bin/rsync ]; then        # Apple Silicon
  RSYNC=/opt/homebrew/bin/rsync
elif [ -x /usr/local/bin/rsync ]; then         # Intel Mac
  RSYNC=/usr/local/bin/rsync
else
  RSYNC=rsync
fi
echo "Using rsync: $RSYNC ($("$RSYNC" --version | head -1))"

# Priority order — edit freely (controls first, then treatments)
FOLDERS=(
  REVO26-C REVO27-C REVO29-C REVO30-C REVO31-C
  REVO26-N4 REVO26-N10 REVO27-N4 REVO27-N10
  REVO29-N4 REVO30-N4 REVO31-N4 REVO31-N8
  REVO26-P4 REVO26-P10 REVO27-P4 REVO27-P10
  REVO29-P4 REV30-P4 REV31-P4 REVO30-P10 REVO31-P10
)

# --- Option A (default): sequential, one folder after another ---
for f in "${FOLDERS[@]}"; do
  echo "=== $(date '+%H:%M:%S') downloading $f ==="
  "$RSYNC" -avh \
    --partial \
    --partial-dir=.rsync-partial \
    --progress \
    "$SRC/$f/" "$DEST/$f/"
done

# --- Option B (tested 2026-09-22): 6 parallel streams, aggregate ~1.5 MB/s ---
# xargs -P N = number of SIMULTANEOUS transfers (not the number of folders).
# Per-stream speed stayed at 250 kB/s with -P 6 (limit is per-connection).
# Uncomment the block below and comment out the loop above to use it.
#
# printf '%s\n' "${FOLDERS[@]}" | xargs -P 6 -I{} \
#   "$RSYNC" -avh \
#     --partial \
#     --partial-dir=.rsync-partial \
#     --progress \
#     "$SRC/{}/" "$DEST/{}/"

echo "=== $(date '+%H:%M:%S') finished — rerun this script to fill any gaps ==="
