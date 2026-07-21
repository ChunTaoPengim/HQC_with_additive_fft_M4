#!/bin/bash
set -euo pipefail

REPO_URL="https://github.com/mupq/pqm4.git"
PQDIR="pqm4"
DEST_BASE="$PQDIR/crypto_kem"

# Pin the exact pqm4 revision used for the results reported in the paper.
# Override with: PQM4_COMMIT=<sha> ./pqm4.sh
PQM4_COMMIT="${PQM4_COMMIT:-a24bb4b}"

VARIANTS=(hqc-1 hqc-3 hqc-5)

if [[ ! -d "$PQDIR/.git" ]]; then
  git clone "$REPO_URL" "$PQDIR"
else
  echo "[info] '$PQDIR' already exists; skipping clone."
fi

# Check out the pinned revision and sync submodules to match it.
pushd "$PQDIR" >/dev/null
CURRENT="$(git rev-parse --short HEAD)"
if [[ "$CURRENT" != "$PQM4_COMMIT"* ]]; then
  echo "[info] Checking out pqm4 $PQM4_COMMIT (was $CURRENT)"
  git fetch --all --tags
  git checkout --detach "$PQM4_COMMIT"
fi
git submodule update --init --recursive
echo "[ok] pqm4 at $(git rev-parse HEAD)"
popd >/dev/null

for v in "${VARIANTS[@]}"; do
  SRC_COMMON="$v/src/common"
  SRC_M4="$v/src/m4"
  DEST_DIR="$DEST_BASE/$v/m4f"

  echo "[info] Processing $v"

  [[ -d "$SRC_COMMON" ]] || { echo "[error] Missing directory: $SRC_COMMON" >&2; exit 1; }
  [[ -d "$SRC_M4"     ]] || { echo "[error] Missing directory: $SRC_M4"     >&2; exit 1; }
  [[ -d "$DEST_BASE"  ]] || { echo "[error] Missing directory: $DEST_BASE (did clone succeed?)" >&2; exit 1; }

  mkdir -p "$DEST_DIR"
  cp -a "$SRC_COMMON/." "$DEST_DIR/"
  cp -a "$SRC_M4/."     "$DEST_DIR/"

  echo "[ok] Copied $SRC_COMMON/* and $SRC_M4/* -> $DEST_DIR/"
done

echo "[done] All variants installed into $DEST_BASE/<variant>/m4f/"

SPEED_SRC="speed.c"
SPEED_DEST="$PQDIR/mupq/crypto_kem/speed.c"

echo "[info] Installing component-level benchmarks"
[[ -f "$SPEED_SRC"  ]] || { echo "[error] Missing file: $SPEED_SRC" >&2; exit 1; }
[[ -f "$SPEED_DEST" ]] || { echo "[error] Missing file: $SPEED_DEST (is the mupq submodule initialised?)" >&2; exit 1; }

cp -a "$SPEED_SRC" "$SPEED_DEST"
echo "[ok] Copied $SPEED_SRC -> $SPEED_DEST"

echo "[done] Setup complete."