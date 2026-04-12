#!/bin/bash
set -euo pipefail
REPO_URL="https://github.com/mupq/pqm4.git"
PQDIR="pqm4"
DEST_BASE="$PQDIR/crypto_kem"

# Variants to process
VARIANTS=(hqc-1 hqc-3 hqc-5)

# Clone once (skip if already cloned)
if [[ ! -d "$PQDIR/.git" ]]; then
  git clone --recursive "$REPO_URL"
else
  echo "[info] '$PQDIR' already exists; skipping clone."
fi

# For each variant, make m4f dir and copy files
for v in "${VARIANTS[@]}"; do
  SRC_COMMON="$v/src/common"
  SRC_M4="$v/src/m4"
  DEST_DIR="$DEST_BASE/$v/m4f"

  echo "[info] Processing $v"
  mkdir -p "$DEST_DIR"

  # Basic checks (helps catch typos early)
  [[ -d "$SRC_COMMON" ]] || { echo "[error] Missing directory: $SRC_COMMON" >&2; exit 1; }
  [[ -d "$SRC_M4"     ]] || { echo "[error] Missing directory: $SRC_M4" >&2; exit 1; }
  [[ -d "$DEST_BASE"  ]] || { echo "[error] Missing directory: $DEST_BASE (did clone succeed?)" >&2; exit 1; }

  # Copy contents of both folders into the same destination
  # (Trailing /. copies "contents" rather than the directory itself.)
  cp -a "$SRC_COMMON/." "$DEST_DIR/"
  cp -a "$SRC_M4/."     "$DEST_DIR/"

  echo "[ok] Copied $SRC_COMMON/* and $SRC_M4/* -> $DEST_DIR/"
done

echo "[done] All variants installed into $DEST_BASE/<variant>/m4f/"