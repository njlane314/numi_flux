#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
#source ./setup_numiana.sh
if ! command -v root-config >/dev/null 2>&1; then
  echo "build.sh: error: root-config not found; ensure ROOT is set up." >&2
  exit 1
fi
make "${1:-all}"
