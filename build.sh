#!/usr/bin/env bash
build_numiana() {
  cd "$(dirname "${BASH_SOURCE[0]}")" || {
    echo "build.sh: error: could not cd to script directory." >&2
    return 1
  }
  if ! command -v root-config >/dev/null 2>&1; then
    echo "build.sh: error: root-config not found; ensure ROOT is set up." >&2
    return 1
  fi
  make "${1:-all}"
}
build_numiana "$@"
