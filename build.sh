#!/usr/bin/env bash
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./setup_numiana.sh "${MODE:-full}"
t="${1:-all}"
if [[ "$t" == clean ]]; then make clean; exit 0; fi
if [[ "$t" == all ]]; then
  make base
  make flugg
  [[ -n "${DK2NU_LIB:-}" ]] && make dk2nu
  exit 0
fi
make "$t"
