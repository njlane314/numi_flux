#!/bin/bash
UBOONE_SETUP=/cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
if [[ -f "$UBOONE_SETUP" ]]; then
  source "$UBOONE_SETUP"
else
  printf 'setup_numiana.sh: warning: missing expected environment setup script: %s\n' "$UBOONE_SETUP" >&2
fi
CFG="${NUMIANA_CONFIG:-$(dirname "${BASH_SOURCE[0]}")/config_numiana.sh}"
[[ -f "$CFG" ]] && source "$CFG"
have_cmd() { command -v "$1" >/dev/null 2>&1; }
maybe_setup() {
  local pkg=$1 ver=$2; shift 2
  if have_cmd setup; then
    if [[ -n "$ver" ]]; then
      setup "$pkg" "$ver" "$@"
    fi
  else
    printf 'setup_numiana.sh: warning: setup command not available; skipping %s\n' "$pkg" >&2
  fi
}
if [[ "${1:-}" == list ]]; then
  if have_cmd ups; then
    if [[ -n "${2:-}" ]]; then ups list -aK+ "$2"; else ups list -aK+ gcc; ups list -aK+ root; ups list -aK+ boost; ups list -aK+ dk2nu; ups list -aK+ ppfx; fi
  else
    printf 'setup_numiana.sh: warning: ups command not available; cannot list products\n' >&2
  fi
  return 0 2>/dev/null || exit 0
fi
if [[ "${1:-}" == show ]]; then
  printf 'QUALS=%s\nGCC=%s\nROOT=%s\nBOOST=%s\nDK2NU_VER=%s\nPPFX_VER=%s\nFLUGG_FHC_DIR=%s\nFLUGG_RHC_DIR=%s\nDK2NU_FHC_DIR=%s\nDK2NU_RHC_DIR=%s\n' \
    "${QUALS:-}" "${GCC:-}" "${ROOT:-}" "${BOOST:-}" "${DK2NU_VER:-}" "${PPFX_VER:-}" "${FLUGG_FHC_DIR:-}" "${FLUGG_RHC_DIR:-}" "${DK2NU_FHC_DIR:-}" "${DK2NU_RHC_DIR:-}"
  return 0 2>/dev/null || exit 0
fi
maybe_setup gcc  "${GCC:-}"
maybe_setup root "${ROOT:-}" ${QUALS:+-q "$QUALS"}
maybe_setup boost "${BOOST:-}" ${QUALS:+-q "$QUALS"}
maybe_setup dk2nu "${DK2NU_VER:-}" ${QUALS:+-q "$QUALS"}
maybe_setup ppfx  "${PPFX_VER:-}" ${QUALS:+-q "$QUALS"}
export NUMIANA_DIR="${NUMIANA_DIR:-$PWD}"
export NUMIANA_INC="$NUMIANA_DIR/include"
export BOOSTROOT="${BOOST_INC:-$BOOSTROOT}"
export DK2NU_INC="${DK2NU:+$DK2NU/include/dk2nu/tree}"
export DK2NU_LIB="${DK2NU:+$DK2NU/lib}"
export FLUGG_FHC_DIR DK2NU_FHC_DIR FLUGG_RHC_DIR DK2NU_RHC_DIR
export LD_LIBRARY_PATH="$NUMIANA_DIR/lib:$NUMIANA_DIR/dict${PPFX_DIR:+:$PPFX_DIR/lib}${PPFX_LIB:+:$PPFX_LIB}:$LD_LIBRARY_PATH"
export LIBRARY_PATH="$LD_LIBRARY_PATH"
