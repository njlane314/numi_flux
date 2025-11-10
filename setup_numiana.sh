#!/bin/bash
set -euo pipefail
: "${NUMIANA_QUALS:=e19:prof}"
: "${NUMIANA_GCC_VER:=}"
: "${NUMIANA_ROOT_VER:=}"
: "${NUMIANA_BOOST_VER:=}"
: "${NUMIANA_DK2NU_VER:=}"
: "${NUMIANA_PPFX_VER:=}"
: "${NUMIANA_WITH_DK2NU:=1}"
NUMIANA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_emph(){ printf '\033[1;36m%s\033[0m\n' "$*"; }
_warn(){ printf '\033[1;33mWARN:\033[0m %s\n' "$*"; }
_err() { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; }
_pick_latest() {
  local prod="$1"; local filt="${2:-}"
  local ver
  ver=$(ups list -aK+ "$prod" 2>/dev/null \
        | { if [[ -n "$filt" ]]; then grep -i -- "$filt" || true; else cat; fi; } \
        | awk '{print $2}' | sort -V | tail -n 1 || true)
  printf '%s' "$ver"
}
numiana_up() {
  _emph "NUMIANA: activating environment in ${NUMIANA_DIR}"
  if ! command -v setup >/dev/null 2>&1; then
    source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
  fi
  local gcc_ver="${NUMIANA_GCC_VER}"
  if [[ -z "${gcc_ver}" ]]; then
    gcc_ver="$(_pick_latest gcc)"
    [[ -z "${gcc_ver}" ]] && { _err "Could not auto-pick a GCC version (ups list). Set NUMIANA_GCC_VER."; return 1; }
  fi
  setup gcc "${gcc_ver}"
  local root_ver="${NUMIANA_ROOT_VER}"
  if [[ -z "${root_ver}" ]]; then
    root_ver="$(_pick_latest root "${NUMIANA_QUALS}")"
    [[ -z "${root_ver}" ]] && { _err "Could not find a ROOT version with quals ${NUMIANA_QUALS}. Set NUMIANA_ROOT_VER."; return 1; }
  fi
  setup root "${root_ver}" -q "${NUMIANA_QUALS}"
  local boost_ver="${NUMIANA_BOOST_VER}"
  if [[ -z "${boost_ver}" ]]; then
    if ups list -aK+ boost 2>/dev/null | grep -q 'v1_66_0'; then
      boost_ver="v1_66_0"
    else
      boost_ver="$(_pick_latest boost "${NUMIANA_QUALS}")"
    fi
  fi
  if [[ -n "${boost_ver}" ]]; then
    setup boost "${boost_ver}" -q "${NUMIANA_QUALS}" || _warn "Boost ${boost_ver} not available; continuing"
  fi
  if [[ "${NUMIANA_WITH_DK2NU}" == "1" ]]; then
    if ups list -aK+ dk2nu >/dev/null 2>&1 || [[ -n "${NUMIANA_DK2NU_VER}" ]]; then
      local dk2nu_ver="${NUMIANA_DK2NU_VER:-$(_pick_latest dk2nu "${NUMIANA_QUALS}")}"
      [[ -n "${dk2nu_ver}" ]] && setup dk2nu "${dk2nu_ver}" -q "${NUMIANA_QUALS}" || _warn "dk2nu not available"
    fi
    if ups list -aK+ ppfx >/dev/null 2>&1 || [[ -n "${NUMIANA_PPFX_VER}" ]]; then
      local ppfx_ver="${NUMIANA_PPFX_VER:-$(_pick_latest ppfx "${NUMIANA_QUALS}")}"
      [[ -n "${ppfx_ver}" ]] && setup ppfx "${ppfx_ver}" -q "${NUMIANA_QUALS}" || _warn "ppfx not available"
    fi
  fi
  export MODE="NUMI"
  export NUMIANA_DIR
  export NUMIANA_INC="${NUMIANA_DIR}/include"
  export BOOSTROOT="${BOOST_INC:-${BOOSTROOT:-}}"
  export DK2NU_INC="${DK2NU}/include/dk2nu/tree"
  export DK2NU_LIB="${DK2NU}/lib"
  export LD_LIBRARY_PATH="${NUMIANA_DIR}/lib:${NUMIANA_DIR}/dict:${PPFX_DIR:-}/lib:${PPFX_LIB:-}:$LD_LIBRARY_PATH"
  export LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
  if ! command -v root-config >/dev/null 2>&1; then
    _err "root-config not in PATH. Did ROOT fail to set up?"
    return 1
  fi
  if ! echo 'int main(){}' | g++ -std=c++17 -x c++ - -o /dev/null >/dev/null 2>&1; then
    _err "Your g++ doesn't support -std=c++17. Ensure UPS gcc is first in PATH."
    return 1
  fi
  _emph "NUMIANA: environment ready"
  numiana_env
}
numiana_down() {
  _emph "NUMIANA: deactivating environment"
  export LD_LIBRARY_PATH="$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' \
    | grep -v -E "${NUMIANA_DIR}/(lib|dict)|${PPFX_DIR:-}/lib|${PPFX_LIB:-}" \
    | paste -sd: -)"
  export LIBRARY_PATH="$LD_LIBRARY_PATH"
  unset MODE NUMIANA_DIR NUMIANA_INC BOOSTROOT DK2NU_INC DK2NU_LIB
}
numiana_env() {
  echo "  QUALS       : ${NUMIANA_QUALS}"
  echo "  GCC         : $(command -v g++ 2>/dev/null) | $(g++ --version | head -1 || true)"
  echo "  ROOT        : ${ROOTSYS:-<unset>} | $(command -v root-config 2>/dev/null || true)"
  echo "  BOOSTROOT   : ${BOOSTROOT:-<unset>}"
  echo "  DK2NU_INC   : ${DK2NU_INC:-<unset>}"
  echo "  DK2NU_LIB   : ${DK2NU_LIB:-<unset>}"
  echo "  NUMIANA_DIR : ${NUMIANA_DIR}"
  echo "  LD_LIBRARY_PATH entries: $(echo "$LD_LIBRARY_PATH" | awk -F: '{print NF}')"
}
if hostname -f 2>/dev/null | grep -q 'uboone'; then
  numiana_up
else
  _warn "Not on a uboone host; run 'numiana_up' manually if desired."
fi
