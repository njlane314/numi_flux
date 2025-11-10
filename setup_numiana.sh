#!/bin/bash
set -e
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
QUALS="${QUALS:-e19:prof}"
setup gcc "${GCC:-v8_2_0}"
setup root "${ROOT:-v6_26_10}" -q "${QUALS}"
if [[ "${1:-full}" == full ]]; then
  setup boost "${BOOST:-v1_66_0}" -q "${QUALS}"
  [[ -n "${DK2NU_VER:-}" ]] && setup dk2nu "${DK2NU_VER}" -q "${QUALS}"
  [[ -n "${PPFX_VER:-}"  ]] && setup ppfx  "${PPFX_VER}"  -q "${QUALS}"
fi
export NUMIANA_DIR="${NUMIANA_DIR:-$PWD}"
export NUMIANA_INC="$NUMIANA_DIR/include"
[[ -n "${BOOST_INC:-}" ]] && export BOOSTROOT="${BOOST_INC}"
[[ -n "${DK2NU:-}" ]] && export DK2NU_INC="${DK2NU}/include/dk2nu/tree" DK2NU_LIB="${DK2NU}/lib"
export LD_LIBRARY_PATH="$NUMIANA_DIR/lib:$NUMIANA_DIR/dict${PPFX_DIR:+:$PPFX_DIR/lib}:$LD_LIBRARY_PATH"
export LIBRARY_PATH="$LD_LIBRARY_PATH:$LIBRARY_PATH"
