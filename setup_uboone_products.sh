#!/bin/bash

# Source this script to set up the UPS products needed by the numi_flux
# build on MicroBooNE gpvm or uboonepro machines.  The defaults are chosen
# to match the e20 toolchain used by uboonecode mcc9.  Override any of the
# NUMIANA_* variables before sourcing the script to pick a different
# product version or qualifier (for example NUMIANA_ROOT_VERSION).

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "[numi_flux] Please source this script instead of executing it." >&2
  echo "Example: source setup_uboone_products.sh" >&2
  exit 1
fi

# Use mcc9 product definitions (provides the e20 compiler stack).
if [[ -z "${SETUP_UBOONE_PRODUCTS:-}" ]]; then
  source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
fi

# Common Fermilab UPS setups (samweb, ifdhc, etc.).
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups

setup sam_web_client
setup ifdhc

# Allow users to override product versions/qualifiers.
: "${NUMIANA_ROOT_VERSION:=v6_28_12}"
: "${NUMIANA_ROOT_QUALS:=e20:p3920:prof}"
: "${NUMIANA_CMAKE_VERSION:=v3_20_0}"
: "${NUMIANA_CMAKE_QUALS:=Linux64bit+3.10-2.17}"
: "${NUMIANA_NLOHMANN_JSON_VERSION:=v3_11_2}"
: "${NUMIANA_NLOHMANN_JSON_QUALS:=e20}"
: "${NUMIANA_TBB_VERSION:=v2021_9_0}"
: "${NUMIANA_TBB_QUALS:=e20}"
: "${NUMIANA_LIBTORCH_VERSION:=v1_13_1b}"
: "${NUMIANA_LIBTORCH_QUALS:=e20:prof}"
: "${NUMIANA_EIGEN_VERSION:=v3_4_0}"
: "${NUMIANA_EIGEN_QUALS:=e20}"
: "${NUMIANA_BOOST_VERSION:=v1_66_0}"
: "${NUMIANA_BOOST_QUALS:=e20:prof}"
: "${NUMIANA_PPFX_VERSION:=v1_72}"
: "${NUMIANA_PPFX_QUALS:=e20:prof}"
: "${NUMIANA_DK2NU_VERSION:=v01_05_05d}"
: "${NUMIANA_DK2NU_QUALS:=e20}"
: "${NUMIANA_JOBSUB_VERSION:=v1_4_2}"
: "${NUMIANA_JOBSUB_QUALS:=}"

setup root        "${NUMIANA_ROOT_VERSION}"        -q "${NUMIANA_ROOT_QUALS}"
setup cmake       "${NUMIANA_CMAKE_VERSION}"       -q "${NUMIANA_CMAKE_QUALS}"
setup nlohmann_json "${NUMIANA_NLOHMANN_JSON_VERSION}" -q "${NUMIANA_NLOHMANN_JSON_QUALS}"
setup tbb         "${NUMIANA_TBB_VERSION}"         -q "${NUMIANA_TBB_QUALS}"
setup libtorch    "${NUMIANA_LIBTORCH_VERSION}"    -q "${NUMIANA_LIBTORCH_QUALS}"
setup eigen       "${NUMIANA_EIGEN_VERSION}"       -q "${NUMIANA_EIGEN_QUALS}"
setup boost       "${NUMIANA_BOOST_VERSION}"       -q "${NUMIANA_BOOST_QUALS}"
setup ppfx        "${NUMIANA_PPFX_VERSION}"        -q "${NUMIANA_PPFX_QUALS}"
setup dk2nu       "${NUMIANA_DK2NU_VERSION}"       -q "${NUMIANA_DK2NU_QUALS}"

if [[ -n "${NUMIANA_JOBSUB_QUALS}" ]]; then
  setup jobsub_client "${NUMIANA_JOBSUB_VERSION}" -q "${NUMIANA_JOBSUB_QUALS}"
else
  setup jobsub_client "${NUMIANA_JOBSUB_VERSION}"
fi

# The UPS boost product defines BOOST_INC but not BOOST_DIR.  Populate the
# latter for consistency with the Makefile expectations.
if [[ -z "${BOOST_DIR:-}" && -n "${BOOST_INC:-}" ]]; then
  export BOOST_DIR="${BOOST_INC%/include}"
fi

# Pull in the repository-specific environment customisation.
source "$(dirname "${BASH_SOURCE[0]}")/setup_numiana.sh"

# Sanity checks so users notice missing dependencies immediately.
missing_vars=()
for v in NUMIANA_DIR NUMIANA_INC PPFX_DIR PPFX_LIB DK2NU BOOST_DIR; do
  if [[ -z "${!v:-}" ]]; then
    missing_vars+=("$v")
  fi
done

if (( ${#missing_vars[@]} )); then
  echo "[numi_flux] Warning: the following variables are unset after setup:" >&2
  printf '  - %s\n' "${missing_vars[@]}" >&2
  echo "You may need to adjust the product versions/qualifiers above." >&2
else
  echo "[numi_flux] UPS environment ready.  NUMIANA_DIR=${NUMIANA_DIR}"
fi
