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
# Defaults tuned to the AlmaLinux 9 MicroBooNE UPS snapshot.
: "${NUMIANA_ROOT_VERSION:=v6_28_12}"
: "${NUMIANA_ROOT_QUALS:=e20:p3915:prof}"
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
: "${NUMIANA_BOOST_VERSION:=v1_82_0}"
: "${NUMIANA_BOOST_QUALS:=e20:prof}"
: "${NUMIANA_PPFX_VERSION:=v02_17_07}"
: "${NUMIANA_PPFX_QUALS:=e20:prof}"
: "${NUMIANA_DK2NU_VERSION:=v01_05_01b}"
# dk2nu has not yet been rebuilt for the e20 toolchain; use the e15 release
# unless the user overrides it explicitly.
: "${NUMIANA_DK2NU_QUALS:=e15:prof}"
: "${NUMIANA_JOBSUB_VERSION:=v1_3_1}"
: "${NUMIANA_JOBSUB_QUALS:=}"

setup_errors=()

setup_product() {
  local product=$1
  local version=$2
  local quals=$3
  local setup_args=("${product}" "${version}")

  if [[ -n "${quals}" ]]; then
    setup_args+=("-q" "${quals}")
  fi

  if ! setup "${setup_args[@]}"; then
    setup_errors+=("${product}")
    echo "[numi_flux] Failed to set up ${product} ${version} ${quals:+(quals: ${quals})}" >&2
    echo "[numi_flux] To inspect available versions run: ups list -aK+ ${product}" >&2
  fi
}

setup_product root        "${NUMIANA_ROOT_VERSION}"        "${NUMIANA_ROOT_QUALS}"
setup_product cmake       "${NUMIANA_CMAKE_VERSION}"       "${NUMIANA_CMAKE_QUALS}"
setup_product nlohmann_json "${NUMIANA_NLOHMANN_JSON_VERSION}" "${NUMIANA_NLOHMANN_JSON_QUALS}"
setup_product tbb         "${NUMIANA_TBB_VERSION}"         "${NUMIANA_TBB_QUALS}"
setup_product libtorch    "${NUMIANA_LIBTORCH_VERSION}"    "${NUMIANA_LIBTORCH_QUALS}"
setup_product eigen       "${NUMIANA_EIGEN_VERSION}"       "${NUMIANA_EIGEN_QUALS}"
setup_product boost       "${NUMIANA_BOOST_VERSION}"       "${NUMIANA_BOOST_QUALS}"
setup_product ppfx        "${NUMIANA_PPFX_VERSION}"        "${NUMIANA_PPFX_QUALS}"
setup_product dk2nu       "${NUMIANA_DK2NU_VERSION}"       "${NUMIANA_DK2NU_QUALS}"

if [[ -n "${NUMIANA_JOBSUB_QUALS}" ]]; then
  setup_product jobsub_client "${NUMIANA_JOBSUB_VERSION}" "${NUMIANA_JOBSUB_QUALS}"
else
  setup_product jobsub_client "${NUMIANA_JOBSUB_VERSION}" ""
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

if (( ${#setup_errors[@]} )); then
  echo "[numi_flux] The UPS setup reported errors for: ${setup_errors[*]}" >&2
fi

if (( ${#missing_vars[@]} )); then
  echo "[numi_flux] Warning: the following variables are unset after setup:" >&2
  printf '  - %s\n' "${missing_vars[@]}" >&2
  echo "You may need to adjust the product versions/qualifiers above." >&2
else
  echo "[numi_flux] UPS environment ready.  NUMIANA_DIR=${NUMIANA_DIR}"
fi
