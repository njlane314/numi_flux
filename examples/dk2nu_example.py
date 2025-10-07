from __future__ import print_function

import io
import os
from pathlib import Path

import ROOT

REPO_DIR = Path(__file__).resolve().parents[1]
LIB_PATH = REPO_DIR / "lib" / "Dk2NuFlux_cc.so"

if not LIB_PATH.exists():
    raise RuntimeError(
        "{} does not exist. Please build the project (e.g. `make dk2nu`) "
        "after sourcing the appropriate setup script.".format(LIB_PATH)
    )

load_result = ROOT.gSystem.Load(str(LIB_PATH))
if load_result != 0:
    raise RuntimeError(
        "Failed to load `{}` (gSystem.Load returned {}). Check that the ROOT "
        "environment matches the build used to produce the library and that "
        "all dependent libraries such as libEG are available on your LD_LIBRARY_PATH.".format(
            LIB_PATH, load_result
        )
    )

ROOT.gROOT.SetBatch(True)

# Default globs for the Geant4 10.4 MicroBooNE dk2nu production.
FHC_PATTERN = "/pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4/me000z200i/run*/files/g4numi_minervame_me000z200i_*.root"
RHC_PATTERN = "/pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4/me000z-200i/run*/files/g4numi_minervame_me000z-200i_*.root"

# Allow overrides so you can point at a cached copy, a subset, or a different
# production without editing this example.
fhc_pattern = os.environ.get("UB_NUMIANA_FHC_PATTERN", FHC_PATTERN)
rhc_pattern = os.environ.get("UB_NUMIANA_RHC_PATTERN", RHC_PATTERN)

filelist_path = os.environ.get("UB_NUMIANA_FILELIST", "dk2nu_g4104_fhc_rhc.files")
outfile = os.environ.get("UB_NUMIANA_OUTFILE", "dk2nu_g4104_fhc_rhc.root")

with io.open(filelist_path, "w", encoding="utf-8") as handle:
    handle.write(u"{}\n".format(fhc_pattern))
    handle.write(u"{}\n".format(rhc_pattern))

print("Wrote horn-polarity patterns to {}".format(filelist_path))
print("  FHC: {}".format(fhc_pattern))
print("  RHC: {}".format(rhc_pattern))
print("Output will be stored in {}".format(outfile))

flux = ROOT.Dk2NuFlux(True, str(filelist_path), outfile)
flux.CalculateFlux()
