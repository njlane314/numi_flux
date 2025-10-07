import os
from pathlib import Path

import ROOT

ROOT.gSystem.Load("../lib/Dk2NuFlux_cc.so")
ROOT.gROOT.SetBatch(True)

# Default globs for the Geant4 10.4 MicroBooNE dk2nu production.
FHC_PATTERN = "/pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4/me000z200i/run*/files/g4numi_minervame_me000z200i_*.root"
RHC_PATTERN = "/pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4/me000z-200i/run*/files/g4numi_minervame_me000z-200i_*.root"

# Allow overrides so you can point at a cached copy, a subset, or a different
# production without editing this example.
fhc_pattern = os.environ.get("UB_NUMIANA_FHC_PATTERN", FHC_PATTERN)
rhc_pattern = os.environ.get("UB_NUMIANA_RHC_PATTERN", RHC_PATTERN)

filelist_path = Path(os.environ.get("UB_NUMIANA_FILELIST", "dk2nu_g4104_fhc_rhc.files"))
outfile = os.environ.get("UB_NUMIANA_OUTFILE", "dk2nu_g4104_fhc_rhc.root")

with filelist_path.open("w", encoding="utf-8") as handle:
    handle.write(f"{fhc_pattern}\n")
    handle.write(f"{rhc_pattern}\n")

print(f"Wrote horn-polarity patterns to {filelist_path}")
print(f"  FHC: {fhc_pattern}")
print(f"  RHC: {rhc_pattern}")
print(f"Output will be stored in {outfile}")

flux = ROOT.Dk2NuFlux(True, str(filelist_path), outfile)
flux.CalculateFlux()
