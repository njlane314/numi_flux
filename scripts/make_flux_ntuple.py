#!/usr/bin/env python3
import os
import sys
import ROOT

ROOT.gROOT.SetBatch(True)

if ROOT.gSystem.Load(os.path.join(os.path.dirname(__file__), "..", "lib", "Dk2NuFlux_cc.so")) != 0:
    if ROOT.gSystem.Load("libDk2NuFlux_cc.so") != 0:
        sys.exit("ERROR: could not load Dk2NuFlux_cc.so. Did you run './build.sh dk2nu'?")


def default_pattern(mode: str) -> str:
    if mode.upper() == "FHC":
        base = os.environ.get(
            "DK2NU_FHC_DIR",
            "/pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4_nothresh/me000z200i/run0/files"
        )
        return os.path.join(base, "g4numi_minervame_me000z200i_*.root")
    else:
        base = os.environ.get(
            "DK2NU_RHC_DIR",
            "/pnfs/uboone/persistent/users/bnayak/flux_files/uboone_beamsim_g4.10.4_nothresh/me000z-200i/run0/files"
        )
        return os.path.join(base, "g4numi_minervame_me000z-200i_*.root")


def main():
    mode = (sys.argv[1] if len(sys.argv) > 1 else "FHC").upper()
    inpat = sys.argv[2] if len(sys.argv) > 2 else default_pattern(mode)
    out = sys.argv[3] if len(sys.argv) > 3 else f"outputs/ubnumi_g4104_nothresh_{mode}.root"

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    print(f"[make_flux_ntuple] mode={mode}\n  in : {inpat}\n  out: {out}")
    ROOT.Dk2NuFlux(inpat, out).CalculateFlux()
    print("[make_flux_ntuple] done.")


if __name__ == "__main__":
    main()
