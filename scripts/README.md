# FHC
python scripts/make_flux_ntuple.py FHC

# RHC
python scripts/make_flux_ntuple.py RHC

# FHC (per-POT, μ-flavour only, 0.25–5 GeV)
root -l -q 'integrate_mu_flux.C("fhc_numi_g4104_ppfx.root",0.25,5.0,8,true)'

# RHC
root -l -q 'integrate_mu_flux.C("rhc_numi_g4104_ppfx.root",0.25,5.0,8,true)'
