// root -l -b -q 'scripts/flux_integrals_minimal.C++'

#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TString.h"

#include <cstdio>
#include <cmath>
#include <limits>

// ------------------------------------------------------------------
// Sum POT by reading the "POT" histogram in each file in the TChain.
// (Same approach as in your plotting macro.)
// ------------------------------------------------------------------
static double sumPOT(TChain &ch) {
  double tot = 0.0;
  if (auto *files = ch.GetListOfFiles()) {
    for (int i = 0, n = files->GetEntriesFast(); i < n; ++i) {
      auto *el = (TChainElement*)files->UncheckedAt(i);
      if (!el) continue;
      TFile f(el->GetTitle(), "READ");
      if (!f.IsOpen()) continue;
      if (auto *hp = dynamic_cast<TH1*>(f.Get("POT"))) tot += hp->Integral();
    }
  }
  return tot;
}

// ------------------------------------------------------------------
// Compute a single integrated flux with selection on PDG and energy.
// This fills a 1-bin histogram so the bin content is the sum of weights.
// ------------------------------------------------------------------
static double integrate_flux(TChain &ch,
                             int pdg,
                             double Emin, double Emax,
                             const char *wexpr,
                             const char *tag)
{
  TString hname = Form("h_int_%s_%d", tag, pdg);
  TH1D h(hname, "", 1, 0., 1.);
  h.Sumw2();
  ch.Draw(Form("0.5>>%s", hname.Data()),
          Form("(%g<=nuE && nuE<%g) * (ntype==%d) * (%s)", Emin, Emax, pdg, wexpr),
          "goff");
  return h.GetBinContent(1);
}

// ------------------------------------------------------------------
// Print a block of integrals for one running mode (FHC or RHC).
// ------------------------------------------------------------------
static void run_mode(const char *file, const char *tag,
                     double Emin, double Emax,
                     bool norm_per_pot, double nominal_pot)
{
  TChain ch("outTree");
  ch.Add(file);

  const bool has_wgt  = ch.GetListOfBranches()->FindObject("wgt");
  const bool has_ppfx_cv  = ch.GetListOfBranches()->FindObject("ppfx_cv") ||
                            ch.GetListOfBranches()->FindObject("wgt_ppfx_cv");
  const bool has_ppfx_vec = ch.GetListOfBranches()->FindObject("wgt_ppfx");

  // Base weight: geometry/acceptance * (optional) PPFX CV
  TString w = has_wgt ? "wgt" : "1";
  if (has_ppfx_cv) {
    // Prefer an explicit CV scalar if present
    w += Form("*%s",
              ch.GetListOfBranches()->FindObject("ppfx_cv") ? "ppfx_cv" : "wgt_ppfx_cv");
  } else if (has_ppfx_vec) {
    // Fall back to the first element of the vector (conventionally CV)
    w += "*wgt_ppfx[0]";
  }

  // Build scaled weight expression
  const double pot_total = sumPOT(ch);
  TString w_scaled;
  if (pot_total > 0.0) {
    if (norm_per_pot) {
      w_scaled.Form("(%s)/(%.*g)", w.Data(), 15, pot_total);             // per POT
    } else {
      w_scaled.Form("(%s)*(%.*g)", w.Data(), 15, nominal_pot/pot_total); // scaled to nominal_pot
    }
  } else {
    std::fprintf(stderr,
      "[flux_integrals_minimal] WARNING: POT<=0 in inputs; proceeding without POT scaling.\n");
    w_scaled = w;
  }

  // Integrate for each species
  const double I_numu  = integrate_flux(ch,  +14, Emin, Emax, w_scaled.Data(), tag);
  const double I_anumu = integrate_flux(ch,  -14, Emin, Emax, w_scaled.Data(), tag);
  const double I_nue   = integrate_flux(ch,  +12, Emin, Emax, w_scaled.Data(), tag);
  const double I_anue  = integrate_flux(ch,  -12, Emin, Emax, w_scaled.Data(), tag);
  const double I_tot   = I_numu + I_anumu + I_nue + I_anue;

  const char *units = norm_per_pot ? "nu / POT / cm^2" : "nu / (6e20 POT) / cm^2";

  // Print nicely
  std::printf("\n================  %s mode  ================\n", tag);
  std::printf("Input file      : %s\n", file);
  std::printf("POT in inputs   : %.6g\n", pot_total);
  std::printf("Energy window   : %.3f – %.3f GeV\n", Emin, Emax);
  std::printf("Normalization   : %s%s\n",
              norm_per_pot ? "per POT" : "scaled to 6e20 POT",
              (has_ppfx_cv || has_ppfx_vec) ? " (includes PPFX CV)" : "");
  std::printf("Integrated flux (units = %s):\n", units);
  std::printf("  nu_mu     : %.6e\n", I_numu);
  std::printf("  nubar_mu  : %.6e\n", I_anumu);
  std::printf("  nu_e      : %.6e\n", I_nue);
  std::printf("  nubar_e   : %.6e\n", I_anue);
  std::printf("  TOTAL     : %.6e\n", I_tot);

  if (I_tot > 0.0) {
    std::printf("Fractions (% of total): "
                "nu_mu %.2f%% | nubar_mu %.2f%% | nu_e %.2f%% | nubar_e %.2f%%\n",
                100.0*I_numu/I_tot, 100.0*I_anumu/I_tot,
                100.0*I_nue/I_tot,  100.0*I_anue/I_tot);
  }
}

// ------------------------------------------------------------------
// Entry point
// ------------------------------------------------------------------
void flux_integrals_minimal() {
  // Avoid attaching created histograms to any TFile directory
  TH1::AddDirectory(kFALSE);

  // ---------- hardcoded inputs ----------
  const char *FHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx_g4_10_4.root";
  const char *RHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc_ppfx_g4_10_4.root";

  const double Emin = 0.25;  // GeV
  const double Emax = 10.0;  // GeV

  // Normalization: per POT (true) or scaled to NOMINAL_POT (false)
  constexpr bool   NORM_PER_POT = false;
  constexpr double NOMINAL_POT  = 6e20;

  run_mode(FHC_FILE, "FHC", Emin, Emax, NORM_PER_POT, NOMINAL_POT);
  run_mode(RHC_FILE, "RHC", Emin, Emax, NORM_PER_POT, NOMINAL_POT);

  std::printf("============================================\n\n");
}




