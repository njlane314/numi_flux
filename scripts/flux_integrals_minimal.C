// root -l -b -q 'scripts/flux_integrals_minimal.C++'

#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

#include <cstdio>
#include <cmath>
#include <limits>
#include <cstring>

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
  // Let ROOT create the 1-bin histogram inline; retrieve it via GetHistogram().
  // This is robust even with TH1::AddDirectory(kFALSE).
  const TString drawspec = Form("0.5>>h_int_tmp_%s_%d(1,0,1)", tag, pdg);
  const TString cut      = Form("(%g<=nuE && nuE<%g) * (ntype==%d) * (%s)",
                                Emin, Emax, pdg, wexpr);
  ch.Draw(drawspec, cut, "goff");
  TH1 *h = ch.GetHistogram();
  const double sumw = h ? h->GetBinContent(1) : 0.0;
  if (h) delete h; // avoid leaking histograms created by Draw
  return sumw;
}

// ------------------------------------------------------------------
// Print a block of integrals for one running mode (FHC or RHC).
// ------------------------------------------------------------------
static void run_mode(const char *file, const char *tag,
                     double Emin, double Emax,
                     bool norm_per_pot, double nominal_pot)
{
  // --- pick a real tree name from the file ---
  TString tname = "outTree";
  {
    TFile f(file, "READ");
    if (!f.IsOpen()) {
      std::fprintf(stderr, "[flux_integrals_minimal] ERROR: cannot open %s\n", file);
      return;
    }
    auto exists_as_tree = [&](const char* nm){
      if (auto *o = f.Get(nm)) return o->InheritsFrom(TTree::Class());
      return false;
    };
    if (!exists_as_tree("outTree")) {
      const char* cand[] = {"dk2nuTree","flux","NuMIFlux","Events","tree"};
      for (auto nm : cand) {
        if (exists_as_tree(nm)) { tname = nm; break; }
      }
      if (tname == "outTree") {
        // fallback: first TTree in the file
        TIter next(f.GetListOfKeys());
        while (TKey *k = (TKey*)next()) {
          if (TString(k->GetClassName()) == "TTree") { tname = k->GetName(); break; }
        }
      }
    }
  }

  TChain ch(tname);
  if (ch.Add(file) <= 0) {
    std::fprintf(stderr, "[flux_integrals_minimal] ERROR: failed to add %s to chain '%s'\n",
                 file, tname.Data());
    return;
  }

  // --- create aliases so your selection terms still work ---
  auto *br = ch.GetListOfBranches();
  auto has = [&](const char* n){ return br && br->FindObject(n); };

  if (!has("nuE")) {
    const char* candE[] = {"Ev","Enu","E","enu"};
    for (auto nm : candE) if (has(nm)) { ch.SetAlias("nuE", nm); break; }
  }
  if (!has("ntype")) {
    const char* candP[] = {"pdg","nuPDG","pdg_nu","inu"};
    for (auto nm : candP) if (has(nm)) { ch.SetAlias("ntype", nm); break; }
  }

  // --- build weight expression with robust PPFX CV handling ---
  TString w = has("wgt") ? "wgt" : "1";
  TString ppfx = "";
  if (has("ppfx_cv"))          ppfx = "ppfx_cv";
  else if (has("wgt_ppfx_cv")) ppfx = "wgt_ppfx_cv";
  else {
    // Decide if wgt_ppfx is a vector/array or a scalar. If it's a scalar,
    // using [0] would zero out the weight selection.
    if (auto *bppfx = ch.GetBranch("wgt_ppfx")) {
      const char *cls = bppfx->GetClassName();
      if (cls && std::strstr(cls, "vector")) {
        ppfx = "wgt_ppfx[0]";   // convention: [0] is CV
      } else {
        ppfx = "wgt_ppfx";      // scalar branch
      }
    }
  }
  if (!ppfx.IsNull()) w += "*" + ppfx;

  const double pot_total = sumPOT(ch);
  TString w_scaled = w;
  if (pot_total > 0.0) {
    if (norm_per_pot) w_scaled = Form("(%s)/%.15g", w.Data(), pot_total);
    else              w_scaled = Form("(%s)*%.15g", w.Data(), nominal_pot/pot_total);
  } else {
    std::fprintf(stderr,
      "[flux_integrals_minimal] WARNING: POT<=0 in inputs; proceeding without POT scaling.\n");
  }

  // --- sanity prints ---
  const Long64_t nEntries = ch.GetEntries();
  std::printf("\n================  %s mode  ================\n", tag);
  std::printf("Input file      : %s\n", file);
  std::printf("Tree name       : %s\n", tname.Data());
  std::printf("Entries in tree : %lld\n", nEntries);
  std::printf("POT in inputs   : %.6g\n", pot_total);
  std::printf("Energy window   : %.3f – %.3f GeV\n", Emin, Emax);
  std::printf("Normalization   : %s%s\n",
              norm_per_pot ? "per POT" : "scaled to 6e20 POT",
              ppfx.IsNull() ? "" : " (includes PPFX CV)");
  if (nEntries == 0) {
    std::fprintf(stderr, "[flux_integrals_minimal] ERROR: chain has 0 entries — wrong tree name.\n");
    return;
  }
  if (!has("nuE") && !ch.GetAlias("nuE")) {
    std::fprintf(stderr, "[flux_integrals_minimal] ERROR: cannot find energy branch (nuE/Ev/Enu/...)\n");
    return;
  }
  if (!has("ntype") && !ch.GetAlias("ntype")) {
    std::fprintf(stderr, "[flux_integrals_minimal] ERROR: cannot find PDG branch (ntype/pdg/...)\n");
    return;
  }

  // --- integrate using your existing helper (still uses nuE & ntype) ---
  const double I_numu  = integrate_flux(ch,  +14, Emin, Emax, w_scaled.Data(), tag);
  const double I_anumu = integrate_flux(ch,  -14, Emin, Emax, w_scaled.Data(), tag);
  const double I_nue   = integrate_flux(ch,  +12, Emin, Emax, w_scaled.Data(), tag);
  const double I_anue  = integrate_flux(ch,  -12, Emin, Emax, w_scaled.Data(), tag);
  const double I_tot   = I_numu + I_anumu + I_nue + I_anue;

  const char *units = norm_per_pot ? "nu / POT / cm^2" : "nu / (6e20 POT) / cm^2";
  std::printf("Integrated flux (units = %s):\n", units);
  std::printf("  nu_mu     : %.6e\n", I_numu);
  std::printf("  nubar_mu  : %.6e\n", I_anumu);
  std::printf("  nu_e      : %.6e\n", I_nue);
  std::printf("  nubar_e   : %.6e\n", I_anue);
  std::printf("  TOTAL     : %.6e\n", I_tot);

  if (I_tot > 0.0) {
    std::printf("Fractions (%% of total): "
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




