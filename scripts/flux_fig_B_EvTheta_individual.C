// ============================================================================
// flux_fig_B_EvTheta_individual.C
//
// Figure B — Energy–angle (Eν, θν) heatmaps per flavor (individual plots).
//
// Why: makes the energy–angle correlation explicit; motivates any acceptance
//       in angle and the analysis’ angle variable.
// Context to cite: 2D analogue of the 1D spectra used to normalise in Sec. 2.2;
//       connects to the integrated fluxes for 0.25–10 GeV (Table 1).
//
// What this macro does:
//   • Opens FHC & RHC files.
//   • For each flavor {numu,numubar,nue,nuebar}, reads
//       <flav>/Detsmear/<flav>_CV_AV_TPC_2D
//   • Draws **one canvas per flavor per mode** with COLZ and **logZ**.
//   • Enforces identical X/Y ranges across ALL flavors and modes,
//     framing Eν on [0.25, 8] GeV (as requested) and using the
//     common θ range intersection across histograms.
//   • Enforces a **shared Z scale** across ALL plots (min>0, max in that ROI).
//   • Extra right margin for the palette; white background; watermark at top‑right.
//
// Usage:
//   root -l -b -q flux_fig_B_EvTheta_individual.C
//   root -l -b -q 'flux_fig_B_EvTheta_individual.C("path/FHC.root","path/RHC.root")'
// Outputs: figB_Etheta_<MODE>_<FLAVOR>.pdf
// ============================================================================

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TLatex.h"
#include "TColor.h"

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <cmath>

namespace CFG {
  // ---- Default input files (your structure) ----
  const char* FILE_FHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* FILE_RHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";

  // ---- Flavors to process ----
  static const std::vector<std::string> FLAVS = {"numu","numubar","nue","nuebar"};

  // ---- Drawing configuration ----
  constexpr bool   USE_LOGZ        = true;
  constexpr int    N_CONTOUR       = 255;
  constexpr double PAD_RIGHT_MARG  = 0.18;     // more room for Z axis (palette)
  constexpr double PAD_LEFT_MARG   = 0.15;
  constexpr double PAD_TOP_MARG    = 0.06;
  constexpr double PAD_BOTTOM_MARG = 0.12;

  // ---- Target display window (requested) ----
  // Frame X to 0.25–8 GeV (even if histogram stops at 5, the extra region appears blank).
  constexpr double FRAME_E_MIN = 0.25;  // GeV
  constexpr double FRAME_E_MAX = 8.0;   // GeV
}

// ------------------------------ Style ---------------------------------------
static void set_minimal_white_style(){
  TStyle* s = new TStyle("FigBStyle","FigBStyle");
  s->SetOptStat(0); s->SetOptTitle(0);
  // Pure white everywhere (no grey)
  s->SetCanvasColor(0); s->SetPadColor(0); s->SetFrameFillColor(0); s->SetFillColor(0);
  s->SetCanvasBorderMode(0); s->SetPadBorderMode(0); s->SetFrameBorderMode(0);
  // Margins
  s->SetPadLeftMargin(CFG::PAD_LEFT_MARG);
  s->SetPadRightMargin(CFG::PAD_RIGHT_MARG);
  s->SetPadTopMargin(CFG::PAD_TOP_MARG);
  s->SetPadBottomMargin(CFG::PAD_BOTTOM_MARG);
  // Fonts & sizes
  s->SetLabelFont(42,"xyz"); s->SetTitleFont(42,"xyz");
  s->SetLabelSize(0.045,"xyz"); s->SetTitleSize(0.045,"xyz");
  s->SetTitleOffset(1.10,"X"); s->SetTitleOffset(1.10,"Y");
  // Palette
  s->SetNumberContours(CFG::N_CONTOUR);
  gROOT->SetStyle("FigBStyle");
  gROOT->ForceStyle();
}

static const char* flav_label(const std::string& f){
  if(f=="numu")    return "#nu_{#mu}";
  if(f=="numubar") return "#bar{#nu}_{#mu}";
  if(f=="nue")     return "#nu_{e}";
  if(f=="nuebar")  return "#bar{#nu}_{e}";
  return f.c_str();
}

// ------------------------------ Helpers -------------------------------------
static TH2D* fetch_Etheta_2D(TFile& f, const std::string& flav){
  std::string path = flav + "/Detsmear/" + flav + "_CV_AV_TPC_2D";
  TH2D* h = dynamic_cast<TH2D*>(f.Get(path.c_str()));
  if(!h){
    ::printf(" [FigB] MISSING: %s\n", path.c_str());
    return nullptr;
  }
  TH2D* c = (TH2D*)h->Clone(("cl_"+path).c_str());
  c->SetDirectory(nullptr);
  return c;
}

struct Ranges {
  double th_lo, th_hi;   // common θ range (intersection)
  double zmin_pos, zmax; // over E∈[0.25,8] GeV and θ∈[th_lo,th_hi]
  bool have_pos;
};

static Ranges compute_common_ranges(const std::vector<TH2D*>& hs){
  Ranges R;
  R.th_lo   = -std::numeric_limits<double>::infinity();
  R.th_hi   =  std::numeric_limits<double>::infinity();
  R.zmin_pos=  std::numeric_limits<double>::infinity();
  R.zmax    =  0.0;
  R.have_pos= false;

  // First, θ intersection
  for(auto* h: hs){
    if(!h) continue;
    R.th_lo = std::max(R.th_lo, h->GetYaxis()->GetXmin());
    R.th_hi = std::min(R.th_hi, h->GetYaxis()->GetXmax());
  }
  if(!(R.th_lo < R.th_hi)){ // fallback: take something from first non-null
    for(auto* h: hs){ if(h){ R.th_lo=h->GetYaxis()->GetXmin(); R.th_hi=h->GetYaxis()->GetXmax(); break; } }
  }

  // Then, Z scale over requested E frame and common θ
  for(auto* h: hs){
    if(!h) continue;
    int ix1 = std::max(1, h->GetXaxis()->FindFixBin(CFG::FRAME_E_MIN + 1e-9));
    int ix2 = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(CFG::FRAME_E_MAX - 1e-9));
    int iy1 = std::max(1, h->GetYaxis()->FindFixBin(R.th_lo + 1e-9));
    int iy2 = std::min(h->GetNbinsY(), h->GetYaxis()->FindFixBin(R.th_hi - 1e-9));
    for(int ix=ix1; ix<=ix2; ++ix){
      for(int iy=iy1; iy<=iy2; ++iy){
        double v = h->GetBinContent(ix,iy);
        if(v>0.0){ R.have_pos=true; if(v<R.zmin_pos) R.zmin_pos=v; }
        if(v>R.zmax) R.zmax=v;
      }
    }
  }
  if(!std::isfinite(R.zmin_pos)) R.zmin_pos = 1e-30;
  if(R.zmax<=0.0) R.zmax = 1.0;
  return R;
}

static void draw_one_plot(const char* mode,
                          const std::string& flav,
                          TH2D* H,
                          const Ranges& R){
  // Canvas (white, no grey), margins per CFG
  TCanvas c(Form("c_%s_%s",mode,flav.c_str()),
            Form("(E_{#nu}, #theta_{#nu}) — %s, %s", mode, flav.c_str()),
            900, 750);
  c.SetFillColor(0); c.SetBorderMode(0);

  gPad->SetLeftMargin(CFG::PAD_LEFT_MARG);
  gPad->SetRightMargin(CFG::PAD_RIGHT_MARG);
  gPad->SetTopMargin(CFG::PAD_TOP_MARG);
  gPad->SetBottomMargin(CFG::PAD_BOTTOM_MARG);
  if(CFG::USE_LOGZ) gPad->SetLogz();

  // Frame to force X=0.25–8 GeV and θ common range across all plots
  TH2D frame("frame","",10, CFG::FRAME_E_MIN, CFG::FRAME_E_MAX,
                          10, R.th_lo, R.th_hi);
  frame.GetXaxis()->SetTitle("E_{#nu} [GeV]");
  frame.GetYaxis()->SetTitle("#theta_{#nu} [rad]");
  frame.GetZaxis()->SetTitle("Flux / 6 #times 10^{20} POT / (GeV#timesrad) / cm^{2}");
  frame.Draw("AXIS"); // white background, no fill

  if(H){
    H->SetContour(CFG::N_CONTOUR);
    if(CFG::USE_LOGZ){
      double zmin = std::max(std::ldexp(R.zmin_pos, -2), 1e-30); // ~ zmin/4
      H->SetMinimum(zmin);
    } else {
      H->SetMinimum(0.0);
    }
    H->SetMaximum(R.zmax);
    H->Draw("COLZ SAME"); // draw on top of the frame (palette on the right)
  } else {
    TLatex miss; miss.SetTextFont(42); miss.SetTextSize(0.046);
    miss.DrawLatexNDC(0.20, 0.88, "missing histogram");
  }

  // Flavor watermark (top-right, safely left of the palette)
  TLatex lab; lab.SetTextFont(42); lab.SetTextSize(0.050); lab.SetTextAlign(33); // right-top
  const double xlab = 1.0 - gPad->GetRightMargin() - 0.02;
  const double ylab = 1.0 - gPad->GetTopMargin()  - 0.02;
  lab.DrawLatexNDC(xlab, ylab, flav_label(flav));

  // Save
  c.Update();
  c.Print(Form("figB_Etheta_%s_%s.pdf", mode, flav.c_str()));
}

// ------------------------------ Driver --------------------------------------
void flux_fig_B_EvTheta_individual(const char* fhc_file,
                                   const char* rhc_file)
{
  set_minimal_white_style();

  TFile fFHC(fhc_file,"READ");
  TFile fRHC(rhc_file,"READ");
  if(fFHC.IsZombie()){ printf("[FigB] Cannot open FHC file: %s\n", fhc_file); return; }
  if(fRHC.IsZombie()){ printf("[FigB] Cannot open RHC file: %s\n", rhc_file); return; }

  // Collect all histograms first to compute global θ and Z ranges (ROI: E in [0.25, 8])
  std::vector<TH2D*> all;
  for(const auto& flav : CFG::FLAVS) all.push_back(fetch_Etheta_2D(fFHC,flav));
  for(const auto& flav : CFG::FLAVS) all.push_back(fetch_Etheta_2D(fRHC,flav));

  Ranges R = compute_common_ranges(all);

  // Clean up the scan clones (we'll refetch fresh ones for drawing to avoid double ownership)
  for(auto* h : all) delete h;
  all.clear();

  // Draw one PDF per flavor per mode
  for(const auto& flav : CFG::FLAVS){
    TH2D* hF = fetch_Etheta_2D(fFHC,flav);
    draw_one_plot("FHC", flav, hF, R);
    delete hF;

    TH2D* hR = fetch_Etheta_2D(fRHC,flav);
    draw_one_plot("RHC", flav, hR, R);
    delete hR;
  }

  fFHC.Close(); fRHC.Close();
}

void flux_fig_B_EvTheta_individual(){
  flux_fig_B_EvTheta_individual(CFG::FILE_FHC, CFG::FILE_RHC);
}
