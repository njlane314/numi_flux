// ============================================================================
// flux_fig_B_EvTheta_minimal.C
//
// Figure B — Energy–angle (E\u03bd, \u03b8\u03bd) heatmaps per flavor, FHC & RHC.
// Why: makes the energy–angle correlation explicit; motivates any acceptance
//       in angle and the analysis\u2019 angle variable.
// Context to cite: These are the 2D analogues of the 1D spectra used to
//       normalise in Sec. 2.2 of your note; they connect directly to the
//       integrated fluxes you tabulate for 0.25\u201310 GeV (Table 1).
//
// What is drawn:
//   For each mode in {FHC, RHC} and each flavor in
//     <flav> \u2208 {numu, numubar, nue, nuebar},
//   draw  <flav>/Detsmear/<flav>_CV_AV_TPC_2D  with COLZ, logz,
//   with *identical X/Y axis ranges across all flavors and both modes*.
//   The Z (color) scale is unified across all 8 heatmaps so colors are comparable.
//
// Defaults:
//   \u2022 X cropped to [0.25, 10] GeV (if that interval lies inside all hist axes);
//     otherwise the code falls back to the common intersection across hist axes.
//   \u2022 \u03b8 range chosen as the common intersection across all hist axes.
//   \u2022 logz on; 255 contours.
//
// I/O (defaults match the file structure you showed):
//   FHC: "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root"
//   RHC: "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root"
//
// Usage:
//   root -l -b -q flux_fig_B_EvTheta_minimal.C
//   root -l -b -q 'flux_fig_B_EvTheta_minimal.C("FHC.root","RHC.root")'
//
// Outputs:
//   uboone_figB_EvTheta_FHC.pdf
//   uboone_figB_EvTheta_RHC.pdf
// ============================================================================

#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TH2.h"
#include "TH2D.h"
#include "TLatex.h"

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <cmath>

namespace CFG {
  // Exact files you showed:
  const char* FILE_FHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* FILE_RHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";

  // Flavors to process
  static const std::vector<std::string> FLAVS = {"numu","numubar","nue","nuebar"};

  // Enforce the analysis window used in the note (Table 1):
  constexpr bool   CROP_X_TO_025_10 = true;
  constexpr double X_CROP_MIN_GEV   = 0.25;
  constexpr double X_CROP_MAX_GEV   = 10.0;

  // Drawing
  constexpr bool LOGZ = true;
  constexpr int  NCONT = 255;
}

// Minimal, legible global style (matches your font/spacing choices)
static void set_global_style(){
  const int f=42;
  TStyle* s=new TStyle("FluxFigStyle","FluxFig Style");
  s->SetTitleFont(f,"X"); s->SetTitleFont(f,"Y"); s->SetTitleFont(f,"Z");
  s->SetTitleSize(0.045,"X"); s->SetTitleSize(0.045,"Y"); s->SetTitleSize(0.045,"Z");
  s->SetLabelFont(f,"X"); s->SetLabelFont(f,"Y"); s->SetLabelFont(f,"Z");
  s->SetLabelSize(0.045,"X"); s->SetLabelSize(0.045,"Y"); s->SetLabelSize(0.045,"Z");
  s->SetLabelOffset(0.005,"X"); s->SetLabelOffset(0.005,"Y"); s->SetLabelOffset(0.005,"Z");
  s->SetTitleOffset(1.10,"X"); s->SetTitleOffset(1.10,"Y");
  s->SetOptStat(0); s->SetOptTitle(0);
  s->SetPadTickX(1); s->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  s->SetPadLeftMargin(0.15); s->SetPadRightMargin(0.14);
  s->SetPadTopMargin(0.07);  s->SetPadBottomMargin(0.12);
  gROOT->SetStyle("FluxFigStyle"); gROOT->ForceStyle();
}

struct HRef {
  std::string flav;
  std::string mode; // "FHC" or "RHC"
  TH2D* hist = nullptr; // owned clone
};

static std::string label_for(const std::string& flav){
  if(flav=="numu")    return "#nu_{#mu}";
  if(flav=="numubar") return "#bar{#nu}_{#mu}";
  if(flav=="nue")     return "#nu_{e}";
  if(flav=="nuebar")  return "#bar{#nu}_{e}";
  return flav.c_str();
}

// Compute common X/Y intersection across all hists; optionally crop X to [0.25,10] GeV.
// Also compute a global positive Z-min and Z-max inside that ROI (for unified color scale).
static void compute_common_ranges(const std::vector<HRef>& H,
                                  double& xlo, double& xhi,
                                  double& ylo, double& yhi,
                                  double& zmin_pos, double& zmax){
  // 1) X/Y intersection of histogram axes
  double xlo_int = -std::numeric_limits<double>::infinity();
  double xhi_int =  std::numeric_limits<double>::infinity();
  double ylo_int = -std::numeric_limits<double>::infinity();
  double yhi_int =  std::numeric_limits<double>::infinity();

  for(const auto& h : H){
    if(!h.hist) continue;
    xlo_int = std::max(xlo_int, h.hist->GetXaxis()->GetXmin());
    xhi_int = std::min(xhi_int, h.hist->GetXaxis()->GetXmax());
    ylo_int = std::max(ylo_int, h.hist->GetYaxis()->GetXmin());
    yhi_int = std::min(yhi_int, h.hist->GetYaxis()->GetXmax());
  }

  // 2) Optionally crop X to [0.25,10] GeV but only if that lies inside the intersection
  if(CFG::CROP_X_TO_025_10){
    double want_lo = CFG::X_CROP_MIN_GEV;
    double want_hi = CFG::X_CROP_MAX_GEV;
    // if desired window is fully within intersection, apply it; otherwise keep intersection
    if(want_lo < xhi_int && want_hi > xlo_int){
      xlo_int = std::max(xlo_int, want_lo);
      xhi_int = std::min(xhi_int, want_hi);
    }
  }

  // Guard: if something went wrong, fallback to a sane default from the first hist
  if(!(xlo_int < xhi_int) || !(ylo_int < yhi_int)){
    for(const auto& h : H){
      if(h.hist){
        xlo_int = h.hist->GetXaxis()->GetXmin();
        xhi_int = h.hist->GetXaxis()->GetXmax();
        ylo_int = h.hist->GetYaxis()->GetXmin();
        yhi_int = h.hist->GetYaxis()->GetXmax();
        break;
      }
    }
  }

  xlo = xlo_int; xhi = xhi_int; ylo = ylo_int; yhi = yhi_int;

  // 3) Global Z range over the chosen ROI (positive min for logz)
  zmin_pos = std::numeric_limits<double>::infinity();
  zmax     = 0.0;

  for(const auto& hr : H){
    if(!hr.hist) continue;
    auto* h = hr.hist;
    int ix1 = std::max(1, h->GetXaxis()->FindFixBin(xlo + 1e-9));
    int ix2 = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xhi - 1e-9));
    int iy1 = std::max(1, h->GetYaxis()->FindFixBin(ylo + 1e-9));
    int iy2 = std::min(h->GetNbinsY(), h->GetYaxis()->FindFixBin(yhi - 1e-9));
    for(int ix=ix1; ix<=ix2; ++ix){
      for(int iy=iy1; iy<=iy2; ++iy){
        double z = h->GetBinContent(ix,iy);
        if(z>0.0 && z<zmin_pos) zmin_pos = z;
        if(z>zmax) zmax = z;
      }
    }
  }

  if(!std::isfinite(zmin_pos)) zmin_pos = 1e-30; // safe tiny positive for logz
  if(!(zmax>0.0)) zmax = 1.0;
  // Optional: pad the range a bit
  // Keep zmin_pos as the smallest positive; inflate zmax to give headroom
  zmax *= 1.05;
}

static TH2D* fetch_clone(TFile& f, const std::string& flav){
  std::string path = flav + "/Detsmear/" + flav + "_CV_AV_TPC_2D";
  TH2D* h = dynamic_cast<TH2D*>( f.Get(path.c_str()) );
  if(!h) {
    ::printf("  [warn] Missing: %s\n", path.c_str());
    return nullptr;
  }
  TH2D* c = (TH2D*)h->Clone(Form("cl_%s", path.c_str()));
  c->SetDirectory(0);
  c->SetContour(CFG::NCONT);
  c->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  c->GetYaxis()->SetTitle("#theta_{#nu} [rad]");
  c->GetZaxis()->SetTitle("Flux / 6 #times 10^{20} POT / (GeV#timesrad) / cm^{2}");
  return c;
}

static void draw_mode_canvas(const char* mode,
                             const std::vector<HRef>& H_mode,
                             double xlo, double xhi, double ylo, double yhi,
                             double zmin, double zmax){
  TCanvas c(Form("c_figB_%s",mode),
            Form("(E_{#nu}, #theta_{#nu}) heatmaps per flavor \u2014 %s",mode),
            1400, 1100);
  c.Divide(2,2);

  // order: numu, numubar, nue, nuebar
  for(size_t i=0;i<H_mode.size();++i){
    c.cd((int)i+1);
    gPad->SetRightMargin(0.14);
    if(CFG::LOGZ) gPad->SetLogz();

    TH2D* H = H_mode[i].hist;
    if(!H){ TLatex t; t.SetTextFont(42); t.SetTextSize(0.05);
            t.DrawLatexNDC(0.2,0.5,"[missing histogram]"); continue; }

    H->GetXaxis()->SetRangeUser(xlo, xhi);
    H->GetYaxis()->SetRangeUser(ylo, yhi);
    H->SetMinimum(zmin);
    H->SetMaximum(zmax);
    H->Draw("COLZ");

    // Flavor label in-pad (top-left)
    TLatex lab;
    lab.SetTextFont(42);
    lab.SetTextSize(0.052);
    lab.DrawLatexNDC(0.16, 0.92, label_for(H_mode[i].flav).c_str());
  }

  c.Update();
  c.Print(Form("uboone_figB_EvTheta_%s.pdf",mode));
}

void flux_fig_B_EvTheta_minimal(const char* fhc_file = CFG::FILE_FHC,
                                const char* rhc_file = CFG::FILE_RHC){
  set_global_style();

  // Open files
  TFile fFHC(fhc_file,"READ");
  TFile fRHC(rhc_file,"READ");
  if(fFHC.IsZombie()){ ::printf("[error] Cannot open FHC file: %s\n", fhc_file); return; }
  if(fRHC.IsZombie()){ ::printf("[error] Cannot open RHC file: %s\n", rhc_file); return; }

  // Gather all histograms (both modes) for uniform ranges
  std::vector<HRef> H_all; H_all.reserve(8);
  std::vector<HRef> H_FHC, H_RHC; H_FHC.reserve(4); H_RHC.reserve(4);

  for(const auto& flav : CFG::FLAVS){
    if(TH2D* h = fetch_clone(fFHC, flav)) { H_FHC.push_back({flav,"FHC",h}); H_all.push_back({flav,"FHC",h}); }
    else                                   { H_FHC.push_back({flav,"FHC",nullptr}); }
  }
  for(const auto& flav : CFG::FLAVS){
    if(TH2D* h = fetch_clone(fRHC, flav)) { H_RHC.push_back({flav,"RHC",h}); H_all.push_back({flav,"RHC",h}); }
    else                                   { H_RHC.push_back({flav,"RHC",nullptr}); }
  }

  if(H_all.empty()){ ::printf("[error] No 2D histograms found.\n"); return; }

  // Compute COMMON ranges across *all* flavors & *both* modes
  double xlo, xhi, ylo, yhi, zmin, zmax;
  compute_common_ranges(H_all, xlo, xhi, ylo, yhi, zmin, zmax);

  ::printf("[info] Using common ranges:\n");
  ::printf("       E_nu:  [%.6g, %.6g] GeV\n", xlo, xhi);
  ::printf("       theta: [%.6g, %.6g] rad\n", ylo, yhi);
  ::printf("       Z(min>0,max): [%.6g, %.6g]\n", zmin, zmax);

  // Draw each mode with the same X/Y/Z ranges
  draw_mode_canvas("FHC", H_FHC, xlo, xhi, ylo, yhi, zmin, zmax);
  draw_mode_canvas("RHC", H_RHC, xlo, xhi, ylo, yhi, zmin, zmax);

  // Cleanup clones
  for(auto& hr : H_all) delete hr.hist;

  fFHC.Close(); fRHC.Close();
}
