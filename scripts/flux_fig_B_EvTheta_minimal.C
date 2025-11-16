// ============================================================================
// flux_fig_B_EvTheta_minimal.C  (updated)
// - Watermark removed
// - Flavor label moved to top-right corner of the plot
// - Canvas width computed so the drawable area is ~square
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
#include "TPaletteAxis.h"

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

  // Layout: room for the palette on the right
  constexpr double PAD_LEFT   = 0.15;
  constexpr double PAD_RIGHT  = 0.22;
  constexpr double PAD_TOP    = 0.07;
  constexpr double PAD_BOTTOM = 0.12;

  // Canvas height in pixels; width will be chosen to square the drawable area
  constexpr int CANVAS_H = 750;
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
  s->SetPadLeftMargin(CFG::PAD_LEFT);
  s->SetPadRightMargin(CFG::PAD_RIGHT);
  s->SetPadTopMargin(CFG::PAD_TOP);
  s->SetPadBottomMargin(CFG::PAD_BOTTOM);

  s->SetCanvasColor(kWhite);
  s->SetPadColor(kWhite);
  s->SetFrameFillColor(kWhite);
  s->SetFrameFillStyle(0);
  s->SetStatColor(kWhite);
  s->SetCanvasBorderMode(0);
  s->SetPadBorderMode(0);
  s->SetFrameBorderMode(0);
  s->SetTitleOffset(1.25,"Z");
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
  return flav;
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
  c->GetZaxis()->CenterTitle(true);
  return c;
}

// One-pad, per-flavor canvas (shared X/Y/Z ranges)
static void draw_single_canvas(const HRef& hr,
                               double xlo, double xhi,
                               double ylo, double yhi,
                               double zmin, double zmax)
{
  // Choose canvas width so the drawable area is ~square:
  const double wfrac = 1.0 - CFG::PAD_LEFT - CFG::PAD_RIGHT;
  const double hfrac = 1.0 - CFG::PAD_TOP  - CFG::PAD_BOTTOM;
  const int W = int( CFG::CANVAS_H * (hfrac / wfrac) + 0.5 );
  const int H = CFG::CANVAS_H;

  TCanvas c(Form("c_single_%s_%s", hr.mode.c_str(), hr.flav.c_str()),
            Form("(E_{#nu}, #theta_{#nu}) — %s, %s",
                 hr.mode.c_str(), hr.flav.c_str()),
            W, H);
  c.cd();
  c.SetFillColor(kWhite);
  c.SetBorderMode(0);
  gPad->SetFillColor(kWhite);
  gPad->SetLeftMargin(CFG::PAD_LEFT);
  gPad->SetRightMargin(CFG::PAD_RIGHT);
  gPad->SetTopMargin(CFG::PAD_TOP);
  gPad->SetBottomMargin(CFG::PAD_BOTTOM);
  gPad->SetBorderMode(0);
  if (CFG::LOGZ) gPad->SetLogz();

  if (!hr.hist) {
    TLatex t; t.SetTextFont(42); t.SetTextSize(0.05);
    t.DrawLatexNDC(0.20, 0.50, "[missing histogram]");
    c.Print(Form("uboone_figB_EvTheta_%s_%s.pdf",
                 hr.mode.c_str(), hr.flav.c_str()));
    return;
  }

  TH2D* H2 = hr.hist;
  H2->GetXaxis()->SetRangeUser(xlo, xhi);
  H2->GetYaxis()->SetRangeUser(ylo, yhi);
  H2->SetMinimum(zmin);
  H2->SetMaximum(zmax);
  H2->Draw("COLZ");

  gPad->Update();
  if (TPaletteAxis* pal =
        dynamic_cast<TPaletteAxis*>(H2->GetListOfFunctions()->FindObject("palette"))) {
    const double x1 = 1.0 - gPad->GetRightMargin() + 0.02;
    const double x2 = x1 + 0.06;
    const double y1 = gPad->GetBottomMargin();
    const double y2 = 1.0 - gPad->GetTopMargin();
    pal->SetX1NDC(x1);
    pal->SetX2NDC(x2);
    pal->SetY1NDC(y1);
    pal->SetY2NDC(y2);
    pal->SetLabelSize(0.035);
    pal->SetTitleSize(0.042);
    pal->SetTitleOffset(1.10);
  }

  // Flavor label — now top-right, inside the plotting area (black)
  TLatex lab; lab.SetTextFont(42); lab.SetTextSize(0.052); lab.SetTextAlign(33);
  const double xr = 1.0 - gPad->GetRightMargin() - 0.01;
  const double yt = 1.0 - gPad->GetTopMargin()   - 0.01;
  lab.DrawLatexNDC(xr, yt, label_for(hr.flav).c_str());
  c.Update();
  c.Print(Form("uboone_figB_EvTheta_%s_%s.pdf",
               hr.mode.c_str(), hr.flav.c_str()));
}

void flux_fig_B_EvTheta_minimal(){
  set_global_style();

  // Open files
  TFile fFHC(CFG::FILE_FHC,"READ");
  TFile fRHC(CFG::FILE_RHC,"READ");
  if(fFHC.IsZombie()){ ::printf("[error] Cannot open FHC file: %s\n", CFG::FILE_FHC); return; }
  if(fRHC.IsZombie()){ ::printf("[error] Cannot open RHC file: %s\n", CFG::FILE_RHC); return; }

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

  // Save one PDF per flavor per mode (identical X/Y/Z across all)
  for (const auto& hr : H_FHC) draw_single_canvas(hr, xlo, xhi, ylo, yhi, zmin, zmax);
  for (const auto& hr : H_RHC) draw_single_canvas(hr, xlo, xhi, ylo, yhi, zmin, zmax);

  // Cleanup clones
  for(auto& hr : H_all) delete hr.hist;

  fFHC.Close(); fRHC.Close();
}
