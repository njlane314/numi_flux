// ============================================================================
// uboone_numu_parentE_angleStyle.C
//
// Draw \u03bd\u03bc energy spectra broken down by parent (\u03c0+, \u03c0-, K+, K-, \u03bc+, \u03bc-, K0L)
// using the same visual style as your angle spectra (split canvas + stacked
// legend, log-y, 1.2 line width, etc.).
// NOTE: All histograms are converted to a flux *density* (divide by bin width).
// ============================================================================

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <map>
#include <string>
#include <vector>

// ------------------------------ CONFIG --------------------------------------

namespace CFG {
  // Your exact files:
  const char* FILE_FHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* FILE_RHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";

  // For legend percentages:
  constexpr double E_FRAC_MIN = 0.060; // GeV (60 MeV)
}

// ------------------------------ STYLE (your angle style) --------------------

static void set_global_style(){
  const int f=42;
  TStyle* s=new TStyle("PlotterStyle","Plotter Style");
  s->SetTitleFont(f,"X"); s->SetTitleFont(f,"Y"); s->SetTitleFont(f,"Z");
  s->SetTitleSize(0.04,"X"); s->SetTitleSize(0.04,"Y"); s->SetTitleSize(0.05,"Z");
  s->SetLabelFont(f,"X"); s->SetLabelFont(f,"Y"); s->SetLabelFont(f,"Z");
  s->SetLabelSize(0.045,"X"); s->SetLabelSize(0.045,"Y"); s->SetLabelSize(0.045,"Z");
  s->SetLabelOffset(0.005,"X"); s->SetLabelOffset(0.005,"Y"); s->SetLabelOffset(0.005,"Z");
  s->SetTitleOffset(1.10,"X"); s->SetTitleOffset(1.10,"Y");
  s->SetOptStat(0); s->SetOptTitle(0);
  s->SetPadTickX(1); s->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  s->SetPadLeftMargin(0.15); s->SetPadRightMargin(0.05);
  s->SetPadTopMargin(0.07); s->SetPadBottomMargin(0.12);
  s->SetMarkerSize(1.0);
  s->SetCanvasColor(0); s->SetPadColor(0); s->SetFrameFillColor(0);
  s->SetCanvasBorderMode(0); s->SetPadBorderMode(0); s->SetStatColor(0); s->SetFrameBorderMode(0);
  s->SetTitleFillColor(0); s->SetTitleBorderSize(0);
  gROOT->SetStyle("PlotterStyle"); gROOT->ForceStyle();
}

static void style_line(TH1* h,int col,int ls, double lw=1.2){
  h->SetLineColor(col); h->SetLineStyle(ls); h->SetLineWidth(lw); h->SetMarkerSize(0);
}

static double integral_in(double xmin, double xmax, const TH1* h, bool width=false){
  int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
  int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
  return width ? h->Integral(bmin, bmax, "width") : h->Integral(bmin, bmax);
}

static void auto_logy_limits_range(TH1* frame, std::initializer_list<TH1*> hs, double xmin, double xmax){
  double mn=std::numeric_limits<double>::infinity(), mx=0.0;
  for(TH1* h:hs){
    int bmin=std::max(1,h->GetXaxis()->FindFixBin(xmin+1e-9));
    int bmax=std::min(h->GetNbinsX(),h->GetXaxis()->FindFixBin(xmax-1e-9));
    for(int b=bmin;b<=bmax;++b){
      double y=h->GetBinContent(b);
      if(y>0&&y<mn) mn=y; if(y>mx) mx=y;
    }
  }
  if(!std::isfinite(mn)) mn=1e-18; if(mx<=0.0) mx=1.0;
  frame->SetMinimum(std::max(1e-30,mn*0.8)); frame->SetMaximum(mx*6.0);
}

static void make_split_canvas(const char* cname, const char* ctitle, double split,
                              TCanvas*& canv, TPad*& p_main, TPad*& p_leg,
                              bool logy=true){
  canv = new TCanvas(cname, ctitle, 1200, 700);
  p_main = new TPad("pad_main","pad_main",0.,0.00,1.,split);
  p_leg  = new TPad("pad_legend","pad_legend",0.,split,1.,1.00);
  p_main->SetTopMargin(0.01); p_main->SetBottomMargin(0.12);
  p_main->SetLeftMargin(0.12); p_main->SetRightMargin(0.05);
  if(logy) p_main->SetLogy();
  p_leg->SetTopMargin(0.05); p_leg->SetBottomMargin(0.01);
  p_leg->SetLeftMargin(0.02); p_leg->SetRightMargin(0.02);
  p_main->Draw(); p_leg->Draw();
}

static TLegend* build_legend_like_stacked(TPad* p_leg,
    const std::vector<std::pair<TH1*, TString>>& items,
    const std::vector<double>& sums, double split, bool show_pct=true){
  p_leg->cd();
  TLegend* L=new TLegend(0.12,0.00,0.95,0.75);
  L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42);
  int n_entries=(int)items.size(); int n_cols=(n_entries>5)?3:2;
  L->SetNColumns(n_cols); L->SetColumnSeparation(0.08); L->SetEntrySeparation(0.00); L->SetMargin(0.25);
  const double s_main=0.045; const double s_leg=s_main*(split/(1.0-split)); L->SetTextSize(s_leg);
  double s_tot=0.0; for(double s:sums) s_tot+=s;
  for(size_t i=0;i<items.size();++i){
    TString lab=items[i].second;
    if(show_pct && s_tot>0 && i<sums.size()) lab = Form("%s (%.2f%%)", lab.Data(), 100.0*sums[i]/s_tot);
    L->AddEntry(items[i].first, lab, "l");
  }
  return L;
}

// Helper: convert a histogram to density (divide by \u0394E of each bin).
static void make_density(TH1* h){
  if(!h) return;
  // ROOT will divide each bin by its width and properly scale errors.
  h->Scale(1.0, "width");
}

// Helper: fetch the best-available \u03bd\u03bc CV energy spectrum (prefer uniform binning)
static TH1D* cv_energy_numu(TFile& f){
  // Prefer uniformly binned version if available, fall back to legacy (often variable-width)
  TH1D* h = (TH1D*)f.Get("numu/Detsmear/numu_CV_AV_TPC_5MeV_bin");
  if(!h)   h = (TH1D*)f.Get("numu/Detsmear/numu_CV_AV_TPC");
  if(!h)   return nullptr;
  h=(TH1D*)h->Clone("cv_numu"); h->SetDirectory(0);
  // Do NOT rebin here: parents may have different native binning; we will convert
  // everything to a density so comparisons are meaningful regardless of \u0394E.
  // If you later want to make the total less "steppy", rebin after cloning and
  // before make_density().
  return h;
}

// ------------------------------ DRAWER --------------------------------------

static void draw_numu_parent_energy(const char* mode, TFile& f){
  // Total \u03bd\u03bc flux (CV)
  TH1D* Htot = cv_energy_numu(f);
  if(!Htot){ printf("[%s] Cannot find numu CV energy histogram.\n",mode); return; }

  // Parent contributions (individual signs)
  auto fetch = [&](const char* path)->TH1D*{
    TH1D* h = (TH1D*)f.Get(path);
    if(!h)  { printf("[%s] Missing: %s\n",mode,path); return nullptr; }
    h=(TH1D*)h->Clone(Form("cl_%s",path)); h->SetDirectory(0); return h;
  };

  TH1D* h_piP = fetch("numu/PI_Plus/Enu_numu_PI_Plus_AV_TPC");
  TH1D* h_piM = fetch("numu/PI_Minus/Enu_numu_PI_Minus_AV_TPC");
  TH1D* h_kP  = fetch("numu/Kaon_Plus/Enu_numu_Kaon_Plus_AV_TPC");
  TH1D* h_kM  = fetch("numu/Kaon_Minus/Enu_numu_Kaon_Minus_AV_TPC");
  TH1D* h_muP = fetch("numu/Mu_Plus/Enu_numu_Mu_Plus_AV_TPC");
  TH1D* h_muM = fetch("numu/Mu_Minus/Enu_numu_Mu_Minus_AV_TPC");
  TH1D* h_KL  = fetch("numu/K0L/Enu_numu_K0L_AV_TPC");

  if(!h_piP||!h_piM||!h_kP||!h_kM||!h_muP||!h_muM||!h_KL){
    printf("[%s] One or more parent histograms are missing \u2014 aborting.\n",mode);
    delete Htot; return;
  }

  // --- Normalize all histograms to *density* (per GeV) so variable-width bins
  //     (e.g. CV total) are comparable to fine-binned parent spectra.
  for(TH1* h : std::vector<TH1*>{Htot,h_piP,h_piM,h_kP,h_kM,h_muP,h_muM,h_KL}){
    make_density(h);
  }

  // Styling: switch to a more vibrant palette and SOLID lines only (no dashed)
  int C_piP = TColor::GetColor("#ff1744"); // bright red
  int C_piM = TColor::GetColor("#ff6d00"); // vivid orange
  int C_kP  = TColor::GetColor("#2979ff"); // bright blue
  int C_kM  = TColor::GetColor("#00c853"); // vivid green
  int C_muP = TColor::GetColor("#7c4dff"); // vivid violet
  int C_muM = TColor::GetColor("#aa00ff"); // vivid magenta
  int C_KL  = TColor::GetColor("#00e5ff"); // bright cyan
  int C_tot = TColor::GetColor("#616161"); // dark gray for total (not black)

  style_line(h_piP, C_piP, 1);
  style_line(h_piM, C_piM, 1);
  style_line(h_kP , C_kP , 1);
  style_line(h_kM , C_kM , 1);
  style_line(h_muP, C_muP, 1);
  style_line(h_muM, C_muM, 1);
  style_line(h_KL , C_KL , 1);
  style_line(Htot , C_tot, 1, /*lw=*/1.8); // keep total slightly thicker

  // Canvas in your stacked/legend style
  const double split=0.85;
  TCanvas* c=nullptr; TPad* p_main=nullptr; TPad* p_leg=nullptr;
  make_split_canvas(Form("c_numu_parentE_%s",mode),
                    Form("#nu_{#mu} parent breakdown \u2014 %s",mode),
                    split,c,p_main,p_leg,/*logy=*/true);

  // Axes & ranges
  p_main->cd();
  // Keep full histogram span for computations, but PLOT only 0–10 GeV
  double xmin_full = Htot->GetXaxis()->GetXmin();
  double xmax_full = Htot->GetXaxis()->GetXmax();
  double xmin_plot = 0.0, xmax_plot = 10.0;
  TH1D* frame=new TH1D(Form("frame_numu_parentE_%s",mode),"",100,xmin_plot,xmax_plot);
  frame->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  frame->GetYaxis()->SetTitle("Flux / 6 #times 10^{20} POT / GeV / cm^{2}");
  auto_logy_limits_range(frame,{Htot,h_piP,h_piM,h_kP,h_kM,h_muP,h_muM,h_KL},xmin_plot,xmax_plot);
  frame->Draw("AXIS");
  // Restrict each histogram’s visible x-range to [0,10] GeV
  for(TH1* h : std::vector<TH1*>{Htot,h_piP,h_piM,h_kP,h_kM,h_muP,h_muM,h_KL}){
    h->GetXaxis()->SetRangeUser(xmin_plot, xmax_plot);
  }

  // Draw (parents first, then total on top)
  h_piP->Draw("HIST SAME"); h_piM->Draw("HIST SAME");
  h_kP ->Draw("HIST SAME"); h_kM ->Draw("HIST SAME");
  h_muP->Draw("HIST SAME"); h_muM->Draw("HIST SAME");
  h_KL ->Draw("HIST SAME");
  Htot ->Draw("HIST SAME");

  // Percentages for E\u03bd > 60 MeV (width=true integrates GeV to match units)
  const double thr = CFG::E_FRAC_MIN;
  double SpiP = integral_in(thr,xmax_full,h_piP,true);
  double SpiM = integral_in(thr,xmax_full,h_piM,true);
  double SkP  = integral_in(thr,xmax_full,h_kP ,true);
  double SkM  = integral_in(thr,xmax_full,h_kM ,true);
  double SmuP = integral_in(thr,xmax_full,h_muP,true);
  double SmuM = integral_in(thr,xmax_full,h_muM,true);
  double SKL  = integral_in(thr,xmax_full,h_KL ,true);

  // Legend (put "Total Flux" last so it has no percentage)
  p_leg->cd();
  auto L = build_legend_like_stacked(
    p_leg,
    {
      {h_piP,"#pi^{+}"}, {h_piM,"#pi^{-}"},
      {h_kP,"K^{+}"},    {h_kM,"K^{-}"},
      {h_muP,"#mu^{+}"}, {h_muM,"#mu^{-}"},
      {h_KL,"K^{0}_{L}"},
      {Htot,"Total Flux"}
    },
    {SpiP,SpiM,SkP,SkM,SmuP,SmuM,SKL}, // no entry for "Total Flux"
    split, /*show_pct=*/true);
  L->Draw();

  c->Update();
  c->Print(Form("uboone_numu_parentE_%s_angleStyle.pdf",mode));

  // cleanup
  delete L; delete frame; delete c;
  delete Htot; delete h_piP; delete h_piM; delete h_kP; delete h_kM;
  delete h_muP; delete h_muM; delete h_KL;
}

// ------------------------------ DRIVER (no-arg) -----------------------------
// Makes BOTH FHC and RHC \u03bd\u03bc parent-breakdown plots with one call.
void uboone_numu_parentE_angleStyle(){
  set_global_style();

  TFile fFHC(CFG::FILE_FHC,"READ");
  TFile fRHC(CFG::FILE_RHC,"READ");
  if(fFHC.IsZombie()){ printf("Cannot open FHC file: %s\n", CFG::FILE_FHC); return; }
  if(fRHC.IsZombie()){ printf("Cannot open RHC file: %s\n", CFG::FILE_RHC); return; }

  // \u03bd\u03bc parent energy spectra (angle-style legend/layout)
  draw_numu_parent_energy("FHC", fFHC);
  draw_numu_parent_energy("RHC", fRHC);

  fFHC.Close();
  fRHC.Close();
}
