// ============================================================================
// flux_figs_minimal_uboone.C
//
// Figures 5.11–5.16, 5.18–5.19 using your file structure and your legend style.
//
// Directory structure (as you showed):
//   <flav> = {numu, numubar, nue, nuebar}
//
//   <flav>/Detsmear/
//     - <flav>_CV_AV_TPC
//     - <flav>_CV_AV_TPC_5MeV_bin
//     - <flav>_CV_AV_TPC_2D     (E vs angle)
//     - Th_<flav>_CV_TPC        (angle 1D)
//
//   <flav>/{PI_Plus, PI_Minus, Kaon_Plus, Kaon_Minus, Mu_Plus, Mu_Minus, K0L}/
//     - Enu_<flav>_<ParentDir>_AV_TPC   (energy 1D)
//     - Th_<flav>_<ParentDir>_AV_TPC    (angle 1D)
//
//   <flav>/Multisims/
//     - <flav>_ppfx_<category>_Uni_<N>_AV_TPC   (many TH1D, ignore *_2D)
//
// What’s drawn:
//  • Fig. 5.11  : CV flux vs angle (4 flavors) — legend in your style
//  • Fig. 5.12  : angle vs parent decay z (2D) — axes titled (legendless)
//  • Fig. 5.13  : parent energy (π/K/μ/KL) — % for Eν>60 MeV, your legend style
//  • Fig. 5.14  : parent angle (π/K/μ/KL) — % from angle-integrated areas
//  • Fig. 5.15  : (Eν,θ) heatmaps per flavor — axes titled/units, logz
//  • Fig. 5.16  : PPFX CV ± 1σ band from summed cov — your legend style
//  • Fig. 5.18  : PPFX fractional covariance (dimensionless)
//  • Fig. 5.19  : PPFX correlation (dimensionless)
//
// Notes:
//  * Line width is 1.2 (matches your style).
//  * All legends use your stacked style (top pad, scaled text, cols, margins).
//  * Energy-axis Y titles include per-bin MeV when applicable; angle uses /rad.
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
#include "TH2.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <map>
#include <set>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <limits>
#include <initializer_list>

// ------------------------------ CONFIG --------------------------------------

namespace CFG {
  // Exact files you showed:
  const char* FILE_FHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* FILE_RHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";

  // Percentages for parent energy in Fig. 5.13 use this threshold (60 MeV):
  constexpr double E_FRAC_MIN = 0.060; // GeV

  // Flavor list (for loops)
  static const std::vector<std::string> FLAVS = {"numu","numubar","nue","nuebar"};
}

// ------------------------------ STYLE (yours verbatim) ----------------------

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

// ------------------------------ HELPERS -------------------------------------

static void style_line(TH1* h,int col,int ls){ h->SetLineColor(col); h->SetLineStyle(ls); h->SetLineWidth(1.2); h->SetMarkerSize(0); }

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
    for(int b=bmin;b<=bmax;++b){ double y=h->GetBinContent(b); if(y>0&&y<mn) mn=y; if(y>mx) mx=y; }
  }
  if(!std::isfinite(mn)) mn=1e-18; if(mx<=0.0) mx=1.0;
  frame->SetMinimum(std::max(1e-30,mn*0.8)); frame->SetMaximum(mx*6.0);
}

static void make_split_canvas(const char* cname, const char* ctitle, double split,
                              TCanvas*& canv, TPad*& p_main, TPad*& p_leg,
                              bool logy=false, bool logz=false){
  canv = new TCanvas(cname, ctitle, 1200, 700);
  p_main = new TPad("pad_main","pad_main",0.,0.00,1.,split);
  p_leg  = new TPad("pad_legend","pad_legend",0.,split,1.,1.00);
  p_main->SetTopMargin(0.01); p_main->SetBottomMargin(0.12);
  p_main->SetLeftMargin(0.12); p_main->SetRightMargin(0.05);
  if(logy) p_main->SetLogy(); if(logz) p_main->SetLogz();
  p_leg->SetTopMargin(0.05); p_leg->SetBottomMargin(0.01);
  p_leg->SetLeftMargin(0.02); p_leg->SetRightMargin(0.02);
  p_main->Draw(); p_leg->Draw();
}

// Generic legend in your stacked style (percentages optional)
static TLegend* build_legend_like_stacked(TPad* p_leg,
    const std::vector<std::pair<TH1*, TString>>& items,
    const std::vector<double>& sums, double split, bool show_pct=true){
  p_leg->cd();
  TLegend* L=new TLegend(0.12,0.00,0.95,0.75);
  L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42);
  int n_entries=(int)items.size(); int n_cols=(n_entries>4)?3:2;
  L->SetNColumns(n_cols); L->SetColumnSeparation(0.08); L->SetEntrySeparation(0.00); L->SetMargin(0.25);
  const double s_main=0.045; const double s_leg=s_main*(split/(1.0-split)); L->SetTextSize(s_leg);
  double s_tot=0.0; for(double s:sums) s_tot+=s;
  for(size_t i=0;i<items.size();++i){
    TString lab=items[i].second;
    if(show_pct && s_tot>0 && i<sums.size()) lab = Form("%s (%.1f%%)", lab.Data(), 100.0*sums[i]/s_tot);
    L->AddEntry(items[i].first, lab, "l");
  }
  return L;
}

// Search for a TH2D under a flavor directory whose name looks like theta vs z
static TH2D* find_theta_vs_z(TFile& f, const std::string& flav){
  auto* d = f.GetDirectory(flav.c_str());
  if(!d) return nullptr;
  std::vector<TDirectory*> stack = {d};
  while(!stack.empty()){
    TDirectory* cur = stack.back(); stack.pop_back();
    TIter it(cur->GetListOfKeys());
    while(auto* k=(TKey*)it()){
      std::string cls = k->GetClassName();
      std::string nm  = k->GetName();
      std::string low(nm.size(),'\0');
      std::transform(nm.begin(),nm.end(),low.begin(),::tolower);
      if(cls=="TDirectoryFile"){
        stack.push_back( (TDirectory*)cur->Get(nm.c_str()) );
      } else if(cls=="TH2D"){
        bool hasTheta = (low.find("th_")!=std::string::npos) || (low.find("theta")!=std::string::npos);
        bool hasZ     = (low.find("zpos")!=std::string::npos) || (low.find("_z")!=std::string::npos) || (low.find("parentz")!=std::string::npos);
        if(hasTheta && hasZ){
          auto* h=(TH2D*)cur->Get(nm.c_str());
          if(h){ auto* c=(TH2D*)h->Clone(Form("cl_%s",nm.c_str())); c->SetDirectory(0); return c; }
        }
      }
    }
  }
  return nullptr;
}

// ------------------------------ FIGURE 5.11 ---------------------------------
// CV flux vs angle (4 flavors), with your legend style
static void fig_5_11_onefile(const char* mode, TFile& f){
  auto* h_numu    = (TH1D*)f.Get("numu/Detsmear/Th_numu_CV_TPC");
  auto* h_numubar = (TH1D*)f.Get("numubar/Detsmear/Th_numubar_CV_TPC");
  auto* h_nue     = (TH1D*)f.Get("nue/Detsmear/Th_nue_CV_TPC");
  auto* h_nuebar  = (TH1D*)f.Get("nuebar/Detsmear/Th_nuebar_CV_TPC");
  if(!h_numu||!h_numubar||!h_nue||!h_nuebar){ printf("[5.11/%s] missing Th_*_CV_TPC\n",mode); return; }
  h_numu   =(TH1D*)h_numu->Clone();    h_numu->SetDirectory(0);
  h_numubar=(TH1D*)h_numubar->Clone(); h_numubar->SetDirectory(0);
  h_nue    =(TH1D*)h_nue->Clone();     h_nue->SetDirectory(0);
  h_nuebar =(TH1D*)h_nuebar->Clone();  h_nuebar->SetDirectory(0);

  int CR=TColor::GetColor("#e41a1c"), CB=TColor::GetColor("#1f78b4");
  style_line(h_numu,CR,1); style_line(h_nue,CR,2); style_line(h_numubar,CB,1); style_line(h_nuebar,CB,3);

  const double split=0.85;
  TCanvas* c=nullptr; TPad* p_main=nullptr; TPad* p_leg=nullptr;
  make_split_canvas(Form("c_511_%s",mode),Form("Flux vs angle — %s",mode),split,c,p_main,p_leg,/*logy=*/true);

  p_main->cd();
  double xmin=h_numu->GetXaxis()->GetXmin(), xmax=h_numu->GetXaxis()->GetXmax();
  TH1D* frame=new TH1D(Form("frame_511_%s",mode),"",100,xmin,xmax);
  frame->GetXaxis()->SetTitle("#theta_{#nu} [rad]");
  frame->GetYaxis()->SetTitle("Flux / 6 #times 10^{20} POT / rad / cm^{2}");
  auto_logy_limits_range(frame,{h_numu,h_numubar,h_nue,h_nuebar},xmin,xmax);
  frame->Draw("AXIS");
  h_numu->Draw("HIST SAME"); h_nue->Draw("HIST SAME"); h_numubar->Draw("HIST SAME"); h_nuebar->Draw("HIST SAME");

  // Legend with angle-integrated fractions (area uses width=true)
  p_leg->cd();
  double s_numu   = integral_in(xmin,xmax,h_numu,true);
  double s_numubar= integral_in(xmin,xmax,h_numubar,true);
  double s_nue    = integral_in(xmin,xmax,h_nue,true);
  double s_nuebar = integral_in(xmin,xmax,h_nuebar,true);
  auto L = build_legend_like_stacked(
    p_leg,
    {{h_numu,"#nu_{#mu}"},{h_numubar,"#bar{#nu}_{#mu}"},{h_nue,"#nu_{e}"},{h_nuebar,"#bar{#nu}_{e}"}},
    {s_numu,s_numubar,s_nue,s_nuebar}, split, /*show_pct=*/true);
  L->Draw();

  c->Update();
  c->Print(Form("uboone_fig5_11_%s_angle.pdf",mode));
  delete frame; delete L; delete c;
}

// ------------------------------ FIGURE 5.12 ---------------------------------
// angle vs parent decay z map (find if present), axes titled, logz
static void fig_5_12_onefile(const char* mode, TFile& f){
  TH2D* H = find_theta_vs_z(f,"numu");
  if(!H){ printf("[5.12/%s] no TH2D(theta,z) found under numu/\n",mode); return; }
  TCanvas c(Form("c_512_%s",mode),Form("#theta vs parent z — %s",mode),900,700);
  c.SetRightMargin(0.14); c.SetLogz();
  H->SetContour(255);
  H->GetXaxis()->SetTitle("#theta_{#nu} [rad]");
  H->GetYaxis()->SetTitle("Parent decay z [m]");
  H->GetZaxis()->SetTitle("Density [arb.]");
  H->Draw("COLZ");
  c.Print(Form("uboone_fig5_12_%s_theta_vs_z.pdf",mode));
  delete H;
}

// --------------------------- FIGURES 5.13 & 5.14 ----------------------------
// Parent energy (E) and angle (θ) with π/K/μ/KL, using your dirs, your legend style

static TH1D* sum2_detach(TH1D* a, TH1D* b){
  if(!a||!b) return nullptr;
  TH1D* c=(TH1D*)a->Clone(); c->SetDirectory(0); c->Add(b); return c;
}

static void fig_5_13_energy_one(const char* mode, TFile& f, const std::string& flav="numu"){
  auto* h_piP   = (TH1D*)f.Get((flav+"/PI_Plus/Enu_"+flav+"_PI_Plus_AV_TPC").c_str());
  auto* h_piM   = (TH1D*)f.Get((flav+"/PI_Minus/Enu_"+flav+"_PI_Minus_AV_TPC").c_str());
  auto* h_kP    = (TH1D*)f.Get((flav+"/Kaon_Plus/Enu_"+flav+"_Kaon_Plus_AV_TPC").c_str());
  auto* h_kM    = (TH1D*)f.Get((flav+"/Kaon_Minus/Enu_"+flav+"_Kaon_Minus_AV_TPC").c_str());
  auto* h_muP   = (TH1D*)f.Get((flav+"/Mu_Plus/Enu_"+flav+"_Mu_Plus_AV_TPC").c_str());
  auto* h_muM   = (TH1D*)f.Get((flav+"/Mu_Minus/Enu_"+flav+"_Mu_Minus_AV_TPC").c_str());
  auto* h_KL    = (TH1D*)f.Get((flav+"/K0L/Enu_"+flav+"_K0L_AV_TPC").c_str());
  if(!h_piP||!h_piM||!h_kP||!h_kM||!h_muP||!h_muM||!h_KL){ printf("[5.13/%s] missing parent energy for %s\n",mode,flav.c_str()); return; }

  TH1D* Hpi = sum2_detach(h_piP,h_piM);
  TH1D* HK  = sum2_detach(h_kP,h_kM);
  TH1D* Hmu = sum2_detach(h_muP,h_muM);
  TH1D* HKL = (TH1D*)h_KL->Clone(); HKL->SetDirectory(0);

  style_line(Hpi, TColor::GetColor("#e41a1c"),1);
  style_line(HK , TColor::GetColor("#377eb8"),1);
  style_line(Hmu, TColor::GetColor("#984ea3"),2);
  style_line(HKL, TColor::GetColor("#4daf4a"),3);

  const double split=0.85;
  TCanvas* c=nullptr; TPad* p_main=nullptr; TPad* p_leg=nullptr;
  make_split_canvas(Form("c_513_%s_%s",mode,flav.c_str()),
                    Form("Parent energy — %s %s",mode,flav.c_str()),
                    split,c,p_main,p_leg,/*logy=*/true);

  p_main->cd();
  double xmin=Hpi->GetXaxis()->GetXmin(), xmax=Hpi->GetXaxis()->GetXmax();
  TH1D* frame=new TH1D(Form("frame_513_%s_%s",mode,flav.c_str()),"",100,xmin,xmax);
  int binwMeV=(int)std::lround(Hpi->GetXaxis()->GetBinWidth(1)*1000.0);
  frame->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  frame->GetYaxis()->SetTitle(Form("Flux / 6 #times 10^{20} POT / %d MeV / cm^{2}",binwMeV));
  auto_logy_limits_range(frame,{Hpi,HK,Hmu,HKL},xmin,xmax);
  frame->Draw("AXIS"); Hpi->Draw("HIST SAME"); HK->Draw("HIST SAME"); Hmu->Draw("HIST SAME"); HKL->Draw("HIST SAME");

  // Percentages with Eν > 60 MeV (width=true for area)
  p_leg->cd();
  const double thr = CFG::E_FRAC_MIN;
  double Spi  = integral_in(thr,xmax,Hpi,true);
  double SK   = integral_in(thr,xmax,HK ,true);
  double Smu  = integral_in(thr,xmax,Hmu,true);
  double SKL  = integral_in(thr,xmax,HKL,true);
  auto L = build_legend_like_stacked(p_leg,
    {{Hpi,"#pi"}, {HK,"K"}, {Hmu,"#mu"}, {HKL,"K^{0}_{L}"}},
    {Spi,SK,Smu,SKL}, split, /*show_pct=*/true);
  L->Draw();

  c->Update(); c->Print(Form("uboone_fig5_13_%s_%s_parentE.pdf",mode,flav.c_str()));
  delete frame; delete L; delete c;
}

static void fig_5_14_angle_one(const char* mode, TFile& f, const std::string& flav="numu"){
  auto* h_piP   = (TH1D*)f.Get((flav+"/PI_Plus/Th_"+flav+"_PI_Plus_AV_TPC").c_str());
  auto* h_piM   = (TH1D*)f.Get((flav+"/PI_Minus/Th_"+flav+"_PI_Minus_AV_TPC").c_str());
  auto* h_kP    = (TH1D*)f.Get((flav+"/Kaon_Plus/Th_"+flav+"_Kaon_Plus_AV_TPC").c_str());
  auto* h_kM    = (TH1D*)f.Get((flav+"/Kaon_Minus/Th_"+flav+"_Kaon_Minus_AV_TPC").c_str());
  auto* h_muP   = (TH1D*)f.Get((flav+"/Mu_Plus/Th_"+flav+"_Mu_Plus_AV_TPC").c_str());
  auto* h_muM   = (TH1D*)f.Get((flav+"/Mu_Minus/Th_"+flav+"_Mu_Minus_AV_TPC").c_str());
  auto* h_KL    = (TH1D*)f.Get((flav+"/K0L/Th_"+flav+"_K0L_AV_TPC").c_str());
  if(!h_piP||!h_piM||!h_kP||!h_kM||!h_muP||!h_muM||!h_KL){ printf("[5.14/%s] missing parent angle for %s\n",mode,flav.c_str()); return; }

  TH1D* Hpi = sum2_detach(h_piP,h_piM);
  TH1D* HK  = sum2_detach(h_kP,h_kM);
  TH1D* Hmu = sum2_detach(h_muP,h_muM);
  TH1D* HKL = (TH1D*)h_KL->Clone(); HKL->SetDirectory(0);

  style_line(Hpi, TColor::GetColor("#e41a1c"),1);
  style_line(HK , TColor::GetColor("#377eb8"),1);
  style_line(Hmu, TColor::GetColor("#984ea3"),2);
  style_line(HKL, TColor::GetColor("#4daf4a"),3);

  const double split=0.85;
  TCanvas* c=nullptr; TPad* p_main=nullptr; TPad* p_leg=nullptr;
  make_split_canvas(Form("c_514_%s_%s",mode,flav.c_str()),
                    Form("Parent angle — %s %s",mode,flav.c_str()),
                    split,c,p_main,p_leg,/*logy=*/true);

  p_main->cd();
  double xmin=Hpi->GetXaxis()->GetXmin(), xmax=Hpi->GetXaxis()->GetXmax();
  TH1D* frame=new TH1D(Form("frame_514_%s_%s",mode,flav.c_str()),"",100,xmin,xmax);
  frame->GetXaxis()->SetTitle("#theta_{#nu} [rad]");
  frame->GetYaxis()->SetTitle("Flux / 6 #times 10^{20} POT / rad / cm^{2}");
  auto_logy_limits_range(frame,{Hpi,HK,Hmu,HKL},xmin,xmax);
  frame->Draw("AXIS"); Hpi->Draw("HIST SAME"); HK->Draw("HIST SAME"); Hmu->Draw("HIST SAME"); HKL->Draw("HIST SAME");

  // Angle‑integrated percentages (width=true integrates over radians)
  p_leg->cd();
  double Spi = integral_in(xmin,xmax,Hpi,true);
  double SK  = integral_in(xmin,xmax,HK ,true);
  double Smu = integral_in(xmin,xmax,Hmu,true);
  double SKL = integral_in(xmin,xmax,HKL,true);
  auto L = build_legend_like_stacked(p_leg,
    {{Hpi,"#pi"}, {HK,"K"}, {Hmu,"#mu"}, {HKL,"K^{0}_{L}"}},
    {Spi,SK,Smu,SKL}, split, /*show_pct=*/true);
  L->Draw();

  c->Update(); c->Print(Form("uboone_fig5_14_%s_%s_parentAngle.pdf",mode,flav.c_str()));
  delete frame; delete L; delete c;
}

// ------------------------------ FIGURE 5.15 ---------------------------------
// E–angle heatmaps (per flavor), axes titled/units, logz
static void fig_5_15_onefile(const char* mode, TFile& f){
  for(const auto& flav : CFG::FLAVS){
    auto* h = (TH2D*)f.Get((flav + "/Detsmear/" + flav + "_CV_AV_TPC_2D").c_str());
    if(!h){ printf("[5.15/%s] missing %s 2D E-θ\n",mode,flav.c_str()); continue; }
    TH2D* H=(TH2D*)h->Clone(Form("cl_515_%s_%s",mode,flav.c_str())); H->SetDirectory(0);
    TCanvas c(Form("c_515_%s_%s",mode,flav.c_str()),Form("E vs #theta — %s %s",flav.c_str(),mode),900,700);
    c.SetRightMargin(0.14); c.SetLogz();
    H->SetContour(255);
    H->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    H->GetYaxis()->SetTitle("#theta_{#nu} [rad]");
    H->GetZaxis()->SetTitle("Flux / 6 #times 10^{20} POT / (GeV#timesrad) / cm^{2}");
    H->Draw("COLZ");
    c.Print(Form("uboone_fig5_15_%s_%s_EvTheta.pdf",mode,flav.c_str()));
    delete H;
  }
}

// -------- PPFX: scan categories & build total hadron-production covariance ----

struct PPFXCat { std::string name; std::map<int,TH1D*> univ; };

static std::map<std::string,PPFXCat> scan_ppfx_categories(TFile& f, const std::string& flav){
  std::map<std::string,PPFXCat> out;
  auto* d = f.GetDirectory( (flav + "/Multisims").c_str() );
  if(!d) return out;

  TIter it(d->GetListOfKeys());
  while(auto* k=(TKey*)it()){
    if(std::string(k->GetClassName())!="TH1D") continue;
    std::string name = k->GetName();
    if(name.find("_2D")!=std::string::npos) continue;
    // Expect: flav_ppfx_<category>_Uni_<N>_AV_TPC
    if(name.find(flav+"_ppfx_")!=0) continue;
    size_t posUni = name.find("_Uni_");
    if(posUni==std::string::npos) continue;
    std::string cat = name.substr( (flav+"_ppfx_").size(), posUni-(flav+"_ppfx_").size() );
    size_t posIdx = posUni + 5;
    size_t posEnd = name.find("_", posIdx);
    int idx = std::stoi(name.substr(posIdx, posEnd-posIdx));

    auto* h=(TH1D*)d->Get(name.c_str());
    if(!h) continue;
    auto* c=(TH1D*)h->Clone(Form("cl_%s",name.c_str())); c->SetDirectory(0);
    if(!out.count(cat)) out[cat] = PPFXCat{cat,{}};
    out[cat].univ[idx] = c;
  }
  for(auto& kv : out) printf("[PPFX] %s: %zu universes in category '%s'\n",
                             flav.c_str(), kv.second.univ.size(), kv.first.c_str());
  return out;
}

static TMatrixD covariance_about_cv(const std::vector<TH1D*>& U, const TH1D* CV){
  const int nb = CV->GetNbinsX();
  TMatrixD C(nb,nb); C.Zero();
  if(U.empty()) return C;
  const int N = (int)U.size();
  for(const auto* h : U){
    for(int i=1;i<=nb;++i){
      double di = h->GetBinContent(i) - CV->GetBinContent(i);
      for(int j=1;j<=nb;++j){
        double dj = h->GetBinContent(j) - CV->GetBinContent(j);
        C(i-1,j-1) += di*dj;
      }
    }
  }
  if(N>1) C *= (1.0/(N-1));
  return C;
}

// ------------------------------ FIGURE 5.16/18/19 ---------------------------
// CV ± 1σ PPFX band (from summed cov), fractional covariance, correlation

static TH1D* cv_energy_for_ppfx(TFile& f, const std::string& flav){
  auto* h = (TH1D*)f.Get((flav + "/Detsmear/" + flav + "_CV_AV_TPC").c_str());
  if(!h)   h = (TH1D*)f.Get((flav + "/Detsmear/" + flav + "_CV_AV_TPC_5MeV_bin").c_str());
  return h;
}

static void fig_5_16_18_19_oneflavor(const char* mode, TFile& f, const std::string& flav){
  TH1D* CV = cv_energy_for_ppfx(f,flav);
  if(!CV){ printf("[PPFX/%s] missing CV energy for %s\n",mode,flav.c_str()); return; }
  CV = (TH1D*)CV->Clone(Form("cv_%s_%s",mode,flav.c_str())); CV->SetDirectory(0);

  auto cats = scan_ppfx_categories(f,flav);
  if(cats.empty()){ printf("[PPFX/%s] no universes for %s\n",mode,flav.c_str()); delete CV; return; }

  const int nb = CV->GetNbinsX();
  TMatrixD Ctot(nb,nb); Ctot.Zero();

  // Sum category covariances about CV
  for(auto& kv : cats){
    std::vector<TH1D*> U; U.reserve(kv.second.univ.size());
    for(auto& ij : kv.second.univ) U.push_back(ij.second);
    TMatrixD Ccat = covariance_about_cv(U, CV);
    Ctot += Ccat;
  }

  // Build ±1σ from diag(Ctot)
  TH1D* hUp=(TH1D*)CV->Clone("ppfx_up");
  TH1D* hDn=(TH1D*)CV->Clone("ppfx_dn");
  TH1D* Sig=(TH1D*)CV->Clone("ppfx_sigma"); Sig->Reset(); Sig->SetDirectory(0);
  for(int i=1;i<=nb;++i){
    double s = std::sqrt(std::max(0.0, Ctot(i-1,i-1)));
    hUp->SetBinContent(i, CV->GetBinContent(i)+s);
    hDn->SetBinContent(i, CV->GetBinContent(i)-s);
    Sig->SetBinContent(i, s);
  }

  // Draw with your legend style (no percentages)
  const double split=0.85;
  TCanvas* c=nullptr; TPad* p_main=nullptr; TPad* p_leg=nullptr;
  make_split_canvas(Form("c_516_%s_%s",mode,flav.c_str()),
                    Form("PPFX CV #pm 1#sigma — %s %s",mode,flav.c_str()),
                    split,c,p_main,p_leg,/*logy=*/true);

  p_main->cd();
  TH1D* frame=(TH1D*)CV->Clone("frame_516"); frame->Reset(); frame->SetDirectory(0);
  int binwMeV=(int)std::lround(CV->GetXaxis()->GetBinWidth(1)*1000.0);
  frame->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  frame->GetYaxis()->SetTitle(Form("Flux / 6 #times 10^{20} POT / %d MeV / cm^{2}",binwMeV));
  double mn=1e99,mx=0; for(int b=1;b<=nb;++b){
    mn=std::min(mn,std::max(1e-30,std::min({CV->GetBinContent(b),hUp->GetBinContent(b),hDn->GetBinContent(b)})));
    mx=std::max(mx,std::max({CV->GetBinContent(b),hUp->GetBinContent(b),hDn->GetBinContent(b)}));
  }
  if(!(mn<1e98)) mn=1e-18; if(mx<=0) mx=1;
  frame->SetMinimum(mn*0.6); frame->SetMaximum(mx*4.0);
  frame->Draw("AXIS");

  style_line(CV, TColor::GetColor("#222222"),1);
  style_line(hUp,TColor::GetColor("#1f78b4"),2);
  style_line(hDn,TColor::GetColor("#1f78b4"),2);
  hUp->SetFillColorAlpha(TColor::GetColor("#1f78b4"),0.20);
  hDn->SetFillColorAlpha(TColor::GetColor("#1f78b4"),0.20);

  hUp->Draw("E3 SAME"); hDn->Draw("E3 SAME"); CV->Draw("HIST SAME");

  p_leg->cd();
  auto L = build_legend_like_stacked(
    p_leg,
    {{CV,"CV"}, {hUp,"CV + 1#sigma (PPFX)"}, {hDn,"CV - 1#sigma (PPFX)"}},
    /*sums=*/{}, split, /*show_pct=*/false);
  L->Draw();

  c->Update(); c->Print(Form("uboone_fig5_16_%s_%s_ppfxBand.pdf",mode,flav.c_str()));

  // Fractional covariance (dimensionless) & correlation
  TH2D Hfrac("ppfx_frac","PPFX fractional covariance;bin i;bin j", nb,0.5,nb+0.5, nb,0.5,nb+0.5);
  TH2D Hcorr("ppfx_corr","PPFX correlation;bin i;bin j", nb,0.5,nb+0.5, nb,0.5,nb+0.5);
  for(int i=1;i<=nb;++i){
    for(int j=1;j<=nb;++j){
      double mi = CV->GetBinContent(i), mj = CV->GetBinContent(j);
      double cij = Ctot(i-1,j-1);
      double frac = (mi!=0.0 && mj!=0.0) ? (cij/(mi*mj)) : 0.0;
      double sii = Ctot(i-1,i-1), sjj = Ctot(j-1,j-1);
      double rho = (sii>0 && sjj>0) ? (cij/std::sqrt(sii*sjj)) : 0.0;
      Hfrac.SetBinContent(i,j, frac);
      Hcorr.SetBinContent(i,j, rho);
    }
  }
  Hfrac.SetContour(255); Hcorr.SetContour(255);

  TCanvas cf(Form("c_518_%s_%s",mode,flav.c_str()),"frac cov",900,800);
  cf.SetRightMargin(0.14); cf.SetLogz(); Hfrac.GetZaxis()->SetTitle("Fractional covariance"); Hfrac.Draw("COLZ");
  cf.Print(Form("uboone_fig5_18_%s_%s_fracCov.pdf",mode,flav.c_str()));

  TCanvas cr(Form("c_519_%s_%s",mode,flav.c_str()),"corr",900,800);
  cr.SetRightMargin(0.14); Hcorr.GetZaxis()->SetTitle("Correlation coefficient"); Hcorr.Draw("COLZ");
  cr.Print(Form("uboone_fig5_19_%s_%s_corr.pdf",mode,flav.c_str()));

  // cleanup
  delete L; delete frame; delete c; delete hUp; delete hDn; delete Sig; delete CV;
  for(auto& kv : cats) for(auto& ij : kv.second.univ) delete ij.second;
}

// ------------------------------ DRIVER --------------------------------------

void flux_figs_minimal_uboone(){
  set_global_style();

  TFile fFHC(CFG::FILE_FHC,"READ");
  TFile fRHC(CFG::FILE_RHC,"READ");
  if(fFHC.IsZombie()){ printf("Cannot open FHC file.\n"); return; }
  if(fRHC.IsZombie()){ printf("Cannot open RHC file.\n"); return; }

  // 5.11 — Flux vs angle (legend in your style)
  fig_5_11_onefile("FHC", fFHC);
  fig_5_11_onefile("RHC", fRHC);

  // 5.12 — θ vs parent-z (legendless, axes titled)
  fig_5_12_onefile("FHC", fFHC);
  fig_5_12_onefile("RHC", fRHC);

  // 5.13 & 5.14 — parent decompositions (π/K/μ/KL) with your legend style
  fig_5_13_energy_one("FHC", fFHC, "numu");
  fig_5_13_energy_one("RHC", fRHC, "numu");
  fig_5_14_angle_one ("FHC", fFHC, "numu");
  fig_5_14_angle_one ("RHC", fRHC, "numu");

  // 5.15 — E–θ heatmaps (all flavors)
  fig_5_15_onefile("FHC", fFHC);
  fig_5_15_onefile("RHC", fRHC);

  // 5.16, 5.18, 5.19 — PPFX CV ± 1σ band + frac-cov + corr (each flavor)
  for(const auto& flav : CFG::FLAVS){
    fig_5_16_18_19_oneflavor("FHC", fFHC, flav);
    fig_5_16_18_19_oneflavor("RHC", fRHC, flav);
  }

  fFHC.Close(); fRHC.Close();
}

