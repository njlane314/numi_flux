#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TString.h"
#include "TColor.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>

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

static TLegend* build_flux_legend_like_stacked(TPad* p_leg, TH1* h_numu, TH1* h_anumu, TH1* h_nue, TH1* h_anue, double split, double s_numu, double s_anumu, double s_nue, double s_anue, double s_tot){
  p_leg->cd();
  TLegend* L=new TLegend(0.12,0.00,0.95,0.75);
  L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42);
  const int n_entries=4; int n_cols=(n_entries>4)?3:2;
  L->SetNColumns(n_cols); L->SetColumnSeparation(0.08); L->SetEntrySeparation(0.00); L->SetMargin(0.25);
  const double s_main=0.045; const double s_leg=s_main*(split/(1.0-split)); L->SetTextSize(s_leg);
  auto pct=[&](double x){ return (s_tot>0?100.0*x/s_tot:0.0); };
  L->AddEntry(h_numu ,Form("#nu_{#mu} (%.1f%%)",      pct(s_numu)) ,"l");
  L->AddEntry(h_anumu,Form("#bar{#nu}_{#mu} (%.1f%%)",pct(s_anumu)),"l");
  L->AddEntry(h_nue  ,Form("#nu_{e} (%.1f%%)",        pct(s_nue))  ,"l");
  L->AddEntry(h_anue ,Form("#bar{#nu}_{e} (%.1f%%)",  pct(s_anue)) ,"l");
  return L;
}

static void style_line(TH1* h,int col,int ls){ h->SetLineColor(col); h->SetLineStyle(ls); h->SetLineWidth(3); h->SetMarkerSize(0); }

static double integral_in(double xmin, double xmax, const TH1* h, bool width=false){
  int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
  int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
  return width ? h->Integral(bmin, bmax, "width") : h->Integral(bmin, bmax);
}

static double integral_and_error(double xmin, double xmax, const TH1* h,
                                 double& err, bool width=true){
  int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
  int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
  err = 0.0;
  return h->IntegralAndError(bmin, bmax, err, width ? "width" : "");
}

static void print_flux_window_integrals(const char* tag,
                                        const TH1* h_numu,  const TH1* h_anumu,
                                        const TH1* h_nue,   const TH1* h_anue,
                                        double xmin, double xmax){
  double e_numu, e_anumu, e_nue, e_anue;
  const double i_numu  = integral_and_error(xmin, xmax, h_numu , e_numu , /*width=*/true);
  const double i_anumu = integral_and_error(xmin, xmax, h_anumu, e_anumu, /*width=*/true);
  const double i_nue   = integral_and_error(xmin, xmax, h_nue  , e_nue  , /*width=*/true);
  const double i_anue  = integral_and_error(xmin, xmax, h_anue , e_anue , /*width=*/true);
  const double i_tot   = i_numu + i_anumu + i_nue + i_anue;

  auto pct = [&](double x){ return (i_tot > 0 ? 100.0*x/i_tot : 0.0); };
  const int binwMeV = (int)std::lround(h_numu->GetXaxis()->GetBinWidth(1) * 1000.0);

  printf("[plot_flux_minimal/%s] Energy-window integrals (%.2f–%.2f GeV)\n", tag, xmin, xmax);
  printf("  (TH1::Integral(...,\"width\") used: area = sum(content × binWidth).\n");
  printf("   Units ≈ # / 6×10^20 POT / cm^2; histogram y-axis is per %d MeV.)\n", binwMeV);

  printf("    nu_mu      : % .6e  ± %.2e   (%.1f%%)\n", i_numu , e_numu , pct(i_numu ));
  printf("    anti-nu_mu : % .6e  ± %.2e   (%.1f%%)\n", i_anumu, e_anumu, pct(i_anumu));
  printf("    nu_e       : % .6e  ± %.2e   (%.1f%%)\n", i_nue  , e_nue  , pct(i_nue  ));
  printf("    anti-nu_e  : % .6e  ± %.2e   (%.1f%%)\n", i_anue , e_anue , pct(i_anue ));
  printf("    TOTAL      : % .6e\n\n", i_tot);
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

static void draw_one(const char* file,const char* tag,const char* out){
  const double Emin=0.0, Emax=10.0, split=0.85;
  TFile f(file,"READ");
  TH1D* a=(TH1D*)f.Get("numu/Detsmear/numu_CV_AV_TPC_5MeV_bin"); if(!a) a=(TH1D*)f.Get("numu/Detsmear/numu_CV_AV_TPC");
  TH1D* b=(TH1D*)f.Get("numubar/Detsmear/numubar_CV_AV_TPC_5MeV_bin"); if(!b) b=(TH1D*)f.Get("numubar/Detsmear/numubar_CV_AV_TPC");
  TH1D* c=(TH1D*)f.Get("nue/Detsmear/nue_CV_AV_TPC_5MeV_bin"); if(!c) c=(TH1D*)f.Get("nue/Detsmear/nue_CV_AV_TPC");
  TH1D* d=(TH1D*)f.Get("nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin"); if(!d) d=(TH1D*)f.Get("nuebar/Detsmear/nuebar_CV_AV_TPC");
  if(!a||!b||!c||!d){ printf("[plot_flux_minimal/%s] missing *_CV_AV_TPC in Detsmear\n",tag); return; }
  a=(TH1D*)a->Clone("h_numu"); b=(TH1D*)b->Clone("h_anumu"); c=(TH1D*)c->Clone("h_nue"); d=(TH1D*)d->Clone("h_anue");
  a->SetDirectory(0); b->SetDirectory(0); c->SetDirectory(0); d->SetDirectory(0); f.Close();

  const double PrintMin = 0.25;
  print_flux_window_integrals(tag, a, b, c, d, PrintMin, Emax);

  int CR=TColor::GetColor("#e41a1c"), CB=TColor::GetColor("#1f78b4");
  style_line(a,CR,1); style_line(c,CR,2); style_line(b,CB,1); style_line(d,CB,3);

  TCanvas canv(Form("c_%s",tag),Form("%s Mode",tag),900,700);
  TPad* p_main=new TPad("pad_main","pad_main",0.,0.00,1.,split);
  TPad* p_leg =new TPad("pad_legend","pad_legend",0.,split,1.,1.00);
  p_main->SetTopMargin(0.01); p_main->SetBottomMargin(0.12); p_main->SetLeftMargin(0.12); p_main->SetRightMargin(0.05); p_main->SetLogy();
  p_leg->SetTopMargin(0.05); p_leg->SetBottomMargin(0.01); p_leg->SetLeftMargin(0.02); p_leg->SetRightMargin(0.02);
  p_main->Draw(); p_leg->Draw();
  p_main->cd();

  TH1D* frame=new TH1D(Form("frame_%s",tag),"",100,Emin,Emax);
  int binwMeV=(int)std::lround(a->GetXaxis()->GetBinWidth(1)*1000.0);
  frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  frame->GetYaxis()->SetTitle(Form("#nu / 6 #times 10^{20} POT / %d MeV / cm^{2}",binwMeV));
  a->GetXaxis()->SetRangeUser(Emin,Emax); b->GetXaxis()->SetRangeUser(Emin,Emax); c->GetXaxis()->SetRangeUser(Emin,Emax); d->GetXaxis()->SetRangeUser(Emin,Emax);
  auto_logy_limits_range(frame,{a,b,c,d},Emin,Emax);
  frame->Draw("AXIS");
  a->Draw("HIST SAME"); c->Draw("HIST SAME"); b->Draw("HIST SAME"); d->Draw("HIST SAME");

  p_leg->cd();
  double s_numu=integral_in(Emin,Emax,a), s_anumu=integral_in(Emin,Emax,b), s_nue=integral_in(Emin,Emax,c), s_anue=integral_in(Emin,Emax,d);
  double s_tot=std::max(1e-300,s_numu+s_anumu+s_nue+s_anue);
  TLegend* L=build_flux_legend_like_stacked(p_leg,a,b,c,d,split,s_numu,s_anumu,s_nue,s_anue,s_tot);
  L->Draw();

  canv.cd(); canv.Update(); canv.Print(out);
  printf("[plot_flux_minimal] %s | bin width ≈ %.3f GeV (%d MeV)\n", tag, a->GetXaxis()->GetBinWidth(1), binwMeV);
  delete frame; delete L; delete p_main; delete p_leg; delete a; delete b; delete c; delete d;
}

void plot_flux_minimal(){
  set_global_style();
  draw_one("/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root","FHC","uboone_flux_FHC.pdf");
  draw_one("/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root","RHC","uboone_flux_RHC.pdf");
}
