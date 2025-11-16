#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TString.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>

static void set_global_style() {
  const int font_style = 42;
  TStyle* style = new TStyle("PlotterStyle", "Plotter Style");
  style->SetTitleFont(font_style, "X");
  style->SetTitleFont(font_style, "Y");
  style->SetTitleFont(font_style, "Z");
  style->SetTitleSize(0.04, "X");
  style->SetTitleSize(0.04, "Y");
  style->SetTitleSize(0.05, "Z");
  style->SetLabelFont(font_style, "X");
  style->SetLabelFont(font_style, "Y");
  style->SetLabelFont(font_style, "Z");
  style->SetLabelSize(0.045, "X");
  style->SetLabelSize(0.045, "Y");
  style->SetLabelSize(0.045, "Z");
  style->SetLabelOffset(0.005, "X");
  style->SetLabelOffset(0.005, "Y");
  style->SetLabelOffset(0.005, "Z");
  style->SetTitleOffset(1.10, "X");
  style->SetTitleOffset(1.10, "Y");
  style->SetOptStat(0);
  style->SetOptTitle(0);
  style->SetPadTickX(1);
  style->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  style->SetPadLeftMargin(0.15);
  style->SetPadRightMargin(0.05);
  style->SetPadTopMargin(0.07);
  style->SetPadBottomMargin(0.12);
  style->SetMarkerSize(1.0);
  style->SetCanvasColor(0);
  style->SetPadColor(0);
  style->SetFrameFillColor(0);
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetStatColor(0);
  style->SetFrameBorderMode(0);
  style->SetTitleFillColor(0);
  style->SetTitleBorderSize(0);
  gROOT->SetStyle("PlotterStyle");
  gROOT->ForceStyle();
}

static TLegend* build_flux_legend_like_stacked(
    TPad* p_leg,
    TH1* h_numu, TH1* h_anumu, TH1* h_nue, TH1* h_anue,
    double split,
    double s_numu, double s_anumu, double s_nue, double s_anue, double s_tot)
{
  p_leg->cd();
  TLegend* L = new TLegend(0.12, 0.00, 0.95, 0.75);
  L->SetBorderSize(0);
  L->SetFillStyle(0);
  L->SetTextFont(42);
  const int n_entries = 4;
  int n_cols = (n_entries > 4) ? 3 : 2;
  L->SetNColumns(n_cols);
  L->SetColumnSeparation(0.08);
  L->SetEntrySeparation(0.00);
  L->SetMargin(0.25);
  const double s_main = 0.045;
  const double s_leg  = s_main * (split / (1.0 - split));
  L->SetTextSize(s_leg);
  auto pct = [&](double x){ return (s_tot>0 ? 100.0*x/s_tot : 0.0); };
  L->AddEntry(h_numu , Form("#nu_{#mu} (%.1f%%)",       pct(s_numu)),  "l");
  L->AddEntry(h_anumu, Form("#bar{#nu}_{#mu} (%.1f%%)", pct(s_anumu)), "l");
  L->AddEntry(h_nue  , Form("#nu_{e} (%.1f%%)",         pct(s_nue)),   "l");
  L->AddEntry(h_anue , Form("#bar{#nu}_{e} (%.1f%%)",   pct(s_anue)),  "l");
  return L;
}

static TH1D* fetch(TFile& f, const char* flav) {
  TH1D* h = (TH1D*)f.Get(Form("%s/Detsmear/%s_CV_AV_TPC_5MeV_bin", flav, flav));
  if (!h) h = (TH1D*)f.Get(Form("%s/Detsmear/%s_CV_AV_TPC", flav, flav));
  if (!h) return nullptr;
  TH1D* c = (TH1D*)h->Clone(Form("h_%s",flav));
  c->SetDirectory(0);
  return c;
}

static void style_line(TH1* h, int col, int ls){
  h->SetLineColor(col);
  h->SetLineStyle(ls);
  h->SetLineWidth(1);
  h->SetMarkerSize(0);
}

static void auto_logy_limits_range(TH1* frame, std::initializer_list<TH1*> hs,
                                   double xmin, double xmax){
  double minpos = std::numeric_limits<double>::infinity();
  double maxval = 0.0;
  for (TH1* h : hs) {
    int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
    int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
    for (int b=bmin; b<=bmax; ++b) {
      const double y = h->GetBinContent(b);
      if (y>0.0 && y<minpos) minpos = y;
      if (y>maxval)          maxval = y;
    }
  }
  if (!std::isfinite(minpos)) minpos = 1e-18;
  if (maxval<=0.0)            maxval = 1.0;
  frame->SetMinimum(std::max(1e-30, minpos*0.8));
  frame->SetMaximum(maxval*6.0);
}

static double integral_in(double xmin, double xmax, const TH1* h){
  int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
  int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
  return h->Integral(bmin, bmax);
}

static void draw_one(const char* file, const char* tag, const char* out) {
  TFile f(file,"READ");
  TH1D* h_numu   = fetch(f,"numu");
  TH1D* h_anumu  = fetch(f,"numubar");
  TH1D* h_nue    = fetch(f,"nue");
  TH1D* h_anue   = fetch(f,"nuebar");
  if (!h_numu || !h_anumu || !h_nue || !h_anue) { std::cerr << "missing spectra in " << file << "\n"; return; }

  style_line(h_numu,  kRed+1, 1);
  style_line(h_nue,   kRed+1, 2);
  style_line(h_anumu, kBlue+2,1);
  style_line(h_anue,  kBlue+2,3);

  const double Emin = 0.0, Emax = 10.0;
  h_numu ->GetXaxis()->SetRangeUser(Emin,Emax);
  h_anumu->GetXaxis()->SetRangeUser(Emin,Emax);
  h_nue  ->GetXaxis()->SetRangeUser(Emin,Emax);
  h_anue ->GetXaxis()->SetRangeUser(Emin,Emax);

  TCanvas c(Form("c_%s",tag), Form("%s Mode",tag), 900, 700);
  const double split = 0.85;
  TPad* p_main = new TPad("pad_main","pad_main", 0., 0.00, 1., split);
  TPad* p_leg  = new TPad("pad_legend","pad_legend", 0., split, 1., 1.00);
  p_main->SetTopMargin(0.01);
  p_main->SetBottomMargin(0.12);
  p_main->SetLeftMargin(0.12);
  p_main->SetRightMargin(0.05);
  p_main->SetLogy();
  p_leg->SetTopMargin(0.05);
  p_leg->SetBottomMargin(0.01);
  p_leg->SetLeftMargin(0.02);
  p_leg->SetRightMargin(0.02);
  p_main->Draw(); p_leg->Draw();
  p_main->cd();

  TH1D* frame = (TH1D*)h_numu->Clone(Form("frame_%s",tag));
  frame->Reset("ICES");
  const int binwMeV = int(std::lround(h_numu->GetXaxis()->GetBinWidth(1)*1000.0));
  frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  TString yTitle = Form("#nu / 6 #times 10^{20} POT / %d MeV / cm^{2}", binwMeV);
  frame->GetYaxis()->SetTitle(yTitle);
  auto_logy_limits_range(frame, {h_numu,h_anumu,h_nue,h_anue}, Emin, Emax);
  frame->Draw("AXIS");

  h_numu ->Draw("HIST SAME");
  h_nue  ->Draw("HIST SAME");
  h_anumu->Draw("HIST SAME");
  h_anue ->Draw("HIST SAME");

  p_leg->cd();
  const double s_numu  = integral_in(Emin,Emax,h_numu);
  const double s_anumu = integral_in(Emin,Emax,h_anumu);
  const double s_nue   = integral_in(Emin,Emax,h_nue);
  const double s_anue  = integral_in(Emin,Emax,h_anue);
  const double s_tot   = std::max(1e-300, s_numu+s_anumu+s_nue+s_anue);
  TLegend* L = build_flux_legend_like_stacked(p_leg, h_numu, h_anumu, h_nue, h_anue, split, s_numu, s_anumu, s_nue, s_anue, s_tot);
  L->Draw();

  c.cd(); c.Update(); c.SaveAs(out);

  delete frame; delete L; delete p_main; delete p_leg;
  delete h_numu; delete h_anumu; delete h_nue; delete h_anue;

  std::cout << "[plot_flux_minimal] " << tag
            << " | bin width = " << (binwMeV/1000.0) << " GeV (" << binwMeV << " MeV)"
            << " | y-axis: " << yTitle << "\n";
}

void plot_flux_minimal() {
  set_global_style();
  draw_one("/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root","FHC","uboone_flux_FHC.pdf");
  draw_one("/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root","RHC","uboone_flux_RHC.pdf");
  std::cout << "Saved: uboone_flux_FHC.pdf and uboone_flux_RHC.pdf\n";
}
