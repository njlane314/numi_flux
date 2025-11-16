#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TString.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdio>

static TH1D* fetch(TFile& f, const char* flav) {
  TString p = Form("%s/Detsmear/%s_CV_AV_TPC_5MeV_bin", flav, flav);
  TH1D* h = (TH1D*)f.Get(p);
  if (!h) { p = Form("%s/Detsmear/%s_CV_AV_TPC", flav, flav); h = (TH1D*)f.Get(p); }
  if (!h) return nullptr;
  h = (TH1D*)h->Clone(Form("h_%s", flav));
  h->SetDirectory(0);
  return h;
}

static void style(TH1D* h, int c, int ls) {
  h->SetLineColor(c);
  h->SetLineStyle(ls);
  h->SetLineWidth(2);
  h->SetMarkerSize(0);
}

static void draw_one(const char* file, const char* out) {
  TFile f(file, "READ");
  TH1D* h_numu   = fetch(f, "numu");
  TH1D* h_anumu  = fetch(f, "numubar");
  TH1D* h_nue    = fetch(f, "nue");
  TH1D* h_anue   = fetch(f, "nuebar");
  if (!h_numu || !h_anumu || !h_nue || !h_anue) { fprintf(stderr, "missing histograms in %s\n", file); return; }

  style(h_numu, kRed+1, 1);
  style(h_nue, kRed+1, 2);
  style(h_anumu, kBlue+2, 1);
  style(h_anue, kBlue+2, 3);

  double Emin = 0.0, Emax = 10.0;
  h_numu->GetXaxis()->SetRangeUser(Emin, Emax);
  h_nue->GetXaxis()->SetRangeUser(Emin, Emax);
  h_anumu->GetXaxis()->SetRangeUser(Emin, Emax);
  h_anue->GetXaxis()->SetRangeUser(Emin, Emax);

  TH1D* frame = (TH1D*)h_numu->Clone("frame"); frame->Reset("ICES");
  int binwMeV = int(std::lround(h_numu->GetXaxis()->GetBinWidth(1) * 1000.0));
  frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  frame->GetYaxis()->SetTitle(Form("#nu / 6 #times 10^{20} POT / %d MeV / cm^{2}", binwMeV));

  double minpos = std::numeric_limits<double>::infinity(), maxval = 0.0;
  auto upd = [&](TH1D* h){ int bmin = std::max(1, h->GetXaxis()->FindFixBin(Emin + 1e-9));
                           int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(Emax - 1e-9));
                           for (int b=bmin; b<=bmax; ++b) { double y=h->GetBinContent(b);
                             if (y>0 && y<minpos) minpos=y; if (y>maxval) maxval=y; } };
  upd(h_numu); upd(h_nue); upd(h_anumu); upd(h_anue);
  if (!std::isfinite(minpos)) minpos = 1e-18;
  if (maxval <= 0.0) maxval = 1.0;

  TCanvas c("c","",900,700);
  c.SetLogy();
  frame->SetMinimum(std::max(1e-30, minpos*0.8));
  frame->SetMaximum(maxval*6.0);
  frame->Draw("AXIS");
  h_numu->Draw("HIST SAME");
  h_nue->Draw("HIST SAME");
  h_anumu->Draw("HIST SAME");
  h_anue->Draw("HIST SAME");

  auto integ = [&](TH1D* h){ int bmin = std::max(1, h->GetXaxis()->FindFixBin(Emin + 1e-9));
                             int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(Emax - 1e-9));
                             return h->Integral(bmin, bmax); };
  double s1 = integ(h_numu), s2 = integ(h_anumu), s3 = integ(h_nue), s4 = integ(h_anue);
  double st = std::max(1e-300, s1+s2+s3+s4);

  TLegend L(0.55,0.72,0.95,0.93);
  L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42); L.SetNColumns(2);
  L.AddEntry(h_numu , Form("#nu_{#mu} (%.1f%%)",       100.0*s1/st), "l");
  L.AddEntry(h_anumu, Form("#bar{#nu}_{#mu} (%.1f%%)", 100.0*s2/st), "l");
  L.AddEntry(h_nue  , Form("#nu_{e} (%.1f%%)",         100.0*s3/st), "l");
  L.AddEntry(h_anue , Form("#bar{#nu}_{e} (%.1f%%)",   100.0*s4/st), "l");
  L.Draw();

  c.Update();
  c.Print(out);

  delete frame; delete h_numu; delete h_anumu; delete h_nue; delete h_anue;
}

void plot_flux_minimal() {
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.12);
  TGaxis::SetMaxDigits(4);
  draw_one("/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root","uboone_flux_FHC.pdf");
  draw_one("/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root","uboone_flux_RHC.pdf");
}
