// root -l -b -q 'scripts/plot_flux_minimal.C++'

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <limits>

void plot_flux_minimal() {
  // ---------------- hardcoded choices ----------------
  const char* FHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx_g4_10_4.root";
  const char* RHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc_ppfx_g4_10_4.root";

  const double Emin = 0.0;
  const double Emax = 10.0;

  // We want 25 MeV/bin:
  const double dE = 0.025;                  // GeV  (25 MeV)

  // fallback y-range (will be overridden by auto scaling below)
  const double ymin_fallback = 1e-8;
  const double ymax_fallback = 1e+6;

  const char* OUT_FHC = "uboone_flux_FHC.pdf";
  const char* OUT_RHC = "uboone_flux_RHC.pdf";

  // ---- Normalization ----
  // We want 6×10^20 POT units:
  constexpr bool   NORM_PER_POT = false;    // false => scale to NOMINAL_POT
  constexpr double NOMINAL_POT  = 6e20;
  // ---------------------------------------------------

  gStyle->SetOptStat(0);
  TH1::AddDirectory(kTRUE);

  const int nbins = std::max(1, int((Emax - Emin)/dE + 0.5));

  auto fresh = [&](const char* name) {
    if (auto* old = dynamic_cast<TH1D*>(gROOT->FindObject(name))) delete old;
    TH1D* h = new TH1D(name, "", nbins, Emin, Emax);
    h->Sumw2();
    return h;
  };

  auto sumPOT = [](TChain& ch) {
    double tot = 0.0;
    if (auto* files = ch.GetListOfFiles()) {
      for (int i = 0, n = files->GetEntriesFast(); i < n; ++i) {
        auto* el = (TChainElement*)files->UncheckedAt(i);
        if (!el) continue;
        TFile f(el->GetTitle(),"READ");
        if (!f.IsOpen()) continue;
        if (auto* hp = dynamic_cast<TH1*>(f.Get("POT"))) tot += hp->Integral();
      }
    }
    return tot;
  };

  auto style = [](TH1* h, int col, int ls){
    h->SetLineColor(col);
    h->SetLineStyle(ls);
    h->SetLineWidth(2);
    h->SetMarkerSize(0);
  };

  auto set_logy_limits = [&](TH1* frame, std::initializer_list<TH1*> hs){
    double minpos = std::numeric_limits<double>::infinity();
    double maxval = 0.0;
    for (TH1* h : hs) {
      for (int b=1; b<=h->GetNbinsX(); ++b) {
        double y = h->GetBinContent(b);
        if (y>0.0 && y<minpos) minpos = y;
        if (y>maxval)          maxval = y;
      }
    }
    if (!std::isfinite(minpos)) minpos = ymin_fallback;
    if (maxval<=0.0)            maxval = ymax_fallback;
    frame->SetMinimum(std::max(1e-30, minpos*0.8));   // must be >0 for log
    frame->SetMaximum(maxval*5.0);
  };

  auto make_and_draw = [&](const char* file, const char* tag, const char* outname){
    // ----- build spectra -----
    TChain ch("outTree"); ch.Add(file);

    const bool has_wgt  = ch.GetListOfBranches()->FindObject("wgt");
    const bool has_ppfx = ch.GetListOfBranches()->FindObject("wgt_ppfx");

    // CV weight: geometry/acceptance * PPFX CV (if present)
    TString w = has_wgt ? "wgt" : "1";
    if (has_ppfx) w += "*wgt_ppfx";

    // POT scaling
    const double pot_total = sumPOT(ch);
    if (pot_total <= 0) {
      std::cerr << "[plot_flux_minimal] WARNING: total POT <= 0; proceeding without POT scaling.\n";
    }
    char w_scaled[256];
    if (pot_total > 0) {
      if (NORM_PER_POT) {
        std::snprintf(w_scaled, sizeof(w_scaled), "(%s)/(%g)", w.Data(), pot_total);
      } else {
        std::snprintf(w_scaled, sizeof(w_scaled), "(%s)*(%g)", w.Data(), NOMINAL_POT / pot_total);
      }
    } else {
      std::snprintf(w_scaled, sizeof(w_scaled), "%s", w.Data());
    }

    TH1D* h_numu  = fresh(Form("h_numu_%s",  tag));
    TH1D* h_anumu = fresh(Form("h_anumu_%s", tag));
    TH1D* h_nue   = fresh(Form("h_nue_%s",   tag));
    TH1D* h_anue  = fresh(Form("h_anue_%s",  tag));

    auto fill = [&](TH1D* h, int pdg){
      ch.Draw(Form("nuE>>%s", h->GetName()),
              Form("(%g<=nuE && nuE<%g) * (ntype==%d) * (%s)", Emin, Emax, pdg, w_scaled),
              "goff");
    };
    fill(h_numu,  14);
    fill(h_anumu,-14);
    fill(h_nue,   12);
    fill(h_anue, -12);

    // ----- style + frame -----
    style(h_numu,  kRed+1, 1);
    style(h_nue,   kRed+1, 2);
    style(h_anumu, kBlue+2,1);
    style(h_anue,  kBlue+2,3);

    TCanvas c(Form("c_%s",tag), Form("%s Mode",tag), 900, 600);
    c.SetLogy();

    TH1D* frame = (TH1D*)h_numu->Clone(Form("frame_%s",tag));
    frame->Reset("ICES");
    frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");

    // Axis label that exactly matches the bin width and 6e20 scaling:
    const double dE_MeV = dE * 1000.0;
    frame->GetYaxis()->SetTitle(Form("#nu / 6#times10^{20} POT / %.0f MeV / cm^{2}", dE_MeV));

    set_logy_limits(frame, {h_numu, h_anumu, h_nue, h_anue});
    frame->Draw("AXIS");

    h_numu ->Draw("HIST SAME");
    h_nue  ->Draw("HIST SAME");
    h_anumu->Draw("HIST SAME");
    h_anue ->Draw("HIST SAME");

    auto integ = [&](const TH1D* h){ return h->Integral(1, nbins); };
    const double s_numu  = integ(h_numu);
    const double s_anumu = integ(h_anumu);
    const double s_nue   = integ(h_nue);
    const double s_anue  = integ(h_anue);
    const double s_tot   = std::max(1e-300, s_numu+s_anumu+s_nue+s_anue);
    auto pct = [&](double x){ return 100.0*x/s_tot; };

    TLegend L(0.62, 0.67, 0.935, 0.92);
    L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextSize(0.038);
    L.AddEntry(h_numu,  Form("#nu_{#mu} (%.1f%%)",       pct(s_numu)),  "l");
    L.AddEntry(h_nue,   Form("#nu_{e} (%.1f%%)",         pct(s_nue)),   "l");
    L.AddEntry(h_anumu, Form("#bar{#nu}_{#mu} (%.1f%%)", pct(s_anumu)), "l");
    L.AddEntry(h_anue,  Form("#bar{#nu}_{e} (%.1f%%)",   pct(s_anue)),  "l");
    L.Draw();

    TLatex t; t.SetNDC(); t.SetTextFont(42);
    t.SetTextSize(0.05); t.SetTextColor(kGray+2);
    t.DrawLatex(0.16, 0.90, Form("%s Mode", tag));
    t.SetTextSize(0.035); t.SetTextColor(kGray+1);
    t.DrawLatex(0.16, 0.84, Form("POT in inputs: %.3g", std::max(0.0, sumPOT(ch))));

    // small printout to confirm bin width actually used
    std::cout << "[plot_flux_minimal] " << tag
              << ": bin width = " << h_numu->GetBinWidth(1) << " GeV ("
              << h_numu->GetBinWidth(1)*1000.0 << " MeV), "
              << (NORM_PER_POT ? "per-POT" : "scaled to 6e20 POT") << "\n";

    c.SaveAs(outname);

    delete frame;
    delete h_numu; delete h_anumu; delete h_nue; delete h_anue;
  };

  make_and_draw(FHC_FILE, "FHC", OUT_FHC);
  make_and_draw(RHC_FILE, "RHC", OUT_RHC);

  std::cout << "Saved: " << OUT_FHC << " and " << OUT_RHC << std::endl;
}
