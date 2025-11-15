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

void plot_flux_minimal() {
  // ---------------- hardcoded choices ----------------
  const char* FHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx_g4_10_4.root";
  const char* RHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc_ppfx_g4_10_4.root";

  const double Emin = 0.0;
  const double Emax = 10.0;     // extend to 8 GeV if you like (bins stay 10 MeV)
  const double ymin = 1e-12;     // y-axis clip
  const double ymax = 1e-3;

  const char* OUT_FHC = "uboone_flux_FHC.pdf";
  const char* OUT_RHC = "uboone_flux_RHC.pdf";

  // Normalization toggle: per POT (true) or scaled to 6e20 POT (false)
  constexpr bool NORM_PER_POT = true;
  constexpr double NOMINAL_POT = 6e20;
  // ---------------------------------------------------

  gStyle->SetOptStat(0);
  TH1::AddDirectory(kTRUE);         // let TTree::Draw find our histograms

  const double dE = 0.01;           // 10 MeV/bin
  const int nbins = std::max(1, int((Emax - Emin)/dE + 0.5));

  // helper: delete stale hist with same name (avoids "Potential memory leak" warning)
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

  auto make_and_draw = [&](const char* file, const char* tag, const char* outname){
    // ----- build spectra -----
    TChain ch("outTree"); ch.Add(file);

    const bool has_wgt  = ch.GetListOfBranches()->FindObject("wgt");
    const bool has_ppfx = ch.GetListOfBranches()->FindObject("wgt_ppfx");

    // (1) Combine geometry/acceptance with PPFX CV
    std::string wexpr = has_wgt ? "wgt" : "1";
    if (has_ppfx) wexpr += "*wgt_ppfx";

    // (2) Scale by POT to match the axis label (per POT) or to 6e20 POT
    const double pot_total = sumPOT(ch);
    if (pot_total <= 0) {
      std::cerr << "[plot_flux_minimal] WARNING: total POT <= 0; proceeding without POT scaling.\n";
    }
    char w_scaled[256];
    if (NORM_PER_POT) {
      std::snprintf(w_scaled, sizeof(w_scaled), "(%s)/(%g)", wexpr.c_str(), std::max(1e-30, pot_total));
    } else {
      std::snprintf(w_scaled, sizeof(w_scaled), "(%s)*(%g)", wexpr.c_str(), NOMINAL_POT / std::max(1e-30, pot_total));
    }

    TH1D* h_numu  = fresh(Form("h_numu_%s",  tag));
    TH1D* h_anumu = fresh(Form("h_anumu_%s", tag));
    TH1D* h_nue   = fresh(Form("h_nue_%s",   tag));
    TH1D* h_anue  = fresh(Form("h_anue_%s",  tag));

    auto fill = [&](TH1D* h, int pdg){
      ch.Draw(Form("nuE>>%s", h->GetName()),
              Form("(%g<=nuE && nuE<%g) * (ntype==%d) * %s", Emin, Emax, pdg, w_scaled),
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
    frame->GetYaxis()->SetTitle(NORM_PER_POT ?
      "#nu / POT / 10 MeV / cm^{2}" :
      "#nu / 6e20 POT / 10 MeV / cm^{2}");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.35);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
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
    t.DrawLatex(0.16, 0.84, Form("POT in inputs: %.3g", pot_total));

    c.SaveAs(outname);

    // cleanup
    delete frame;
    delete h_numu; delete h_anumu; delete h_nue; delete h_anue;
  };

  // Make the two separate plots
  make_and_draw(FHC_FILE, "FHC", OUT_FHC);
  make_and_draw(RHC_FILE, "RHC", OUT_RHC);

  std::cout << "Saved: " << OUT_FHC << " and " << OUT_RHC << std::endl;
}
