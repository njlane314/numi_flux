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
#include "TGaxis.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>
#include <limits>

// -----------------------------
// Global style (verbatim settings from your Plotter::set_global_style)
// -----------------------------
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

void plot_flux_minimal() {
  set_global_style();

  // ---------------- hardcoded choices ----------------
  const char* FHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx_g4_10_4.root";
  const char* RHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc_ppfx_g4_10_4.root";

  const double Emin = 0.0;
  const double Emax = 10.0;

  // 25 MeV/bin to match requested units
  const double dE = 0.025;                  // GeV  (25 MeV)

  const char* OUT_FHC = "uboone_flux_FHC.pdf";
  const char* OUT_RHC = "uboone_flux_RHC.pdf";

  // ---- Normalization ----
  // We want 6×10^20 POT units:
  constexpr bool   NORM_PER_POT = false;    // false => scale to NOMINAL_POT
  constexpr double NOMINAL_POT  = 6e20;
  // ---------------------------------------------------

  TH1::AddDirectory(kTRUE);

  const int nbins = std::max(1, int((Emax - Emin)/dE + 0.5));
  const double binwGeV = (Emax - Emin) / nbins;
  const int    binwMeV = (int)std::lround(binwGeV * 1000.0);

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

  auto style_line = [](TH1* h, int col, int ls){
    h->SetLineColor(col);
    h->SetLineStyle(ls);
    h->SetLineWidth(2);
    h->SetMarkerSize(0);
  };

  auto auto_logy_limits = [&](TH1* frame, std::initializer_list<TH1*> hs){
    double minpos = std::numeric_limits<double>::infinity();
    double maxval = 0.0;
    for (TH1* h : hs) {
      for (int b=1; b<=h->GetNbinsX(); ++b) {
        const double y = h->GetBinContent(b);
        if (y>0.0 && y<minpos) minpos = y;
        if (y>maxval)          maxval = y;
      }
    }
    if (!std::isfinite(minpos)) minpos = 1e-18;
    if (maxval<=0.0)            maxval = 1.0;
    frame->SetMinimum(std::max(1e-30, minpos*0.8));   // must be >0 for log
    frame->SetMaximum(maxval*6.0);
  };

  auto make_and_draw = [&](const char* file, const char* tag, const char* outname){
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

    // Build spectra
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

    // Lines & colors (ν in red, ν̄ in blue; e-flavors dashed/dotted)
    style_line(h_numu,  kRed+1, 1);
    style_line(h_nue,   kRed+1, 2);
    style_line(h_anumu, kBlue+2,1);
    style_line(h_anue,  kBlue+2,3);

    // -----------------------------
    // Canvas with "legend on top" pad like rarexsec style
    // -----------------------------
    TCanvas c(Form("c_%s",tag), Form("%s Mode",tag), 900, 700);

    const double split = 0.83; // bottom main pad up to here, legend pad on top
    TPad* p_main = new TPad("pad_main","pad_main", 0., 0.00, 1., split);
    TPad* p_leg  = new TPad("pad_legend","pad_legend", 0., split, 1., 1.00);

    // match margins to your style
    p_main->SetTopMargin(0.02);
    p_main->SetBottomMargin(0.12);
    p_main->SetLeftMargin(0.15);
    p_main->SetRightMargin(0.05);
    p_main->SetLogy();

    p_leg->SetTopMargin(0.15);
    p_leg->SetBottomMargin(0.05);
    p_leg->SetLeftMargin(0.02);
    p_leg->SetRightMargin(0.02);

    p_main->Draw();
    p_leg->Draw();

    // ---- draw main pad ----
    p_main->cd();

    // Frame with exact axis titles
    TH1D* frame = (TH1D*)h_numu->Clone(Form("frame_%s",tag));
    frame->Reset("ICES");
    frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");

    TString potLabel = NORM_PER_POT ? "POT" : "6 #times 10^{20} POT";
    TString yTitle = Form("#nu / %s / %d MeV / cm^{2}", potLabel.Data(), binwMeV);
    frame->GetYaxis()->SetTitle(yTitle);

    auto_logy_limits(frame, {h_numu, h_anumu, h_nue, h_anue});
    frame->Draw("AXIS");

    h_numu ->Draw("HIST SAME");
    h_nue  ->Draw("HIST SAME");
    h_anumu->Draw("HIST SAME");
    h_anue ->Draw("HIST SAME");

    // ---- legend pad on top ----
    p_leg->cd();

    // percentages for labels
    auto integ = [&](const TH1D* h){ return h->Integral(1, nbins); };
    const double s_numu  = integ(h_numu);
    const double s_anumu = integ(h_anumu);
    const double s_nue   = integ(h_nue);
    const double s_anue  = integ(h_anue);
    const double s_tot   = std::max(1e-300, s_numu+s_anumu+s_nue+s_anue);
    auto pct = [&](double x){ return 100.0*x/s_tot; };

    // Optional header on the left
    TLatex hdr; hdr.SetNDC(); hdr.SetTextFont(62); hdr.SetTextSize(0.070);
    hdr.DrawLatex(0.12, 0.86, Form("%s Mode", tag));
    TLatex pot; pot.SetNDC(); pot.SetTextFont(42); pot.SetTextSize(0.045); pot.SetTextColor(kGray+2);
    pot.DrawLatex(0.12, 0.66, Form("POT in inputs: %.3g", pot_total));

    // Legend block (right/top), 2 columns, row-wise pairing
    // Narrower box + margin so line samples are a sensible length
    TLegend* L = new TLegend(0.48, 0.20, 0.97, 0.92);
    L->SetBorderSize(0);
    L->SetFillStyle(0);
    L->SetTextFont(42);
    L->SetTextSize(0.055);
    L->SetNColumns(2);
    L->SetColumnSeparation(0.07);
    L->SetMargin(0.18);  // space between sample and text

    // *** ORDER MATTERS ***
    // Add entries in [νμ, ν̄μ, νe, ν̄e] so each row pairs ν with ν̄
    L->AddEntry(h_numu , Form("#nu_{#mu} (%.1f%%)",       pct(s_numu)),  "l");
    L->AddEntry(h_anumu, Form("#bar{#nu}_{#mu} (%.1f%%)", pct(s_anumu)), "l");
    L->AddEntry(h_nue  , Form("#nu_{e} (%.1f%%)",         pct(s_nue)),   "l");
    L->AddEntry(h_anue , Form("#bar{#nu}_{e} (%.1f%%)",   pct(s_anue)),  "l");

    L->Draw();

    // ---- save ----
    c.cd();
    c.Update();
    c.SaveAs(outname);

    // cleanup
    delete frame;
    delete L;
    delete p_main;
    delete p_leg;
    delete h_numu; delete h_anumu; delete h_nue; delete h_anue;

    std::cout << "[plot_flux_minimal] " << tag
              << ": bin width = " << binwGeV << " GeV (" << binwMeV << " MeV), "
              << (NORM_PER_POT ? "per-POT" : "scaled to 6e20 POT")
              << " | y-axis: " << yTitle << "\n";
  };

  make_and_draw(FHC_FILE, "FHC", OUT_FHC);
  make_and_draw(RHC_FILE, "RHC", OUT_RHC);

  std::cout << "Saved: " << OUT_FHC << " and " << OUT_RHC << std::endl;
}
