// Compile & run example:
// root -l -b -q 'scripts/plot_flux_perpot.C++("/path/to/dk2nu_fhc.root", "/path/to/dk2nu_rhc.root", 0.0, 5.0, "flux_FHC_RHC.pdf")'

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TChainElement.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

namespace {
  struct FluxSet {
    TH1D* h_numu   = nullptr; //  νμ
    TH1D* h_anumu  = nullptr; //  ν̄μ
    TH1D* h_nue    = nullptr; //  νe
    TH1D* h_anue   = nullptr; //  ν̄e
    double totPOT  = 0.0;
    int nbins      = 0;
  };

  // Return sum of POT histograms in all files on the chain (if present)
  double SumPOT(TChain& ch) {
    double totPOT = 0.0;
    if (auto* files = ch.GetListOfFiles()) {
      const int nfiles = files->GetEntriesFast();
      for (int i = 0; i < nfiles; ++i) {
        auto* el = static_cast<TChainElement*>(files->UncheckedAt(i));
        if (!el) continue;
        TFile f(el->GetTitle(), "READ");
        if (!f.IsOpen()) continue;
        if (auto* hp = dynamic_cast<TH1*>(f.Get("POT"))) totPOT += hp->Integral();
      }
    }
    return totPOT;
  }

  // Fill all four flavor histograms for one running mode (FHC or RHC)
  FluxSet BuildFluxSet(const char* inpat, double Emin, double Emax) {
    FluxSet fs;
    TChain ch("outTree");
    ch.Add(inpat);

    // Choose weight: PPFX if available, else CV only
    const bool has_ppfx = ch.GetListOfBranches()->FindObject("wgt_ppfx");
    const char* w = has_ppfx ? "wgt_ppfx" : "wgt";

    // --- 10 MeV binning by construction ---
    const double dE_GeV = 0.01; // 10 MeV
    const int nbins = std::max(1, int(std::floor((Emax - Emin)/dE_GeV + 0.5)));
    fs.nbins = nbins;

    auto make = [&](const char* name){
      auto* h = new TH1D(name,"",nbins,Emin,Emax);
      h->SetDirectory(nullptr);
      h->Sumw2();
      return h;
    };

    fs.h_numu  = make("h_numu");
    fs.h_anumu = make("h_anumu");
    fs.h_nue   = make("h_nue");
    fs.h_anue  = make("h_anue");

    // Helper to fill one histogram with a PDG selection
    auto fill = [&](TH1D* h, int pdg){
      ch.Draw(Form("nuE>>%s", h->GetName()),
              Form("(%g<=nuE && nuE<%g) * (ntype==%d) * %s",Emin,Emax,pdg,w),
              "goff");
    };

    fill(fs.h_numu,   14);
    fill(fs.h_anumu, -14);
    fill(fs.h_nue,    12);
    fill(fs.h_anue,  -12);

    fs.totPOT = SumPOT(ch);
    return fs;
  }

  // Draw one panel (either FHC or RHC)
  void DrawPanel(TPad* pad, const FluxSet& fs, const char* modeLabel, bool drawXAxis) {
    pad->cd();
    pad->SetLogy();
    pad->SetGrid(0,0);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.05);
    pad->SetTopMargin(0.08);
    pad->SetBottomMargin(drawXAxis ? 0.16 : 0.04);

    // Style: red for ν (neutrinos), blue for ν̄ (antineutrinos)
    const int col_nu  = kRed+1;
    const int col_anu = kBlue+2;

    auto style = [](TH1* h, int col, int ls){
      h->SetLineColor(col);
      h->SetLineStyle(ls);
      h->SetLineWidth(2);
      h->SetMarkerSize(0);
    };
    // solid = μ, dashed/dotted = e
    style(fs.h_numu,  col_nu,  1);
    style(fs.h_anumu, col_anu, 1);
    style(fs.h_nue,   col_nu,  2);
    style(fs.h_anue,  col_anu, 3);

    // Axis cosmetics on a "frame" copy (use νμ to define axes)
    TH1D* frame = (TH1D*)fs.h_numu->Clone("frame");
    frame->Reset("ICES");
    frame->SetDirectory(nullptr);
    frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    frame->GetYaxis()->SetTitle("#nu / POT / 10 MeV / cm^{2}");
    frame->GetXaxis()->SetLabelSize(drawXAxis ? 0.045 : 0.0);
    frame->GetXaxis()->SetTitleSize(drawXAxis ? 0.05 : 0.0);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1.35);

    // Y range
    const double ymax = 1.5 * std::max({fs.h_numu->GetMaximum(), fs.h_anumu->GetMaximum(),
                                        fs.h_nue->GetMaximum(),  fs.h_anue->GetMaximum()});
    frame->SetMinimum(1e-15);
    frame->SetMaximum(std::max(1e-12, ymax));
    frame->Draw("AXIS");

    fs.h_numu->Draw("HIST SAME");
    fs.h_anumu->Draw("HIST SAME");
    fs.h_nue->Draw("HIST SAME");
    fs.h_anue->Draw("HIST SAME");

    // Fractions in the plotted energy window
    auto integ = [&](const TH1D* h){ return h->Integral(1, fs.nbins); }; // constant bin width ⇒ proportional
    const double s_numu  = integ(fs.h_numu);
    const double s_anumu = integ(fs.h_anumu);
    const double s_nue   = integ(fs.h_nue);
    const double s_anue  = integ(fs.h_anue);
    const double s_tot   = std::max(1e-300, s_numu + s_anumu + s_nue + s_anue);

    auto pct = [&](double x){ return 100.0 * x / s_tot; };

    // Legend
    TLegend* L = new TLegend(0.62, 0.67, 0.935, 0.92);
    L->SetBorderSize(0);
    L->SetFillStyle(0);
    L->SetTextSize(0.038);
    L->AddEntry(fs.h_numu,  Form("#nu_{#mu} (%.1f%%)", pct(s_numu)),  "l");
    L->AddEntry(fs.h_nue,   Form("#nu_{e} (%.1f%%)",  pct(s_nue)),   "l");
    L->AddEntry(fs.h_anumu, Form("#bar{#nu}_{#mu} (%.1f%%)", pct(s_anumu)), "l");
    L->AddEntry(fs.h_anue,  Form("#bar{#nu}_{e} (%.1f%%)",  pct(s_anue)),  "l");
    L->Draw();

    // Label like the original figure
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.05);
    t.SetTextFont(42);
    t.SetTextColor(kGray+2);
    t.DrawLatex(0.16, 0.90, modeLabel);

    // Optional: POT readback (comment out if not needed)
    if (fs.totPOT > 0) {
      t.SetTextSize(0.035);
      t.SetTextColor(kGray+1);
      t.DrawLatex(0.16, 0.84, Form("POT in inputs: %.3g", fs.totPOT));
    }
  }
} // namespace

void plot_flux_perpot(const char* fhc_inpat,
                      const char* rhc_inpat,
                      double Emin = 0.0, double Emax = 5.0,
                      const char* outname = "flux_FHC_RHC.pdf")
{
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  // Build flux histograms for each running mode
  FluxSet FHC = BuildFluxSet(fhc_inpat, Emin, Emax);
  FluxSet RHC = BuildFluxSet(rhc_inpat, Emin, Emax);

  // Canvas with two stacked pads
  TCanvas c("c","Flux per POT",900,1100);
  c.Divide(1,2,0,0);

  // Top (FHC) – hide x labels to mimic the reference figure
  DrawPanel((TPad*)c.cd(1), FHC, "FHC Mode", /*drawXAxis=*/false);
  // Bottom (RHC)
  DrawPanel((TPad*)c.cd(2), RHC, "RHC Mode", /*drawXAxis=*/true);

  c.SaveAs(outname);
  std::cout << "Saved: " << outname << std::endl;
}
