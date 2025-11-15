// root -l -b -q 'scripts/plot_flux_extras.C++'
//
// Produces add-on Central Flux Model plots in the same style as plot_flux_minimal.C++
//   (1) Wrong-sign fraction vs E (FHC & RHC)
//   (2) Flavor ratios vs E:  (nu_e/nu_mu)_FHC and (anue/anumu)_RHC
//   (3) Cumulative (integrated) flux vs E (CDF) for nu_mu, FHC & RHC
//   (4) RHC/FHC shape ratios for nu_mu and for anumu
//   (5) [Optional] Parentage (pi/K/mu) stacked vs E per mode — requires correct outTree.nuE
//
// NOTES
// - This script uses the pre-filled flux histograms in your files (already scaled to 6e20 POT per cm^2).
// - Parentage stacks use outTree branches (ptype, ntype, wgt, wgt_ppfx, nuE). Ensure nuE = MicroBooNE energy.
//   If your current tree has the SHS-energy bug, either patch the producer or comment out make_parentage_stacks().
//
// Author: you, adapted by ChatGPT
// ----------------------------------------------------------------------------

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TColor.h"
#include "TGaxis.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

// ---------------------------------------------------------------------------------
// Global style (copied from your plot_flux_minimal.C++ with identical settings)
// ---------------------------------------------------------------------------------
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

// Simple styling helpers
static void style_line(TH1* h, int col, int ls = 1, int lw = 2) {
  h->SetLineColor(col);
  h->SetLineStyle(ls);
  h->SetLineWidth(lw);
  h->SetMarkerSize(0);
}

static void style_fill(TH1* h, int col) {
  h->SetFillColor(col);
  h->SetLineColor(TColor::GetColor(0,0,0));
  h->SetLineWidth(1);
}

// Divide hnum by hden bin-by-bin handling zeros safely
static TH1D* safe_ratio(const TH1D* hnum, const TH1D* hden, const char* name) {
  TH1D* r = (TH1D*)hnum->Clone(name);
  r->Reset("ICES");
  r->Sumw2();
  const int nb = hnum->GetNbinsX();
  for (int b = 1; b <= nb; ++b) {
    const double a = hnum->GetBinContent(b);
    const double ea = hnum->GetBinError(b);
    const double d = hden->GetBinContent(b);
    const double ed = hden->GetBinError(b);
    double y=0, e=0;
    if (d > 0) {
      y = a/d;
      // error propagation: (a/d)
      const double rel2 = (ea>0? (ea*ea) : 0)/(d*d) + (a>0? (a*a*ed*ed)/(d*d*d*d) : 0);
      e = std::sqrt(rel2);
    } else {
      y = 0.0; e = 0.0;
    }
    r->SetBinContent(b, y);
    r->SetBinError(b, e);
  }
  return r;
}

// CDF builder: cumulative integral / total integral
static TH1D* make_cdf(const TH1D* h, const char* name) {
  TH1D* c = (TH1D*)h->Clone(name);
  c->Reset("ICES");
  double tot = h->Integral(1, h->GetNbinsX());
  if (tot <= 0) tot = 1.0;
  double cum = 0.0;
  for (int b=1; b<=h->GetNbinsX(); ++b) {
    cum += h->GetBinContent(b);
    c->SetBinContent(b, cum/tot);
  }
  return c;
}

// Fetch the four flavor flux histograms from a file (already scaled to 6e20 POT / cm^2 by your producer)
struct FluxH {
  TH1D* numu  = nullptr;
  TH1D* anumu = nullptr;
  TH1D* nue   = nullptr;
  TH1D* anue  = nullptr;
};

static FluxH load_flux_hists(const char* fpath, const char* tag) {
  FluxH H;
  TFile* F = TFile::Open(fpath, "READ");
  if (!F || F->IsZombie()) { std::cerr << "[load_flux_hists] Cannot open " << fpath << "\n"; return H; }
  H.numu  = dynamic_cast<TH1D*>(F->Get("numuFluxHisto"));
  H.anumu = dynamic_cast<TH1D*>(F->Get("anumuFluxHisto"));
  H.nue   = dynamic_cast<TH1D*>(F->Get("nueFluxHisto"));
  H.anue  = dynamic_cast<TH1D*>(F->Get("anueFluxHisto"));
  if (!H.numu || !H.anumu || !H.nue || !H.anue) {
    std::cerr << "[load_flux_hists] Missing flux histos in " << fpath << "\n";
  } else {
    // Make them unique in the global directory (clone+set directory to gROOT)
    H.numu  = (TH1D*)H.numu ->Clone(Form("numu_%s",  tag));  H.numu ->SetDirectory(gROOT);
    H.anumu = (TH1D*)H.anumu->Clone(Form("anumu_%s", tag));  H.anumu->SetDirectory(gROOT);
    H.nue   = (TH1D*)H.nue  ->Clone(Form("nue_%s",   tag));  H.nue  ->SetDirectory(gROOT);
    H.anue  = (TH1D*)H.anue ->Clone(Form("anue_%s",  tag));  H.anue ->SetDirectory(gROOT);
  }
  F->Close();
  delete F;
  return H;
}

// ------------------- Plotters -------------------

static void make_wrong_sign_fraction(const FluxH& FHC, const FluxH& RHC, const char* outpdf) {
  // FHC: WS = anumu / (numu + anumu)
  TH1D* denF = (TH1D*)FHC.numu->Clone("denF");  denF->Add(FHC.anumu);
  TH1D* wsF  = safe_ratio(FHC.anumu, denF, "wsF");
  // RHC: WS = numu / (anumu + numu)
  TH1D* denR = (TH1D*)RHC.anumu->Clone("denR"); denR->Add(RHC.numu);
  TH1D* wsR  = safe_ratio(RHC.numu, denR, "wsR");

  style_line(wsF, kRed+1, 1, 2);
  style_line(wsR, kBlue+2, 2, 2);

  wsF->SetTitle("");
  wsF->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  wsF->GetYaxis()->SetTitle("Wrong-sign fraction");
  wsF->SetMinimum(0.0);
  wsF->SetMaximum(1.0);

  TCanvas c("c_ws","Wrong-sign fraction",900,700);
  wsF->Draw("hist");
  wsR->Draw("hist same");

  TLegend L(0.55,0.18,0.88,0.34);
  L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42); L.SetTextSize(0.04);
  L.AddEntry(wsF, "FHC: #bar{#nu}_{#mu}/(#nu_{#mu}+#bar{#nu}_{#mu})", "l");
  L.AddEntry(wsR, "RHC: #nu_{#mu}/(#nu_{#mu}+#bar{#nu}_{#mu})", "l");
  L.Draw();

  c.SaveAs(outpdf);
}

static void make_flavor_ratios(const FluxH& FHC, const FluxH& RHC, const char* outpdf) {
  // FHC: nue/numu ;  RHC: anue/anumu
  TH1D* rF = safe_ratio(FHC.nue , FHC.numu , "rF_nue_over_numu");
  TH1D* rR = safe_ratio(RHC.anue, RHC.anumu, "rR_anue_over_anumu");

  style_line(rF, kRed+1, 1, 2);
  style_line(rR, kBlue+2, 2, 2);

  rF->SetTitle("");
  rF->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  rF->GetYaxis()->SetTitle("Flavour ratio");
  rF->SetMinimum(0.0);
  rF->SetMaximum(std::max(0.15, std::max(rF->GetMaximum(), rR->GetMaximum())*1.2));

  TCanvas c("c_fr","Flavour ratios",900,700);
  rF->Draw("hist");
  rR->Draw("hist same");

  TLegend L(0.55,0.18,0.88,0.34);
  L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42); L.SetTextSize(0.04);
  L.AddEntry(rF, "FHC: #nu_{e}/#nu_{#mu}", "l");
  L.AddEntry(rR, "RHC: #bar{#nu}_{e}/#bar{#nu}_{#mu}", "l");
  L.Draw();

  c.SaveAs(outpdf);
}

static void make_cdf_numu(const FluxH& FHC, const FluxH& RHC, const char* outpdf) {
  TH1D* cF = make_cdf(FHC.numu , "cdf_FHC_numu");
  TH1D* cR = make_cdf(RHC.numu , "cdf_RHC_numu");

  style_line(cF, kRed+1, 1, 2);
  style_line(cR, kBlue+2, 2, 2);

  cF->SetTitle("");
  cF->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  cF->GetYaxis()->SetTitle("CDF: #int_{E_{min}}^{E} #phi(E') dE' / #int #phi dE");
  cF->SetMinimum(0.0);
  cF->SetMaximum(1.05);

  TCanvas c("c_cdf","Cumulative flux (nu_{#mu})",900,700);
  cF->Draw("hist");
  cR->Draw("hist same");

  TLegend L(0.55,0.18,0.88,0.34);
  L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42); L.SetTextSize(0.04);
  L.AddEntry(cF, "FHC: #nu_{#mu}", "l");
  L.AddEntry(cR, "RHC: #nu_{#mu}", "l");
  L.Draw();

  c.SaveAs(outpdf);
}

static void make_rhc_over_fhc(const FluxH& FHC, const FluxH& RHC, const char* outpdf) {
  TH1D* r_numu  = safe_ratio(RHC.numu , FHC.numu , "r_numu_RHC_over_FHC");
  TH1D* r_anumu = safe_ratio(RHC.anumu, FHC.anumu, "r_anumu_RHC_over_FHC");

  style_line(r_numu , kRed+1 , 1, 2);
  style_line(r_anumu, kBlue+2, 2, 2);

  r_numu->SetTitle("");
  r_numu->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  r_numu->GetYaxis()->SetTitle("RHC / FHC");
  r_numu->SetMinimum(0.0);
  r_numu->SetMaximum(std::max(2.0, std::max(r_numu->GetMaximum(), r_anumu->GetMaximum())*1.2));

  TCanvas c("c_rhc_fhc","RHC/FHC shape ratio",900,700);
  r_numu->Draw("hist");
  r_anumu->Draw("hist same");

  TLegend L(0.55,0.18,0.88,0.34);
  L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42); L.SetTextSize(0.04);
  L.AddEntry(r_numu , "#nu_{#mu}",  "l");
  L.AddEntry(r_anumu, "#bar{#nu}_{#mu}", "l");
  L.Draw();

  c.SaveAs(outpdf);
}

// ---------------- Parentage stacks (optional; requires correct outTree.nuE) ----------------
struct ParentCat { const char* label; std::function<bool(int)> pass; int color; };

static void make_parentage_stack(const char* file, int ntype_sel, const char* tag, const char* outpdf,
                                 double Emin = 0.0, double Emax = 10.0, int nbins = 400, bool apply_ppfx = true)
{
  // Build from outTree
  TChain ch("outTree"); ch.Add(file);
  float E=0, w=0, cv=1; int ntype=0, ptype=0;
  bool has_ppfx = ch.GetListOfBranches()->FindObject("wgt_ppfx");
  ch.SetBranchAddress("nuE", &E);
  ch.SetBranchAddress("wgt", &w);
  if (has_ppfx) ch.SetBranchAddress("wgt_ppfx", &cv);
  ch.SetBranchAddress("ntype", &ntype);
  ch.SetBranchAddress("ptype", &ptype);

  std::vector<ParentCat> cats = {
    {"#pi^{#pm}", [](int pdg){ return std::abs(pdg)==211; }, kOrange+1},
    {"K^{#pm}",   [](int pdg){ return std::abs(pdg)==321; }, kGreen+2 },
    {"K^{0}",     [](int pdg){ return pdg==130 || pdg==310 || std::abs(pdg)==311; }, kAzure+1},
    {"#mu^{#pm}", [](int pdg){ return std::abs(pdg)==13;  }, kViolet+1}
  };

  std::vector<TH1D*> H;
  H.reserve(cats.size());
  for (size_t i=0;i<cats.size();++i) {
    auto* h = new TH1D(Form("h_%s_%s", cats[i].label, tag), "", nbins, Emin, Emax);
    h->Sumw2();
    style_fill(h, cats[i].color);
    H.push_back(h);
  }

  Long64_t N = ch.GetEntries();
  for (Long64_t i=0;i<N;++i) {
    ch.GetEntry(i);
    if (ntype != ntype_sel) continue;
    double ww = apply_ppfx && has_ppfx ? double(w)*double(cv) : double(w);
    for (size_t k=0;k<cats.size();++k) {
      if (cats[k].pass(ptype)) H[k]->Fill(E, ww);
    }
  }

  // Stack
  THStack st("st",""); 
  for (auto* h : H) st.Add(h);

  // Frame (axes)
  TH1D* frame = (TH1D*)H[0]->Clone(Form("frame_parent_%s", tag));
  frame->Reset("ICES");
  frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
  frame->GetYaxis()->SetTitle("#nu / (6 #times 10^{20} POT) / cm^{2} / bin");
  frame->SetMinimum(1e-9);
  frame->SetMaximum( std::max(1.0, 4.0*std::max_element(H.begin(),H.end(),
                          [](TH1D* a, TH1D* b){ return a->GetMaximum()<b->GetMaximum(); })[0]->GetMaximum()) );

  TCanvas c(Form("c_parent_%s", tag), Form("Parentage %s", tag), 900, 700);
  c.SetLogy();
  frame->Draw("AXIS");
  st.Draw("hist same");

  TLegend L(0.65,0.64,0.88,0.88);
  L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42); L.SetTextSize(0.035);
  for (size_t k=0;k<cats.size();++k) L.AddEntry(H[k], cats[k].label, "f");
  L.Draw();

  c.SaveAs(outpdf);

  delete frame;
  for (auto* h : H) delete h;
}

// ------------------- Main driver -------------------
void plot_flux_extras() {
  set_global_style();

  // ---------------- user paths (match your minimal script) ----------------
  const char* FHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx_g4_10_4.root";
  const char* RHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc_ppfx_g4_10_4.root";

  // Load pre-filled, POT-scaled flux histograms
  auto FHC = load_flux_hists(FHC_FILE, "FHC");
  auto RHC = load_flux_hists(RHC_FILE, "RHC");
  if (!FHC.numu || !RHC.numu) {
    std::cerr << "[plot_flux_extras] Missing flux histograms. Abort.\n";
    return;
  }

  // Bin width for y-axis annotations, if you want to add it to titles
  const int binwMeV = int( std::lround( (FHC.numu->GetXaxis()->GetBinWidth(1))*1000.0 ) );

  // ---------------- make plots ----------------
  make_wrong_sign_fraction(FHC, RHC, "uboone_flux_wrong_sign_fraction.pdf");
  make_flavor_ratios     (FHC, RHC, "uboone_flux_flavour_ratios.pdf");
  make_cdf_numu          (FHC, RHC, "uboone_flux_cdf_numu.pdf");
  make_rhc_over_fhc      (FHC, RHC, "uboone_flux_RHC_over_FHC.pdf");

  // ---------------- optional: parentage stacks (requires correct outTree.nuE) ----------------
  // Uncomment if your outTree.nuE is MicroBooNE energy (post-fix).
  // make_parentage_stack(FHC_FILE, +14, "FHC_numu", "uboone_flux_parentage_FHC_numu.pdf", 0.0, 10.0, 400, true);
  // make_parentage_stack(RHC_FILE, -14, "RHC_anumu", "uboone_flux_parentage_RHC_anumu.pdf", 0.0, 10.0, 400, true);

  std::cout << "[plot_flux_extras] Wrote:\n"
            << "  uboone_flux_wrong_sign_fraction.pdf\n"
            << "  uboone_flux_flavour_ratios.pdf\n"
            << "  uboone_flux_cdf_numu.pdf\n"
            << "  uboone_flux_RHC_over_FHC.pdf\n"
            << "  (parentage stacks available if enabled)\n";
}
