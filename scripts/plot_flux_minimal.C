#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

// -----------------------------
// Global style (your original)
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

// --------------------------------------------------------------------
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

// ---------- helpers ----------
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
  frame->SetMinimum(std::max(1e-30, minpos*0.8)); // must be >0 for log
  frame->SetMaximum(maxval*6.0);
}

static double integral_in(double xmin, double xmax, const TH1* h){
  int bmin = std::max(1, h->GetXaxis()->FindFixBin(xmin + 1e-9));
  int bmax = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(xmax - 1e-9));
  return h->Integral(bmin, bmax);
}

// Try to read POT whether it's a TH1 or a TTree
static double read_pot_from_file(const char* path) {
  TFile f(path, "READ");
  if (!f.IsOpen()) return 0.0;

  if (auto* hp = dynamic_cast<TH1*>(f.Get("POT"))) {
    return hp->Integral();
  }
  if (auto* tp = dynamic_cast<TTree*>(f.Get("POT"))) {
    double tot = 0.0;
    Long64_t n = tp->GetEntries();
    for (Long64_t i=0;i<n;++i) {
      tp->GetEntry(i);
      if (auto* leaves = tp->GetListOfLeaves()) {
        // heuristic: sum numeric leaves (typical POT trees have a single leaf)
        for (int j=0, m=leaves->GetEntriesFast(); j<m; ++j) {
          if (auto* L = dynamic_cast<TLeaf*>(leaves->UncheckedAt(j))) {
            tot += L->GetValue();
            break; // first leaf is usually the POT
          }
        }
      }
    }
    return tot;
  }
  return 0.0;
}

static bool file_has_outtree(const char* path) {
  TFile f(path, "READ");
  if (!f.IsOpen()) return false;
  return (dynamic_cast<TTree*>(f.Get("outTree")) != nullptr);
}

// find the first matching TH1(1D) in a directory by simple patterns
static TH1D* find_hist_pass(TDirectory* d,
                            const std::vector<TString>& must,
                            const std::vector<TString>& reject) {
  if (!d) return nullptr;
  TIter next(d->GetListOfKeys());
  while (TKey* k = (TKey*)next()) {
    TString cls = k->GetClassName();
    if (!cls.BeginsWith("TH1")) continue;   // skip TH2/TH3 by class name prefix
    TString nm  = k->GetName();

    bool ok = true;
    for (auto& s : must)   { if (!nm.Contains(s)) { ok=false; break; } }
    for (auto& r : reject) { if ( nm.Contains(r)) { ok=false; break; } }
    if (!ok) continue;

    TH1* h = (TH1*)k->ReadObj();
    if (!h) continue;
    // If it's a TH1 but not TH1D, clone as TH1D
    TH1D* h1d = dynamic_cast<TH1D*>(h);
    if (!h1d) h1d = (TH1D*)h->Clone();
    h1d->SetDirectory(0);
    return h1d;
  }
  return nullptr;
}

// Pick a “best” CV 1D flux histogram for a flavor directory
static TH1D* pick_flux_hist_for_flavor(TFile& f, const char* flavor, TString& picked_from) {
  picked_from = "";
  TDirectory* base = (TDirectory*)f.Get(flavor);
  if (!base) return nullptr;

  const char* subdirs[] = {"OtherPlots", "Detsmear", "Multisims", ""};
  for (const char* sub : subdirs) {
    TDirectory* D = (sub && *sub) ? base->GetDirectory(sub) : base;
    if (!D) continue;

    // 1) prefer CV + AV_TPC (case-sensitive tries)
    TH1D* h = nullptr;
    h = find_hist_pass(D, {"CV","AV_TPC"}, {"2D"});
    if (!h) h = find_hist_pass(D, {"Cv","AV_TPC"}, {"2D"});
    if (!h) h = find_hist_pass(D, {"cv","AV_TPC"}, {"2D"});
    if (h) { picked_from = Form("%s/%s", flavor, (sub && *sub)?sub:""); return h; }

    // 2) fall back to Uni_0 + AV_TPC (typical central universe)
    h = find_hist_pass(D, {"Uni_0","AV_TPC"}, {"2D"});
    if (h) { picked_from = Form("%s/%s", flavor, (sub && *sub)?sub:""); return h; }

    // 3) any AV_TPC 1D
    h = find_hist_pass(D, {"AV_TPC"}, {"2D"});
    if (h) { picked_from = Form("%s/%s", flavor, (sub && *sub)?sub:""); return h; }

    // 4) any TH1 1D
    h = find_hist_pass(D, {}, {"2D"});
    if (h) { picked_from = Form("%s/%s", flavor, (sub && *sub)?sub:""); return h; }
  }
  return nullptr;
}

// ---------------- hardcoded choices ----------------
static constexpr double Emin = 0.0;
static constexpr double Emax = 10.0;

// Normalization (kept from your original)
static constexpr bool   NORM_PER_POT = false;   // false => scale to NOMINAL_POT
static constexpr double NOMINAL_POT  = 6e20;

// For histogram-only files, many are already scaled to 6e20 POT.
// Leave this OFF unless you know your histograms are per-POT and need scaling.
static constexpr bool   SCALE_HISTS_BY_POT = false;

// ---------------- drawing from TREE (your original logic with robust POT) -----------
static void make_and_draw_from_tree(const char* file, const char* tag, const char* outname){

  TChain ch("outTree"); ch.Add(file);
  const bool has_wgt  = ch.GetListOfBranches()->FindObject("wgt");
  const bool has_ppfx = ch.GetListOfBranches()->FindObject("wgt_ppfx");

  TString w = has_wgt ? "wgt" : "1";
  if (has_ppfx) w += "*wgt_ppfx";

  // Sum POT from file(s) (accepts POT as TH1 or TTree)
  double pot_total = 0.0;
  if (auto* files = ch.GetListOfFiles()) {
    for (int i = 0, n = files->GetEntriesFast(); i < n; ++i) {
      auto* el = (TChainElement*)files->UncheckedAt(i);
      if (!el) continue;
      pot_total += read_pot_from_file(el->GetTitle());
    }
  }

  if (pot_total <= 0) {
    std::cerr << "[plot_flux_minimal/" << tag << "] WARNING: total POT <= 0; proceeding without POT scaling.\n";
  }

  char w_scaled[256];
  if (pot_total > 0) {
    if (NORM_PER_POT) std::snprintf(w_scaled, sizeof(w_scaled), "(%s)/(%g)", w.Data(), pot_total);
    else              std::snprintf(w_scaled, sizeof(w_scaled), "(%s)*(%g)", w.Data(), NOMINAL_POT / pot_total);
  } else {
    std::snprintf(w_scaled, sizeof(w_scaled), "%s", w.Data());
  }

  auto fresh = [&](const char* name) {
    if (auto* old = dynamic_cast<TH1D*>(gROOT->FindObject(name))) delete old;
    TH1D* h = new TH1D(name, "", int((Emax-Emin)/0.025 + 0.5), Emin, Emax); // 25 MeV default
    h->Sumw2();
    return h;
  };

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

  style_line(h_numu,  kRed+1, 1);
  style_line(h_nue,   kRed+1, 2);
  style_line(h_anumu, kBlue+2,1);
  style_line(h_anue,  kBlue+2,3);

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

  // Axis frame
  TH1D* frame = (TH1D*)h_numu->Clone(Form("frame_%s",tag));
  frame->Reset("ICES");
  frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");

  const int binwMeV = int(std::lround(h_numu->GetXaxis()->GetBinWidth(1)*1000.0));
  TString potLabel = NORM_PER_POT ? "POT" : "6 #times 10^{20} POT";
  TString yTitle = Form("#nu / %s / %d MeV / cm^{2}", potLabel.Data(), binwMeV);
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

  TLegend* L = build_flux_legend_like_stacked(p_leg, h_numu, h_anumu, h_nue, h_anue,
                                              split, s_numu, s_anumu, s_nue, s_anue, s_tot);
  L->Draw();

  c.cd(); c.Update(); c.SaveAs(outname);

  delete frame; delete L; delete p_main; delete p_leg;
  delete h_numu; delete h_anumu; delete h_nue; delete h_anue;

  std::cout << "[plot_flux_minimal] " << tag
            << " | bin width = " << (binwMeV/1000.0) << " GeV (" << binwMeV << " MeV)"
            << " | y-axis: " << yTitle << "\n";
}

// ---------------- drawing from HISTOGRAMS ------------------------------
static void make_and_draw_from_hists(const char* file, const char* tag, const char* outname){

  TFile f(file,"READ");
  if (!f.IsOpen()) { std::cerr << "ERROR: cannot open " << file << "\n"; return; }

  TString pickedFrom[4];
  TH1D* h_numu  = pick_flux_hist_for_flavor(f, "numu",    pickedFrom[0]);
  TH1D* h_anumu = pick_flux_hist_for_flavor(f, "numubar", pickedFrom[1]);
  TH1D* h_nue   = pick_flux_hist_for_flavor(f, "nue",     pickedFrom[2]);
  TH1D* h_anue  = pick_flux_hist_for_flavor(f, "nuebar",  pickedFrom[3]);

  if (!h_numu || !h_anumu || !h_nue || !h_anue) {
    std::cerr << "[plot_flux_minimal/" << tag << "] ERROR: could not locate 1D flux hists for all flavors.\n";
    return;
  }

  std::cout << "[pick/" << tag << "] numu   <= " << pickedFrom[0] << " : " << h_numu->GetName()  << "\n";
  std::cout << "[pick/" << tag << "] numubar<= " << pickedFrom[1] << " : " << h_anumu->GetName() << "\n";
  std::cout << "[pick/" << tag << "] nue    <= " << pickedFrom[2] << " : " << h_nue->GetName()   << "\n";
  std::cout << "[pick/" << tag << "] nuebar <= " << pickedFrom[3] << " : " << h_anue->GetName()  << "\n";

  // Optional POT scaling for histogram mode (OFF by default)
  if (SCALE_HISTS_BY_POT) {
    double pot = read_pot_from_file(file);
    if (pot>0) {
      const double sf = NORM_PER_POT ? (1.0/pot) : (NOMINAL_POT/pot);
      h_numu->Scale(sf); h_anumu->Scale(sf); h_nue->Scale(sf); h_anue->Scale(sf);
      std::cout << "[scale/" << tag << "] Applied POT scaling factor " << sf
                << " using POT=" << pot << "\n";
    } else {
      std::cerr << "[scale/" << tag << "] WARNING: POT not found or <=0; not scaling.\n";
    }
  }

  style_line(h_numu,  kRed+1, 1);
  style_line(h_nue,   kRed+1, 2);
  style_line(h_anumu, kBlue+2,1);
  style_line(h_anue,  kBlue+2,3);

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

  // Use a simple frame over [Emin,Emax]
  TH1D* frame = new TH1D(Form("frame_%s",tag), "", 100, Emin, Emax);
  frame->GetXaxis()->SetTitle("Neutrino Energy [GeV]");

  const int binwMeV = int(std::lround(h_numu->GetXaxis()->GetBinWidth(1)*1000.0));
  TString potLabel = NORM_PER_POT ? "POT" : "6 #times 10^{20} POT";
  TString yTitle = Form("#nu / %s / %d MeV / cm^{2}", potLabel.Data(), binwMeV);
  frame->GetYaxis()->SetTitle(yTitle);

  auto_logy_limits_range(frame, {h_numu,h_anumu,h_nue,h_anue}, Emin, Emax);
  frame->Draw("AXIS");

  // Draw, restricting x-range for display
  h_numu ->GetXaxis()->SetRangeUser(Emin, Emax);
  h_nue  ->GetXaxis()->SetRangeUser(Emin, Emax);
  h_anumu->GetXaxis()->SetRangeUser(Emin, Emax);
  h_anue ->GetXaxis()->SetRangeUser(Emin, Emax);

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

  TLegend* L = build_flux_legend_like_stacked(p_leg, h_numu, h_anumu, h_nue, h_anue,
                                              split, s_numu, s_anumu, s_nue, s_anue, s_tot);
  L->Draw();

  c.cd(); c.Update(); c.SaveAs(outname);

  delete frame; delete L; delete p_main; delete p_leg;
  delete h_numu; delete h_anumu; delete h_nue; delete h_anue;

  std::cout << "[plot_flux_minimal] " << tag
            << " | bin width ≈ " << (binwMeV/1000.0) << " GeV (" << binwMeV << " MeV)"
            << " | y-axis: " << yTitle << "\n";
}

// ---------------- orchestrator ----------------
void plot_flux_minimal() {
  set_global_style();

  const char* FHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* RHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";

  const char* OUT_FHC = "uboone_flux_FHC.pdf";
  const char* OUT_RHC = "uboone_flux_RHC.pdf";

  auto run_one = [&](const char* file, const char* tag, const char* out) {
    if (file_has_outtree(file)) {
      std::cout << "[mode/" << tag << "] Detected outTree: using tree-based method.\n";
      make_and_draw_from_tree(file, tag, out);
    } else {
      std::cout << "[mode/" << tag << "] No outTree found: using histogram-based method.\n";
      make_and_draw_from_hists(file, tag, out);
    }
  };

  run_one(FHC_FILE, "FHC", OUT_FHC);
  run_one(RHC_FILE, "RHC", OUT_RHC);

  std::cout << "Saved: " << OUT_FHC << " and " << OUT_RHC << std::endl;
}
