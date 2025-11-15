// root -l -b -q 'scripts/plot_flux_systematics.C++'
// scripts/plot_flux_systematics.C++
//
// Generates main flux-systematics figures + a joint (modes×flavours) correlation heatmap.
//  - Uses per-event weights: CV  = wgt * wgt_ppfx
//                            univ = wgt * wgt_ppfxunivs[u]
//  - Scales to 6e20 POT by reading/adding the 'POT' histogram from each file.
//  - Implements CV-anchored covariance (Sec. 3.3, Eq. (3.4)) and joint covariance (Sec. 3.5, Eq. (3.9)).
//    See the analysis note §3 for definitions.  [Lane et al., Technical Note v0.1]  <-- see PDF in your repo.
//
// Outputs (PDF):
//   fluxsys_FHC_numu_band.pdf
//   fluxsys_RHC_anumu_band.pdf
//   fluxsys_frac_FHC_numu.pdf
//   fluxsys_frac_RHC_anumu.pdf
//   fluxsys_WS_fraction_FHC.pdf
//   fluxsys_WS_fraction_RHC.pdf
//   fluxsys_flavour_ratio_FHC.pdf
//   fluxsys_flavour_ratio_RHC.pdf
//   fluxsys_joint_correlation_blocks.pdf
//
// Author: you (skeleton by ChatGPT)

// ---------------- ROOT includes ----------------
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdio>
#include <iostream>

// ---------------- Global style (matching your minimal script) ----------------
static void set_global_style() {
  const int font_style = 42;
  TStyle* style = new TStyle("PlotterStyle", "Plotter Style");
  style->SetTitleFont(font_style, "XYZ");
  style->SetTitleSize(0.04, "X");
  style->SetTitleSize(0.04, "Y");
  style->SetTitleSize(0.05, "Z");
  style->SetLabelFont(font_style, "XYZ");
  style->SetLabelSize(0.045, "XYZ");
  style->SetLabelOffset(0.005, "XY");
  style->SetTitleOffset(1.10, "XY");
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

static void style_line(TH1* h, int col, int ls=1, int lw=2) {
  h->SetLineColor(col);
  h->SetLineStyle(ls);
  h->SetLineWidth(lw);
  h->SetMarkerSize(0);
}

static void style_band(TGraphAsymmErrors* g, int col, int fill_alpha=100) {
  g->SetFillColorAlpha(col, fill_alpha/255.0);
  g->SetLineColor(col);
  g->SetLineWidth(1);
}

// ---------------- Utilities ----------------
static double sumPOT_in_chain(TChain& ch) {
  double tot = 0.0;
  if (auto* files = ch.GetListOfFiles()) {
    for (int i = 0, n = files->GetEntriesFast(); i < n; ++i) {
      auto* el = (TChainElement*)files->UncheckedAt(i);
      if (!el) continue;
      TFile f(el->GetTitle(), "READ");
      if (!f.IsOpen()) continue;
      if (auto* hp = dynamic_cast<TH1*>(f.Get("POT"))) {
        // producer used a 1-bin hist with binContent = AccumulatedPOT
        tot += hp->GetBinContent(1);
      }
    }
  }
  return tot;
}

static TH1D* fresh_hist(const char* name, int nbins, double Emin, double Emax) {
  if (auto* old = dynamic_cast<TH1D*>(gROOT->FindObject(name))) delete old;
  TH1D* h = new TH1D(name, "", nbins, Emin, Emax);
  h->Sumw2();
  return h;
}

// simple quantile helper (p in [0,1])
static double quantile(std::vector<double>& v, double p) {
  if (v.empty()) return 0.0;
  const size_t n = v.size();
  const double fp = std::clamp(p, 0.0, 1.0) * (n - 1);
  size_t i = (size_t)std::floor(fp);
  size_t j = std::min(i + 1, n - 1);
  std::nth_element(v.begin(), v.begin()+i, v.end());
  double vi = v[i];
  if (j==i) return vi;
  double w = fp - i;
  auto vj_it = v.begin()+j; std::nth_element(v.begin(), vj_it, v.end());
  return (1.0-w)*vi + w*(*vj_it);
}

// ---------- Core container for one mode (FHC or RHC) ----------
struct ModeSpectra {
  // fine binning for detailed plots
  int nbins; double Emin, Emax;
  // CV spectra per flavour: 0:nu_mu, 1:anumu, 2:nue, 3:anue
  TH1D* cv[4] {nullptr,nullptr,nullptr,nullptr};
  // universes: vector-of-hists per flavour
  std::vector<TH1D*> univ[4];
  // number of universes found
  int nU = 0;
};

// Build spectra for a mode by looping the outTree once, filling CV and universes
static ModeSpectra build_mode_spectra(const char* file, const char* tag,
                                      int nbins, double Emin, double Emax) {
  ModeSpectra M; M.nbins=nbins; M.Emin=Emin; M.Emax=Emax;

  TChain ch("outTree"); ch.Add(file);

  // Branches
  float nuE=0.f, wgt=1.f, wgt_ppfx=1.f;
  int ntype=0;
  std::vector<float>* univs=nullptr;

  ch.SetBranchAddress("nuE", &nuE);
  ch.SetBranchAddress("wgt", &wgt);
  if (ch.GetListOfBranches()->FindObject("wgt_ppfx"))
    ch.SetBranchAddress("wgt_ppfx", &wgt_ppfx);
  if (ch.GetListOfBranches()->FindObject("wgt_ppfxunivs"))
    ch.SetBranchAddress("wgt_ppfxunivs", &univs);
  ch.SetBranchAddress("ntype", &ntype);

  // POT scaling to 6e20 POT
  const double NOMINAL_POT = 6e20;
  const double pot_total = sumPOT_in_chain(ch);
  const double pot_scale = (pot_total>0 ? NOMINAL_POT/pot_total : 1.0);
  if (pot_total<=0) std::cerr << "[build_mode_spectra] WARNING: POT<=0 for " << tag << "\n";

  // Make CV hists
  const char* flav_name[4] = {"numu","anumu","nue","anue"};
  const int   flav_pdg [4] = { 14,   -14,   12,  -12};
  for (int f=0; f<4; ++f)
    M.cv[f] = fresh_hist(Form("h_%s_%s_CV", flav_name[f], tag), nbins, Emin, Emax);

  // Discover nU
  Long64_t N = ch.GetEntries();
  int discovered_nU = -1;
  for (Long64_t i=0; i<std::min<Long64_t>(N,1000); ++i) {
    ch.GetEntry(i);
    if (univs && !univs->empty()) { discovered_nU = (int)univs->size(); break; }
  }
  if (discovered_nU<=0) {
    std::cerr << "[build_mode_spectra] ERROR: no PPFX universes found for " << tag
              << " (wgt_ppfxunivs missing or empty). Aborting.\n";
    return M;
  }
  M.nU = discovered_nU;
  for (int f=0; f<4; ++f) {
    M.univ[f].reserve(M.nU);
    for (int u=0; u<M.nU; ++u) {
      M.univ[f].push_back( fresh_hist(Form("h_%s_%s_U%03d", flav_name[f], tag, u), nbins, Emin, Emax) );
    }
  }

  // Event loop
  for (Long64_t i=0; i<N; ++i) {
    ch.GetEntry(i);
    if (!(nuE>=Emin && nuE<Emax)) continue;

    int fidx = -1;
    if      (ntype==14)  fidx = 0;
    else if (ntype==-14) fidx = 1;
    else if (ntype==12)  fidx = 2;
    else if (ntype==-12) fidx = 3;
    else continue;

    const double w_cv = double(wgt) * double(wgt_ppfx) * pot_scale;
    M.cv[fidx]->Fill(nuE, w_cv);

    if (!univs) continue;
    const int thisU = std::min<int>(M.nU, (int)univs->size());
    for (int u=0; u<thisU; ++u) {
      const double w_u = double(wgt) * double((*univs)[u]) * pot_scale;
      M.univ[fidx][u]->Fill(nuE, w_u);
    }
  }

  return M;
}

// ---------------- Bands & ratios ----------------
static TGraphAsymmErrors* band_from_universes(const std::vector<TH1D*>& U, const TH1D* CV,
                                              double qlo=0.16, double qhi=0.84) {
  const int nb = CV->GetNbinsX();
  auto* g = new TGraphAsymmErrors(nb);
  for (int b=1; b<=nb; ++b) {
    const double x  = CV->GetXaxis()->GetBinCenter(b);
    const double ex = CV->GetXaxis()->GetBinWidth(b)/2.0;
    const double c  = CV->GetBinContent(b);

    std::vector<double> vals; vals.reserve(U.size());
    for (auto* h: U) vals.push_back(h->GetBinContent(b));
    std::sort(vals.begin(), vals.end());
    const double lo = quantile(vals, qlo);
    const double hi = quantile(vals, qhi);

    const double y  = c;
    const double eyL= std::max(0.0, c - lo);
    const double eyH= std::max(0.0, hi - c);

    g->SetPoint(b-1, x, y);
    g->SetPointError(b-1, ex, ex, eyL, eyH);
  }
  return g;
}

static TH1D* fractional_uncertainty(const std::vector<TH1D*>& U, const TH1D* CV) {
  const int nb = CV->GetNbinsX();
  TH1D* h = (TH1D*)CV->Clone(Form("%s_frac", CV->GetName())); h->Reset("ICES");
  for (int b=1; b<=nb; ++b) {
    const double c = CV->GetBinContent(b);
    if (c<=0) { h->SetBinContent(b, 0.0); continue; }
    // use std-dev across universes, CV-anchored
    double mu = 0, mu2 = 0; int m=0;
    for (auto* hu: U) { double v = hu->GetBinContent(b); mu += v; mu2 += v*v; ++m; }
    if (m>0) {
      mu  /= m; mu2 /= m;
      double s = std::sqrt(std::max(0.0, mu2 - mu*mu));
      h->SetBinContent(b, s / c);
    }
  }
  return h;
}

// ratio per-universe → band
static TGraphAsymmErrors* ratio_band(const std::vector<TH1D*>& NumU, const std::vector<TH1D*>& DenU,
                                     const TH1D* NumCV, const TH1D* DenCV,
                                     double qlo=0.16, double qhi=0.84) {
  const int nb = NumCV->GetNbinsX();
  auto* g = new TGraphAsymmErrors(nb);
  for (int b=1; b<=nb; ++b) {
    const double x  = NumCV->GetXaxis()->GetBinCenter(b);
    const double ex = NumCV->GetXaxis()->GetBinWidth(b)/2.0;
    const double ncv= NumCV->GetBinContent(b);
    const double dcv= DenCV->GetBinContent(b);
    const double y  = (dcv>0 ? ncv/dcv : 0.0);

    std::vector<double> vals; vals.reserve(NumU.size());
    for (size_t u=0; u<NumU.size(); ++u) {
      const double n = NumU[u]->GetBinContent(b);
      const double d = DenU[u]->GetBinContent(b);
      vals.push_back(d>0 ? n/d : 0.0);
    }
    std::sort(vals.begin(), vals.end());
    const double lo = quantile(vals, qlo);
    const double hi = quantile(vals, qhi);

    const double eyL= std::max(0.0, y - lo);
    const double eyH= std::max(0.0, hi - y);
    g->SetPoint(b-1, x, y);
    g->SetPointError(b-1, ex, ex, eyL, eyH);
  }
  return g;
}

// ---------------- Joint correlation heat map ----------------
// Build coarse-binned spectra for readability in the heat map (fewer bins)
static ModeSpectra build_mode_spectra_coarse(const char* file, const char* tag,
                                             int nbins, double Emin, double Emax) {
  return build_mode_spectra(file, tag, nbins, Emin, Emax);
}

// Stack [FHC 4 flavours | RHC 4 flavours] and form CV-anchored joint correlation
static TH2D* make_joint_correlation(const ModeSpectra& F, const ModeSpectra& R,
                                    const char* name="corr_joint") {
  // Ensure same binning & same number of universes
  const int nb = F.nbins;
  const int nbR = R.nbins;
  const int nU = std::min(F.nU, R.nU);
  if (nU<=0 || nb!=nbR) { std::cerr << "[make_joint_correlation] incompatible inputs.\n"; return nullptr; }

  const int NS = 8; // samples: FHC {numu, anumu, nue, anue}, RHC {numu, anumu, nue, anue}
  const int NbinsTot = NS * nb;

  // Build CV stacked vector
  std::vector<double> T0; T0.reserve(NbinsTot);
  auto append_hist = [&](TH1D* h){ for (int b=1; b<=nb; ++b) T0.push_back(h->GetBinContent(b)); };

  append_hist(F.cv[0]); append_hist(F.cv[1]); append_hist(F.cv[2]); append_hist(F.cv[3]);
  append_hist(R.cv[0]); append_hist(R.cv[1]); append_hist(R.cv[2]); append_hist(R.cv[3]);

  // Covariance (CV-anchored): V = (1/U) Σ (T0 - Tu)(T0 - Tu)^T  (Eq. 3.4 / 3.9)
  std::vector<double> V(NbinsTot*NbinsTot, 0.0);

  std::vector<double> Tu; Tu.resize(NbinsTot);
  for (int u=0; u<nU; ++u) {
    int idx=0;
    for (int s=0; s<4; ++s) for (int b=1; b<=nb; ++b) Tu[idx++] = F.univ[s][u]->GetBinContent(b);
    for (int s=0; s<4; ++s) for (int b=1; b<=nb; ++b) Tu[idx++] = R.univ[s][u]->GetBinContent(b);

    // subtract CV
    for (int i=0; i<NbinsTot; ++i) Tu[i] = T0[i] - Tu[i];

    // outer product add
    for (int i=0; i<NbinsTot; ++i) {
      const double vi = Tu[i];
      if (vi==0.0) continue;
      double* row = &V[i*NbinsTot];
      for (int j=0; j<NbinsTot; ++j) row[j] += vi * Tu[j];
    }
  }
  // average
  const double invU = (nU>0 ? 1.0/nU : 1.0);
  for (double& x : V) x *= invU;

  // Correlation
  std::vector<double> D(NbinsTot, 1.0);
  for (int i=0; i<NbinsTot; ++i) {
    double vii = V[i*NbinsTot + i];
    D[i] = (vii>0 ? 1.0/std::sqrt(vii) : 0.0);
  }

  TH2D* H = new TH2D(name, "", NbinsTot, -0.5, NbinsTot-0.5, NbinsTot, -0.5, NbinsTot-0.5);
  for (int i=0; i<NbinsTot; ++i) {
    for (int j=0; j<NbinsTot; ++j) {
      double cij = (D[i]>0 && D[j]>0 ? V[i*NbinsTot + j]*D[i]*D[j] : 0.0);
      H->SetBinContent(i+1, j+1, cij);
    }
  }
  return H;
}

// Draw the joint correlation with visual block guides
static void draw_joint_correlation(TH2D* H, int nb_per_sample, const char* outpdf) {
  if (!H) return;
  TCanvas c("c_joint","Joint correlation", 1000, 900);
  H->GetXaxis()->SetTitle("Stacked bins: FHC [#nu_{#mu}, #bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e}] | RHC [#nu_{#mu}, #bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e}]");
  H->GetYaxis()->SetTitle("Same ordering");
  H->SetStats(0);
  H->SetContour(100);
  H->Draw("COLZ");

  // block guides
  const int NS=8;
  for (int k=1; k<NS; ++k) {
    const double x = k*nb_per_sample - 0.5;
    TLine* vx = new TLine(x, -0.5, x, NS*nb_per_sample-0.5);
    TLine* hy = new TLine(-0.5, x, NS*nb_per_sample-0.5, x);
    vx->SetLineColor(kGray+2); vx->SetLineStyle(3); vx->Draw();
    hy->SetLineColor(kGray+2); hy->SetLineStyle(3); hy->Draw();
  }

  c.SaveAs(outpdf);
}

// ---------------- Main driver ----------------
void plot_flux_systematics() {
  set_global_style();

  // ---------- user paths ----------
  const char* FHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx_g4_10_4.root";
  const char* RHC_FILE = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_rhc_ppfx_g4_10_4.root";

  // Energy binning for main figures (match central model plots, e.g. 25 MeV bins)
  const double Emin = 0.0, Emax = 10.0;
  const int    nbins = 400; // 0.025 GeV per bin

  // Coarser binning for the big joint heatmap (visual clarity)
  const int    nbins_heat = 50; // 0.2 GeV per bin

  // ---------- Build spectra (fine) ----------
  auto F = build_mode_spectra(FHC_FILE, "FHC", nbins, Emin, Emax);
  auto R = build_mode_spectra(RHC_FILE, "RHC", nbins, Emin, Emax);
  if (F.nU<=0 || R.nU<=0) { std::cerr << "[plot_flux_systematics] Missing universes; abort.\n"; return; }

  // ---------------- 1) Flux with ±1σ bands (RS components) ----------------
  // FHC RS = nu_mu; RHC RS = anumu
  {
    // FHC νμ
    auto* gF = band_from_universes(F.univ[0], F.cv[0], 0.16, 0.84);
    style_band(gF, kRed+1, 120);
    style_line(F.cv[0], kBlack, 1, 2);

    TCanvas c("c1","FHC nu_mu + band",900,700);
    F.cv[0]->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    F.cv[0]->GetYaxis()->SetTitle("#nu / (6 #times 10^{20} POT) / cm^{2} / bin");
    F.cv[0]->SetMinimum(1e-12); F.cv[0]->SetMaximum(F.cv[0]->GetMaximum()*3.0);
    c.SetLogy();
    F.cv[0]->Draw("hist");
    gF->Draw("e2 same");
    F.cv[0]->Draw("hist same");

    TLegend L(0.55,0.75,0.87,0.88); L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42);
    L.AddEntry(F.cv[0],  "FHC: #nu_{#mu} (CV)", "l");
    L.AddEntry(gF,       "PPFX 68% band",      "f");
    L.Draw();
    c.SaveAs("fluxsys_FHC_numu_band.pdf");
  }
  {
    // RHC ν̄μ
    auto* gR = band_from_universes(R.univ[1], R.cv[1], 0.16, 0.84);
    style_band(gR, kBlue+2, 120);
    style_line(R.cv[1], kBlack, 1, 2);

    TCanvas c("c2","RHC anumu + band",900,700);
    R.cv[1]->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    R.cv[1]->GetYaxis()->SetTitle("#nu / (6 #times 10^{20} POT) / cm^{2} / bin");
    R.cv[1]->SetMinimum(1e-12); R.cv[1]->SetMaximum(R.cv[1]->GetMaximum()*3.0);
    c.SetLogy();
    R.cv[1]->Draw("hist");
    gR->Draw("e2 same");
    R.cv[1]->Draw("hist same");

    TLegend L(0.55,0.75,0.87,0.88); L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42);
    L.AddEntry(R.cv[1],  "RHC: #bar{#nu}_{#mu} (CV)", "l");
    L.AddEntry(gR,       "PPFX 68% band",            "f");
    L.Draw();
    c.SaveAs("fluxsys_RHC_anumu_band.pdf");
  }

  // ---------------- 2) Fractional uncertainty (#sigma/#phi) ----------------
  {
    auto* fF = fractional_uncertainty(F.univ[0], F.cv[0]);
    style_line(fF, kRed+1, 1, 2);
    TCanvas c("c3","FHC frac",900,700);
    fF->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    fF->GetYaxis()->SetTitle("Fractional uncertainty #sigma/ #phi");
    fF->SetMinimum(0.0);
    fF->SetMaximum(std::max(0.5, 1.2*fF->GetMaximum()));
    fF->Draw("hist");
    c.SaveAs("fluxsys_frac_FHC_numu.pdf");
  }
  {
    auto* fR = fractional_uncertainty(R.univ[1], R.cv[1]);
    style_line(fR, kBlue+2, 1, 2);
    TCanvas c("c4","RHC frac",900,700);
    fR->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    fR->GetYaxis()->SetTitle("Fractional uncertainty #sigma/ #phi");
    fR->SetMinimum(0.0);
    fR->SetMaximum(std::max(0.5, 1.2*fR->GetMaximum()));
    fR->Draw("hist");
    c.SaveAs("fluxsys_frac_RHC_anumu.pdf");
  }

  // ---------------- 3) Wrong-sign fraction with band ----------------
  // FHC: anumu/(numu+anumu),  RHC: numu/(anumu+numu)
  {
    // Build denominators per universe
    std::vector<TH1D*> F_den, R_den;
    F_den.reserve(F.nU); R_den.reserve(R.nU);
    for (int u=0; u<F.nU; ++u) {
      auto* h = (TH1D*)F.univ[0][u]->Clone(Form("F_den_%03d",u)); // numu
      h->Add(F.univ[1][u]); // + anumu
      F_den.push_back(h);
    }
    for (int u=0; u<R.nU; ++u) {
      auto* h = (TH1D*)R.univ[1][u]->Clone(Form("R_den_%03d",u)); // anumu
      h->Add(R.univ[0][u]); // + numu
      R_den.push_back(h);
    }
    auto* F_denCV = (TH1D*)F.cv[0]->Clone("F_denCV"); F_denCV->Add(F.cv[1]);
    auto* R_denCV = (TH1D*)R.cv[1]->Clone("R_denCV"); R_denCV->Add(R.cv[0]);

    // FHC WS
    auto* gF = ratio_band(/*NumU*/F.univ[1], /*DenU*/F_den, /*NumCV*/F.cv[1], /*DenCV*/F_denCV);
    style_band(gF, kRed+1, 120);
    TCanvas c("c5","FHC WS",900,700);
    gF->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    gF->GetYaxis()->SetTitle("Wrong-sign fraction");
    gF->SetMinimum(0.0);
    gF->SetMaximum(1.0);
    gF->Draw("A e2");
    TLegend LF(0.55,0.18,0.88,0.34); LF.SetBorderSize(0); LF.SetFillStyle(0); LF.SetTextFont(42);
    LF.AddEntry(gF, "FHC: #bar{#nu}_{#mu}/(#nu_{#mu}+#bar{#nu}_{#mu})", "f");
    LF.Draw();
    c.SaveAs("fluxsys_WS_fraction_FHC.pdf");

    // RHC WS
    auto* gR = ratio_band(/*NumU*/R.univ[0], /*DenU*/R_den, /*NumCV*/R.cv[0], /*DenCV*/R_denCV);
    style_band(gR, kBlue+2, 120);
    TCanvas c2("c6","RHC WS",900,700);
    gR->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    gR->GetYaxis()->SetTitle("Wrong-sign fraction");
    gR->SetMinimum(0.0);
    gR->SetMaximum(1.0);
    gR->Draw("A e2");
    TLegend LR(0.55,0.18,0.88,0.34); LR.SetBorderSize(0); LR.SetFillStyle(0); LR.SetTextFont(42);
    LR.AddEntry(gR, "RHC: #nu_{#mu}/(#nu_{#mu}+#bar{#nu}_{#mu})", "f");
    LR.Draw();
    c2.SaveAs("fluxsys_WS_fraction_RHC.pdf");

    // cleanup
    for (auto* h: F_den) delete h; delete F_denCV;
    for (auto* h: R_den) delete h; delete R_denCV;
  }

  // ---------------- 4) Flavour-contamination ratio with band ----------------
  // FHC: nue/numu,  RHC: anue/anumu
  {
    auto* gF = ratio_band(F.univ[2], F.univ[0], F.cv[2], F.cv[0]); // FHC nue/numu
    style_band(gF, kRed+1, 120);
    TCanvas c("c7","FHC flavour ratio",900,700);
    gF->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    gF->GetYaxis()->SetTitle("#nu_{e}/#nu_{#mu}");
    gF->SetMinimum(0.0);
    gF->SetMaximum( std::max(0.2, 1.2*gF->GetHistogram()->GetMaximum()) );
    gF->Draw("A e2");
    TLegend L(0.55,0.18,0.88,0.34); L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42);
    L.AddEntry(gF, "FHC: #nu_{e}/#nu_{#mu}", "f");
    L.Draw();
    c.SaveAs("fluxsys_flavour_ratio_FHC.pdf");
  }
  {
    auto* gR = ratio_band(R.univ[3], R.univ[1], R.cv[3], R.cv[1]); // RHC anue/anumu
    style_band(gR, kBlue+2, 120);
    TCanvas c("c8","RHC flavour ratio",900,700);
    gR->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    gR->GetYaxis()->SetTitle("#bar{#nu}_{e}/#bar{#nu}_{#mu}");
    gR->SetMinimum(0.0);
    gR->SetMaximum( std::max(0.2, 1.2*gR->GetHistogram()->GetMaximum()) );
    gR->Draw("A e2");
    TLegend L(0.55,0.18,0.88,0.34); L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42);
    L.AddEntry(gR, "RHC: #bar{#nu}_{e}/#bar{#nu}_{#mu}", "f");
    L.Draw();
    c.SaveAs("fluxsys_flavour_ratio_RHC.pdf");
  }

  // ---------------- 5) Joint correlation heat map (coarse bins) ----------------
  auto Fc = build_mode_spectra_coarse(FHC_FILE, "FHCc", nbins_heat, Emin, Emax);
  auto Rc = build_mode_spectra_coarse(RHC_FILE, "RHCc", nbins_heat, Emin, Emax);
  auto* Corr = make_joint_correlation(Fc, Rc, "corr_joint");
  draw_joint_correlation(Corr, nbins_heat, "fluxsys_joint_correlation_blocks.pdf");

  std::cout << "[plot_flux_systematics] Wrote:\n"
            << "  fluxsys_FHC_numu_band.pdf\n"
            << "  fluxsys_RHC_anumu_band.pdf\n"
            << "  fluxsys_frac_FHC_numu.pdf\n"
            << "  fluxsys_frac_RHC_anumu.pdf\n"
            << "  fluxsys_WS_fraction_FHC.pdf\n"
            << "  fluxsys_WS_fraction_RHC.pdf\n"
            << "  fluxsys_flavour_ratio_FHC.pdf\n"
            << "  fluxsys_flavour_ratio_RHC.pdf\n"
            << "  fluxsys_joint_correlation_blocks.pdf\n";
}
