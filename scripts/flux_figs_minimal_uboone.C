// ============================================================================
// flux_figs_minimal_uboone.C — Figures 5.11–5.16, 5.18–5.19 from your file layout
// This version is wired to your structure shown in the ls() outputs:
//   • <flav>/Detsmear/{flav_CV_AV_TPC, Th_flav_CV_TPC, flav_CV_AV_TPC_2D}
//   • <flav>/{PI_Plus, PI_Minus, Kaon_Plus, Kaon_Minus, Mu_Plus, Mu_Minus, K0L}/
//       Enu_flav_<ParentDir>_AV_TPC  and  Th_flav_<ParentDir>_AV_TPC
//   • <flav>/Multisims/flav_ppfx_<category>_Uni_<N>_AV_TPC  (and *_2D we ignore)
// ============================================================================

#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TMatrixD.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <cstdio>
#include <cmath>

namespace CFG {
  // Use your exact files
  const char* FILE_FHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* FILE_RHC = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";
  const char* OUT_PREFIX = "uboone_";

  // For parent fractions on energy plots (avoid μDAR dominance):
  constexpr double E_FRAC_MIN = 0.06; // GeV

  // Flavors to process
  static const std::vector<std::string> FLAVS = {"numu","numubar","nue","nuebar"};

  // Parent directories present in your file
  static const std::map<std::string, std::vector<std::string>> PARENT_DIRS = {
    {"pi", {"PI_Plus","PI_Minus"}},
    {"K",  {"Kaon_Plus","Kaon_Minus"}},
    {"mu", {"Mu_Plus","Mu_Minus"}},
    {"KL", {"K0L"}} // single dir
  };

  // Styling
  int col_pi() { return TColor::GetColor("#e41a1c"); }
  int col_K () { return TColor::GetColor("#377eb8"); }
  int col_mu() { return TColor::GetColor("#984ea3"); }
  int col_KL() { return TColor::GetColor("#4daf4a"); }

  void set_style(){
    TStyle* s = new TStyle("FluxStyle","Flux Style");
    s->SetOptStat(0);
    s->SetTitleFont(42,"XYZ"); s->SetLabelFont(42,"XYZ");
    s->SetTitleSize(0.045,"XYZ"); s->SetLabelSize(0.045,"XYZ");
    s->SetPadLeftMargin(0.14); s->SetPadRightMargin(0.12);
    s->SetPadTopMargin(0.07);  s->SetPadBottomMargin(0.12);
    s->SetPadTickX(1); s->SetPadTickY(1);
    TGaxis::SetMaxDigits(4);
    gROOT->SetStyle("FluxStyle"); gROOT->ForceStyle();
  }
}

// ----------------------- small helpers ---------------------------------------
static TH1D* get_clone_1d(TFile& f, const std::string& path){
  if(auto* o = f.Get(path.c_str())){
    if(o->InheritsFrom(TH1D::Class())){
      TH1D* h = (TH1D*)o->Clone(); h->SetDirectory(nullptr);
      printf("[found] %s\n", path.c_str());
      return h;
    }
  }
  return nullptr;
}
static TH2D* get_clone_2d(TFile& f, const std::string& path){
  if(auto* o = f.Get(path.c_str())){
    if(o->InheritsFrom(TH2D::Class())){
      TH2D* h = (TH2D*)o->Clone(); h->SetDirectory(nullptr);
      printf("[found] %s\n", path.c_str());
      return h;
    }
  }
  return nullptr;
}
static void style_line(TH1* h, int c, int ls=1){ h->SetLineColor(c); h->SetLineWidth(3); h->SetLineStyle(ls); h->SetMarkerSize(0); }
static void style_heat(TH2* h){ h->SetContour(120); }

// integral over [xmin,xmax] with bin-width weighting (true area)
static double integral_area(const TH1* h, double xmin, double xmax){
  if(!h) return 0.0; const TAxis* ax=h->GetXaxis(); if(!ax) return 0.0;
  int bmin=std::max(1, ax->FindFixBin(xmin+1e-12));
  int bmax=std::min(h->GetNbinsX(), ax->FindFixBin(xmax-1e-12));
  double s=0.0;
  for(int b=bmin;b<=bmax;++b){
    double w=ax->GetBinUpEdge(b)-ax->GetBinLowEdge(b);
    s += h->GetBinContent(b)*w;
  }
  return s;
}

// ----------------------- 5.11 — Flux vs angle -------------------------------
static TH1D* angle_cv(TFile& f, const std::string& flav){
  return get_clone_1d(f, flav + "/Detsmear/Th_" + flav + "_CV_TPC");
}
static void fig_5_11_onefile(const char* mode, TFile& f){
  auto* h_nu   = angle_cv(f,"numu");
  auto* h_nub  = angle_cv(f,"numubar");
  auto* h_e    = angle_cv(f,"nue");
  auto* h_eb   = angle_cv(f,"nuebar");
  if(!h_nu||!h_nub||!h_e||!h_eb){
    printf("[5.11/%s] missing one or more Th_*_CV_TPC; skipping.\n", mode);
    delete h_nu; delete h_nub; delete h_e; delete h_eb; return;
  }

  style_line(h_nu,  TColor::GetColor("#e41a1c"),1);
  style_line(h_e,   TColor::GetColor("#e41a1c"),2);
  style_line(h_nub, TColor::GetColor("#1f78b4"),1);
  style_line(h_eb,  TColor::GetColor("#1f78b4"),3);

  double mn=1e99, mx=0;
  for(TH1* h : std::vector<TH1*>{h_nu,h_nub,h_e,h_eb}){
    for(int b=1;b<=h->GetNbinsX();++b){ double y=h->GetBinContent(b); if(y>0 && y<mn) mn=y; if(y>mx) mx=y; }
  }
  if(!(mn<1e98)) mn=1e-18; if(mx<=0) mx=1;

  TCanvas c(Form("c_511_%s",mode),Form("Flux vs angle (%s)",mode),900,700);
  c.SetLogy();
  h_nu->GetXaxis()->SetTitle("#theta_{#nu} [rad]");
  h_nu->GetYaxis()->SetTitle("Flux / 6#times10^{20} POT / bin / cm^{2}");
  h_nu->SetMinimum(mn*0.6); h_nu->SetMaximum(mx*4.0);
  h_nu->Draw("HIST"); h_e->Draw("HIST SAME"); h_nub->Draw("HIST SAME"); h_eb->Draw("HIST SAME");

  TLegend L(0.60,0.70,0.92,0.90); L.SetBorderSize(0); L.SetFillStyle(0); L.SetTextFont(42);
  L.AddEntry(h_nu,  "#nu_{#mu}", "l");
  L.AddEntry(h_nub, "#bar{#nu}_{#mu}", "l");
  L.AddEntry(h_e,   "#nu_{e}", "l");
  L.AddEntry(h_eb,  "#bar{#nu}_{e}", "l");
  L.Draw();
  c.Print(Form("%sfig5_11_%s_angle.pdf", CFG::OUT_PREFIX, mode));

  delete h_nu; delete h_nub; delete h_e; delete h_eb;
}

// ------------- 5.12 — Angle vs parent decay z (search & draw if present) ----
static TH2D* find_theta_vs_z(TFile& f, const std::string& flav){
  // Look for a 2D with "Th"/"theta" and "z"/"zpos" in the name, under flav/
  TH2D* found=nullptr;
  auto* d = f.GetDirectory(flav.c_str());
  if(!d) return nullptr;
  std::vector<TDirectory*> stack = {d};
  while(!stack.empty()){
    TDirectory* cur = stack.back(); stack.pop_back();
    TIter it(cur->GetListOfKeys());
    while(auto* k=(TKey*)it()){
      if(std::string(k->GetClassName())=="TDirectoryFile"){
        stack.push_back( (TDirectory*)cur->Get(k->GetName()) );
      } else if(std::string(k->GetClassName())=="TH2D"){
        std::string name = k->GetName();
        std::string low; low.resize(name.size());
        std::transform(name.begin(),name.end(),low.begin(),::tolower);
        if( (low.find("th_")!=std::string::npos || low.find("theta")!=std::string::npos)
         && (low.find("z")!=std::string::npos   || low.find("zpos")!=std::string::npos) ){
          found = (TH2D*)cur->Get(k->GetName());
          if(found){ found=(TH2D*)found->Clone(); found->SetDirectory(nullptr); printf("[found θ–z] %s/%s\n", cur->GetPath(), k->GetName()); return found; }
        }
      }
    }
  }
  return nullptr;
}
static void fig_5_12_onefile(const char* mode, TFile& f){
  // Try numu first (representative flavor)
  TH2D* h = find_theta_vs_z(f,"numu");
  if(!h){ printf("[5.12/%s] Did not find a θ–z TH2D under numu/; skipping.\n", mode); return; }
  style_heat(h);
  TCanvas c(Form("c_512_%s",mode),Form("#theta vs z (%s)",mode),900,700);
  c.SetRightMargin(0.14);
  if(h->GetXaxis()->GetXmax() > h->GetYaxis()->GetXmax()){
    // If axes are swapped, draw as-is (we can't know without title). Just draw.
  }
  h->Draw("COLZ");
  c.Print(Form("%sfig5_12_%s_theta_vs_z.pdf", CFG::OUT_PREFIX, mode));
  delete h;
}

// -------- 5.13 (energy) & 5.14 (angle) — Parent decomposition ----------------
static TH1D* parent_1d_energy(TFile& f, const std::string& flav, const std::string& parentDir){
  return get_clone_1d(f, flav + "/" + parentDir + "/Enu_" + flav + "_" + parentDir + "_AV_TPC");
}
static TH1D* parent_1d_angle (TFile& f, const std::string& flav, const std::string& parentDir){
  return get_clone_1d(f, flav + "/" + parentDir + "/Th_"  + flav + "_" + parentDir + "_AV_TPC");
}
static TH1D* sum_hists(const std::vector<TH1D*>& hs, const char* name){
  if(hs.empty()) return nullptr;
  for(auto* h:hs) if(!h) return nullptr;
  TH1D* out=(TH1D*)hs.front()->Clone(name); out->SetDirectory(nullptr);
  for(size_t i=1;i<hs.size();++i) out->Add(hs[i]);
  return out;
}

static void fig_5_13_energy_one(const char* mode, TFile& f, const std::string& flav="numu"){
  // Build parent sums exactly from your dirs
  auto Hpi  = std::vector<TH1D*>{ parent_1d_energy(f,flav,"PI_Plus"),   parent_1d_energy(f,flav,"PI_Minus") };
  auto HK   = std::vector<TH1D*>{ parent_1d_energy(f,flav,"Kaon_Plus"), parent_1d_energy(f,flav,"Kaon_Minus") };
  auto Hmu  = std::vector<TH1D*>{ parent_1d_energy(f,flav,"Mu_Plus"),   parent_1d_energy(f,flav,"Mu_Minus") };
  auto HKL  = std::vector<TH1D*>{ parent_1d_energy(f,flav,"K0L") };

  TH1D* h_pi = sum_hists(Hpi,  "pi_E");  TH1D* h_K  = sum_hists(HK,  "K_E");
  TH1D* h_mu = sum_hists(Hmu,  "mu_E");  TH1D* h_KL = sum_hists(HKL, "KL_E");
  if(!h_pi||!h_K||!h_mu||!h_KL){ printf("[5.13/%s] missing parent energy hists (%s)\n", mode, flav.c_str()); return; }

  style_line(h_pi, CFG::col_pi(),1);
  style_line(h_K , CFG::col_K (),1);
  style_line(h_mu, CFG::col_mu(),2);
  style_line(h_KL, CFG::col_KL(),3);

  // Fractions for E>60 MeV
  double e0=CFG::E_FRAC_MIN, e1=h_pi->GetXaxis()->GetXmax();
  double Ipi=integral_area(h_pi,e0,e1), IK=integral_area(h_K,e0,e1),
         Imu=integral_area(h_mu,e0,e1), IKL=integral_area(h_KL,e0,e1);
  double Itot = std::max(1e-300, Ipi+IK+Imu+IKL);

  TCanvas c(Form("c_513_%s_%s",mode,flav.c_str()),Form("Parent E (%s %s)",mode,flav.c_str()),900,700);
  c.SetLogy();
  TH1D* frame=(TH1D*)h_pi->Clone("frame"); frame->Reset(); frame->SetDirectory(nullptr);
  frame->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  frame->GetYaxis()->SetTitle("Flux / 6#times10^{20} POT / GeV / cm^{2}");
  double mn=1e99,mx=0;
  for(TH1* h: {h_pi,h_K,h_mu,h_KL}) for(int b=1;b<=h->GetNbinsX();++b){ double y=h->GetBinContent(b); if(y>0&&y<mn) mn=y; if(y>mx) mx=y; }
  if(!(mn<1e98)) mn=1e-18; if(mx<=0) mx=1;
  frame->SetMinimum(mn*0.6); frame->SetMaximum(mx*4.0);
  frame->Draw("AXIS"); h_pi->Draw("HIST SAME"); h_K->Draw("HIST SAME"); h_mu->Draw("HIST SAME"); h_KL->Draw("HIST SAME");

  TLegend L(0.53,0.68,0.92,0.90); L.SetBorderSize(0); L.SetFillStyle(0);
  L.AddEntry(h_pi, Form("#pi (%.1f%%)",100.*Ipi/Itot), "l");
  L.AddEntry(h_K,  Form("K (%.1f%%)",  100.*IK /Itot), "l");
  L.AddEntry(h_mu, Form("#mu (%.1f%%)",100.*Imu/Itot), "l");
  L.AddEntry(h_KL, Form("K^{0}_{L} (%.1f%%)",100.*IKL/Itot), "l");
  L.Draw();
  c.Print(Form("%sfig5_13_%s_%s_parentE.pdf", CFG::OUT_PREFIX, mode, flav.c_str()));

  delete h_pi; delete h_K; delete h_mu; delete h_KL; delete frame;
}

static void fig_5_14_angle_one(const char* mode, TFile& f, const std::string& flav="numu"){
  auto Hpi  = std::vector<TH1D*>{ parent_1d_angle(f,flav,"PI_Plus"),   parent_1d_angle(f,flav,"PI_Minus") };
  auto HK   = std::vector<TH1D*>{ parent_1d_angle(f,flav,"Kaon_Plus"), parent_1d_angle(f,flav,"Kaon_Minus") };
  auto Hmu  = std::vector<TH1D*>{ parent_1d_angle(f,flav,"Mu_Plus"),   parent_1d_angle(f,flav,"Mu_Minus") };
  auto HKL  = std::vector<TH1D*>{ parent_1d_angle(f,flav,"K0L") };

  TH1D* h_pi = sum_hists(Hpi,  "pi_A");  TH1D* h_K  = sum_hists(HK,  "K_A");
  TH1D* h_mu = sum_hists(Hmu,  "mu_A");  TH1D* h_KL = sum_hists(HKL, "KL_A");
  if(!h_pi||!h_K||!h_mu||!h_KL){ printf("[5.14/%s] missing parent angle hists (%s)\n", mode, flav.c_str()); return; }

  style_line(h_pi, CFG::col_pi(),1);
  style_line(h_K , CFG::col_K (),1);
  style_line(h_mu, CFG::col_mu(),2);
  style_line(h_KL, CFG::col_KL(),3);

  TCanvas c(Form("c_514_%s_%s",mode,flav.c_str()),Form("Parent #theta (%s %s)",mode,flav.c_str()),900,700);
  c.SetLogy();
  TH1D* frame=(TH1D*)h_pi->Clone("frame"); frame->Reset(); frame->SetDirectory(nullptr);
  frame->GetXaxis()->SetTitle("#theta_{#nu} [rad]");
  frame->GetYaxis()->SetTitle("Flux / 6#times10^{20} POT / rad / cm^{2}");
  double mn=1e99,mx=0;
  for(TH1* h: {h_pi,h_K,h_mu,h_KL}) for(int b=1;b<=h->GetNbinsX();++b){ double y=h->GetBinContent(b); if(y>0&&y<mn) mn=y; if(y>mx) mx=y; }
  if(!(mn<1e98)) mn=1e-18; if(mx<=0) mx=1;
  frame->SetMinimum(mn*0.6); frame->SetMaximum(mx*4.0);
  frame->Draw("AXIS"); h_pi->Draw("HIST SAME"); h_K->Draw("HIST SAME"); h_mu->Draw("HIST SAME"); h_KL->Draw("HIST SAME");

  TLegend L(0.53,0.68,0.92,0.90); L.SetBorderSize(0); L.SetFillStyle(0);
  L.AddEntry(h_pi, "#pi", "l"); L.AddEntry(h_K, "K", "l");
  L.AddEntry(h_mu, "#mu", "l"); L.AddEntry(h_KL, "K^{0}_{L}", "l");
  L.Draw();
  c.Print(Form("%sfig5_14_%s_%s_parentAngle.pdf", CFG::OUT_PREFIX, mode, flav.c_str()));

  delete h_pi; delete h_K; delete h_mu; delete h_KL; delete frame;
}

// -------------------------- 5.15 — E–angle heatmaps --------------------------
static TH2D* eangle_cv(TFile& f, const std::string& flav){
  // exactly matches: flav/Detsmear/flav_CV_AV_TPC_2D
  return get_clone_2d(f, flav + "/Detsmear/" + flav + "_CV_AV_TPC_2D");
}
static void fig_5_15_onefile(const char* mode, TFile& f){
  for(const auto& flav : CFG::FLAVS){
    if(auto* h=eangle_cv(f,flav)){
      style_heat(h);
      TCanvas c(Form("c_515_%s_%s",mode,flav.c_str()),Form("E vs #theta (%s %s)",flav.c_str(),mode),900,700);
      c.SetRightMargin(0.14);
      h->Draw("COLZ");
      c.Print(Form("%sfig5_15_%s_%s_EvTheta.pdf", CFG::OUT_PREFIX, mode, flav.c_str()));
      delete h;
    }else{
      printf("[5.15/%s] missing 2D E-θ for %s\n", mode, flav.c_str());
    }
  }
}

// ---- 5.16 + 5.18–5.19 — PPFX universes → total hadron-production cov --------
// We build covariance about the CV spectrum, per PPFX *category* under Multisims,
// then sum covariances across categories to get the total hadron-production cov.
static TH1D* cv_energy_for_ppfx(TFile& f, const std::string& flav){
  // Use the non-5MeV version to match most Multisims hist binnings
  if(auto* h = get_clone_1d(f, flav + "/Detsmear/" + flav + "_CV_AV_TPC")) return h;
  return get_clone_1d(f, flav + "/Detsmear/" + flav + "_CV_AV_TPC_5MeV_bin");
}

struct PPFXCat { std::string name; std::map<int,TH1D*> univ; };
static std::map<std::string,PPFXCat> scan_ppfx_categories(TFile& f, const std::string& flav){
  std::map<std::string,PPFXCat> out;
  auto* d = f.GetDirectory( (flav + "/Multisims").c_str() );
  if(!d) return out;
  TIter it(d->GetListOfKeys());
  while(auto* k=(TKey*)it()){
    if(std::string(k->GetClassName())!="TH1D") continue;
    std::string name = k->GetName();
    // Accept only energy-like hist (exclude 2D)
    if(name.find("_2D")!=std::string::npos) continue;
    // Expect "flav_ppfx_<category>_Uni_<N>_AV_TPC"
    if(name.find(flav+"_ppfx_")!=0) continue;
    size_t posUni = name.find("_Uni_");
    if(posUni==std::string::npos) continue;
    // Extract category between "flav_ppfx_" and "_Uni_"
    std::string cat = name.substr( (flav+"_ppfx_").size(), posUni-(flav+"_ppfx_").size() );
    // Parse universe index after "_Uni_"
    size_t posIdx = posUni + 5;
    size_t posEnd = name.find("_", posIdx);
    int idx = std::stoi(name.substr(posIdx, posEnd-posIdx));

    TH1D* h = (TH1D*)d->Get(name.c_str());
    if(!h) continue;
    h = (TH1D*)h->Clone(); h->SetDirectory(nullptr);

    if(!out.count(cat)) out[cat] = PPFXCat{cat,{}};
    out[cat].univ[idx] = h;
  }
  // Report what we found
  for(auto& kv : out){
    printf("[PPFX] %s: found %zu universes in category '%s'\n",
           flav.c_str(), kv.second.univ.size(), kv.first.c_str());
  }
  return out;
}

static TMatrixD covariance_about_cv(const std::vector<TH1D*>& U, const TH1D* CV){
  const int nb = CV->GetNbinsX();
  TMatrixD C(nb,nb); C.Zero();
  if(U.empty()) return C;
  for(const auto* h : U){
    for(int i=1;i<=nb;++i){
      const double di = h->GetBinContent(i) - CV->GetBinContent(i);
      for(int j=1;j<=nb;++j){
        const double dj = h->GetBinContent(j) - CV->GetBinContent(j);
        C(i-1,j-1) += di*dj;
      }
    }
  }
  const double norm = (U.size()>1)? (1.0/(U.size()-1)) : 0.0;
  C *= norm;
  return C;
}

static void fig_5_16_18_19_oneflavor(const char* mode, TFile& f, const std::string& flav){
  TH1D* CV = cv_energy_for_ppfx(f, flav);
  if(!CV){ printf("[PPFX/%s] missing CV energy for %s\n", mode, flav.c_str()); return; }

  // Scan categories
  auto cats = scan_ppfx_categories(f, flav);
  if(cats.empty()){ printf("[PPFX/%s] no PPFX universes under %s/Multisims\n", mode, flav.c_str()); delete CV; return; }

  // Sum covariances across categories
  const int nb = CV->GetNbinsX();
  TMatrixD Ctot(nb,nb); Ctot.Zero();

  // (Optional) show per-category stats in stdout
  for(auto& kv : cats){
    // Collect in ascending universe index
    std::vector<TH1D*> U; U.reserve(kv.second.univ.size());
    for(auto& ij : kv.second.univ) U.push_back(ij.second);
    TMatrixD Ccat = covariance_about_cv(U, CV);
    Ctot += Ccat;
  }

  // Build ±1σ band on CV from diagonal of total covariance (Fig. 5.16)
  TH1D* up=(TH1D*)CV->Clone("ppfx_up"); TH1D* dn=(TH1D*)CV->Clone("ppfx_dn");
  for(int i=1;i<=nb;++i){
    double s = std::sqrt(std::max(0.0, Ctot(i-1,i-1)));
    up->SetBinContent(i, CV->GetBinContent(i) + s);
    dn->SetBinContent(i, CV->GetBinContent(i) - s);
  }

  TCanvas c(Form("c_516_%s_%s",mode,flav.c_str()),Form("PPFX band (%s %s)",flav.c_str(),mode),900,700);
  c.SetLogy();
  CV->SetLineColor(kBlack); CV->SetLineWidth(3); CV->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  CV->GetYaxis()->SetTitle("Flux / 6#times10^{20} POT / GeV / cm^{2}");
  CV->Draw("HIST"); up->SetLineColor(kGray+2); up->SetLineStyle(2); up->Draw("HIST SAME");
  dn->SetLineColor(kGray+2); dn->SetLineStyle(2); dn->Draw("HIST SAME");
  TLegend L(0.60,0.75,0.92,0.90); L.SetBorderSize(0); L.SetFillStyle(0);
  L.AddEntry(CV,"CV","l"); L.AddEntry(up,"CV #pm 1#sigma (PPFX)","l"); L.Draw();
  c.Print(Form("%sfig5_16_%s_%s_ppfx_band.pdf", CFG::OUT_PREFIX, mode, flav.c_str()));

  // Fractional covariance & correlation (Figs. 5.18–5.19)
  TH2D Hfrac("ppfx_frac",";bin i;bin j", nb,0.5,nb+0.5, nb,0.5,nb+0.5);
  TH2D Hcorr("ppfx_corr",";bin i;bin j", nb,0.5,nb+0.5, nb,0.5,nb+0.5);
  for(int i=1;i<=nb;++i){
    for(int j=1;j<=nb;++j){
      const double mi = CV->GetBinContent(i), mj = CV->GetBinContent(j);
      const double cij = Ctot(i-1,j-1);
      const double denom = (mi!=0.0 && mj!=0.0) ? mi*mj : 1.0;
      Hfrac.SetBinContent(i,j, cij/denom);
      const double sii = Ctot(i-1,i-1), sjj = Ctot(j-1,j-1);
      const double rho = (sii>0 && sjj>0) ? (cij/std::sqrt(sii*sjj)) : 0.0;
      Hcorr.SetBinContent(i,j, rho);
    }
  }
  style_heat(&Hfrac); style_heat(&Hcorr);
  TCanvas cf(Form("c_518_%s_%s",mode,flav.c_str()),"frac cov",900,800); cf.SetRightMargin(0.14); Hfrac.Draw("COLZ"); cf.Print(Form("%sfig5_18_%s_%s_fracCov.pdf", CFG::OUT_PREFIX, mode, flav.c_str()));
  TCanvas cr(Form("c_519_%s_%s",mode,flav.c_str()),"corr",    900,800); cr.SetRightMargin(0.14); Hcorr.Draw("COLZ"); cr.Print(Form("%sfig5_19_%s_%s_corr.pdf",    CFG::OUT_PREFIX, mode, flav.c_str()));

  // cleanup
  delete up; delete dn;
  for(auto& kv : cats) for(auto& ij : kv.second.univ) delete ij.second;
  delete CV;
}

// ------------------------------ DRIVER ---------------------------------------
void flux_figs_minimal_uboone(){
  CFG::set_style();

  TFile fFHC(CFG::FILE_FHC,"READ");
  TFile fRHC(CFG::FILE_RHC,"READ");
  if(fFHC.IsZombie()){ printf("Cannot open FHC file.\n"); return; }
  if(fRHC.IsZombie()){ printf("Cannot open RHC file.\n"); return; }

  // 5.11 — CV flux vs angle (from Th_*_CV_TPC)
  fig_5_11_onefile("FHC", fFHC);
  fig_5_11_onefile("RHC", fRHC);

  // 5.12 — θ vs z (search under numu/, esp. OtherPlots/)
  fig_5_12_onefile("FHC", fFHC);
  fig_5_12_onefile("RHC", fRHC);

  // 5.13 & 5.14 — parent decomposition (use PI_±, Kaon_±, Mu_±, K0L dirs)
  fig_5_13_energy_one("FHC", fFHC, "numu");
  fig_5_13_energy_one("RHC", fRHC, "numu");
  fig_5_14_angle_one ("FHC", fFHC, "numu");
  fig_5_14_angle_one ("RHC", fRHC, "numu");

  // 5.15 — E–θ heatmaps (use flav/Detsmear/flav_CV_AV_TPC_2D)
  fig_5_15_onefile("FHC", fFHC);
  fig_5_15_onefile("RHC", fRHC);

  // 5.16, 5.18–5.19 — PPFX hadron-production ensemble (from Multisims/)
  for(const auto& flav : CFG::FLAVS){
    fig_5_16_18_19_oneflavor("FHC", fFHC, flav);
    fig_5_16_18_19_oneflavor("RHC", fRHC, flav);
  }

  fFHC.Close(); fRHC.Close();
}
