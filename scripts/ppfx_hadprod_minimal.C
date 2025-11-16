// ============================================================================
// ppfx_hadprod_minimal.C
//
// Figure 3: Hadron-production uncertainties (PPFX universes)
//   E) PPFX ensemble band & covariance (per flavor)
//   F) Cross-flavor fractional correlation (FHC & RHC side-by-side)
//
// What it does
//   - Scans */Multisims/ for universes matching:
//       <flav>_ppfx_<category>_Uni_<N>_AV_TPC             (energy 1D)
//       (skips *_2D)
//   - Draws all universes (spaghetti), overlays CV, ensemble mean, and ±1σ band
//     where σ is from the **summed per-category covariance about the CV**.
//   - Builds per-flavor covariances and a **joint** covariance across flavors
//     (stacked vector across samples; same method used in fit scaffolding).
//   - Exports per-flavor and cross-flavor: covariance, fractional covariance,
//     and correlation heatmaps, plus a ROOT file with TMatrixD/TVectorD objects.
//
// References (visual & methodology context):
//   NuMI_Flux_Reweighting_Scheme_co… (see Figs. 6–8)
//
// Usage:
//   root -l -b -q 'ppfx_hadprod_minimal.C()'
//   root -l -b -q 'ppfx_hadprod_minimal.C("/path/FHC.root","/path/RHC.root",50)'
// ============================================================================

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdio>

// ------------------------------ CONFIG ----------------------------------------
namespace CFG {
  // Default files (replace if needed when calling the macro)
  const char* FILE_FHC_DEF = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* FILE_RHC_DEF = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";

  // Flavors (order also used for stacking in the joint covariance)
  static const std::vector<std::string> FLAVS = {"numu","numubar","nue","nuebar"};

  // Default target energy bin width for cross-flavor (MeV); set 0 to skip rebin
  int TARGET_BIN_MEV = 50;
}

// ------------------------------ STYLE -----------------------------------------
static void set_global_style(){
  const int f=42;
  TStyle* s=new TStyle("PlotterStyle","Plotter Style");
  s->SetTitleFont(f,"X"); s->SetTitleFont(f,"Y"); s->SetTitleFont(f,"Z");
  s->SetTitleSize(0.04,"X"); s->SetTitleSize(0.04,"Y"); s->SetTitleSize(0.05,"Z");
  s->SetLabelFont(f,"X"); s->SetLabelFont(f,"Y"); s->SetLabelFont(f,"Z");
  s->SetLabelSize(0.045,"X"); s->SetLabelSize(0.045,"Y"); s->SetLabelSize(0.045,"Z");
  s->SetLabelOffset(0.005,"X"); s->SetLabelOffset(0.005,"Y"); s->SetLabelOffset(0.005,"Z");
  s->SetTitleOffset(1.10,"X"); s->SetTitleOffset(1.10,"Y");
  s->SetOptStat(0); s->SetOptTitle(0);
  s->SetPadTickX(1); s->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  s->SetPadLeftMargin(0.15); s->SetPadRightMargin(0.12);
  s->SetPadTopMargin(0.07);  s->SetPadBottomMargin(0.12);
  s->SetMarkerSize(1.0);
  s->SetCanvasColor(0); s->SetPadColor(0); s->SetFrameFillColor(0);
  s->SetCanvasBorderMode(0); s->SetPadBorderMode(0); s->SetStatColor(0); s->SetFrameBorderMode(0);
  s->SetTitleFillColor(0); s->SetTitleBorderSize(0);
  gROOT->SetStyle("PlotterStyle"); gROOT->ForceStyle();
}

static void style_line(TH1* h,int col,int ls=1,double lw=1.2){
  if(!h) return;
  h->SetLineColor(col); h->SetLineStyle(ls); h->SetLineWidth(lw); h->SetMarkerSize(0);
}

static void make_split_canvas(const char* cname,const char* ctitle,double split,
                              TCanvas*& c, TPad*& p_main, TPad*& p_leg,
                              bool logy=false){
  c = new TCanvas(cname, ctitle, 1200, 800);
  p_main = new TPad("pad_main","pad_main",0.,0.00,1.,split);
  p_leg  = new TPad("pad_legend","pad_legend",0.,split,1.,1.00);
  p_main->SetTopMargin(0.01); p_main->SetBottomMargin(0.12);
  p_main->SetLeftMargin(0.12); p_main->SetRightMargin(0.12);
  if(logy) p_main->SetLogy();
  p_leg->SetTopMargin(0.05); p_leg->SetBottomMargin(0.01);
  p_leg->SetLeftMargin(0.02); p_leg->SetRightMargin(0.02);
  p_main->Draw(); p_leg->Draw();
}

static TLegend* build_legend_toppad(TPad* p_leg,
                                    const std::vector<std::pair<TH1*,TString>>& items,
                                    double split){
  p_leg->cd();
  TLegend* L=new TLegend(0.12,0.02,0.95,0.78);
  L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextFont(42);
  int n_entries=(int)items.size(); int n_cols=(n_entries>4)?3:2;
  L->SetNColumns(n_cols); L->SetColumnSeparation(0.08);
  L->SetEntrySeparation(0.00); L->SetMargin(0.25);
  const double s_main=0.045; const double s_leg=s_main*(split/(1.0-split));
  L->SetTextSize(s_leg);
  for(auto& it: items) L->AddEntry(it.first, it.second, (it.first && it.first->GetFillColorAlpha()>0)? "f":"l");
  return L;
}

static bool same_binning(const TH1* a, const TH1* b){
  if(!a || !b) return false;
  if(a->GetNbinsX()!=b->GetNbinsX()) return false;
  for(int i=1;i<=a->GetNbinsX()+1;++i){
    double ea=a->GetXaxis()->GetBinLowEdge(i), eb=b->GetXaxis()->GetBinLowEdge(i);
    if(std::abs(ea-eb)>1e-9) return false;
  }
  return true;
}

// Create a symmetric blue-white-red palette and lock z range to [-1,1]
static void set_diverging_palette_11(){
  const Int_t NRGBs = 3;
  Double_t stops[NRGBs] = {0.00, 0.50, 1.00};
  Double_t red[NRGBs]   = {0.231, 1.00, 0.706};
  Double_t green[NRGBs] = {0.298, 1.00, 0.016};
  Double_t blue[NRGBs]  = {0.753, 1.00, 0.150};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
}

// ------------------------------ IO helpers ------------------------------------
static TH1D* get_cv_energy(TFile& f, const std::string& flav){
  // Prefer the 5 MeV-binned CV if present; otherwise the default
  TH1D* h = (TH1D*)f.Get((flav + "/Detsmear/" + flav + "_CV_AV_TPC_5MeV_bin").c_str());
  if(!h) h = (TH1D*)f.Get((flav + "/Detsmear/" + flav + "_CV_AV_TPC").c_str());
  if(!h) return nullptr;
  TH1D* c = (TH1D*)h->Clone(Form("cv_%s",flav.c_str()));
  c->SetDirectory(0);
  return c;
}

struct PPFXCat { std::string name; std::map<int,TH1D*> univ; };

// Scan */Multisims for TH1D universes (ignore any *_2D)
static std::map<std::string,PPFXCat> scan_ppfx_categories(TFile& f, const std::string& flav){
  std::map<std::string,PPFXCat> out;
  auto* d = f.GetDirectory( (flav + "/Multisims").c_str() );
  if(!d) return out;

  TIter it(d->GetListOfKeys());
  while(auto* k=(TKey*)it()){
    if(std::string(k->GetClassName())!="TH1D") continue;
    std::string name = k->GetName();
    if(name.find("_2D")!=std::string::npos) continue;
    // Expected: <flav>_ppfx_<category>_Uni_<N>_AV_TPC
    if(name.find(flav+"_ppfx_")!=0) continue;
    size_t posUni = name.find("_Uni_");
    if(posUni==std::string::npos) continue;

    std::string cat = name.substr( (flav+"_ppfx_").size(), posUni-(flav+"_ppfx_").size() );
    size_t posIdx = posUni + 5;
    size_t posEnd = name.find("_", posIdx);
    int idx = std::stoi(name.substr(posIdx, posEnd-posIdx));

    auto* h=(TH1D*)d->Get(name.c_str());
    if(!h) continue;
    auto* c=(TH1D*)h->Clone(Form("cl_%s",name.c_str()));
    c->SetDirectory(0);

    if(!out.count(cat)) out[cat] = PPFXCat{cat,{}};
    out[cat].univ[idx] = c;
  }
  return out;
}

// Sample covariance about the CV
static TMatrixD covariance_about_cv(const std::vector<TH1D*>& U, const TH1D* CV){
  const int nb = CV->GetNbinsX();
  TMatrixD C(nb,nb); C.Zero();
  if(U.empty()) return C;
  const int N = (int)U.size();
  for(const auto* h : U){
    for(int i=1;i<=nb;++i){
      const double di = h->GetBinContent(i) - CV->GetBinContent(i);
      for(int j=1;j<=nb;++j){
        const double dj = h->GetBinContent(j) - CV->GetBinContent(j);
        C(i-1,j-1) += di*dj;
      }
    }
  }
  if(N>1) C *= (1.0/(N-1));
  return C;
}

// Rebin to fixed width (MeV) if a 5 MeV-binned CV is available; otherwise return clone.
static TH1D* rebin_to_common_or_clone(const TH1D* h_in, int target_bin_mev){
  if(!h_in || target_bin_mev<=0) return (TH1D*)h_in->Clone(Form("%s_reb",h_in->GetName()));
  const double w0 = h_in->GetXaxis()->GetBinWidth(1);
  // If already integer multiple of target width, prefer simple Rebin
  const int target_factor = (int)std::lround((target_bin_mev/1000.0)/w0);
  if(target_factor>0 && std::abs(target_factor*w0 - (target_bin_mev/1000.0))<1e-9){
    TH1D* h = (TH1D*)h_in->Clone(Form("%s_reb",h_in->GetName()));
    h->SetDirectory(0);
    TH1D* out = (TH1D*)h->Rebin(target_factor, Form("%s_r%02d",h->GetName(),target_factor));
    out->SetDirectory(0);
    delete h;
    return out;
  }
  // Otherwise, return a safe clone (we’ll only use this if all flavors match already)
  TH1D* h = (TH1D*)h_in->Clone(Form("%s_reb",h_in->GetName()));
  h->SetDirectory(0);
  return h;
}

// ------------------------------ E: per-flavor band ----------------------------
static void draw_ppfx_ensemble_and_band(const char* mode, TFile& f, const std::string& flav){
  TH1D* CV = get_cv_energy(f, flav);
  if(!CV){ printf("[E/%s] Missing CV for %s\n",mode,flav.c_str()); return; }

  // Gather universes by category
  auto cats = scan_ppfx_categories(f,flav);
  if(cats.empty()){ printf("[E/%s] No PPFX universes for %s\n",mode,flav.c_str()); delete CV; return; }

  // Build total covariance (sum of per-category covariances about CV)
  const int nb = CV->GetNbinsX();
  TMatrixD Ctot(nb,nb); Ctot.Zero();
  std::vector<TH1D*> allU; allU.reserve(1024);

  for(auto& kv : cats){
    std::vector<TH1D*> U; U.reserve(kv.second.univ.size());
    for(auto& ij : kv.second.univ){ U.push_back(ij.second); allU.push_back(ij.second); }
    TMatrixD Ccat = covariance_about_cv(U, CV);
    Ctot += Ccat;
  }

  // Ensemble mean (across all universes from all categories; for visualization only)
  TH1D* Hmean = (TH1D*)CV->Clone(Form("ppfx_mean_%s_%s",mode,flav.c_str())); Hmean->Reset(); Hmean->SetDirectory(0);
  if(!allU.empty()){
    for(auto* u : allU) Hmean->Add(u);
    Hmean->Scale(1.0/double(allU.size()));
  }

  // ±1σ band from diag(Ctot)
  TH1D* hUp=(TH1D*)CV->Clone("ppfx_up"); TH1D* hDn=(TH1D*)CV->Clone("ppfx_dn");
  for(int i=1;i<=nb;++i){
    const double s = std::sqrt(std::max(0.0, Ctot(i-1,i-1)));
    hUp->SetBinContent(i, CV->GetBinContent(i)+s);
    hDn->SetBinContent(i, CV->GetBinContent(i)-s);
  }

  // Canvas (legend up top)
  const double split=0.83;
  TCanvas* c=nullptr; TPad* p_main=nullptr; TPad* p_leg=nullptr;
  make_split_canvas(Form("c_E_%s_%s",mode,flav.c_str()),
                    Form("PPFX universes & band — %s %s",mode,flav.c_str()),
                    split,c,p_main,p_leg,/*logy=*/true);

  p_main->cd();
  // Frame y-range: include CV, band, and faint universes
  double xmin=CV->GetXaxis()->GetXmin(), xmax=CV->GetXaxis()->GetXmax();
  double ymin=1e30,ymax=0;
  auto update_range=[&](const TH1D* h){
    if(!h) return;
    for(int b=1;b<=h->GetNbinsX();++b){
      const double x=h->GetXaxis()->GetBinCenter(b);
      if(x<xmin||x>xmax) continue;
      double y=h->GetBinContent(b);
      if(y>0 && y<ymin) ymin=y;
      if(y>ymax) ymax=y;
    }
  };
  update_range(CV); update_range(hUp); update_range(hDn);
  // (avoid scanning all universes for speed; the band already brackets them)

  if(!(ymin<1e29)) ymin=1e-18;
  if(ymax<=0) ymax=1;

  TH1D* frame=(TH1D*)CV->Clone("frame_E"); frame->Reset(); frame->SetDirectory(0);
  frame->GetXaxis()->SetTitle("E_{#nu} [GeV]");
  int binwMeV=(int)std::lround(CV->GetXaxis()->GetBinWidth(1)*1000.0);
  frame->GetYaxis()->SetTitle(Form("Flux / 6 #times 10^{20} POT / %d MeV / cm^{2}",binwMeV));
  frame->SetMinimum(ymin*0.6); frame->SetMaximum(ymax*4.0);
  frame->Draw("AXIS");

  // Spaghetti universes (faint)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,08,00)
  const int ucol = TColor::GetColor("#1f78b4");
  for(auto* u : allU){ style_line(u,ucol,1,1.0); u->SetLineColorAlpha(ucol,0.12); u->Draw("HIST SAME"); }
#else
  for(auto* u : allU){ style_line(u,38,1,1.0); u->Draw("HIST SAME"); }
#endif

  // Band + CV + mean
  style_line(CV, TColor::GetColor("#222222"),1,1.6);
  style_line(Hmean, TColor::GetColor("#4d4d4d"),2,1.4);

  style_line(hUp,TColor::GetColor("#1f78b4"),1,1.0);
  style_line(hDn,TColor::GetColor("#1f78b4"),1,1.0);
  hUp->SetFillColorAlpha(TColor::GetColor("#1f78b4"),0.25);
  hDn->SetFillColorAlpha(TColor::GetColor("#1f78b4"),0.25);
  hUp->Draw("E3 SAME");  // draws filled band between up & dn when paired
  hDn->Draw("E3 SAME");
  CV->Draw("HIST SAME");
  Hmean->Draw("HIST SAME");

  p_leg->cd();
  auto L = build_legend_toppad(p_leg,
    {{CV,"CV"},
     {Hmean,"Mean of universes"},
     {hUp,"CV #pm 1#sigma (from #Sigma_{ppfx})"},
     {(TH1D*)nullptr,"Universes (faint)"}},
    split);
  L->Draw();

  c->Update();
  c->Print(Form("ppfx_E_ensemble_band_%s_%s.pdf",mode,flav.c_str()));
  c->Print(Form("ppfx_E_ensemble_band_%s_%s.png",mode,flav.c_str()));

  // --- Per-flavor fractional covariance & correlation (energy) ---
  TH2D Hfrac(Form("ppfx_frac_%s_%s",mode,flav.c_str()),
             "PPFX fractional covariance;bin i;bin j", nb,0.5,nb+0.5, nb,0.5,nb+0.5);
  TH2D Hcorr(Form("ppfx_corr_%s_%s",mode,flav.c_str()),
             "PPFX correlation;bin i;bin j", nb,0.5,nb+0.5, nb,0.5,nb+0.5);
  for(int i=1;i<=nb;++i){
    for(int j=1;j<=nb;++j){
      const double mi = CV->GetBinContent(i), mj = CV->GetBinContent(j);
      const double cij = Ctot(i-1,j-1);
      const double sii = Ctot(i-1,i-1), sjj = Ctot(j-1,j-1);
      const double frac = (mi!=0 && mj!=0)? (cij/(mi*mj)) : 0.0;
      const double rho  = (sii>0 && sjj>0)? (cij/std::sqrt(sii*sjj)) : 0.0;
      Hfrac.SetBinContent(i,j, frac);
      Hcorr.SetBinContent(i,j, rho);
    }
  }

  TCanvas cf(Form("c_E_frac_%s_%s",mode,flav.c_str()),"frac",900,800);
  cf.SetRightMargin(0.14); cf.SetLogz();
  Hfrac.SetContour(255);
  Hfrac.GetZaxis()->SetTitle("Fractional covariance");
  Hfrac.Draw("COLZ");
  cf.Print(Form("ppfx_E_fracCov_%s_%s.pdf",mode,flav.c_str()));

  TCanvas cr(Form("c_E_corr_%s_%s",mode,flav.c_str()),"corr",900,800);
  cr.SetRightMargin(0.14);
  set_diverging_palette_11();
  Hcorr.SetMinimum(-1.0); Hcorr.SetMaximum( 1.0);
  Hcorr.GetZaxis()->SetTitle("Correlation coefficient #rho");
  Hcorr.Draw("COLZ");
  cr.Print(Form("ppfx_E_corr_%s_%s.pdf",mode,flav.c_str()));

  // --- Persist per-flavor outputs for propagation
  TFile fout("ppfx_hadprod_outputs.root","UPDATE");
  fout.cd();
  Ctot.Write(Form("C_ppfx_%s_%s",mode,flav.c_str()), TObject::kOverwrite);
  TVectorD cvv(nb); for(int i=0;i<nb;++i) cvv[i]=CV->GetBinContent(i);
  cvv.Write(Form("CV_%s_%s",mode,flav.c_str()), TObject::kOverwrite);
  Hfrac.Write("", TObject::kOverwrite);
  Hcorr.Write("", TObject::kOverwrite);
  fout.Close();

  // cleanup (keep cats’ cloned universes alive until end of function)
  delete L; delete frame; delete c; delete hUp; delete hDn; delete Hmean; delete CV;
}

// ------------------------------ F: cross-flavor joint -------------------------
struct JointPack {
  // Stacked CV and covariance over flavors (in order CFG::FLAVS)
  TVectorD cv_stacked;
  TMatrixD C_joint;
  int nb; // bins per flavor
};

static bool build_joint_covariance(TFile& f, const char* mode,
                                   int target_bin_mev,
                                   JointPack& JP){
  // 1) Load CVs (rebinned if a uniform 5 MeV CV exists)
  std::vector<TH1D*> cvs; cvs.reserve(CFG::FLAVS.size());
  for(const auto& flav : CFG::FLAVS){
    TH1D* cv = get_cv_energy(f,flav);
    if(!cv){ printf("[F/%s] Missing CV for %s\n",mode,flav.c_str()); return false; }
    TH1D* cv_reb = rebin_to_common_or_clone(cv, target_bin_mev);
    cvs.push_back(cv_reb);
    delete cv; // cloned copy kept
  }

  // Ensure common binning across flavors
  for(size_t i=1;i<cvs.size();++i){
    if(!same_binning(cvs[0],cvs[i])){
      printf("[F/%s] CV binnings are not common across flavors; using native if they match, else abort.\n",mode);
      // try native match (target_bin_mev=0 would have done this); if mismatch, abort
      for(auto* h : cvs) delete h;
      return false;
    }
  }
  const int nb = cvs[0]->GetNbinsX();
  const int NF = (int)CFG::FLAVS.size();
  const int Ntot = nb*NF;

  // 2) Scan universes per flavor and category
  std::map<std::string, std::map<std::string,PPFXCat>> cats_by_flav; // flav -> cat -> PPFXCat
  for(const auto& flav : CFG::FLAVS){
    cats_by_flav[flav] = scan_ppfx_categories(f,flav);
    if(cats_by_flav[flav].empty()){
      printf("[F/%s] No PPFX universes for %s\n",mode,flav.c_str());
      for(auto* h : cvs) delete h;
      return false;
    }
  }

  // 3) For each category, find the set of universe indices common to ALL flavors
  TMatrixD C(Ntot,Ntot); C.Zero();
  for(const auto& catName : cats_by_flav[CFG::FLAVS[0]]){
    const std::string& cat = catName.first;
    // Skip if any flavor misses this category
    bool all_have = true;
    for(const auto& flav : CFG::FLAVS) if(!cats_by_flav[flav].count(cat)) { all_have=false; break; }
    if(!all_have) continue;

    // Intersection of indices
    std::set<int> idx = {};
    bool first=true;
    for(const auto& flav : CFG::FLAVS){
      std::set<int> s;
      for(const auto& ij : cats_by_flav[flav][cat].univ) s.insert(ij.first);
      if(first){ idx = s; first=false; }
      else {
        std::set<int> inter;
        std::set_intersection(idx.begin(),idx.end(), s.begin(),s.end(),
                              std::inserter(inter,inter.begin()));
        idx.swap(inter);
      }
    }
    if(idx.empty()) continue;

    // Accumulate sample covariance for this category (about the stacked CV)
    const int Nu = (int)idx.size();
    if(Nu<2) continue;
    // Prebuild stacked CV
    TVectorD vCV(Ntot); vCV.Zero();
    for(int fli=0; fli<NF; ++fli)
      for(int b=1;b<=nb;++b)
        vCV[fli*nb + (b-1)] = cvs[fli]->GetBinContent(b);

    // Loop universes
    for(int ui : idx){
      TVectorD v(Ntot); v.Zero();
      bool ok=true;
      for(int fli=0; fli<NF && ok; ++fli){
        const auto& flav = CFG::FLAVS[fli];
        auto itU = cats_by_flav[flav][cat].univ.find(ui);
        if(itU==cats_by_flav[flav][cat].univ.end()){ ok=false; break; }
        TH1D* u = itU->second;
        if(!same_binning(cvs[fli],u)){ ok=false; break; }
        for(int b=1;b<=nb;++b) v[fli*nb + (b-1)] = u->GetBinContent(b);
      }
      if(!ok) continue;
      // delta = (v - vCV)
      for(int i=0;i<Ntot;++i){
        const double di = v[i]-vCV[i];
        for(int j=0;j<Ntot;++j){
          const double dj = v[j]-vCV[j];
          C(i,j) += di*dj;
        }
      }
    }
    C *= (1.0/(Nu-1));
  }

  // 4) Fill joint pack
  JP.nb = nb;
  JP.cv_stacked.ResizeTo(Ntot);
  JP.C_joint.ResizeTo(Ntot,Ntot);
  // cv
  for(int fli=0; fli<NF; ++fli)
    for(int b=1;b<=nb;++b)
      JP.cv_stacked[fli*nb + (b-1)] = cvs[fli]->GetBinContent(b);
  // C
  JP.C_joint = C;

  for(auto* h : cvs) delete h;
  return true;
}

static void draw_cross_flavor_corr_side_by_side(const char* modeF, TFile& fF,
                                                const char* modeR, TFile& fR,
                                                int target_bin_mev){
  JointPack JF, JR;
  if(!build_joint_covariance(fF, modeF, target_bin_mev, JF)){ printf("[F] Could not build joint (%s)\n",modeF); return; }
  if(!build_joint_covariance(fR, modeR, target_bin_mev, JR)){ printf("[F] Could not build joint (%s)\n",modeR); return; }

  const int NF = (int)CFG::FLAVS.size();
  const int nb = JF.nb; // assume same nb in F and R (from same target rebin)
  const int Ntot = NF*nb;

  // Build correlation matrices
  TH2D HF(Form("Hcorr_%s",modeF),"Cross-flavor correlation;stack index i;stack index j",
          Ntot,0.5,Ntot+0.5,Ntot,0.5,Ntot+0.5);
  TH2D HR(Form("Hcorr_%s",modeR),"Cross-flavor correlation;stack index i;stack index j",
          Ntot,0.5,Ntot+0.5,Ntot,0.5,Ntot+0.5);

  auto fill_corr=[&](const TMatrixD& C, TH2D& H){
    for(int i=0;i<Ntot;++i){
      for(int j=0;j<Ntot;++j){
        const double sii = C(i,i), sjj = C(j,j);
        const double cij = C(i,j);
        const double rho = (sii>0 && sjj>0)? (cij/std::sqrt(sii*sjj)) : 0.0;
        H.SetBinContent(i+1,j+1, rho);
      }
    }
  };
  fill_corr(JF.C_joint, HF);
  fill_corr(JR.C_joint, HR);

  // Visual
  set_diverging_palette_11();
  HF.SetMinimum(-1.0); HF.SetMaximum( 1.0);
  HR.SetMinimum(-1.0); HR.SetMaximum( 1.0);

  TCanvas c("c_F_corr_dual","Cross-flavor correlation (FHC vs RHC)",1600,750);
  c.Divide(2,1);
  c.cd(1); gPad->SetRightMargin(0.14); HF.GetZaxis()->SetTitle("#rho"); HF.Draw("COLZ");
  // Flavor block separators
  for(int k=1;k<NF;++k){ int cut = k*nb; TLine* L1=new TLine(cut+0.5,0.5,cut+0.5,Ntot+0.5); L1->SetLineStyle(2); L1->Draw();
                         TLine* L2=new TLine(0.5,cut+0.5,Ntot+0.5,cut+0.5); L2->SetLineStyle(2); L2->Draw(); }

  c.cd(2); gPad->SetRightMargin(0.14); HR.GetZaxis()->SetTitle("#rho"); HR.Draw("COLZ");
  for(int k=1;k<NF;++k){ int cut = k*nb; TLine* R1=new TLine(cut+0.5,0.5,cut+0.5,Ntot+0.5); R1->SetLineStyle(2); R1->Draw();
                         TLine* R2=new TLine(0.5,cut+0.5,Ntot+0.5,cut+0.5); R2->SetLineStyle(2); R2->Draw(); }

  // Flavor labels along diagonal (minimal but handy)
  TLatex tx; tx.SetTextFont(42); tx.SetTextSize(0.03); tx.SetNDC(false);
  for(int fli=0; fli<NF; ++fli){
    const double cbin = fli*nb + nb*0.5 + 0.5;
    c.cd(1); tx.DrawLatex(cbin, Ntot+0.5 + 0.5, CFG::FLAVS[fli].c_str());
    c.cd(2); tx.DrawLatex(cbin, Ntot+0.5 + 0.5, CFG::FLAVS[fli].c_str());
  }

  c.Update();
  c.Print("ppfx_F_corr_crossFlavor_FHC_RHC.pdf");
  c.Print("ppfx_F_corr_crossFlavor_FHC_RHC.png");

  // Persist joint (for propagation)
  TFile fout("ppfx_hadprod_outputs.root","UPDATE");
  fout.cd();
  JF.C_joint.Write("C_ppfx_joint_FHC", TObject::kOverwrite);
  JF.cv_stacked.Write("CV_joint_FHC", TObject::kOverwrite);
  JR.C_joint.Write("C_ppfx_joint_RHC", TObject::kOverwrite);
  JR.cv_stacked.Write("CV_joint_RHC", TObject::kOverwrite);
  fout.Close();
}

// ------------------------------ DRIVER ----------------------------------------
void ppfx_hadprod_minimal(const char* fileFHC = CFG::FILE_FHC_DEF,
                          const char* fileRHC = CFG::FILE_RHC_DEF,
                          int target_bin_mev = CFG::TARGET_BIN_MEV)
{
  set_global_style();
  CFG::TARGET_BIN_MEV = target_bin_mev;

  TFile fFHC(fileFHC,"READ");
  TFile fRHC(fileRHC,"READ");
  if(fFHC.IsZombie()){ printf("Cannot open FHC file: %s\n",fileFHC); return; }
  if(fRHC.IsZombie()){ printf("Cannot open RHC file: %s\n",fileRHC); return; }

  // --- E) per-flavor ensemble band & covariance
  for(const auto& flav : CFG::FLAVS){
    draw_ppfx_ensemble_and_band("FHC", fFHC, flav);
    draw_ppfx_ensemble_and_band("RHC", fRHC, flav);
  }

  // --- F) cross-flavor fractional correlation (side-by-side FHC vs RHC)
  draw_cross_flavor_corr_side_by_side("FHC", fFHC, "RHC", fRHC, CFG::TARGET_BIN_MEV);

  fFHC.Close(); fRHC.Close();
}
