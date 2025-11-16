#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <cstdio>

namespace CFG {
  const char* FILE_FHC_DEF = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_FHC.root";
  const char* FILE_RHC_DEF = "/exp/uboone/data/users/bnayak/ppfx/flugg_studies/NuMIFlux_dk2nu_RHC.root";
  const char* ENV_FHC = "PPFX_FHC_FILE";
  const char* ENV_RHC = "PPFX_RHC_FILE";
  static const std::vector<std::string> FLAVS = {"numu","numubar","nue","nuebar"};
}

static std::string pick_file_path(const char* env_name, const char* def_path, const char* label){
  const char* env_val = gSystem ? gSystem->Getenv(env_name) : nullptr;
  if(env_val && env_val[0]){
    printf("Using %s file override from %s: %s\n", label, env_name, env_val);
    return env_val;
  }
  printf("Using default %s file: %s\n", label, def_path);
  return def_path;
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

static TH1D* get_cv_energy(TFile& f, const std::string& flav){
  TH1D* h = (TH1D*)f.Get((flav + "/Detsmear/" + flav + "_CV_AV_TPC_5MeV_bin").c_str());
  if(!h) h = (TH1D*)f.Get((flav + "/Detsmear/" + flav + "_CV_AV_TPC").c_str());
  if(!h) return nullptr;
  TH1D* c = (TH1D*)h->Clone((std::string("cv_")+flav).c_str());
  c->SetDirectory(0);
  return c;
}

struct PPFXCat { std::map<int,TH1D*> univ; };

static std::map<std::string,PPFXCat> scan_ppfx_categories(TFile& f, const std::string& flav){
  std::map<std::string,PPFXCat> out;
  auto* d = f.GetDirectory( (flav + "/Multisims").c_str() );
  if(!d) return out;
  TIter it(d->GetListOfKeys());
  while(auto* k=(TKey*)it()){
    if(std::string(k->GetClassName())!="TH1D") continue;
    std::string name = k->GetName();
    if(name.find("_2D")!=std::string::npos) continue;
    if(name.find(flav+"_ppfx_")!=0) continue;
    size_t posUni = name.find("_Uni_");
    if(posUni==std::string::npos) continue;
    std::string cat = name.substr( (flav+"_ppfx_").size(), posUni-(flav+"_ppfx_").size() );
    size_t posIdx = posUni + 5;
    size_t posEnd = name.find("_", posIdx);
    int idx = std::stoi(name.substr(posIdx, posEnd-posIdx));
    auto* h=(TH1D*)d->Get(name.c_str());
    if(!h) continue;
    auto* c=(TH1D*)h->Clone((std::string("cl_")+name).c_str());
    c->SetDirectory(0);
    out[cat].univ[idx] = c;
  }
  return out;
}

struct JointPack {
  TVectorD cv_stacked;
  TMatrixD C_joint;
  int nb=0;
};

static bool build_joint_covariance(TFile& f, const char* mode, JointPack& JP){
  std::vector<TH1D*> cvs; cvs.reserve(CFG::FLAVS.size());
  for(const auto& flav : CFG::FLAVS){
    TH1D* cv = get_cv_energy(f,flav);
    if(!cv){ printf("[%s] Missing CV for %s\n",mode,flav.c_str()); for(auto* h:cvs) delete h; return false; }
    cvs.push_back(cv);
  }
  for(size_t i=1;i<cvs.size();++i){
    if(!same_binning(cvs[0],cvs[i])){
      printf("[%s] CV binnings differ across flavors\n",mode);
      for(auto* h:cvs) delete h;
      return false;
    }
  }
  const int nb = cvs[0]->GetNbinsX();
  const int NF = (int)CFG::FLAVS.size();
  const int Ntot = nb*NF;

  std::map<std::string, std::map<std::string,PPFXCat>> cats_by_flav;
  std::set<std::string> all_categories;
  for(const auto& flav : CFG::FLAVS){
    cats_by_flav[flav] = scan_ppfx_categories(f,flav);
    if(cats_by_flav[flav].empty()){
      printf("[%s] No PPFX universes for %s\n",mode,flav.c_str());
      for(auto* h:cvs) delete h;
      return false;
    }
    for(const auto& kv : cats_by_flav[flav]) all_categories.insert(kv.first);
  }

  TVectorD vCV(Ntot); vCV.Zero();
  for(int fli=0; fli<NF; ++fli)
    for(int b=1;b<=nb;++b)
      vCV[fli*nb + (b-1)] = cvs[fli]->GetBinContent(b);

  TMatrixD C(Ntot,Ntot); C.Zero();
  int cats_considered=0, cats_used=0, total_univ_used=0;

  for(const auto& cat : all_categories){
    ++cats_considered;

    std::set<int> idx_union;
    for(const auto& flav : CFG::FLAVS){
      auto itf = cats_by_flav[flav].find(cat);
      if(itf==cats_by_flav[flav].end()) continue;
      for(const auto& ij : itf->second.univ) idx_union.insert(ij.first);
    }
    if(idx_union.size()<2) continue;

    TMatrixD Cc(Ntot,Ntot); Cc.Zero();
    int Nu = 0;
    for(int ui : idx_union){
      ++Nu;
      TVectorD v(Ntot); v.Zero();
      bool ok=true, has_variation=false;
      for(int fli=0; fli<NF && ok; ++fli){
        const auto& flav = CFG::FLAVS[fli];
        auto itf = cats_by_flav[flav].find(cat);
        TH1D* src = nullptr;
        if(itf!=cats_by_flav[flav].end()){
          auto itu = itf->second.univ.find(ui);
          if(itu!=itf->second.univ.end()) src = itu->second;
        }
        if(src){
          if(!same_binning(cvs[fli],src)){ ok=false; break; }
          for(int b=1;b<=nb;++b) v[fli*nb + (b-1)] = src->GetBinContent(b);
          has_variation = true;
        } else {
          for(int b=1;b<=nb;++b) v[fli*nb + (b-1)] = cvs[fli]->GetBinContent(b);
        }
      }
      if(!ok || !has_variation) { Nu--; continue; }
      for(int i=0;i<Ntot;++i){
        const double di = v[i]-vCV[i];
        for(int j=0;j<Ntot;++j){
          const double dj = v[j]-vCV[j];
          Cc(i,j) += di*dj;
        }
      }
    }
    if(Nu>1){
      Cc *= (1.0/(Nu-1));
      C += Cc;
      ++cats_used;
      total_univ_used += Nu;
    }
  }

  printf("[%s] categories: considered=%d, used=%d, universes=%d\n", mode, cats_considered, cats_used, total_univ_used);

  JP.nb = nb;
  JP.cv_stacked.ResizeTo(Ntot);
  JP.C_joint.ResizeTo(Ntot,Ntot);
  for(int fli=0; fli<NF; ++fli)
    for(int b=1;b<=nb;++b)
      JP.cv_stacked[fli*nb + (b-1)] = cvs[fli]->GetBinContent(b);
  JP.C_joint = C;

  for(auto* h : cvs) delete h;
  return (cats_used>0);
}

static TH2D* make_corr_hist(const TMatrixD& C, int nb, const char* name, const char* title){
  const int N = C.GetNrows();
  TH2D* H = new TH2D(name, title, N, 0.5, N+0.5, N, 0.5, N+0.5);
  for(int i=0;i<N;++i){
    for(int j=0;j<N;++j){
      const double sii = C(i,i), sjj = C(j,j);
      const double cij = C(i,j);
      const double rho = (sii>0 && sjj>0)? (cij/std::sqrt(sii*sjj)) : 0.0;
      H->SetBinContent(i+1,j+1, rho);
    }
  }
  H->GetZaxis()->SetTitle("#rho");
  H->SetMinimum(-1.0);
  H->SetMaximum( 1.0);
  return H;
}

void ppfx_hadprod_minimal(){
  gStyle->SetOptStat(0);
  const std::string fhc_path = pick_file_path(CFG::ENV_FHC, CFG::FILE_FHC_DEF, "FHC");
  const std::string rhc_path = pick_file_path(CFG::ENV_RHC, CFG::FILE_RHC_DEF, "RHC");
  TFile fFHC(fhc_path.c_str(),"READ");
  TFile fRHC(rhc_path.c_str(),"READ");
  if(fFHC.IsZombie()){ printf("Cannot open FHC file: %s\n",fhc_path.c_str()); return; }
  if(fRHC.IsZombie()){ printf("Cannot open RHC file: %s\n",rhc_path.c_str()); return; }
  JointPack JF, JR;
  if(!build_joint_covariance(fFHC, "FHC", JF)) { printf("[FHC] Joint build failed.\n"); return; }
  if(!build_joint_covariance(fRHC, "RHC", JR)) { printf("[RHC] Joint build failed.\n"); return; }
  TH2D* HF = make_corr_hist(JF.C_joint, JF.nb, "Hcorr_FHC", "PPFX cross-flavor correlation (FHC);stack index i;stack index j");
  TH2D* HR = make_corr_hist(JR.C_joint, JR.nb, "Hcorr_RHC", "PPFX cross-flavor correlation (RHC);stack index i;stack index j");
  TCanvas c("c_corr_dual","PPFX cross-flavor correlation — FHC vs RHC",1400,700);
  c.Divide(2,1);
  c.cd(1); gPad->SetRightMargin(0.13); HF->Draw("COLZ");
  c.cd(2); gPad->SetRightMargin(0.13); HR->Draw("COLZ");
  c.Update();
  c.Print("ppfx_crossFlavor_corr_FHC_RHC.png");
}
