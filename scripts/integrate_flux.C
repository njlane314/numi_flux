#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
using namespace ROOT;
using namespace ROOT::RDF;

void integrate_mu_flux(std::string pattern, double Emin=0.25, double Emax=5.0, int nthreads=1, bool use_ppfx_cv=true) {
  if (nthreads>1) ROOT::EnableImplicitMT(nthreads);
  TChain ch("outTree"); ch.Add(pattern.c_str());
  RDataFrame df(ch);
  auto sel = [Emin,Emax](int t, float E){ return (t==14 || t==-14) && E>=Emin && E<Emax; };
  auto d0  = df.Filter(sel, {"ntype","nuE"});
  auto d   = use_ppfx_cv ? d0.Define("w", [](float a,float b){ return (double)a*(double)b; }, {"wgt","wgt_ppfx"})
                         : d0.Define("w", [](float a){ return (double)a; }, {"wgt"});
  auto sum = d.Sum<double>("w");
  auto s14 = d.Filter([](int t){return t==14;},{"ntype"}).Sum<double>("w");
  auto s14b= d.Filter([](int t){return t==-14;},{"ntype"}).Sum<double>("w");
  std::cout.setf(std::ios::scientific); std::cout.precision(6);
  std::cout<<"E=["<<Emin<<","<<Emax<<"] GeV, per-POT mu-flavour integrated flux"<<(use_ppfx_cv?" (CV)":" (raw)")<<": "<<*sum<<" nu/cm^2/POT\\n";
  std::cout<<"  numu: "<<*s14<<"  numubar: "<<*s14b<<"\\n";
}