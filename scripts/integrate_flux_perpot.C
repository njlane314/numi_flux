#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TChainElement.h"
#include <iostream>

void integrate_flux_perpot(const char* inpat,
                           const char* mode="FHC",
                           double Emin=0.25, double Emax=5.0)
{
  TChain ch("outTree");
  ch.Add(inpat);

  // choose weight column present in the file
  bool has_ppfx = ch.GetListOfBranches()->FindObject("wgt_ppfx");
  const char* w = has_ppfx ? "wgt_ppfx" : "wgt";  // CV only

  // histogram to integrate over [Emin,Emax]
  TH1D h("h","",2000,0.,5.0);
  ch.Draw("nuE>>h",
          Form("(%g<=nuE && nuE<%g) * ((ntype==14)||(ntype==-14)) * %s",Emin,Emax,w),
          "goff");

  // sum contents in window (bins already restricted by the cut)
  double sum = h.Integral(1, h.GetNbinsX());  // integrated μ-flavour flux in [Emin,Emax]

  // total POT in these inputs (for reference; result below is already per POT)
  double totPOT = 0.;
  TList* files = ch.GetListOfFiles();
  for (int i=0; i<files->GetEntries(); ++i) {
    TFile f(((TChainElement*)files->At(i))->GetTitle());
    if (auto hp = (TH1*)f.Get("POT")) totPOT += hp->Integral();
  }

  std::cout << "[" << mode << "] ∫_"
            << Emin << "→" << Emax << " GeV φ_μ(E) dE = "
            << std::scientific << sum << "  (per POT)"
            << "   (POT in files: " << totPOT << ")\n";
}
