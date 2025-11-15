// root -l -b -q 'scripts/validate_uboone_flux.C++("/exp/uboone/data/users/bnayak/ppfx/flugg_studies/comparisons/dk2nu_fhc_ppfx_g4_10_4.root")'

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TKey.h"
#include "TStyle.h"
#include "TString.h"
#include <string>
#include <iostream>
#include <cstdio>

static std::string build_weight(TTree* t) {
  const bool has_wgt  = t->GetListOfBranches()->FindObject("wgt");
  const bool has_ppfx = t->GetListOfBranches()->FindObject("wgt_ppfx");
  std::string w = has_wgt ? "wgt" : "1";
  if (has_ppfx) w += "*wgt_ppfx";
  return w;
}

void validate_uboone_flux(const char* file="dk2nu_fhc_ppfx_g4_10_4.root") {
  gStyle->SetOptStat(0);

  TFile f(file, "READ");
  if (!f.IsOpen()) { std::cerr << "Could not open " << file << "\n"; return; }

  TTree* t = (TTree*)f.Get("outTree");
  if (!t) { std::cerr << "No outTree in " << file << "\n"; return; }

  TH1* hPOT = (TH1*)f.Get("POT");
  const double pot = hPOT ? hPOT->Integral() : 0.0;

  const std::string w = build_weight(t);

  // Preferred test: compare to pre-made histo (already POT-scaled to nominal POT)
  TH1D* href = (TH1D*)f.Get("numuFluxHisto");
  if (href) {
    // Build tree-based νμ spectrum with the correct total weight
    TH1D htr("htr","#nu_{#mu};E_{#nu} [GeV];arb.", href->GetNbinsX(),
             href->GetXaxis()->GetXmin(), href->GetXaxis()->GetXmax());
    htr.Sumw2();
    t->Draw("nuE>>htr", Form("(%s)*(ntype==14)", w.c_str()), "goff");

    // Scale the tree-based spectrum to the nominal POT used by producer hist
    if (pot > 0) htr.Scale(6e20 / pot);  // producer titles use 6e20 POT

    // Draw overlay + ratio
    TCanvas c("c_validate","Validate uB (NuMI off-axis)",900,800);
    c.Divide(1,2);
    c.cd(1); gPad->SetLogy();
    href->SetLineWidth(2);
    href->SetLineColor(kBlack);
    href->SetTitle("#nu_{#mu} at MicroBooNE;E_{#nu} [GeV];Flux / 6e20 POT / bin");
    href->Draw("hist");
    htr.SetLineColor(kRed+1);
    htr.SetLineWidth(2);
    htr.Draw("hist same");
    TLegend L(0.60,0.70,0.89,0.89);
    L.SetBorderSize(0); L.SetFillStyle(0); L.AddEntry(href,"pre-made (producer)","l");
    L.AddEntry(&htr,"tree-drawn (nuE, wgt*wgt_ppfx)","l");
    L.Draw();

    c.cd(2);
    TH1D ratio("ratio","tree / pre-made;E_{#nu} [GeV];ratio", href->GetNbinsX(),
               href->GetXaxis()->GetXmin(), href->GetXaxis()->GetXmax());
    ratio.Sumw2();
    ratio.Divide(&htr, href, 1.0, 1.0, "B");
    ratio.SetMinimum(0.8); ratio.SetMaximum(1.2);
    ratio.Draw("hist");
    c.SaveAs("validate_uB_overlay_ratio.pdf");

    std::cout << "[OK] Wrote validate_uB_overlay_ratio.pdf\n";
    return;
  }

  // Fallback tests if no pre-made histos available:
  // A) uB (nuE) vs MINERvA (E_minerva) shapes
  TH1D huB("huB","uB (NuMI off-axis);E_{#nu} [GeV];arb.",200,0,10);
  TH1D hMIN("hMIN","MINERvA (on-axis);E_{#nu} [GeV];arb.",200,0,10);
  huB.Sumw2(); hMIN.Sumw2();

  const bool has_Eminerva = t->GetListOfBranches()->FindObject("E_minerva");
  const bool has_wMIN     = t->GetListOfBranches()->FindObject("wgt_minerva");

  if (has_Eminerva && has_wMIN) {
    t->Draw("nuE>>huB"     , Form("(%s)*(ntype==14)", w.c_str()), "goff");
    t->Draw("E_minerva>>hMIN", "(wgt_minerva)*(ntype==14)", "goff");

    // Normalize shapes for a fair visual comparison
    if (huB.Integral() > 0)  huB.Scale(1.0/huB.Integral());
    if (hMIN.Integral() > 0) hMIN.Scale(1.0/hMIN.Integral());

    TCanvas c2("c_validate2","uB vs MINERvA shapes",900,600);
    c2.SetLogy();
    huB.SetLineColor(kRed+1); huB.SetLineWidth(2);
    hMIN.SetLineColor(kBlue+2); hMIN.SetLineWidth(2);
    huB.Draw("hist"); hMIN.Draw("hist same");
    TLegend L2(0.55,0.70,0.88,0.88);
    L2.SetBorderSize(0); L2.SetFillStyle(0);
    L2.AddEntry(&huB ,"uB (nuE, wgt*wgt_ppfx)","l");
    L2.AddEntry(&hMIN,"MINERvA (E_minerva, wgt_minerva)","l");
    L2.Draw();
    c2.SaveAs("validate_uB_vs_MINERvA.pdf");
    std::cout << "[OK] Wrote validate_uB_vs_MINERvA.pdf\n";
  }

  // B) Decay z distributions (uB vs MINERvA weighting)
  TH1D zuB ("zuB" ,"Decay z @ uB;z_{decay} [cm];weighted entries",240,-1000,30000);
  TH1D zMIN("zMIN","Decay z @ MINERvA;z_{decay} [cm];weighted entries",240,-1000,30000);
  zuB.Sumw2(); zMIN.Sumw2();
  t->Draw("nvz>>zuB" , Form("(%s)", w.c_str()), "goff");
  if (has_wMIN) t->Draw("nvz>>zMIN", "wgt_minerva", "goff");
  TCanvas c3("c_validate3","Decay z checks",900,600);
  c3.SetLogy();
  zuB.SetLineColor(kRed+1); zuB.SetLineWidth(2);
  zMIN.SetLineColor(kBlue+2); zMIN.SetLineWidth(2);
  zuB.Draw("hist"); if (has_wMIN) zMIN.Draw("hist same");
  TLegend L3(0.55,0.70,0.88,0.88);
  L3.SetBorderSize(0); L3.SetFillStyle(0);
  L3.AddEntry(&zuB ,"uB weighting","l");
  if (has_wMIN) L3.AddEntry(&zMIN,"MINERvA weighting","l");
  L3.Draw();
  c3.SaveAs("validate_decay_z.pdf");
  std::cout << "[OK] Wrote validate_decay_z.pdf\n";
}
