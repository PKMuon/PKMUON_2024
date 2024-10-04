#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TVector3.h>
#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../include/Object.hh"

using namespace std;

TCanvas *draw_scatter(const char *infile = "../../build/root_file/muppoca_50GeV.root",
    const char *outfile = "../../build/root_file/muppoca_50GeV.pdf")
{
  // 输入
  TFile *file_in = TFile::Open(infile);
  TTree *tree = (TTree *)file_in->Get("tree");
  Double_t cos_theta, cos_theta_smeared;
  TClonesArray *Scatters = NULL;
  tree->SetBranchAddress("CosThetaEdep", &cos_theta);
  tree->SetBranchAddress("CosThetaSmeared", &cos_theta_smeared);
  tree->SetBranchAddress("Scatters", &Scatters);

  TCanvas *c1 = new TCanvas("c1", "c1");
  TH1D *theta_scatter[2] = {
    new TH1D("background_theta_scatter", "", 100, 0, 0.005),
    new TH1D("signal_theta_scatter", "", 100, 0, 0.005),
  };
  TH1D *theta[2] = {
    new TH1D("background_theta", "", 100, 0, 0.005),
    new TH1D("signal_theta", "", 100, 0, 0.005),
  };
  TH1D *theta_smeared[2] = {
    new TH1D("background_theta_smeared", "", 100, 0, 0.005),
    new TH1D("signal_theta_smeared", "", 100, 0, 0.005),
  };

  Long64_t nentry = tree->GetEntries();
  for(Long64_t i = 0; i < nentry; i++) {
    tree->GetEntry(i);
    if(i % 10000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2) << (i / (double)nentry) * 100 << "%" << endl;
    }
    bool is_signal = Scatters->GetEntries();
    Double_t cos_theta_scatter = 1.0;
    if(is_signal) {
      auto scatter = (Scatter *)Scatters->UncheckedAt(0);
      TVector3 p_mup[2];
      for(size_t j = 0, k = 0; j < 3; ++j) {
        if(scatter->Pid[j] == -13) p_mup[k++].SetXYZ(scatter->Px[j], scatter->Py[j], scatter->Pz[j]);
      }
      cos_theta_scatter = p_mup[0].Dot(p_mup[1]) / (p_mup[0].Mag() * p_mup[1].Mag());
    }
    theta[is_signal]->Fill(acos(cos_theta));
    theta_smeared[is_signal]->Fill(acos(cos_theta_smeared));
    theta_scatter[is_signal]->Fill(acos(cos_theta_scatter));
  }

  vector<TH1D *> hists = { theta[0], theta[1], theta_smeared[0], theta_smeared[1], theta_scatter[1] };
  for(size_t j = 0; j < hists.size(); ++j) {
    hists[j]->SetLineColor(2 + j);
    hists[j]->SetStats(kFALSE);
    if(j == 0) {
      hists[j]->SetXTitle("#theta");
      hists[j]->SetYTitle("Events");
      hists[j]->Draw("hist");
    } else {
      hists[j]->Draw("hist same");
    }
  }
  c1->BuildLegend(0.7, 0.75, 0.9, 0.9)->Draw();
  c1->SetLogy();
  c1->Draw();
  c1->SaveAs(outfile);

  for(size_t j = 0; j < hists.size(); ++j) {
    if(j == 0) {
      hists[j]->DrawNormalized("hist");
    } else {
      hists[j]->DrawNormalized("hist same");
    }
  }
  c1->BuildLegend(0.7, 0.75, 0.9, 0.9)->Draw();
  c1->SetLogy(kFALSE);
  c1->Draw();
  string filename = outfile;
  filename = filename.substr(0, filename.length() - 4) + "_norm.pdf";
  c1->SaveAs(filename.c_str());

  return c1;
}
