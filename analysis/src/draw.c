#include <math.h>
#include <fstream>
using namespace std;

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"

Int_t draw(TString SimROOTfilename="new_Simdata.root")
{
  //---读取用于重建图像的数据------------------------------
  TFile *f_sim = TFile::Open(SimROOTfilename, "read");
  TTree *Trec = (TTree *)f_sim->Get("Trec");
  Int_t i, j, k;

  Float_t x, y, z, xt, yt, zt, dthetaX, dthetaY, mp, nom_ang, ang, S, L, Ep=1., d;
  Double_t x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
  Trec->SetBranchAddress("x", &x);
  Trec->SetBranchAddress("y", &y);
  Trec->SetBranchAddress("z", &z);
  Trec->SetBranchAddress("x1", &x1);
  Trec->SetBranchAddress("y1", &y1);
  Trec->SetBranchAddress("z1", &z1);
  Trec->SetBranchAddress("x2", &x2);
  Trec->SetBranchAddress("y2", &y2);
  Trec->SetBranchAddress("z2", &z2);
  Trec->SetBranchAddress("x3", &x3);
  Trec->SetBranchAddress("y3", &y3);
  Trec->SetBranchAddress("z3", &z3);
  Trec->SetBranchAddress("x4", &x4);
  Trec->SetBranchAddress("y4", &y4);
  Trec->SetBranchAddress("z4", &z4);
  Trec->SetBranchAddress("ang", &ang);
  Trec->SetBranchAddress("d", &d);
  //--------------------------------------------------
  Double_t x1s, y1s, z1s, x2s, y2s, z2s, x3s, y3s, z3s, x4s, y4s, z4s;
  Double_t xs, ys, zs, angs, ds;
  Int_t nentries_sim = Trec->GetEntries(); 

  int count = 0;
  
  //Draw Y:X//
  TCanvas *c2 = new TCanvas("c2", "SIM Y:X", 800, 1600);
  TH2D *sim0 = new TH2D("sim0", "#0 SIM Y:X", 150, -150, 150, 150, -150, 150);
  TH2D *sim1 = new TH2D("sim1", "#1 SIM Y:X", 150, -150, 150, 150, -150, 150);
  TH2D *sim2 = new TH2D("sim2", "#2 SIM Y:X", 150, -150, 150, 150, -150, 150);
  TH2D *sim3 = new TH2D("sim3", "#3 SIM Y:X", 150, -150, 150, 150, -150, 150);

  c2->Divide(2, 2);
  c2->cd(1);
  Trec->Draw("y1:x1>>sim0");
  c2->cd(2);
  Trec->Draw("y2:x2>>sim1");
  c2->cd(3);
  Trec->Draw("y3:x3>>sim2");
  c2->cd(4);
  Trec->Draw("y4:x4>>sim3");

  return nentries_sim;
}
