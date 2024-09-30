#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "../../include/Object.hh"

using namespace std;

TCanvas *draw_sim(const char *infile = "../../build/root_file/CryMuPoca.root",
    const char *outfile = "../../build/root_file/CryMuPoca.pdf", const char *selection = NULL,
    Double_t min_theta = 0.001)
{
  const Int_t Nx = 40, Ny = 40, Nz = 40;  // 控制像元大小
  const Double_t Xdown = -200, Xup = 200, Ydown = -200, Yup = 200, Zdown = -360, Zup = 360;

  // 输入
  TFile *file_in = TFile::Open(infile);
  TTree *tree = (TTree *)file_in->Get("tree");
  TTreeFormula *formula = selection ? new TTreeFormula("tree_cut", selection, tree) : NULL;
  Double_t d, x, y, z, cos_theta;
  tree->SetBranchAddress("DPoCAEdep", &d);
  tree->SetBranchAddress("XPoCAEdep", &x);
  tree->SetBranchAddress("YPoCAEdep", &y);
  tree->SetBranchAddress("ZPoCAEdep", &z);
  tree->SetBranchAddress("CosThetaEdep", &cos_theta);

  Double_t vox = (Xup - Xdown) / Nx, voy = (Yup - Ydown) / Ny, voz = (Zup - Zdown) / Nz;
  Double_t sig[Nx][Ny][Nz] = { 0 };
  Int_t count[Nx][Ny][Nz] = { 0 };

  Long64_t nentry = tree->GetEntries();
  for(Long64_t i = 0; i < nentry; i++) {
    tree->GetEntry(i);
    if(i % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2) << (i / (double)nentry) * 100 << "%" << endl;
    }
    if(acos(cos_theta) < min_theta) continue;
    if(formula && !formula->EvalInstance()) continue;

    // 给 PoCA 点以及路径上的像元赋值
    Int_t u = (x - Xdown) / vox;
    Int_t v = (y - Ydown) / voy;
    Int_t w = (z - Zdown) / voz;
    if(u < 0 || u >= Nx || v < 0 || v >= Ny || w < 0 || w >= Nz) continue;
    sig[u][v][w] += acos(cos_theta) * 1000;  // mrad
    count[u][v][w]++;
  }

  TH2D *hxy = new TH2D("hxy", "xy", Nx, Xdown, Xup, Ny, Ydown, Yup);
  TH2D *hyz = new TH2D("hyz", "yz", Ny, Ydown, Yup, Nz, Zdown, Zup);
  TH2D *hxz = new TH2D("hxz", "xz", Nx, Xdown, Xup, Nz, Zdown, Zup);
  TH3D *hxyz = new TH3D("hxyz", "xyz", Nx, Xdown, Xup, Ny, Ydown, Yup, Nz, Zdown, Zup);

  // 计算每个像元的平均散射强度
  for(Int_t i = 0; i < Nx; i++)
    for(Int_t j = 0; j < Ny; j++)
      for(Int_t k = 0; k < Nz; k++) {
        if(count[i][j][k]) sig[i][j][k] = sig[i][j][k] / count[i][j][k];

        Double_t xt = i * vox + Xdown + vox / 2;
        Double_t yt = j * voy + Ydown + voy / 2;
        Double_t zt = k * voz + Zdown + voz / 2;
        hxy->Fill(xt, yt, count[i][j][k]);
        hyz->Fill(yt, zt, count[i][j][k]);
        hxz->Fill(xt, zt, count[i][j][k]);
        hxyz->Fill(xt, yt, zt, count[i][j][k]);

        // [TODO] Plot sig.
      }

  TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1600);
  c1->Divide(2, 2);

  c1->cd(1);
  hxy->Draw("colz");
  hxy->SetTitle(("XY Events (ang >= " + to_string(min_theta) + " rad)").c_str());
  hxy->SetXTitle("x [mm]");
  hxy->SetYTitle("y [mm]");
  hxy->SetStats(0);

  c1->cd(2);
  hyz->Draw("colz");
  hyz->SetTitle(("YZ Events (ang >= " + to_string(min_theta) + " rad)").c_str());
  hyz->SetXTitle("y [mm]");
  hyz->SetYTitle("z [mm]");
  hyz->SetStats(0);

  c1->cd(3);
  hxz->Draw("colz");
  hxz->SetTitle(("XZ Events (ang >= " + to_string(min_theta) + " rad)").c_str());
  hxz->SetXTitle("x [mm]");
  hxz->SetYTitle("z [mm]");
  hxz->SetStats(0);

  c1->cd(4);
  hxyz->Draw("colz");
  hxyz->SetTitle(("XYZ Events (ang >= " + to_string(min_theta) + " rad)").c_str());
  hxyz->SetXTitle("x [mm]");
  hxyz->SetYTitle("y [mm]");
  hxyz->SetZTitle("z [mm]");
  hxyz->SetStats(0);

  c1->SaveAs(outfile);
  delete formula;
  return c1;
}
