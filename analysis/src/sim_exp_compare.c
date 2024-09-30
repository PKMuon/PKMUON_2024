#include <math.h>
#include <fstream>
using namespace std;

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"

Int_t sim_exp_compare(TString SimROOTfilename="new_Simdata.root")
{
  //---读取用于重建图像的数据------------------------------
  Double_t x1s,y1s,z1s,x2s,y2s,z2s,x3s,y3s,z3s,x4s,y4s,z4s;
  Double_t xs,ys,zs,angs,ds;
  TFile *f_sim = TFile::Open (SimROOTfilename,"read");
  TTree *Trec = (TTree *)f_sim->Get("Trec");
  Int_t nentries_sim = Trec->GetEntries(); 
  Trec->SetBranchAddress("x",&xs); 	  
  Trec->SetBranchAddress("y",&ys); 
  Trec->SetBranchAddress("z",&zs); 
  Trec->SetBranchAddress("x1",&x1s);
  Trec->SetBranchAddress("y1",&y1s);
  Trec->SetBranchAddress("z1",&z1s);
  Trec->SetBranchAddress("x2",&x2s);
  Trec->SetBranchAddress("y2",&y2s);
  Trec->SetBranchAddress("z2",&z2s);
  Trec->SetBranchAddress("x3",&x3s);
  Trec->SetBranchAddress("y3",&y3s); 
  Trec->SetBranchAddress("z3",&z3s);   
  Trec->SetBranchAddress("x4",&x4s);
  Trec->SetBranchAddress("y4",&y4s);
  Trec->SetBranchAddress("z4",&z4s);
  Trec->SetBranchAddress("ang",&angs);
  Trec->SetBranchAddress("d",&ds);
  
  TCanvas *canvas = new TCanvas("canvas", "cos(ang)", 800, 600);
  TH1F *hcos_sim = new TH1F("hcos_sim", "distribution of scattering angle between muon and air(selected)", 100, 0., 1.);

  for(Int_t i=0; i<nentries_sim; i++){
    Trec->GetEntry(i);
    if(zs>200||zs<-200) continue;
    if(ds<-10||ds>10) continue;
    Double_t cosangle_sim=cos(angs);
    hcos_sim->Fill(cosangle_sim);
  }

  Double_t maxBin_sim = hcos_sim->GetMaximum();

  // 归一化直方图
  hcos_sim->Scale(1.0/maxBin_sim);
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  hcos_sim->SetMinimum(1e-6);
  hcos_sim->SetMaximum(2);
  hcos_sim->GetYaxis()->SetTitle("relative events");
  hcos_sim->GetXaxis()->SetTitle("cos #theta");
  hcos_sim->SetLineColor(kBlue);
  hcos_sim->Draw("HIST");

  TLegend *legend = new TLegend(0.15, 0.70, 0.5, 0.9);
  legend->SetX1NDC(0.15);
  legend->SetY1NDC(0.70);
  legend->SetX2NDC(0.50);
  legend->SetY2NDC(0.9);
  legend->AddEntry(hcos_sim, "Simulation", "l");
  legend->Draw();

  return nentries_sim;
}
