#include "dca.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"

void drawDCA(){
   TFile *f =  new TFile("mu1GeV_new.root");
   TTree *T1 =  (TTree*) f->Get("T1");

//Declaration of leaves types
   Double_t        GunEng;
   Double_t        totEdep;
   Double_t        cosTh;
   Double_t        X0;
   Double_t        Y0;
   Double_t        Z0;
   Double_t        px0;
   Double_t        py0;
   Double_t        pz0;
   vector<double>  *vX=0;
   vector<double>  *vY=0;
   vector<double>  *vZ=0;
   vector<double>  *vEdep=0;
   vector<int>     *vTrkID=0;
   vector<int>     *vReadoutTrkid=0;
   vector<double>  *vReadoutEdep=0;
   vector<double>  *vReadoutPosX=0;
   vector<double>  *vReadoutPosY=0;
   vector<double>  *vReadoutPosZ=0;
   vector<double>  *vPx=0;
   vector<double>  *vPy=0;
   vector<double>  *vPz=0;
   Double_t        rec_x[4];
   Double_t        rec_y[4];
   Double_t        rec_z[4];
   Double_t        rec_x_smear[4];
   Double_t        rec_y_smear[4];
   Double_t        angle;
   Double_t        angle_smear;
   Double_t        dca_x;
   Double_t        dca_y;
   Double_t        dca_z;
   Double_t        dca;
   Double_t        dca_x_smear;
   Double_t        dca_y_smear;
   Double_t        dca_z_smear;
   Double_t        dca_smear;
   Double_t        angle_gem12;

   // Set branch addresses.
   T1->SetBranchAddress("GunEng",&GunEng);
   T1->SetBranchAddress("totEdep",&totEdep);
   T1->SetBranchAddress("cosTh",&cosTh);
   T1->SetBranchAddress("X0",&X0);
   T1->SetBranchAddress("Y0",&Y0);
   T1->SetBranchAddress("Z0",&Z0);
   T1->SetBranchAddress("px0",&px0);
   T1->SetBranchAddress("py0",&py0);
   T1->SetBranchAddress("pz0",&pz0);
   T1->SetBranchAddress("vX",&vX);
   T1->SetBranchAddress("vY",&vY);
   T1->SetBranchAddress("vZ",&vZ);
   T1->SetBranchAddress("vEdep",&vEdep);
   T1->SetBranchAddress("vTrkID",&vTrkID);
   T1->SetBranchAddress("vReadoutTrkid",&vReadoutTrkid);
   T1->SetBranchAddress("vReadoutEdep",&vReadoutEdep);
   T1->SetBranchAddress("vReadoutPosX",&vReadoutPosX);
   T1->SetBranchAddress("vReadoutPosY",&vReadoutPosY);
   T1->SetBranchAddress("vReadoutPosZ",&vReadoutPosZ);
   T1->SetBranchAddress("vPx",&vPx);
   T1->SetBranchAddress("vPy",&vPy);
   T1->SetBranchAddress("vPz",&vPz);
   T1->SetBranchAddress("rec_x",rec_x);
   T1->SetBranchAddress("rec_y",rec_y);
   T1->SetBranchAddress("rec_z",rec_z);
   T1->SetBranchAddress("rec_x_smear",rec_x_smear);
   T1->SetBranchAddress("rec_y_smear",rec_y_smear);
   T1->SetBranchAddress("angle",&angle);
   T1->SetBranchAddress("angle_smear",&angle_smear);
   T1->SetBranchAddress("dca_x",&dca_x);
   T1->SetBranchAddress("dca_y",&dca_y);
   T1->SetBranchAddress("dca_z",&dca_z);
   T1->SetBranchAddress("dca",&dca);
   T1->SetBranchAddress("dca_x_smear",&dca_x_smear);
   T1->SetBranchAddress("dca_y_smear",&dca_y_smear);
   T1->SetBranchAddress("dca_z_smear",&dca_z_smear);
   T1->SetBranchAddress("dca_smear",&dca_smear);
   T1->SetBranchAddress("angle_gem12",&angle_gem12);

   Long64_t nentries = T1->GetEntries();
   //TH1D *hDCA_x = new TH1D("hDCA_x",";DCA_{x};Counts/100",200,-10000,10000);
   TH1D *hDCA_x = new TH1D("hDCA_x",";DCA_{x};Counts/100",50,-300,300);
   //TH1D *hDCA_xsmear = new TH1D("hDCA_xsmear",";DCA_{x};Counts/100",200,-10000,10000);
   TH1D *hDCA_xsmear = new TH1D("hDCA_xsmear",";DCA_{x};Counts/100",50,-300,300);
   TRandom3 rand;
   Double_t sigma = 0.1;
   Long64_t nbytes = 0;
  for (Long64_t i=0; i<nentries;i++) {
     nbytes += T1->GetEntry(i);
     struct Line l1 = {{rec_x[0],rec_y[0],rec_z[0]}, {rec_x[1],rec_y[1],rec_z[1]}};
     struct Line l2 = {{rec_x[2],rec_y[2],rec_z[2]}, {rec_x[3],rec_y[3],rec_z[3]}};
     line_perpendicular(l1,l2);
     if(fabs(pub_dca_z)<1000 && fabs(pub_dca_y)<1000)
       hDCA_x->Fill(pub_dca_x);
     struct Line sl1 = {{rec_x[0]+rand.Gaus(0,sigma),rec_y[0]+rand.Gaus(0,sigma),rec_z[0]}, {rec_x[1]+rand.Gaus(0,sigma),rec_y[1]+rand.Gaus(0,sigma),rec_z[1]}};
     struct Line sl2 = {{rec_x[2]+rand.Gaus(0,sigma),rec_y[2]+rand.Gaus(0,sigma),rec_z[2]}, {rec_x[3]+rand.Gaus(0,sigma),rec_y[3]+rand.Gaus(0,sigma),rec_z[3]}};
     line_perpendicular(sl1,sl2);
     if(fabs(pub_dca_z)<1000 && fabs(pub_dca_y)<1000)
       hDCA_xsmear->Fill(pub_dca_x);
  }
  hDCA_x->Draw();
  //hDCA_x->SetMaximum((hDCA_xsmear->GetMaximum()));
  hDCA_xsmear->Scale((hDCA_x->GetMaximum())/(hDCA_xsmear->GetMaximum()));
  hDCA_xsmear->Draw("same,hist");
  hDCA_xsmear->SetLineColor(kRed);
}
