#include "boot.h"
void SaveHist(){

setStyle();

TString var="costheta_smear";
TString xtitle="Cos#theta";
double xmin=-1;
double xmax= 1;
double xbins=100;

TCut cut = "fabs(poca_x_smear)<500 && fabs(poca_y_smear)<500 && fabs(poca_z_smear)<500";

TFile* f_const_0p005 = new TFile("../root/DMmuon_const_0p005.root","");
TFile* f_const_0p05 = new TFile("../root/DMmuon_const_0p05.root","");
TFile* f_const_0p5 = new TFile("../root/DMmuon_const_0p5.root","");
TFile* f_const_0p2 = new TFile("../root/DMmuon_const_0p2.root","");
TFile* f_const_0p1 = new TFile("../root/DMmuon_const_0p1.root","");
TFile* f_const_1 = new TFile("../root/DMmuon_const_1.root","");
TFile* f_const_10 = new TFile("../root/DMmuon_const_10.root","");
TFile* f_const_100 = new TFile("../root/DMmuon_const_100.root","");

TFile* f_maxwell_0p005 = new TFile("../root/DMmuon_maxwell_0p005.root","");
TFile* f_maxwell_0p05 = new TFile("../root/DMmuon_maxwell_0p05.root","");
TFile* f_maxwell_0p5 = new TFile("../root/DMmuon_maxwell_0p5.root","");
TFile* f_maxwell_0p2 = new TFile("../root/DMmuon_maxwell_0p2.root","");
TFile* f_maxwell_0p1 = new TFile("../root/DMmuon_maxwell_0p1.root","");
TFile* f_maxwell_1 = new TFile("../root/DMmuon_maxwell_1.root","");
TFile* f_maxwell_10 = new TFile("../root/DMmuon_maxwell_10.root","");
TFile* f_maxwell_100 = new TFile("../root/DMmuon_maxwell_100.root","");

TChain* t_air_new = new TChain("T1");
//t_air_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job_air/analysis/root/muCRY_air_1.root");
t_air_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job_air/analysis/root/muCRY_air_*.root");

TChain* t_vac_new = new TChain("T1");
//t_vac_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job_vacuum/analysis/root/muCRY_vac_1.root");
t_vac_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job_vacuum/analysis/root/muCRY_vac_*.root");

TTree* t_const_0p005 = (TTree*)f_const_0p005->Get("T1");
TTree* t_const_0p05 = (TTree*)f_const_0p05->Get("T1");
TTree* t_const_0p5 = (TTree*)f_const_0p5->Get("T1");
TTree* t_const_0p2 = (TTree*)f_const_0p2->Get("T1");
TTree* t_const_0p1 = (TTree*)f_const_0p1->Get("T1");
TTree* t_const_1 = (TTree*)f_const_1->Get("T1");
TTree* t_const_10 = (TTree*)f_const_10->Get("T1");
TTree* t_const_100 = (TTree*)f_const_100->Get("T1");

TTree* t_maxwell_0p005 = (TTree*)f_maxwell_0p005->Get("T1");
TTree* t_maxwell_0p05 = (TTree*)f_maxwell_0p05->Get("T1");
TTree* t_maxwell_0p5 = (TTree*)f_maxwell_0p5->Get("T1");
TTree* t_maxwell_0p2 = (TTree*)f_maxwell_0p2->Get("T1");
TTree* t_maxwell_0p1 = (TTree*)f_maxwell_0p1->Get("T1");
TTree* t_maxwell_1 = (TTree*)f_maxwell_1->Get("T1");
TTree* t_maxwell_10 = (TTree*)f_maxwell_10->Get("T1");
TTree* t_maxwell_100 = (TTree*)f_maxwell_100->Get("T1");

TH1D *h_const_0p005 = new TH1D("h_const_0p005","h_const_0p005",xbins,xmin,xmax);
TH1D *h_const_0p05 = new TH1D("h_const_0p05","h_const_0p05",xbins,xmin,xmax);
TH1D *h_const_0p5 = new TH1D("h_const_0p5","h_const_0p5",xbins,xmin,xmax);
TH1D *h_const_0p2 = new TH1D("h_const_0p2","h_const_0p2",xbins,xmin,xmax);
TH1D *h_const_0p1 = new TH1D("h_const_0p1","h_const_0p1",xbins,xmin,xmax);
TH1D *h_const_1 = new TH1D("h_const_1","h_const_1",xbins,xmin,xmax);
TH1D *h_const_10 = new TH1D("h_const_10","h_const_10",xbins,xmin,xmax);
TH1D *h_const_100 = new TH1D("h_const_100","h_const_100",xbins,xmin,xmax);

TH1D *h_maxwell_0p005 = new TH1D("h_maxwell_0p005","h_maxwell_0p005",xbins,xmin,xmax);
TH1D *h_maxwell_0p05 = new TH1D("h_maxwell_0p05","h_maxwell_0p05",xbins,xmin,xmax);
TH1D *h_maxwell_0p5 = new TH1D("h_maxwell_0p5","h_maxwell_0p5",xbins,xmin,xmax);
TH1D *h_maxwell_0p2 = new TH1D("h_maxwell_0p2","h_maxwell_0p2",xbins,xmin,xmax);
TH1D *h_maxwell_0p1 = new TH1D("h_maxwell_0p1","h_maxwell_0p1",xbins,xmin,xmax);
TH1D *h_maxwell_1 = new TH1D("h_maxwell_1","h_maxwell_1",xbins,xmin,xmax);
TH1D *h_maxwell_10 = new TH1D("h_maxwell_10","h_maxwell_10",xbins,xmin,xmax);
TH1D *h_maxwell_100 = new TH1D("h_maxwell_100","h_maxwell_100",xbins,xmin,xmax);

TH1D *h_air = new TH1D("h_air","h_air",xbins,xmin,xmax);
TH1D *h_vac = new TH1D("h_vac","h_vac",xbins,xmin,xmax);

TString a("Events/"); char b[20];  sprintf(b, "(%.5f",(xmax-xmin)/xbins); TString c(")");
TString ytitle = a + b + c;

//project
t_const_0p005->Project("h_const_0p005",var,cut);
t_const_0p05->Project("h_const_0p05",var,cut);
t_const_0p5->Project("h_const_0p5",var,cut);
t_const_0p2->Project("h_const_0p2",var,cut);
t_const_0p1->Project("h_const_0p1",var,cut);
t_const_1->Project("h_const_1",var,cut);
t_const_10->Project("h_const_10",var,cut);
t_const_100->Project("h_const_100",var,cut);

t_maxwell_0p005->Project("h_maxwell_0p005",var,cut);
t_maxwell_0p05->Project("h_maxwell_0p05",var,cut);
t_maxwell_0p5->Project("h_maxwell_0p5",var,cut);
t_maxwell_0p2->Project("h_maxwell_0p2",var,cut);
t_maxwell_0p1->Project("h_maxwell_0p1",var,cut);
t_maxwell_1->Project("h_maxwell_1",var,cut);
t_maxwell_10->Project("h_maxwell_10",var,cut);
t_maxwell_100->Project("h_maxwell_100",var,cut);

//scale
h_const_0p005->Scale(0.1);
h_const_0p05->Scale(0.001);
h_const_0p5->Scale(0.001);
h_const_0p2->Scale(0.001);
h_const_0p1->Scale(0.001);
h_const_1->Scale(0.001);
h_const_10->Scale(0.001);
h_const_100->Scale(0.001);

h_maxwell_0p005->Scale(0.1);
h_maxwell_0p05->Scale(0.001);
h_maxwell_0p5->Scale(0.001);
h_maxwell_0p2->Scale(0.001);
h_maxwell_0p1->Scale(0.001);
h_maxwell_1->Scale(0.001);
h_maxwell_10->Scale(0.001);
h_maxwell_100->Scale(0.001);

t_air_new->Project("h_air",var,cut);
t_vac_new->Project("h_vac",var,cut);

h_const_0p005->SetLineColor(2);
h_const_0p05->SetLineColor(3);
h_const_0p5->SetLineColor(4);
h_const_0p2->SetLineColor(5);
h_const_0p1->SetLineColor(6);
h_const_1->SetLineColor(7);
h_const_10->SetLineColor(8);
h_const_100->SetLineColor(9);

h_maxwell_0p005->SetLineColor(2);
h_maxwell_0p05->SetLineColor(3);
h_maxwell_0p5->SetLineColor(4);
h_maxwell_0p2->SetLineColor(5);
h_maxwell_0p1->SetLineColor(6);
h_maxwell_1->SetLineColor(7);
h_maxwell_10->SetLineColor(8);
h_maxwell_100->SetLineColor(9);

h_air->SetLineColor(1);
h_vac->SetLineColor(1);

cout<<"Constant"<<endl;
cout<<"MDM=5MeV, Nevent = "<<h_const_0p005->Integral()<<", eff = "<<h_const_0p005->Integral()/1000000.<<endl;
cout<<"MDM=50MeV, Nevent = "<<h_const_0p05->Integral()<<", eff = "<<h_const_0p05->Integral()/1000000.<<endl;
cout<<"MDM=100MeV, Nevent = "<<h_const_0p1->Integral()<<", eff = "<<h_const_0p1->Integral()/1000000.<<endl;
cout<<"MDM=200MeV, Nevent = "<<h_const_0p2->Integral()<<", eff = "<<h_const_0p2->Integral()/1000000.<<endl;
cout<<"MDM=500MeV, Nevent = "<<h_const_0p5->Integral()<<", eff = "<<h_const_0p5->Integral()/1000000.<<endl;
cout<<"MDM=1GeV, Nevent = "<<h_const_1->Integral()<<", eff = "<<h_const_1->Integral()/1000000.<<endl;
cout<<"MDM=10GeV, Nevent = "<<h_const_10->Integral()<<", eff = "<<h_const_10->Integral()/1000000.<<endl;
cout<<"MDM=100GeV, Nevent = "<<h_const_100->Integral()<<", eff = "<<h_const_100->Integral()/1000000.<<endl;

cout<<"Maxwell"<<endl;
cout<<"MDM=5MeV, Nevent = "<<h_maxwell_0p005->Integral()<<", eff = "<<h_maxwell_0p005->Integral()/1000000.<<endl;
cout<<"MDM=50MeV, Nevent = "<<h_maxwell_0p05->Integral()<<", eff = "<<h_maxwell_0p05->Integral()/1000000.<<endl;
cout<<"MDM=100MeV, Nevent = "<<h_maxwell_0p1->Integral()<<", eff = "<<h_maxwell_0p1->Integral()/1000000.<<endl;
cout<<"MDM=200MeV, Nevent = "<<h_maxwell_0p2->Integral()<<", eff = "<<h_maxwell_0p2->Integral()/1000000.<<endl;
cout<<"MDM=500MeV, Nevent = "<<h_maxwell_0p5->Integral()<<", eff = "<<h_maxwell_0p5->Integral()/1000000.<<endl;
cout<<"MDM=1GeV, Nevent = "<<h_maxwell_1->Integral()<<", eff = "<<h_maxwell_1->Integral()/1000000.<<endl;
cout<<"MDM=10GeV, Nevent = "<<h_maxwell_10->Integral()<<", eff = "<<h_maxwell_10->Integral()/1000000.<<endl;
cout<<"MDM=100GeV, Nevent = "<<h_maxwell_100->Integral()<<", eff = "<<h_maxwell_100->Integral()/1000000.<<endl;

cout<<"Air, Nevent = "<<h_air->Integral()<<endl;
cout<<"Vac, Nevent = "<<h_vac->Integral()<<endl;

/*
h_const_0p005->Scale((h_air->Integral())/(h_const_0p005->Integral()));
h_const_0p05->Scale((h_air->Integral())/(h_const_0p05->Integral()));
h_const_0p5->Scale((h_air->Integral())/(h_const_0p5->Integral()));
//h_const_0p4->Scale((h_air->Integral())/(h_const_0p4->Integral()));
//h_const_0p3->Scale((h_air->Integral())/(h_const_0p3->Integral()));
h_const_0p2->Scale((h_air->Integral())/(h_const_0p2->Integral()));
h_const_0p1->Scale((h_air->Integral())/(h_const_0p1->Integral()));
h_const_1->Scale((h_air->Integral())/(h_const_1->Integral()));
h_const_10->Scale((h_air->Integral())/(h_const_10->Integral()));
h_const_100->Scale((h_air->Integral())/(h_const_100->Integral()));
//h_const_1000->Scale((h_air->Integral())/(h_const_1000->Integral()));
*/

//TFile* fout = new TFile("Histtest.root","recreate");
TFile* fout = new TFile("Hist.root","recreate");
fout->cd();
h_const_0p005->Write();
h_const_0p05->Write();
h_const_0p5->Write();
h_const_0p2->Write();
h_const_0p1->Write();
h_const_1->Write();
h_const_10->Write();
h_const_100->Write();

h_maxwell_0p005->Write();
h_maxwell_0p05->Write();
h_maxwell_0p5->Write();
h_maxwell_0p2->Write();
h_maxwell_0p1->Write();
h_maxwell_1->Write();
h_maxwell_10->Write();
h_maxwell_100->Write();

h_air->Write();
h_vac->Write();

delete h_const_0p005;
delete h_const_0p05;
delete h_const_0p5;
delete h_const_0p2;
delete h_const_0p1;
delete h_const_1;
delete h_const_10;
delete h_const_100;

delete h_maxwell_0p005;
delete h_maxwell_0p05;
delete h_maxwell_0p5;
delete h_maxwell_0p2;
delete h_maxwell_0p1;
delete h_maxwell_1;
delete h_maxwell_10;
delete h_maxwell_100;

delete h_air;
delete h_vac;

}
