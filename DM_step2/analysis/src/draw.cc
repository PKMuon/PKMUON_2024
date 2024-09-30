#include "boot.h"
void draw(){

setStyle();
TCanvas *can = new TCanvas("can","can",800,600);setPad();
can->SetLogy();

TFile* f_const_0p005 = new TFile("../root/DMmuon_const_0p005.root","");
TFile* f_const_0p05 = new TFile("../root/DMmuon_const_0p05.root","");
TFile* f_const_0p5 = new TFile("../root/DMmuon_const_0p5.root","");
//TFile* f_const_0p4 = new TFile("../root/DMmuon_const_0p4.root","");
//TFile* f_const_0p3 = new TFile("../root/DMmuon_const_0p3.root","");
TFile* f_const_0p2 = new TFile("../root/DMmuon_const_0p2.root","");
TFile* f_const_0p1 = new TFile("../root/DMmuon_const_0p1.root","");
TFile* f_const_1 = new TFile("../root/DMmuon_const_1.root","");
TFile* f_const_10 = new TFile("../root/DMmuon_const_10.root","");
TFile* f_const_100 = new TFile("../root/DMmuon_const_100.root","");
//TFile* f_const_1000 = new TFile("../root/DMmuon_const_1000.root","");

/*
TFile* f_const_0p05 = new TFile("../root_0905/DMmuon_const_0p05.root","");
TFile* f_const_0p5 = new TFile("../root_0905/DMmuon_const_0p5.root","");
TFile* f_const_1 = new TFile("../root_0905/DMmuon_const_1.root","");
TFile* f_const_10 = new TFile("../root_0905/DMmuon_const_10.root","");
TFile* f_const_100 = new TFile("../root_0905/DMmuon_const_100.root","");
*/

//TFile* f_vac = new TFile("/Users/yuxd/Desktop/PKMUON_G4sim/analysis/root/muCRY_vac.root","");
TFile* f_air = new TFile("/Users/yuxd/Desktop/PKMUON_G4sim/analysis/root/muCRY_air.root","");

TTree* t_const_0p005 = (TTree*)f_const_0p005->Get("T1");
TTree* t_const_0p05 = (TTree*)f_const_0p05->Get("T1");
TTree* t_const_0p5 = (TTree*)f_const_0p5->Get("T1");
//TTree* t_const_0p4 = (TTree*)f_const_0p4->Get("T1");
//TTree* t_const_0p3 = (TTree*)f_const_0p3->Get("T1");
TTree* t_const_0p2 = (TTree*)f_const_0p2->Get("T1");
TTree* t_const_0p1 = (TTree*)f_const_0p1->Get("T1");
TTree* t_const_1 = (TTree*)f_const_1->Get("T1");
TTree* t_const_10 = (TTree*)f_const_10->Get("T1");
TTree* t_const_100 = (TTree*)f_const_100->Get("T1");
//TTree* t_const_1000 = (TTree*)f_const_1000->Get("T1");
//TTree* t_vac = (TTree*)f_vac->Get("T1");
TTree* t_air = (TTree*)f_air->Get("T1");

double xmin = -1, xmax = 1, xbins = 100;
//double xmin = 0.8, xmax = 1, xbins = 100;
TH1F *h_const_0p005 = new TH1F("h_const_0p005","h_const_0p005",xbins,xmin,xmax);
TH1F *h_const_0p05 = new TH1F("h_const_0p05","h_const_0p05",xbins,xmin,xmax);
TH1F *h_const_0p5 = new TH1F("h_const_0p5","h_const_0p5",xbins,xmin,xmax);
//TH1F *h_const_0p4 = new TH1F("h_const_0p4","h_const_0p4",xbins,xmin,xmax);
//TH1F *h_const_0p3 = new TH1F("h_const_0p3","h_const_0p3",xbins,xmin,xmax);
TH1F *h_const_0p2 = new TH1F("h_const_0p2","h_const_0p2",xbins,xmin,xmax);
TH1F *h_const_0p1 = new TH1F("h_const_0p1","h_const_0p1",xbins,xmin,xmax);
TH1F *h_const_1 = new TH1F("h_const_1","h_const_1",xbins,xmin,xmax);
TH1F *h_const_10 = new TH1F("h_const_10","h_const_10",xbins,xmin,xmax);
TH1F *h_const_100 = new TH1F("h_const_100","h_const_100",xbins,xmin,xmax);
//TH1F *h_const_1000 = new TH1F("h_const_1000","h_const_1000",xbins,xmin,xmax);
//TH1F *h_vac = new TH1F("h_vac","h_vac",xbins,xmin,xmax);
TH1F *h_air = new TH1F("h_air","h_air",xbins,xmin,xmax);

TString a("Events/"); char b[20];  sprintf(b, "(%.5f",(xmax-xmin)/xbins); TString c(")");
TString ytitle = a + b + c;

//TString var = "costheta";
TString var = "costheta_smear";
TString xtitle = "Cos#theta";
//TString plot_name = "../figure/costheta.eps";
TString plot_name = "../figure/costheta_smear.eps";
TCut cut = "1";
//TCut cut = "costheta_smear<0.8";
//project
t_const_0p005->Project("h_const_0p005",var,cut);
t_const_0p05->Project("h_const_0p05",var,cut);
t_const_0p5->Project("h_const_0p5",var,cut);
//t_const_0p4->Project("h_const_0p4",var,cut);
//t_const_0p3->Project("h_const_0p3",var,cut);
t_const_0p2->Project("h_const_0p2",var,cut);
t_const_0p1->Project("h_const_0p1",var,cut);
t_const_1->Project("h_const_1",var,cut);
t_const_10->Project("h_const_10",var,cut);
t_const_100->Project("h_const_100",var,cut);
//t_const_1000->Project("h_const_1000",var,cut);
//t_vac->Project("h_vac",var,cut);
t_air->Project("h_air",var,cut);

h_const_0p005->SetLineColor(2);
h_const_0p05->SetLineColor(3);
h_const_0p5->SetLineColor(4);
//h_const_0p4->SetLineColor(5);
//h_const_0p3->SetLineColor(6);
h_const_0p2->SetLineColor(7);
h_const_0p1->SetLineColor(8);
h_const_1->SetLineColor(9);
h_const_10->SetLineColor(11);
h_const_100->SetLineColor(12);
//h_const_1000->SetLineColor(13);
//h_vac->SetLineColor(14);
h_air->SetLineColor(14);

cout<<"MDM=5MeV, Nevent = "<<h_const_0p005->Integral()<<", eff = "<<h_const_0p005->Integral()/1000000.<<endl;
cout<<"MDM=50MeV, Nevent = "<<h_const_0p05->Integral()<<", eff = "<<h_const_0p05->Integral()/1000000.<<endl;
cout<<"MDM=100MeV, Nevent = "<<h_const_0p1->Integral()<<", eff = "<<h_const_0p1->Integral()/1000000.<<endl;
cout<<"MDM=200MeV, Nevent = "<<h_const_0p2->Integral()<<", eff = "<<h_const_0p2->Integral()/1000000.<<endl;
//cout<<"MDM=300MeV, Nevent = "<<h_const_0p3->Integral()<<", eff = "<<h_const_0p3->Integral()/1000000.<<endl;
//cout<<"MDM=400MeV, Nevent = "<<h_const_0p4->Integral()<<", eff = "<<h_const_0p4->Integral()/1000000.<<endl;
cout<<"MDM=500MeV, Nevent = "<<h_const_0p5->Integral()<<", eff = "<<h_const_0p5->Integral()/1000000.<<endl;
cout<<"MDM=1GeV, Nevent = "<<h_const_1->Integral()<<", eff = "<<h_const_1->Integral()/1000000.<<endl;
cout<<"MDM=10GeV, Nevent = "<<h_const_10->Integral()<<", eff = "<<h_const_10->Integral()/1000000.<<endl;
cout<<"MDM=100GeV, Nevent = "<<h_const_100->Integral()<<", eff = "<<h_const_100->Integral()/1000000.<<endl;
//cout<<"MDM=1TeV, Nevent = "<<h_const_1000->Integral()<<", eff = "<<h_const_1000->Integral()/1000000.<<endl;
cout<<"Air, Nevent = "<<h_air->Integral()<<endl;

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

//h_const_0p005->SetMaximum(1.1*(h_air->GetMaximum()));
h_const_0p005->SetMaximum(2*(h_air->GetMaximum()));
h_const_0p005->Draw("hist");
h_const_0p005->GetXaxis()->SetTitle(xtitle);
h_const_0p005->GetYaxis()->SetTitle(ytitle);
h_const_0p05->Draw("hist,same");
h_const_0p5->Draw("hist,same");
//h_const_0p4->Draw("hist,same");
//h_const_0p3->Draw("hist,same");
h_const_0p2->Draw("hist,same");
h_const_0p1->Draw("hist,same");
h_const_1->Draw("hist,same");
h_const_10->Draw("hist,same");
h_const_100->Draw("hist,same");
//h_const_1000->Draw("hist,same");
//h_vac->Draw("hist,same");
h_air->Draw("hist,same");

TLegend *legend = new TLegend(0.25,0.60,0.45,0.95,NULL,"brNDC");
legend->AddEntry(h_const_0p005,"Const M_{DM}=0.005 GeV","F");
legend->AddEntry(h_const_0p05,"Const M_{DM}=0.05 GeV","F");
legend->AddEntry(h_const_0p1,"Const M_{DM}=0.1 GeV","F");
legend->AddEntry(h_const_0p2,"Const M_{DM}=0.2 GeV","F");
//legend->AddEntry(h_const_0p3,"Const M_{DM}=0.3 GeV","F");
//legend->AddEntry(h_const_0p4,"Const M_{DM}=0.4 GeV","F");
legend->AddEntry(h_const_0p5,"Const M_{DM}=0.5 GeV","F");
legend->AddEntry(h_const_1,"Const M_{DM}=1 GeV","F");
legend->AddEntry(h_const_10,"Const M_{DM}=10 GeV","F");
legend->AddEntry(h_const_100,"Const M_{DM}=100 GeV","F");
//legend->AddEntry(h_const_1000,"Const M_{DM}=1000 GeV","F");
//legend->AddEntry(h_vac,"Vacuum","F");
legend->AddEntry(h_air,"Air","F");
setLegend(legend);
legend->Draw();

can->Update();
can->Print(plot_name);

}
