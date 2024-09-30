#include "boot.h"
void draw1(TString var, TString xtitle, TString plot_name, double xmin, double xmax, double xbins, TCut cut, int legend_type, TChain* t_air, TTree* t_const_0p005, TTree* t_const_0p05, TTree* t_const_0p5, TTree* t_const_0p2, TTree* t_const_0p1, TTree* t_const_1, TTree* t_const_10, TTree* t_const_100);
void draw(){

const int n=7;

int root_type[n]={2,2,2,2,2,1,1};
TString var[n]={"costheta_smear","dca_smear","poca_x_smear","poca_y_smear","poca_z_smear","GemTrkEdep[2]","GemTrkEdep[3]"};
TString xtitle[n]={"Cos#theta","DCA (mm)","X_{PoCA} (mm)","Y_{PoCA} (mm)","Z_{PoCA} (mm)","Edep2 (MeV)","Edep3 (MeV)"};
TString plot_name[n]={
"../figure/costheta_smear.eps",
"../figure/DCA_smear.eps",
"../figure/POCA_X_smear.eps",
"../figure/POCA_Y_smear.eps",
"../figure/POCA_Z_smear.eps",
"../figure/Edep2.eps",
"../figure/Edep3.eps"
};
int legend_type[n]={1,4,3,3,3,4,4};
double xmin[n]={-1,   0,-1*1000,-1*1000,-1*1000,0,0};
double xmax[n]={ 1,1000, 1*1000, 1*1000, 1*1000,2,2};
double xbins = 100;
//double xbins = 1000;

TCut cut = "1";
//TCut cut = "costheta_smear<0.8";

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

TChain* t_air_raw = new TChain("T1");
//t_air_raw->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job/root/muCRY_0.root");
t_air_raw->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job/root/muCRY_*.root");

TChain* t_air_new = new TChain("T1");
//t_air_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job/analysis/root/muCRY_air_1.root");
t_air_new->Add("/home/pku/yuxd/bond/PKMUON_G4sim/analysis/job/analysis/root/muCRY_air_*.root");

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

for(int i=2; i<5; i++){
//for(int i=0; i<n; i++){
	cout<<"Drawing "<<var[i]<<endl;
	if(root_type[i]==1)
	  draw1(var[i], xtitle[i], plot_name[i], xmin[i], xmax[i], xbins, cut, legend_type[i], t_air_raw, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);
	else if(root_type[i]==2)
	  draw1(var[i], xtitle[i], plot_name[i], xmin[i], xmax[i], xbins, cut, legend_type[i], t_air_new, t_const_0p005, t_const_0p05, t_const_0p5, t_const_0p2, t_const_0p1, t_const_1, t_const_10, t_const_100);

}

} 

void draw1(TString var, TString xtitle, TString plot_name, double xmin, double xmax, double xbins, TCut cut, int legend_type, TChain* t_air, TTree* t_const_0p005, TTree* t_const_0p05, TTree* t_const_0p5, TTree* t_const_0p2, TTree* t_const_0p1, TTree* t_const_1, TTree* t_const_10, TTree* t_const_100){


setStyle();
TCanvas *can = new TCanvas("can","can",800,600);setPad();
can->SetLogy();
//if(var==Form("dca_smear")) can->SetLogx();

TH1D *h_const_0p005 = new TH1D("h_const_0p005","h_const_0p005",xbins,xmin,xmax);
TH1D *h_const_0p05 = new TH1D("h_const_0p05","h_const_0p05",xbins,xmin,xmax);
TH1D *h_const_0p5 = new TH1D("h_const_0p5","h_const_0p5",xbins,xmin,xmax);
//TH1D *h_const_0p4 = new TH1D("h_const_0p4","h_const_0p4",xbins,xmin,xmax);
//TH1D *h_const_0p3 = new TH1D("h_const_0p3","h_const_0p3",xbins,xmin,xmax);
TH1D *h_const_0p2 = new TH1D("h_const_0p2","h_const_0p2",xbins,xmin,xmax);
TH1D *h_const_0p1 = new TH1D("h_const_0p1","h_const_0p1",xbins,xmin,xmax);
TH1D *h_const_1 = new TH1D("h_const_1","h_const_1",xbins,xmin,xmax);
TH1D *h_const_10 = new TH1D("h_const_10","h_const_10",xbins,xmin,xmax);
TH1D *h_const_100 = new TH1D("h_const_100","h_const_100",xbins,xmin,xmax);
//TH1D *h_const_1000 = new TH1D("h_const_1000","h_const_1000",xbins,xmin,xmax);
TH1D *h_air = new TH1D("h_air","h_air",xbins,xmin,xmax);

TString a("Events/"); char b[20];  sprintf(b, "(%.5f",(xmax-xmin)/xbins); TString c(")");
TString ytitle = a + b + c;

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
h_const_0p2->SetLineColor(5);
h_const_0p1->SetLineColor(6);
h_const_1->SetLineColor(7);
h_const_10->SetLineColor(8);
h_const_100->SetLineColor(9);
//h_const_1000->SetLineColor(13);
h_air->SetLineColor(1);

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

h_const_0p005->SetMaximum(10*(h_air->GetMaximum()));
h_const_0p005->SetMinimum(0.1);
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

TLegend *legend;
if(legend_type==1) legend = new TLegend(0.18,0.57,0.28,0.93,NULL,"brNDC");
if(legend_type==2) legend = new TLegend(0.23,0.27,0.38,0.63,NULL,"brNDC");
if(legend_type==3) legend = new TLegend(0.43,0.27,0.58,0.63,NULL,"brNDC");
if(legend_type==4) legend = new TLegend(0.70,0.57,0.85,0.93,NULL,"brNDC");
legend->AddEntry(h_const_0p005,"M_{DM}=0.005 GeV","F");
legend->AddEntry(h_const_0p05,"M_{DM}=0.05 GeV","F");
legend->AddEntry(h_const_0p1,"M_{DM}=0.1 GeV","F");
legend->AddEntry(h_const_0p2,"M_{DM}=0.2 GeV","F");
//legend->AddEntry(h_const_0p3,"M_{DM}=0.3 GeV","F");
//legend->AddEntry(h_const_0p4,"M_{DM}=0.4 GeV","F");
legend->AddEntry(h_const_0p5,"M_{DM}=0.5 GeV","F");
legend->AddEntry(h_const_1,"M_{DM}=1 GeV","F");
legend->AddEntry(h_const_10,"M_{DM}=10 GeV","F");
legend->AddEntry(h_const_100,"M_{DM}=100 GeV","F");
//legend->AddEntry(h_const_1000,"M_{DM}=1000 GeV","F");
//legend->AddEntry(h_vac,"Vacuum","F");
legend->AddEntry(h_air,"Air","F");
setLegend(legend);
legend->Draw();

can->Update();
can->Print(plot_name);

/*
for (int i = 1; i <= h_air->GetNbinsX(); ++i) {
    double binLowEdge = h_air->GetBinLowEdge(i);
    double binContent = h_air->GetBinContent(i);
    cout << "Bin " << i << ": Low Edge = " << binLowEdge << ", Content = " << binContent << endl;
}
*/

delete h_const_0p005;
delete h_const_0p05;
delete h_const_0p5;
delete h_const_0p2;
delete h_const_0p1;
delete h_const_1;
delete h_const_10;
delete h_const_100;
delete h_air;
delete can;

}
