#include "boot.h"
void draw_costheta(){

TFile* fhist = new TFile("./Hist.root");

setStyle();
TCanvas *can = new TCanvas("can","can",800,600);setPad();
gStyle->SetTitleYOffset(1.0);     // Y title offset
gStyle->SetTitleXOffset(0.8);     // X title offset
gPad->SetBottomMargin(0.15);       // pad bottom margin for writing X title
gPad->SetRightMargin(0.02);       // pad right margin  for writing Y title
gPad->SetTopMargin(0.02);         // pad top margin for writing X title
can->SetLogy();

TH1D *h_const_0p005 = (TH1D*)fhist->Get("h_const_0p005");
TH1D *h_const_0p05 = (TH1D*)fhist->Get("h_const_0p05");
TH1D *h_const_0p1 = (TH1D*)fhist->Get("h_const_0p1");
TH1D *h_const_0p2 = (TH1D*)fhist->Get("h_const_0p2");
TH1D *h_const_0p5 = (TH1D*)fhist->Get("h_const_0p5");
TH1D *h_const_1 = (TH1D*)fhist->Get("h_const_1");
TH1D *h_const_10 = (TH1D*)fhist->Get("h_const_10");
TH1D *h_const_100 = (TH1D*)fhist->Get("h_const_100");

TH1D *h_maxwell_0p005 = (TH1D*)fhist->Get("h_maxwell_0p005");
TH1D *h_maxwell_0p05 = (TH1D*)fhist->Get("h_maxwell_0p05");
TH1D *h_maxwell_0p1 = (TH1D*)fhist->Get("h_maxwell_0p1");
TH1D *h_maxwell_0p2 = (TH1D*)fhist->Get("h_maxwell_0p2");
TH1D *h_maxwell_0p5 = (TH1D*)fhist->Get("h_maxwell_0p5");
TH1D *h_maxwell_1 = (TH1D*)fhist->Get("h_maxwell_1");
TH1D *h_maxwell_10 = (TH1D*)fhist->Get("h_maxwell_10");
TH1D *h_maxwell_100 = (TH1D*)fhist->Get("h_maxwell_100");

TH1D *h_air = (TH1D*)fhist->Get("h_air");
TH1D *h_vac = (TH1D*)fhist->Get("h_vac");

double xmin=-1, xmax=1, xbins=100;
TString a("Events/"); char b[20];  sprintf(b, "(%.2f",(xmax-xmin)/xbins); TString c(")");
TString ytitle = a + b + c;
TString xtitle="cos#it{#theta}";

/*
h_const_0p005->SetLineColor(2);
h_const_0p05->SetLineColor(3);
h_const_0p5->SetLineColor(4);
h_const_0p2->SetLineColor(5);
h_const_0p1->SetLineColor(6);
h_const_1->SetLineColor(7);
h_const_10->SetLineColor(8);
h_const_100->SetLineColor(9);
h_air->SetLineColor(1);
*/

h_maxwell_0p005->Scale((h_air->Integral())/(h_maxwell_0p005->Integral()));
h_maxwell_0p05->Scale((h_air->Integral())/(h_maxwell_0p05->Integral()));
h_maxwell_0p5->Scale((h_air->Integral())/(h_maxwell_0p5->Integral()));
h_maxwell_0p2->Scale((h_air->Integral())/(h_maxwell_0p2->Integral()));
h_maxwell_0p1->Scale((h_air->Integral())/(h_maxwell_0p1->Integral()));
h_maxwell_1->Scale((h_air->Integral())/(h_maxwell_1->Integral()));
h_maxwell_10->Scale((h_air->Integral())/(h_maxwell_10->Integral()));
h_maxwell_100->Scale((h_air->Integral())/(h_maxwell_100->Integral()));

h_maxwell_0p005->SetMaximum(10*(h_air->GetMaximum()));
h_maxwell_0p005->SetMinimum(0.1);
h_maxwell_0p005->Draw("hist");
h_maxwell_0p005->GetXaxis()->SetTitle(xtitle);
h_maxwell_0p005->GetYaxis()->SetTitle(ytitle);
h_maxwell_0p05->Draw("hist,same");
h_maxwell_0p5->Draw("hist,same");
h_maxwell_0p2->Draw("hist,same");
h_maxwell_0p1->Draw("hist,same");
h_maxwell_1->Draw("hist,same");
h_maxwell_10->Draw("hist,same");
h_maxwell_100->Draw("hist,same");
h_air->Draw("hist,same");

TLegend *legend = new TLegend(0.18,0.57,0.28,0.93,NULL,"brNDC");
legend->AddEntry(h_maxwell_0p005,"#it{M}_{DM}=0.005 GeV","F");
legend->AddEntry(h_maxwell_0p05,"#it{M}_{DM}=0.05 GeV","F");
legend->AddEntry(h_maxwell_0p1,"#it{M}_{DM}=0.1 GeV","F");
legend->AddEntry(h_maxwell_0p2,"#it{M}_{DM}=0.2 GeV","F");
legend->AddEntry(h_maxwell_0p5,"#it{M}_{DM}=0.5 GeV","F");
legend->AddEntry(h_maxwell_1,"#it{M}_{DM}=1 GeV","F");
legend->AddEntry(h_maxwell_10,"#it{M}_{DM}=10 GeV","F");
legend->AddEntry(h_maxwell_100,"#it{M}_{DM}=100 GeV","F");
legend->AddEntry(h_air,"Air","F");
setLegend(legend);
legend->Draw();

gPad->RedrawAxis();

can->Update();
can->Print("./figure/costheta.eps");

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

delete can;

}
