#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooCurve.h"
#include "RooHist.h"
using namespace RooFit ;

double gmean=0, gmean_err=0;
double gsigma=0, gsigma_err=0;
int gstatus=0;

void fit(){

        TFile *infile = new TFile("mu1GeV_new.root","");
        TTree *intree = (TTree*)infile->Get("T1");

        TFile *tmpfile = new TFile("tmp.root","RECREATE");
        //TCut cut = "fabs(dca_x)<500 && fabs(dca_y)<500 && fabs(dca_z)<500";
        TCut cut = "fabs(dca_x)<500 && fabs(dca_y)<500 && fabs(dca_z)<500 && angle>0.5";
        //TCut cut = "fabs(dca_x)<500 && fabs(dca_y)<500 && fabs(dca_z)<500 && angle<0.5";
        TTree *intree2 = intree->CopyTree(cut);

        double xmin=-500, xmax=500, xbins=100;
	RooRealVar *dca_z=new RooRealVar("dca_z","dca_z",xmin,xmax) ;
        dca_z->setRange("fitRange", -100, 100);
        RooDataSet *data = new RooDataSet("data","data",intree2,RooArgList(*dca_z));

	TString a("Events/"); char b[20];  sprintf(b, "(%.1f",(xmax-xmin)/xbins); TString d(" mm)");
	TString ytitle = a + b + d;

	// Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
	RooRealVar *mean1=new RooRealVar("mean1","mean of gauss1", -1.8533e+01, -100, 100);
	RooRealVar *mean2=new RooRealVar("mean2","mean of gauss2", -1.0000e+01, -100, 100);
	RooRealVar *sig1frac=new RooRealVar("sig1frac", "fraction", 2.0565e-01, 0, 1.0);
	RooRealVar *sigma1=new RooRealVar("sigma1","width of gauss1", 2.0000e+01, 0, 100);
	RooRealVar *sigma2=new RooRealVar("sigma2","width of gauss2", 3.0000e+01, 0, 100);

	RooGaussian *gaus1=new RooGaussian("gaus1","Signal component 1",*dca_z,*mean1,*sigma1) ;
	RooGaussian *gaus2=new RooGaussian("gaus2","Signal component 2",*dca_z,*mean2,*sigma2) ;
	RooAddModel  *sig=new RooAddModel("sig","",RooArgList(*gaus1,*gaus2),*sig1frac) ;//pdf!!

//      RooBifurGauss *gaus3=new RooBifurGauss("gaus3","BifurGauss",*dca_z,*mean1,*sigma1,*sigma2);

	//RooRealVar nsig("nsig","nsig", data->numEntries(), 0, 2.0*(data->numEntries()));
	RooRealVar nsig("nsig","nsig", 0.8*(intree2->GetEntries()), 0, 1.2*(intree2->GetEntries()));

	// Sum the composite signal and background 
	RooAddPdf  *model=new RooAddPdf("model","",RooArgList(*sig),RooArgList(nsig)) ;//create pdf
//	RooAddPdf  *model=new RooAddPdf("model","",RooArgList(*gaus1),RooArgList(nsig)) ;//create pdf
//	RooAddPdf  *model=new RooAddPdf("model","",RooArgList(*gaus3),RooArgList(nsig)) ;//create pdf

	// Fit model to data!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	RooFitResult *m_fitres = model->fitTo(*data,RooFit::Range("fitRange"),Extended(kTRUE),Save(kTRUE)) ;
	// Print structure of composite p.d.f.
	model->Print("t") ;
	m_fitres->Print("v");

        //calculate Double Gaussian mean and sigma
        //considering correlation term in error matrix
        RooFormulaVar mean("mean","sig1frac*mean1+(1-sig1frac)*mean2",RooArgSet(*sig1frac,*mean1,*mean2));
        RooFormulaVar sigma("sigma","sqrt(sig1frac*sigma1*sigma1+(1.0-sig1frac)*sigma2*sigma2+sig1frac*(1.0-sig1frac)*(mean1-mean2)*(mean1-mean2))",RooArgSet(*sig1frac,*mean1,*mean2,*sigma1,*sigma2));
        const TMatrixDSym& covMatrix = m_fitres->covarianceMatrix();
        Double_t meanValue = mean.getVal();
        Double_t sigmaValue = sigma.getVal();
        Double_t meanError = mean.getPropagatedError(*m_fitres);
        Double_t sigmaError = sigma.getPropagatedError(*m_fitres);
        //cout<<"meanValue = "<<meanValue<<endl;
        //cout<<"meanError = "<<meanError<<endl;
        //cout<<"sigmaValue = "<<sigmaValue<<endl;
        //cout<<"sigmaError = "<<sigmaError<<endl;
        //cout<<"arithmetic mean = "<<mean_ar<<" , sigma = "<<rms_ar<<endl;

		// --- Draw ---
		RooPlot* frame = dca_z->frame(xmin,xmax,xbins);
		data->plotOn(frame, RooFit::Name("Data"), Binning(xbins));
		model->plotOn(frame, RooFit::Name("fit"), LineColor(kRed));
		data->plotOn(frame, RooFit::Name("Data"), Binning(xbins));

	        //frame->SetTitle("Fit Z-vertex");
	        //frame->SetTitle(Form("Run %d",run));
	        frame->SetXTitle("PoCA_z (mm)");
		frame->GetYaxis()->SetTitle(ytitle);

		TPaveText *pt = new TPaveText(0.65,0.70,0.90,0.85,"BRNDC");//create pavetest
		pt->SetBorderSize(1);
		pt->SetFillColor(10);
		pt->SetTextAlign(12);
		pt->SetTextFont(22);
		pt->SetTextSize(0.035);

		TString Par0V = Form("%.2f", meanValue );
		TString Par0E = Form("%.2f", meanError );
		TString Par0  = "Mean = " + Par0V + " #pm " + Par0E ;
		TString Par1V = Form("%.2f", sigmaValue );
		TString Par1E = Form("%.2f", sigmaError );
		TString Par1  = "Sigma = " + Par1V + " #pm " + Par1E ;

		TText *text;//test!!
		text = pt->AddText(Par0);
		text = pt->AddText(Par1);

        TCanvas* Canvas = new  TCanvas( "Canvas", "can", 800, 600 );
        gPad->SetLogy(1);
        Canvas->cd();
        frame->Draw();
        pt->Draw();
        Canvas->Update();
        Canvas->Print("./figure/fit.eps");

        system("rm tmp.root");
}

