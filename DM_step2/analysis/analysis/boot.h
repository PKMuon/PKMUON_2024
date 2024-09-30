#ifndef TSTYLE_H
#define TSTYLE_H

#include "TStyle.h"
#include "TPad.h"
#include "TArrow.h"
#include "TLine.h"
#include "TLegend.h"
#include "iostream"
#include <algorithm>
#include <iomanip>

using namespace std;

void setStyle(){

		TStyle *myStyle = new TStyle("myStyle","my Style");
		myStyle->SetPalette(1);
		//myStyle->SetFillColor(0);
		myStyle->SetOptFit(0);             //1= display fitted result, 0 = no
		myStyle->SetOptTitle(0);             //1= display fitted result, 0 = no
		myStyle->SetOptStat(0);            //1= display the statistics, 0 = no
		myStyle->SetTitleFont(22);         // Font for X, Y, and Top title
		myStyle->SetNdivisions(505,"X");
		myStyle->SetNdivisions(505,"Y");

		myStyle->SetPadLeftMargin(0.2);         // pad left margin  for writing Y title
		myStyle->SetPadBottomMargin(0.22);       // pad bottom margin for writing X title
		myStyle->SetPadRightMargin(0.04);       // 0.04 pad right margin  for writing Y title
		myStyle->SetPadTopMargin(0.06);         // 0.04 pad top margin for writing X title

		myStyle->SetMarkerStyle(20);
		myStyle->SetMarkerSize(0.8);
		myStyle->SetStatFont(22);          // Font for statistics
		myStyle->SetLabelFont(22,"X");     // Font for X label
		myStyle->SetLabelFont(22,"Y");     // Font for Y label
		myStyle->SetTitleFont(22,"X");     // Font for X title 
		myStyle->SetTitleFont(22,"Y");     // Font for X title 
		myStyle->SetLabelSize(0.04,"X");   // Size for X label
		myStyle->SetLabelSize(0.04,"Y");   // Size for Y label
		myStyle->SetLabelOffset(0.015,"X");   // Size for X label
		myStyle->SetLabelOffset(0.015,"Y");   // Size for Y label
		myStyle->SetTitleXOffset(0.8);     // X title offset
		myStyle->SetTitleYOffset(1.2);     // Y title offset
		myStyle->SetTitleXSize(0.06);      // X title size
		myStyle->SetTitleYSize(0.06);      // Y title size
		myStyle->SetHistLineWidth(2);      // Histogram line width
		myStyle->SetFrameLineWidth(2);      // Frame line width
		myStyle->SetLineWidth(2);      // line width
		myStyle->SetStatX(0.95);           // Centroid X position of statistics box
		myStyle->SetStatY(0.95);           // Centroid Y position of statistics box
		myStyle->SetTitleX(0.2);          // Centroid X position of title box
		myStyle->SetTitleY(0.2);          // Centroid Y position of title box
		myStyle->SetTextSize(0.04);
		myStyle->SetFillStyle(1001);
		gROOT->SetStyle("myStyle");
}

void setPad()
{
		//gPad->SetGrid();
		gPad->SetBorderMode(0);
		gPad->SetBorderSize(0);
		gPad->SetFrameFillColor(0);
		gPad->SetFrameBorderMode(0);
		gPad->SetFillColor(0);
		gPad->SetLeftMargin(0.15);         // pad left margin  for writing Y title
		gPad->SetBottomMargin(0.22);       // pad bottom margin for writing X title
		gPad->SetRightMargin(0.08);       // pad right margin  for writing Y title
		gPad->SetTopMargin(0.06);         // pad top margin for writing X title
		//gPad->SetLogy(1);  
}

void setLegend(TLegend *legend)
{
	legend->SetBorderSize(0);
	legend->SetFillColor(0);
	legend->SetTextFont(22);
	legend->SetTextSize(0.04);
	legend->SetFillStyle(1001);
}

void setArrow(TArrow *arrow, Int_t fillcolor = 4, Int_t linecolor = 4)
{
	arrow->SetFillColor(fillcolor);
	arrow->SetLineColor(linecolor);
}

void setLine(TLine *line, Int_t color = 1, Int_t width = 2)
{
	line->SetLineColor(1);
	line->SetLineWidth(2);
}

void PrintProgressBar(int iteration, int total, char progress_char = '=')
{
    if(((iteration + 1) % max(1, (int)(total * 0.001)) != 0) )
        return;

    double progress = (1 + iteration)/(double)total;
    int n_progress_char = floor(50*progress);


    cout << "\r[";
    for (int i = 0; i < 50; ++i) 
    {
        if(i < n_progress_char)
            cout<<progress_char;
        else
            cout<<" ";
    }
    cout << "]";

    std::ios_base::fmtflags f( cout.flags() );
    cout << setprecision(1) << std::fixed;
    cout<<" "<< progress*100. << " %" << flush;
    cout.flags( f );


}

#endif
