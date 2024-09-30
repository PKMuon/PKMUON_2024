#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>

#include "TRandom.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using namespace std;

void analysis(){
TRandom *rand = new TRandom();
double sigma=100*0.001;//100um

TFile *f = new TFile("/home/pku/yuxd/bond/PKMUON_G4sim/DM_step1/build/root_file/muCRY.root","");
TTree *t = (TTree*)f->Get("T1");

double GemTrkX[3], GemTrkY[3], GemTrkZ[3];
double GemTrkPx[3], GemTrkPy[3], GemTrkPz[3], GemTrkE[3];

t->SetBranchAddress("GemTrkX",&GemTrkX);
t->SetBranchAddress("GemTrkY",&GemTrkY);
t->SetBranchAddress("GemTrkZ",&GemTrkZ);
t->SetBranchAddress("GemTrkPx",&GemTrkPx);
t->SetBranchAddress("GemTrkPy",&GemTrkPy);
t->SetBranchAddress("GemTrkPz",&GemTrkPz);
t->SetBranchAddress("GemTrkE",&GemTrkE);

double Eout, Pxout, Pyout, Pzout, Xout, Yout, Zout;

struct timeval start;
struct timeval end;
unsigned long timer;

gettimeofday(&start, NULL); // 计时开始

ofstream outfile("UpperGEMmuon.dat");
outfile.setf(ios::fixed,ios::floatfield); outfile.precision(6);

int number = 0;
int nevent = t->GetEntries();
for(int ievent=0; ievent<nevent; ievent++){
        if (ievent % (int)(nevent / 10) == 0) cout << "Processing progress: " << ievent / (int)(nevent / 10) << "0%" << endl;
        t->GetEntry(ievent);

	Eout=GemTrkE[2];
	Pxout=GemTrkPx[2];
	Pyout=GemTrkPy[2];
	Pzout=GemTrkPz[2];
	Xout=GemTrkX[2];
	Yout=GemTrkY[2];
	Zout=GemTrkZ[2];

        outfile<<ievent<<" "<<Xout<<" "<<Yout<<" "<<Zout<<" "<<Pxout<<" "<<Pyout<<" "<<Pzout<<" "<<Eout<<endl;

}

gettimeofday(&end, NULL); // 计时结束
timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
printf("time = %ld us\n", timer);

f->Close();

}

