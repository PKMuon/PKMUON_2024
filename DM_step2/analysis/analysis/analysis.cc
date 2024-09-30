#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>

#include "PoCA.h"

using namespace std;
extern double Z1;
extern double Z2;
extern double Z3;
extern double Z4;
double gem_z_sim[4] = {Z1,Z2,Z3,Z4};

void readtxt(ifstream &infile, int *eventid, double *xout, double *yout, double *pxout, double *pyout, double *pzout, double *eout);

void analysis(TString inrootfile1 , TString intxtfile , TString inrootfile2 , TString outrootfile){
TRandom *rand = new TRandom();
rand->SetSeed(2023);
double sigma=100*0.001;//100um

TFile *f1 = new TFile(inrootfile1,"");
TTree *t1 = (TTree*)f1->Get("T1");

double GemTrkPx1[3], GemTrkPy1[3], GemTrkPz1[3], GemTrkE1[3];
double GemTrkX1[3], GemTrkY1[3], GemTrkZ1[3];
double GemTrkEdep1[3];

t1->SetBranchAddress("GemTrkPx",&GemTrkPx1);
t1->SetBranchAddress("GemTrkPy",&GemTrkPy1);
t1->SetBranchAddress("GemTrkPz",&GemTrkPz1);
t1->SetBranchAddress("GemTrkE",&GemTrkE1);
t1->SetBranchAddress("GemTrkEdep",&GemTrkEdep1);
t1->SetBranchAddress("GemTrkX",&GemTrkX1);
t1->SetBranchAddress("GemTrkY",&GemTrkY1);
t1->SetBranchAddress("GemTrkZ",&GemTrkZ1);
//t1->Print();

ifstream infile(intxtfile);
TFile *f2 = new TFile(inrootfile2,"");

TTree *t2 = (TTree*)f2->Get("T1");

double GemTrkPx2[2], GemTrkPy2[2], GemTrkPz2[2], GemTrkE2[2];
double GemTrkX2[2], GemTrkY2[2], GemTrkZ2[2];
double GemTrkEdep2[2];
int EvtID2;

t2->SetBranchAddress("GemTrkX",&GemTrkX2);
t2->SetBranchAddress("GemTrkY",&GemTrkY2);
t2->SetBranchAddress("GemTrkZ",&GemTrkZ2);
t2->SetBranchAddress("GemTrkPx",&GemTrkPx2);
t2->SetBranchAddress("GemTrkPy",&GemTrkPy2);
t2->SetBranchAddress("GemTrkPz",&GemTrkPz2);
t2->SetBranchAddress("GemTrkE",&GemTrkE2);
t2->SetBranchAddress("GemTrkEdep",&GemTrkEdep2);
t2->SetBranchAddress("EvtID",&EvtID2);
//t2->Print();

TFile * fn = new TFile(outrootfile,"recreate");
TTree * tn = new TTree("T1","tree");
//tn = t2->CloneTree(0);

double m_GemTrkPx[4], m_GemTrkPy[4], m_GemTrkPz[4], m_GemTrkE[4];
double m_GemTrkX[4], m_GemTrkY[4], m_GemTrkZ[4];
double m_GemTrkX_smear[4], m_GemTrkY_smear[4];
double m_GemTrkEdep[4];

double angle, angle_smear;
double costheta, costheta_smear;
double poca_x, poca_y, poca_z, dca;
int poca_status;
double poca_x_smear, poca_y_smear, poca_z_smear, dca_smear;
int poca_status_smear;
tn->Branch("GemTrkE", m_GemTrkE, "GemTrkE[4]/D");
tn->Branch("GemTrkPx", m_GemTrkPx, "GemTrkPx[4]/D");
tn->Branch("GemTrkPy", m_GemTrkPy, "GemTrkPy[4]/D");
tn->Branch("GemTrkPz", m_GemTrkPz, "GemTrkPz[4]/D");
tn->Branch("GemTrkX", m_GemTrkX, "GemTrkX[4]/D");
tn->Branch("GemTrkY", m_GemTrkY, "GemTrkY[4]/D");
tn->Branch("GemTrkZ", m_GemTrkZ, "GemTrkZ[4]/D");
tn->Branch("GemTrkX_smear", m_GemTrkX_smear, "GemTrkX_smear[4]/D");
tn->Branch("GemTrkY_smear", m_GemTrkY_smear, "GemTrkY_smear[4]/D");
tn->Branch("GemTrkEdep", m_GemTrkEdep, "GemTrkEdep[4]/D");
tn->Branch("angle", &angle, "angle/D");
tn->Branch("angle_smear", &angle_smear, "angle_smear/D");
tn->Branch("costheta", &costheta, "costheta/D");
tn->Branch("costheta_smear", &costheta_smear, "costheta_smear/D");
tn->Branch("poca_x", &poca_x, "poca_x/D");
tn->Branch("poca_y", &poca_y, "poca_y/D");
tn->Branch("poca_z", &poca_z, "poca_z/D");
tn->Branch("dca", &dca, "dca/D");
tn->Branch("poca_status", &poca_status, "poca_status/I");
tn->Branch("poca_x_smear", &poca_x_smear, "poca_x_smear/D");
tn->Branch("poca_y_smear", &poca_y_smear, "poca_y_smear/D");
tn->Branch("poca_z_smear", &poca_z_smear, "poca_z_smear/D");
tn->Branch("dca_smear", &dca_smear, "dca_smear/D");
tn->Branch("poca_status_smear", &poca_status_smear, "poca_status_smear/I");

struct timeval start;
struct timeval end;
unsigned long timer;

gettimeofday(&start, NULL); // 计时开始

if(infile.is_open()){
int nevent = t2->GetEntries();
for(int ievent=0; ievent<nevent; ievent++){
        if (ievent % (int)(nevent / 10) == 0) cout << "Processing progress: " << ievent / (int)(nevent / 10) << "0%" << endl;

	//event matching
        t2->GetEntry(ievent);
	//cout<<"ievent = "<<ievent<<", EvtID in file1 = "<<EvtID2<<endl;
	t1->GetEntry(EvtID2);

	//momentum incoming to the detection region
	//obtain from rootfile1
	TVector3 vin(GemTrkPx1[2],GemTrkPy1[2],GemTrkPz1[2]);

	//momentum outcoming to the detection region
	//obtain from datfile
	int eventid=-1; double xout, yout;
        double pxout, pyout, pzout, eout;

	//while loop, matching rootfile2 and txtfile
	while (eventid!=EvtID2) {
		readtxt(infile, &eventid, &xout, &yout, &pxout, &pyout, &pzout, &eout);
		//cout<<"eventid = "<<eventid<<", EvtID in file1 = "<<EvtID2<<endl;
	}
        TVector3 vout(pxout,pyout,pzout);
	//cout<<"xout = "<<xout<<" , yout = "<<yout<<endl;

	//cout<<"vin = "<<vin.Px()<<","<<vin.Py()<<","<<vin.Pz()<<endl;
	//cout<<"vout = "<<vout.Px()<<","<<vout.Py()<<","<<vout.Pz()<<endl;
        //cout<<"costheta(vin,vout) = "<<cos(vin.Angle(vout))<<endl;

	//igem = 0,1 from file1
        for(int igem=0; igem<2; igem++){
		m_GemTrkPx[igem]=GemTrkPx1[igem];
		m_GemTrkPy[igem]=GemTrkPy1[igem];
		m_GemTrkPz[igem]=GemTrkPz1[igem];
		m_GemTrkE[igem]=GemTrkE1[igem];
		m_GemTrkEdep[igem]=GemTrkEdep1[igem];
		m_GemTrkX[igem]=GemTrkX1[igem];
		m_GemTrkY[igem]=GemTrkY1[igem];
		m_GemTrkZ[igem]=gem_z_sim[igem];
        }
	//igem = 2,3 from file2
	//notice igem-2
	//but!!! gem_z_sim do not need to igem-2
        for(int igem=2; igem<4; igem++){
		m_GemTrkPx[igem]=GemTrkPx2[igem-2];
		m_GemTrkPy[igem]=GemTrkPy2[igem-2];
		m_GemTrkPz[igem]=GemTrkPz2[igem-2];
		m_GemTrkE[igem]=GemTrkE2[igem-2];
		m_GemTrkEdep[igem]=GemTrkEdep2[igem-2];
		m_GemTrkX[igem]=GemTrkX2[igem-2];
		m_GemTrkY[igem]=GemTrkY2[igem-2];
		m_GemTrkZ[igem]=gem_z_sim[igem];
        }

        for(int igem=0; igem<4; igem++){
        m_GemTrkX_smear[igem] = m_GemTrkX[igem]+rand->Gaus(0,sigma);
        m_GemTrkY_smear[igem] = m_GemTrkY[igem]+rand->Gaus(0,sigma);
        }

        // angle <incoming vector, outcoming vector>
        TVector3 *Pos1 = new TVector3(m_GemTrkX[0],m_GemTrkY[0],m_GemTrkZ[0]);
        TVector3 *Pos2 = new TVector3(m_GemTrkX[1],m_GemTrkY[1],m_GemTrkZ[1]);
        TVector3 *Pos3 = new TVector3(m_GemTrkX[2],m_GemTrkY[2],m_GemTrkZ[2]);
        TVector3 *Pos4 = new TVector3(m_GemTrkX[3],m_GemTrkY[3],m_GemTrkZ[3]);
        TVector3 Veci = *Pos2 - *Pos1;
        TVector3 Veco = *Pos4 - *Pos3;
        angle = Veci.Angle(Veco) * 180 / M_PI;
        costheta = cos(Veci.Angle(Veco));
        //cout<<"Veci = "<<Veci.X()<<","<<Veci.Y()<<","<<Veci.Z()<<endl;
        //cout<<"Veco = "<<Veco.X()<<","<<Veco.Y()<<","<<Veco.Z()<<endl;
        //cout<<"costheta(position) = "<<costheta<<endl;

        TVector3 *Pos1_smear = new TVector3(m_GemTrkX_smear[0],m_GemTrkY_smear[0],m_GemTrkZ[0]);
        TVector3 *Pos2_smear = new TVector3(m_GemTrkX_smear[1],m_GemTrkY_smear[1],m_GemTrkZ[1]);
        TVector3 *Pos3_smear = new TVector3(m_GemTrkX_smear[2],m_GemTrkY_smear[2],m_GemTrkZ[2]);
        TVector3 *Pos4_smear = new TVector3(m_GemTrkX_smear[3],m_GemTrkY_smear[3],m_GemTrkZ[3]);
        TVector3 Veci_smear = *Pos2_smear - *Pos1_smear;
        TVector3 Veco_smear = *Pos4_smear - *Pos3_smear;
        angle_smear = Veci_smear.Angle(Veco_smear) * 180 / M_PI;
        costheta_smear = cos(Veci_smear.Angle(Veco_smear));

        // calculate PoCA, DCA

        if(CheckPoCAStatus({m_GemTrkX[0],m_GemTrkY[0],m_GemTrkZ[0]},
                           {m_GemTrkX[1],m_GemTrkY[1],m_GemTrkZ[1]},
                           {m_GemTrkX[2],m_GemTrkY[2],m_GemTrkZ[2]},
                           {m_GemTrkX[3],m_GemTrkY[3],m_GemTrkZ[3]}))
        {
                poca_status = 1;
                V3 poca = GetPoCAPoint({m_GemTrkX[0],m_GemTrkY[0],m_GemTrkZ[0]},
                                       {m_GemTrkX[1],m_GemTrkY[1],m_GemTrkZ[1]},
                                       {m_GemTrkX[2],m_GemTrkY[2],m_GemTrkZ[2]},
                                       {m_GemTrkX[3],m_GemTrkY[3],m_GemTrkZ[3]});
                poca_x = poca.x;
                poca_y = poca.y;
                poca_z = poca.z;
                dca = GetDCA({m_GemTrkX[0],m_GemTrkY[0],m_GemTrkZ[0]},
                             {m_GemTrkX[1],m_GemTrkY[1],m_GemTrkZ[1]},
                             {m_GemTrkX[2],m_GemTrkY[2],m_GemTrkZ[2]},
                             {m_GemTrkX[3],m_GemTrkY[3],m_GemTrkZ[3]});
        }
        else poca_status = 0;


        V3 poca_smear = GetPoCAPoint({m_GemTrkX_smear[0],m_GemTrkY_smear[0],m_GemTrkZ[0]},
                                     {m_GemTrkX_smear[1],m_GemTrkY_smear[1],m_GemTrkZ[1]},
                                     {m_GemTrkX_smear[2],m_GemTrkY_smear[2],m_GemTrkZ[2]},
                                     {m_GemTrkX_smear[3],m_GemTrkY_smear[3],m_GemTrkZ[3]});
        poca_x_smear = poca_smear.x;
        poca_y_smear = poca_smear.y;
        poca_z_smear = poca_smear.z;
        dca_smear = GetDCA({m_GemTrkX_smear[0],m_GemTrkY_smear[0],m_GemTrkZ[0]},
                           {m_GemTrkX_smear[1],m_GemTrkY_smear[1],m_GemTrkZ[1]},
                           {m_GemTrkX_smear[2],m_GemTrkY_smear[2],m_GemTrkZ[2]},
                           {m_GemTrkX_smear[3],m_GemTrkY_smear[3],m_GemTrkZ[3]});

        if(CheckPoCAStatus({m_GemTrkX_smear[0],m_GemTrkY_smear[0],m_GemTrkZ[0]},
                           {m_GemTrkX_smear[1],m_GemTrkY_smear[1],m_GemTrkZ[1]},
                           {m_GemTrkX_smear[2],m_GemTrkY_smear[2],m_GemTrkZ[2]},
                           {m_GemTrkX_smear[3],m_GemTrkY_smear[3],m_GemTrkZ[3]})) poca_status_smear = 1;
        else poca_status_smear = 0;

        tn->Fill();
}
}

gettimeofday(&end, NULL); // 计时结束
timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
printf("time = %ld us\n", timer);

fn->cd();
fn->Write();
fn->Close();
f1->Close();
f2->Close();

}

void readtxt(ifstream &infile, int *eventid, double *xout, double *yout, double *pxout, double *pyout, double *pzout, double *eout){
        infile>>*eventid; infile>>*xout; infile>>*yout;
        infile>>*pxout; infile>>*pyout; infile>>*pzout; infile>>*eout;
}
