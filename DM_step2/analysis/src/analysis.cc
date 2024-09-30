//#pragma cling add_include_path("/opt/homebrew/opt/boost/include")
//#pragma cling add_library_path("/opt/homebrew/opt/boost/lib")
//#pragma cling load("libboost_filesystem.dylib")
//#pragma cling load("libboost_system.dylib")

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>

#include "../include/PoCA.h"

using namespace std;
extern double Z1;
extern double Z2;
extern double Z3;
extern double Z4;

void analysis(TString inrootfile1 , TString intxtfile , TString inrootfile2 , TString outrootfile){
TRandom *rand = new TRandom();
double sigma=100*0.001;//100um

TFile *f1 = new TFile(inrootfile1,"");
TTree *t1 = (TTree*)f1->Get("T1");

std::vector<double> *vReadoutPosX1=0; //Edep X in readout bar
std::vector<double> *vReadoutPosY1=0; //Edep Y in readout bar
std::vector<double> *vReadoutPosZ1=0; //Edep Z in readout bar
std::vector<double> *vReadoutEdep1=0; //Edep Energy in readout bar
std::vector<double> *vReadoutE1=0; //Total Energy in readout bar
std::vector<int> *vReadoutTrkid1=0; //Trk id in readout bar
std::vector<int> *vReadoutTrkparentid1=0; //Trk parent id in readout bar
std::vector<double> *vPx1=0;
std::vector<double> *vPy1=0;
std::vector<double> *vPz1=0;

std::vector<double> *vPbPosX1=0;
std::vector<double> *vPbPosY1=0;
std::vector<double> *vPbPosZ1=0;
std::vector<int> *vPbTrkid1=0;

t1->SetBranchAddress("vReadoutPosX",&vReadoutPosX1);
t1->SetBranchAddress("vReadoutPosY",&vReadoutPosY1);
t1->SetBranchAddress("vReadoutPosZ",&vReadoutPosZ1);
t1->SetBranchAddress("vReadoutEdep",&vReadoutEdep1);
t1->SetBranchAddress("vReadoutE",&vReadoutE1);
t1->SetBranchAddress("vReadoutTrkid",&vReadoutTrkid1);
t1->SetBranchAddress("vReadoutTrkparentid",&vReadoutTrkparentid1);
t1->SetBranchAddress("vPx",&vPx1);
t1->SetBranchAddress("vPy",&vPy1);
t1->SetBranchAddress("vPz",&vPz1);
t1->SetBranchAddress("vPbPosX",&vPbPosX1);
t1->SetBranchAddress("vPbPosY",&vPbPosY1);
t1->SetBranchAddress("vPbPosZ",&vPbPosZ1);
t1->SetBranchAddress("vPbTrkid",&vPbTrkid1);
//t1->Print();

ifstream infile(intxtfile);
TFile *f2 = new TFile(inrootfile2,"");

TTree *t2 = (TTree*)f2->Get("T1");

std::vector<double> *vReadoutPosX2=0; //Edep X in readout bar
std::vector<double> *vReadoutPosY2=0; //Edep Y in readout bar
std::vector<double> *vReadoutPosZ2=0; //Edep Z in readout bar
std::vector<double> *vReadoutEdep2=0; //Edep Energy in readout bar
std::vector<double> *vReadoutE2=0; //Total Energy in readout bar
std::vector<int> *vReadoutTrkid2=0; //Trk id in readout bar
std::vector<int> *vReadoutTrkparentid2=0; //Trk parent id in readout bar
std::vector<double> *vPx2=0;
std::vector<double> *vPy2=0;
std::vector<double> *vPz2=0;

std::vector<double> *vPbPosX2=0;
std::vector<double> *vPbPosY2=0;
std::vector<double> *vPbPosZ2=0;
std::vector<int> *vPbTrkid2=0;

t2->SetBranchAddress("vReadoutPosX",&vReadoutPosX2);
t2->SetBranchAddress("vReadoutPosY",&vReadoutPosY2);
t2->SetBranchAddress("vReadoutPosZ",&vReadoutPosZ2);
t2->SetBranchAddress("vReadoutEdep",&vReadoutEdep2);
t2->SetBranchAddress("vReadoutE",&vReadoutE2);
t2->SetBranchAddress("vReadoutTrkid",&vReadoutTrkid2);
t2->SetBranchAddress("vReadoutTrkparentid",&vReadoutTrkparentid2);
t2->SetBranchAddress("vPx",&vPx2);
t2->SetBranchAddress("vPy",&vPy2);
t2->SetBranchAddress("vPz",&vPz2);
t2->SetBranchAddress("vPbPosX",&vPbPosX2);
t2->SetBranchAddress("vPbPosY",&vPbPosY2);
t2->SetBranchAddress("vPbPosZ",&vPbPosZ2);
t2->SetBranchAddress("vPbTrkid",&vPbTrkid2);
//t2->Print();

TFile * fn = new TFile(outrootfile,"recreate");
TTree * tn = new TTree("T1","tree");
tn = t2->CloneTree(0);

double muonE[4];
double rec_x[4], rec_y[4], rec_z[4];
double px[4], py[4], pz[4];
double rec_x_smear[4], rec_y_smear[4];
double angle, angle_smear;
double costheta, costheta_smear;
double poca_x, poca_y, poca_z, dca;
int poca_status;
double poca_x_smear, poca_y_smear, poca_z_smear, dca_smear;
int poca_status_smear;
double angle_gem12;
tn->Branch("muonE", muonE, "muonE[4]/D");
tn->Branch("rec_x", rec_x, "rec_x[4]/D");
tn->Branch("rec_y", rec_y, "rec_y[4]/D");
tn->Branch("rec_z", rec_z, "rec_z[4]/D");
tn->Branch("rec_x_smear", rec_x_smear, "rec_x_smear[4]/D");
tn->Branch("rec_y_smear", rec_y_smear, "rec_y_smear[4]/D");
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
tn->Branch("angle_gem12", &angle_gem12, "angle_gem12/D");

double gem_z_sim[4] = {Z1,Z2,Z3,Z4};
double gem_z_sim2[2] = {0.,13.34};

struct timeval start;
struct timeval end;
unsigned long timer;

gettimeofday(&start, NULL); // 计时开始

if(infile.is_open()){
int nevent = t2->GetEntries();
for(int ievent=0; ievent<nevent; ievent++){
        if (ievent % (int)(nevent / 10) == 0) cout << "Processing progress: " << ievent / (int)(nevent / 10) << "0%" << endl;
        t2->GetEntry(ievent);

        // initialization
        for(int igem=0; igem<4; igem++){
                rec_x[igem]=0;
                rec_y[igem]=0;
                rec_z[igem]=0;
                px[igem]=0;
                py[igem]=0;
                pz[igem]=0;
        }
        double Etot[4]={0,0,0,0};
        bool gemstatus[4]={false};

        TVector3 *pmom[4];

        // match file2 and file1 eventid
        // Warning!!! you must put this part of code in front of others. Otherwise, the 'continue' in the other part of code will make the event miss match
        int eventid; double tmp;
        infile>>eventid;
        //for(int itmp=0; itmp<6; itmp++) infile>>tmp;
        for(int itmp=0; itmp<2; itmp++) infile>>tmp;
        double pxout, pyout, pzout, eout;
        infile>>pxout; infile>>pyout; infile>>pzout; infile>>eout;
        TVector3 vout(pxout,pyout,pzout);
        t1->GetEntry(eventid);

        double nhit2 = vReadoutTrkid2->size();
//cout<<"nhit2 = "<<nhit2<<endl;
        if(nhit2 < 2) continue;
//cout<<"nhit2 = "<<nhit2<<endl;
        for(int itrk=0; itrk<nhit2; itrk++){
                // muon track id
                // trk id = 1 <==> parent id = 0
//cout<<__LINE__<<endl;
//cout<<"vReadoutTrkid2[itrk] = "<<(*vReadoutTrkid2)[itrk]<<endl;
//cout<<"vReadoutTrkparentid2[itrk] = "<<(*vReadoutTrkparentid2)[itrk]<<endl;
                if((*vReadoutTrkid2)[itrk]!=1) continue;
//cout<<__LINE__<<endl;

                // gem layer number
                int igem=-1;
                if( fabs((*vReadoutPosZ2)[itrk]-gem_z_sim2[0])<0.1 ) igem=2;
                if( fabs((*vReadoutPosZ2)[itrk]-gem_z_sim2[1])<0.1 ) igem=3;
                if(igem==-1) continue;
//cout<<"igem = "<<igem<<endl;
//cout<<"(*vReadoutPosZ2)[itrk] = "<<(*vReadoutPosZ2)[itrk]<<endl;

                // energy weight position
                muonE[igem]=(*vReadoutE2)[itrk];
                Etot[igem]+=(*vReadoutEdep2)[itrk];
                rec_x[igem]+=(*vReadoutEdep2)[itrk]*(*vReadoutPosX2)[itrk];
                rec_y[igem]+=(*vReadoutEdep2)[itrk]*(*vReadoutPosY2)[itrk];
                //rec_z[igem]+=(*vReadoutEdep2)[itrk]*(*vReadoutPosZ2)[itrk];
                rec_z[igem]+=(*vReadoutEdep2)[itrk]*gem_z_sim[igem];
                px[igem]=(*vPx2)[itrk];
                py[igem]=(*vPy2)[itrk];
                pz[igem]=(*vPz2)[itrk];
                pmom[igem] = new TVector3((*vPx2)[itrk],(*vPy2)[itrk],(*vPz2)[itrk]);
                gemstatus[igem]=true;
        }

        double nhit1 = vReadoutTrkid1->size();
//cout<<"nhit1 = "<<nhit1<<endl;
        if(nhit1 < 2) continue;
        //cout<<"nhit = "<<nhit<<endl;
        for(int itrk=0; itrk<nhit1; itrk++){
                // muon track id
                // trk id = 1 <==> parent id = 0
//cout<<__LINE__<<endl;
//cout<<"vReadoutTrkid1[itrk] = "<<(*vReadoutTrkid1)[itrk]<<endl;
                if((*vReadoutTrkid1)[itrk]!=1) continue;

                // gem layer number
                int igem=-1;
//cout<<"(*vReadoutPosZ1)[itrk] = "<<(*vReadoutPosZ1)[itrk]<<endl;
                if( fabs((*vReadoutPosZ1)[itrk]-gem_z_sim2[0])<0.1 ) igem=0;
                if( fabs((*vReadoutPosZ1)[itrk]-gem_z_sim2[1])<0.1 ) igem=1;
                if( fabs((*vReadoutPosZ1)[itrk]-13.44)<0.01 ) continue;
                if(igem==-1) continue;
//cout<<"igem = "<<igem<<endl;

                // energy weight position
                muonE[igem]=(*vReadoutE1)[itrk];
                Etot[igem]+=(*vReadoutEdep1)[itrk];
                rec_x[igem]+=(*vReadoutEdep1)[itrk]*(*vReadoutPosX1)[itrk];
                rec_y[igem]+=(*vReadoutEdep1)[itrk]*(*vReadoutPosY1)[itrk];
                //rec_z[igem]+=(*vReadoutEdep1)[itrk]*(*vReadoutPosZ1)[itrk];
                rec_z[igem]+=(*vReadoutEdep1)[itrk]*gem_z_sim[igem];
                px[igem]=(*vPx1)[itrk];
                py[igem]=(*vPy1)[itrk];
                pz[igem]=(*vPz1)[itrk];
                pmom[igem] = new TVector3((*vPx1)[itrk],(*vPy1)[itrk],(*vPz1)[itrk]);
                gemstatus[igem]=true;
        }

        bool hasFalse=false;
        for (bool element : gemstatus) {
            if (!element) {
                hasFalse = true;
                break;
            }
        }
        if (hasFalse) {
            //cout<<"Don't have signals in all GEM detector"<<endl;
            continue;
        }
        //cout<<__LINE__<<endl;
        //TVector3 vin(px[1],py[1],pz[1]);
        //cout<<"costheta(vin,vout) = "<<cos(vin.Angle(vout))<<endl;
        //cout<<"vout = "<<pxout<<","<<pyout<<","<<pzout<<endl;
        //cout<<"vin = "<<px[1]<<","<<py[1]<<","<<pz[1]<<endl;

        for(int igem=0; igem<4; igem++){
        rec_x[igem] = rec_x[igem]/Etot[igem];
        rec_y[igem] = rec_y[igem]/Etot[igem];
        rec_z[igem] = rec_z[igem]/Etot[igem];
        rec_x_smear[igem] = rec_x[igem]+rand->Gaus(0,sigma);
        rec_y_smear[igem] = rec_y[igem]+rand->Gaus(0,sigma);
        }
        //if(fabs(rec_x[0])>1/sqrt(2)/2*1000 || fabs(rec_y[0])>1/sqrt(2)/2*1000) continue;


        for(int igem=0; igem<4; igem++){
        //cout<<"gem layer "<<igem<<": x="<<rec_x[igem]<<" , y="<<rec_y[igem]<<" , z="<<rec_z[igem]<<endl;
        //cout<<"gem layer "<<igem<<": muon E ="<<muonE[igem]<<endl;
        }

        // angle <incoming vector, outcoming vector>
        TVector3 *Pos1 = new TVector3(rec_x[0],rec_y[0],rec_z[0]);
        TVector3 *Pos2 = new TVector3(rec_x[1],rec_y[1],rec_z[1]);
        TVector3 *Pos3 = new TVector3(rec_x[2],rec_y[2],rec_z[2]);
        TVector3 *Pos4 = new TVector3(rec_x[3],rec_y[3],rec_z[3]);
        TVector3 Veci = *Pos2 - *Pos1;
        TVector3 Veco = *Pos4 - *Pos3;
        angle = Veci.Angle(Veco) * 180 / M_PI;
        costheta = cos(Veci.Angle(Veco));
        //cout<<"costheta(position) = "<<costheta<<endl;
        //cout<<"Veci = "<<Veci.X()<<","<<Veci.Y()<<","<<Veci.Z()<<endl;
        //cout<<"Veco = "<<Veco.X()<<","<<Veco.Y()<<","<<Veco.Z()<<endl;
        //cout<<"gem0 pxpypz = "<<px[0]<<","<<py[0]<<","<<pz[0]<<endl;
        //cout<<"gem1 pxpypz = "<<px[1]<<","<<py[1]<<","<<pz[1]<<endl;
        //cout<<"gem2 pxpypz = "<<px[2]<<","<<py[2]<<","<<pz[2]<<endl;
        //cout<<"gem3 pxpypz = "<<px[3]<<","<<py[3]<<","<<pz[3]<<endl;
        //cout<<"angle = "<<angle/3.1415926*180.<<endl;
        //angle = cal_ang(rec_x[0],rec_y[0],rec_z[0],rec_x[1],rec_y[1],rec_z[1],rec_x[2],rec_y[2],rec_z[2],rec_x[3],rec_y[3],rec_z[3]); // old function

        TVector3 *Pos1_smear = new TVector3(rec_x_smear[0],rec_y_smear[0],rec_z[0]);
        TVector3 *Pos2_smear = new TVector3(rec_x_smear[1],rec_y_smear[1],rec_z[1]);
        TVector3 *Pos3_smear = new TVector3(rec_x_smear[2],rec_y_smear[2],rec_z[2]);
        TVector3 *Pos4_smear = new TVector3(rec_x_smear[3],rec_y_smear[3],rec_z[3]);
        TVector3 Veci_smear = *Pos2_smear - *Pos1_smear;
        TVector3 Veco_smear = *Pos4_smear - *Pos3_smear;
        angle_smear = Veci_smear.Angle(Veco_smear) * 180 / M_PI;
        costheta_smear = cos(Veci_smear.Angle(Veco_smear));
        //angle_smear = cal_ang(rec_x_smear[0],rec_y_smear[0],rec_z[0],rec_x_smear[1],rec_y_smear[1],rec_z[1],rec_x_smear[2],rec_y_smear[2],rec_z[2],rec_x_smear[3],rec_y_smear[3],rec_z[3]);
        angle_gem12 = pmom[0]->Angle(*pmom[1]);


        // calculate PoCA, DCA

        if(CheckPoCAStatus({rec_x[0],rec_y[0],rec_z[0]},
                           {rec_x[1],rec_y[1],rec_z[1]},
                           {rec_x[2],rec_y[2],rec_z[2]},
                           {rec_x[3],rec_y[3],rec_z[3]}))
        {
                poca_status = 1;
                V3 poca = GetPoCAPoint({rec_x[0],rec_y[0],rec_z[0]},
                                       {rec_x[1],rec_y[1],rec_z[1]},
                                       {rec_x[2],rec_y[2],rec_z[2]},
                                       {rec_x[3],rec_y[3],rec_z[3]});
                poca_x = poca.x;
                poca_y = poca.y;
                poca_z = poca.z;
                dca = GetDCA({rec_x[0],rec_y[0],rec_z[0]},
                             {rec_x[1],rec_y[1],rec_z[1]},
                             {rec_x[2],rec_y[2],rec_z[2]},
                             {rec_x[3],rec_y[3],rec_z[3]});
        }
        else poca_status = 0;


        // old function
        /*
        struct Line l1 = {{rec_x[0],rec_y[0],rec_z[0]}, {rec_x[1],rec_y[1],rec_z[1]}};
        struct Line l2 = {{rec_x[2],rec_y[2],rec_z[2]}, {rec_x[3],rec_y[3],rec_z[3]}};
        struct Line lper = line_perpendicular(l1,l2);
        poca_x = (lper.p1.x + lper.p2.x) / 2;
        poca_y = (lper.p1.y + lper.p2.y) / 2;
        poca_z = (lper.p1.z + lper.p2.z) / 2;
        dca = sqrt((lper.p1.x-lper.p2.x)*(lper.p1.x-lper.p2.x)+(lper.p1.y-lper.p2.y)*(lper.p1.y-lper.p2.y)+(lper.p1.z-lper.p2.z)*(lper.p1.z-lper.p2.z));
        poca_status = lper.status;
        */

        V3 poca_smear = GetPoCAPoint({rec_x_smear[0],rec_y_smear[0],rec_z[0]},
                                     {rec_x_smear[1],rec_y_smear[1],rec_z[1]},
                                     {rec_x_smear[2],rec_y_smear[2],rec_z[2]},
                                     {rec_x_smear[3],rec_y_smear[3],rec_z[3]});
        poca_x_smear = poca_smear.x;
        poca_y_smear = poca_smear.y;
        poca_z_smear = poca_smear.z;
        dca_smear = GetDCA({rec_x_smear[0],rec_y_smear[0],rec_z[0]},
                           {rec_x_smear[1],rec_y_smear[1],rec_z[1]},
                           {rec_x_smear[2],rec_y_smear[2],rec_z[2]},
                           {rec_x_smear[3],rec_y_smear[3],rec_z[3]});

        if(CheckPoCAStatus({rec_x_smear[0],rec_y_smear[0],rec_z[0]},
                           {rec_x_smear[1],rec_y_smear[1],rec_z[1]},
                           {rec_x_smear[2],rec_y_smear[2],rec_z[2]},
                           {rec_x_smear[3],rec_y_smear[3],rec_z[3]})) poca_status_smear = 1;
        else poca_status_smear = 0;

        //cout<<__LINE__<<endl;
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

