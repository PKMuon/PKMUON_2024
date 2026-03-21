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

void draw_map(){

//TFile *f = new TFile("../root/mu1GeV_box.root","");
//TFile *f = new TFile("../root/mu1GeV_pbfew.root","");
//TFile *f = new TFile("../root/mu1GeV_pku.root","");
//TFile *f = new TFile("../root/muCRY_pku.root","");
//TFile *f = new TFile("../root/mu1GeV_2pb.root","");
//TFile *f = new TFile("../root/mu1GeV_pbbox.root","");
//TFile *f = new TFile("../root/muCRY_pbbox.root","");
//TFile *f = new TFile("../root/muCRY_all.root","");
//TFile *f = new TFile("../root/mu1GeV_all.root","");
//TFile *f = new TFile("../root/muCRY_vac.root","");
TFile *f = new TFile("../root/muCRY_air.root","");
TTree *t = (TTree*)f->Get("T1");

double muonE[4];
double rec_x[4], rec_y[4], rec_z[4];
double rec_x_smear[4], rec_y_smear[4];
double angle, angle_smear;
double poca_x, poca_y, poca_z, dca;
int poca_status;
double poca_x_smear, poca_y_smear, poca_z_smear, dca_smear;
int poca_status_smear;

t->SetBranchAddress("muonE",&muonE);
t->SetBranchAddress("rec_x",&rec_x);
t->SetBranchAddress("rec_y",&rec_y);
t->SetBranchAddress("rec_z",&rec_z);
t->SetBranchAddress("rec_x_smear",&rec_x_smear);
t->SetBranchAddress("rec_y_smear",&rec_y_smear);
t->SetBranchAddress("angle",&angle);
t->SetBranchAddress("angle_smear",&angle_smear);
t->SetBranchAddress("poca_x",&poca_x);
t->SetBranchAddress("poca_y",&poca_y);
t->SetBranchAddress("poca_z",&poca_z);
t->SetBranchAddress("dca",&dca);
t->SetBranchAddress("poca_status",&poca_status);
t->SetBranchAddress("poca_x_smear",&poca_x_smear);
t->SetBranchAddress("poca_y_smear",&poca_y_smear);
t->SetBranchAddress("poca_z_smear",&poca_z_smear);
t->SetBranchAddress("dca_smear",&dca_smear);
t->SetBranchAddress("poca_status_smear",&poca_status_smear);

double gem_z_sim[4] = {Z1,Z2,Z3,Z4};

vector<int> count(SPACE_N);
vector<double> sig(SPACE_N);
vector<double> lambda(SPACE_N);

struct timeval start;
struct timeval end;
unsigned long timer;

gettimeofday(&start, NULL); // 计时开始

int point_count = 0;
int bad_point_count = 0;
//int nevent = 0.1*(t->GetEntries());
int nevent = t->GetEntries();
cout<<"nevent = "<<nevent<<endl;
for(int ievent=0; ievent<nevent; ievent++){
//for(int ievent=4*nevent; ievent<5*nevent; ievent++){
        if (ievent % (int)(nevent / 10) == 0) cout << "Processing progress: " << ievent / (int)(nevent / 10) << "0%" << endl;
        t->GetEntry(ievent);

        ++point_count;
        if(poca_status == 0) continue;

        // construct muon path
        if (!CollectPath({rec_x[0], rec_y[0], Z1}, // p1
                         {rec_x[1], rec_y[1], Z2}, // p2
                         {rec_x[2], rec_y[2], Z3}, // p3
                         {rec_x[3], rec_y[3], Z4}, // p4
                         //1000, count, sig))
                         //muonE[1], count, sig))
                         //muonE[2], count, sig))
                         3000, count, sig))
        {
            ++bad_point_count;
        }

}
gettimeofday(&end, NULL); // 计时结束
timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
printf("time = %ld us\n", timer);

if (point_count)
{
    cerr << "Summary: good/total ="
         << point_count - bad_point_count << '/' << point_count << '\n';

    double lambda0_Pb = 0.00367;

    for (int i = 0; i < SPACE_N; ++i)
        lambda[i] = count[i] == 0 ? 0 : sig[i] / count[i] / lambda0_Pb;

    DumpStep(lambda, "PoCA Algorithm", "../figure/poca_cache", 0, 0); // 注意要加第五个参数 flag，否则无法生成二维图像

    vector<double> med_lambda(SPACE_N);
    median_filter(lambda, med_lambda);
    DumpStep(med_lambda, "PoCA Algorithm", "../figure/poca_median_cache", 0, 0); // 注意要加第五个参数 flag，否则无法生成二维图像
}

f->Close();

}

