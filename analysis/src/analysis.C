#include <TRandom.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>

using namespace std;

void analysis(const char *infile = "../../build/root_file/CryMu_1.root",
    const char *outfile = "../../build/root_file/CryMuAna_1.root")
{
  TRandom *rand = new TRandom();

  // Input file and tree.
  TFile *file_in = new TFile(infile);
  TTree *tree_in = (TTree *)file_in->Get("T1");
  double RpcTrkPx[16], RpcTrkPy[16], RpcTrkPz[16], RpcTrkE[16];
  double RpcTrkEdep[16];
  double RpcTrkX[16], RpcTrkY[16], RpcTrkZ[16];
  bool RpcStatus;
  tree_in->SetBranchAddress("RpcTrkPx",   &RpcTrkPx);
  tree_in->SetBranchAddress("RpcTrkPy",   &RpcTrkPy);
  tree_in->SetBranchAddress("RpcTrkPz",   &RpcTrkPz);
  tree_in->SetBranchAddress("RpcTrkE",    &RpcTrkE);
  tree_in->SetBranchAddress("RpcTrkEdep", &RpcTrkEdep);
  tree_in->SetBranchAddress("RpcTrkX",    &RpcTrkX);
  tree_in->SetBranchAddress("RpcTrkY",    &RpcTrkY);
  tree_in->SetBranchAddress("RpcTrkZ",    &RpcTrkZ);
  tree_in->SetBranchAddress("RpcStatus",  &RpcStatus);

  // Output file and tree.
  TFile *file_out = new TFile(outfile, "RECREATE");
  TTree *tree_out = new TTree("T1", "Analysis Out Tree");
  tree_out = tree_in->CloneTree(0);  // copy 0 entries
  double RPCX[4], RPCY[4], RPCZ[4];
  double RPCX_smeared[4], RPCY_smeared[4];
  double costheta, costheta_smeared;
  double RpcTrkE_last;
  double RPCTrkPz_last;
  tree_out->Branch("RPCX", RPCX, "RPCX[4]/D");
  tree_out->Branch("RPCY", RPCY, "RPCY[4]/D");
  tree_out->Branch("RPCZ", RPCZ, "RPCZ[4]/D");
  tree_out->Branch("RPCX_smeared", RPCX_smeared, "RPCX_smeared[4]/D");
  tree_out->Branch("RPCY_smeared", RPCY_smeared, "RPCY_smeared[4]/D");
  tree_out->Branch("costheta", &costheta);
  tree_out->Branch("costheta_smeared", &costheta_smeared);
  tree_out->Branch("RpcTrkE_last", &RpcTrkE_last);
  tree_out->Branch("RPCTrkPz_last", &RPCTrkPz_last);

  struct timeval start, end;
  gettimeofday(&start, NULL);
  Long64_t nvalid = 0, ncos8 = 0;

  Long64_t nentry = tree_in->GetEntries();
  for(Long64_t ientry = 0; ientry < nentry; ientry++) {
    if(ientry % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2)
           << (ientry / (double)nentry) * 100 << "%" << endl;
    }
    tree_in->GetEntry(ientry);
    if(!RpcStatus) continue;
    nvalid++;

    double validRpcTrkX[8], validRpcTrkY[8], validRpcTrkZ[8];
    for(int i = 0; i < 8; i++) {
      validRpcTrkX[i] = (RpcTrkX[2*i] + RpcTrkX[2*i + 1]) / 2.0;
      validRpcTrkY[i] = (RpcTrkY[2*i] + RpcTrkY[2*i + 1]) / 2.0;
      validRpcTrkZ[i] = (RpcTrkZ[2*i] + RpcTrkZ[2*i + 1]) / 2.0;
    }

    // Compute RPC[XYZ].
    for(int i = 0; i < 4; ++i) {
      RPCX[i] = validRpcTrkX[2*i + 1];  // X-readout
      RPCY[i] = validRpcTrkY[2*i];      // Y-readout
      RPCZ[i] = (validRpcTrkZ[2*i] + validRpcTrkZ[2*i + 1]) / 2.0;
    }

    // Compute RPC[XY]_smeared.
    for(int i = 0; i < 4; i++) {
      double radius, phi, deltaphi, newphi;
      double sigma = 0.057;  // mm

      radius = hypot(RPCX[i], RPCY[i]);  // mm
      if(radius == 0) {
        RPCX_smeared[i] = RPCX[i];
        RPCY_smeared[i] = RPCY[i];
      } else {
        phi = atan2(RPCY[i], RPCX[i]);
        deltaphi = rand->Gaus(0, sigma) / radius;
        newphi = phi + deltaphi;
        RPCX_smeared[i] = radius * cos(newphi);
        RPCY_smeared[i] = radius * sin(newphi);
      }
    }

    // Compute costheta and costheta_smeared.
    TVector3 vec_in  = TVector3(RPCX[1], RPCY[1], RPCZ[1]) - TVector3(RPCX[0], RPCY[0], RPCZ[0]);
    TVector3 vec_out = TVector3(RPCX[3], RPCY[3], RPCZ[3]) - TVector3(RPCX[2], RPCY[2], RPCZ[2]);
    costheta = cos(vec_in.Angle(vec_out));
    vec_in  = TVector3(RPCX_smeared[1], RPCY_smeared[1], RPCZ[1]) - TVector3(RPCX_smeared[0], RPCY_smeared[0], RPCZ[0]);
    vec_out = TVector3(RPCX_smeared[3], RPCY_smeared[3], RPCZ[3]) - TVector3(RPCX_smeared[2], RPCY_smeared[2], RPCZ[2]);
    costheta_smeared = cos(vec_in.Angle(vec_out));

    if(costheta_smeared <= 0.8) ncos8++;
    tree_out->Fill();

    // Compute RpcTrkE_last.
     RpcTrkE_last = RpcTrkE[15];
    
    // Compute RpcTrkPz_last.
    RPCTrkPz_last = RpcTrkPz[15];

  }

  gettimeofday(&end, NULL);
  time_t time = 1000000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec);
  printf("time = %lf s\n", time / 1e6);
  cout << "Event count: " << nentry << endl;
  cout << "Event valid: " << nvalid << endl;
  cout << "cos < 0.8: " << ncos8 << endl;
  double quotient = static_cast<double>(nvalid) / nentry;
  double efficiency = quotient * 100;
  cout << "Efficiency: " << fixed << setprecision(2) << efficiency << "%" << endl;

  file_out->cd();
  file_out->Write();
  file_out->Close();
  file_in->Close();

  delete rand;
}
