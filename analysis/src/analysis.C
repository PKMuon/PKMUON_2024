#include <TClonesArray.h>
#include <TFile.h>
#include <TRandom.h>
#include <TTree.h>
#include <sys/time.h>

#include <iomanip>
#include <iostream>
#include <vector>

#include "../../include/Object.hh"

using namespace std;

static Double_t GetCosTheta(Double_t x1, Double_t y1, double_t z1, Double_t x2, Double_t y2, double_t z2)
{
  Double_t dot = x1 * x2 + y1 * y2 + z1 * z2;
  Double_t square1 = x1 * x1 + y1 * y1 + z1 * z1;
  Double_t square2 = x2 * x2 + y2 * y2 + z2 * z2;
  return dot / sqrt(square1 * square2);
}

static Double_t GetCosTheta(const vector<Double_t> &X, const vector<Double_t> &Y, const vector<Double_t> &Z)
{
  size_t n = X.size();
  assert(n >= 4 && n % 2 == 0);
  size_t i1 = 0, i2 = n / 2 - 1, i3 = n / 2, i4 = n - 1;
  Double_t x1 = X[i2] - X[i1];
  Double_t y1 = Y[i2] - Y[i1];
  Double_t z1 = Z[i2] - Z[i1];
  Double_t x2 = X[i4] - X[i3];
  Double_t y2 = Y[i4] - Y[i3];
  Double_t z2 = Z[i4] - Z[i3];
  return GetCosTheta(x1, y1, z1, x2, y2, z2);
}

void analysis(const char *infile = "../../build/root_file/CryMu.root",
    const char *outfile = "../../build/root_file/CryMuAna.root")
{
  TRandom *rand = new TRandom();

  // Input file and tree.
  TFile *file_in = TFile::Open(infile);
  TTree *tree_in = (TTree *)file_in->Get("tree");
  TClonesArray *Edeps = NULL;
  tree_in->SetBranchAddress("Edeps", &Edeps);
  TTree *params_in = (TTree *)file_in->Get("params");
  TClonesArray *Params = NULL, *Processes = NULL;
  params_in->SetBranchAddress("Params", &Params);
  params_in->SetBranchAddress("Processes", &Processes);
  params_in->GetEntry(0);
  auto params = (::Params *)Params->At(0);
  size_t nlayer = params->LayerZ.size() / 2;

  // Output file and tree.
  TFile *file_out = TFile::Open(outfile, "RECREATE");
  TTree *tree_out = new TTree("tree", "tree");
  //tree_out = tree_in->CloneTree(0);  // copy 0 entries
  vector<Double_t> XEdep(nlayer), YEdep(nlayer), ZEdep(nlayer);
  vector<Double_t> XSmeared(nlayer), YSmeared(nlayer);
  Double_t CosThetaEdep, CosThetaSmeared;
  tree_out->Branch("XEdep", &XEdep);
  tree_out->Branch("YEdep", &YEdep);
  tree_out->Branch("ZEdep", &ZEdep);
  tree_out->Branch("XSmeared", &XSmeared);
  tree_out->Branch("YSmeared", &YSmeared);
  tree_out->Branch("CosThetaEdep", &CosThetaEdep);
  tree_out->Branch("CosThetaSmeared", &CosThetaSmeared);

  // Temporaries.
  vector<Double_t> X2(nlayer * 2), Y2(nlayer * 2), Z2(nlayer * 2), E2(nlayer * 2);
  Long64_t nvalid = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);

  Long64_t nentry = tree_in->GetEntries();
  for(Long64_t ientry = 0; ientry < nentry; ientry++) {
    if(ientry % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2) << (ientry / (double)nentry) * 100 << "%" << endl;
    }
    tree_in->GetEntry(ientry);

    // Simulate detector response.
    E2.assign(E2.size(), 0);
    X2.assign(X2.size(), 0);
    Y2.assign(Y2.size(), 0);
    Z2.assign(Z2.size(), 0);
    Int_t nedep = Edeps->GetEntries();
    for(Int_t iedep = 0; iedep < nedep; ++iedep) {
      auto edep = (Edep *)Edeps->UncheckedAt(iedep);
      assert((size_t)edep->Id < E2.size());
      assert(edep->Process < Processes->GetEntries());
      string process = edep->Process >= 0 ? ((Process *)Processes->UncheckedAt(edep->Process))->Name : "";
      cout << "Processing Edep: id=" << edep->Id << " pid=" << edep->Pid << " process=" << process << endl;
      E2[edep->Id] += edep->Value;
      X2[edep->Id] += edep->Value * edep->X;
      Y2[edep->Id] += edep->Value * edep->Y;
      Z2[edep->Id] += edep->Value * params->LayerZ[edep->Id];
    }
    bool valid = true;
    for(size_t l = 0; l < nlayer * 2; ++l) {
      if(!(E2[l] > 0)) {
        valid = false;
        break;
      }
      X2[l] /= E2[l], Y2[l] /= E2[l], Z2[l] /= E2[l];
    }
    if(!valid) continue;
    nvalid++;

    // Simulate readout system.
    for(size_t l = 0; l < nlayer; ++l) {
      XEdep[l] = X2[2 * l + 1];                      // XEdep-readout
      YEdep[l] = Y2[2 * l];                          // YEdep-readout
      ZEdep[l] = (Z2[2 * l] + Z2[2 * l + 1]) / 2.0;  // ZEdep-constant
    }

    // Simulate detector resolution.
    for(size_t l = 0; l < nlayer; l++) {
      Double_t radius, phi, deltaphi, newphi;
      Double_t sigma = 0.057;              // mm
      radius = hypot(XEdep[l], YEdep[l]);  // mm
      if(radius == 0) {
        XSmeared[l] = XEdep[l];
        YSmeared[l] = YEdep[l];
      } else {
        phi = atan2(YEdep[l], XEdep[l]);
        deltaphi = rand->Gaus(0, sigma) / radius;
        newphi = phi + deltaphi;
        XSmeared[l] = radius * cos(newphi);
        YSmeared[l] = radius * sin(newphi);
      }
    }

    // Compute scatter angle.
    CosThetaEdep = GetCosTheta(XEdep, YEdep, ZEdep);
    CosThetaSmeared = GetCosTheta(XSmeared, YSmeared, ZEdep);

    tree_out->Fill();
  }

  gettimeofday(&end, NULL);
  time_t time = 1000000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec);
  printf("time = %lf s\n", time / 1e6);
  cout << "Event count: " << nentry << endl;
  cout << "Event valid: " << nvalid << endl;
  double quotient = static_cast<double>(nvalid) / nentry;
  double efficiency = quotient * 100;
  cout << "Efficiency: " << fixed << setprecision(2) << efficiency << "%" << endl;

  file_out->cd();
  file_out->Write(NULL, TObject::kOverwrite);
  file_out->Close();
  file_in->Close();

  delete rand;
}
