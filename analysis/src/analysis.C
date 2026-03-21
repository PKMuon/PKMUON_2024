#include <TClonesArray.h>
#include <TFile.h>
#include <TRandom.h>
#include <TTree.h>
#include <sys/time.h>
#include <TMath.h>

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

static Double_t GetDistance(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

void analysis(const char *infile = "../../build/root_file/CryMu.root",
    const char *outfile = "../../build/root_file/CryMuAna.root")
{
  TRandom *rand = new TRandom();

  TFile *file_in = TFile::Open(infile);
  if (!file_in || !file_in->IsOpen()) {
    cerr << "Failed to open input file: " << infile << endl;
    return;
  }

  TTree *tree_in = (TTree *)file_in->Get("tree");
  TClonesArray *Edeps = NULL;
  //TClonesArray *Tracks = NULL;
  TClonesArray *Event = NULL;
  tree_in->SetBranchAddress("Edeps", &Edeps);
  //tree_in->SetBranchAddress("Tracks", &Tracks);
  tree_in->SetBranchAddress("Event", &Event);
  tree_in->GetEntry(0);
  TTree *params_in = (TTree *)file_in->Get("params");
  TClonesArray *Params = NULL, *Processes = NULL;
  params_in->SetBranchAddress("Params", &Params);
  params_in->SetBranchAddress("Processes", &Processes);
  params_in->GetEntry(0);
  auto params = (::Params *)Params->At(0);
  size_t nlayer = params->LayerZ.size();

  TFile *file_out = TFile::Open(outfile, "RECREATE");
  TTree *tree_out = new TTree("tree", "tree");
  vector<Double_t> XEdep(nlayer), YEdep(nlayer), ZEdep(nlayer);
  vector<Double_t> XSmeared(nlayer), YSmeared(nlayer);
  Double_t CosThetaEdep, CosThetaSmeared;
  Double_t Ang, AngSmeared;
  Double_t E;
  Int_t SourcePid;
  Double_t SourcePx, SourcePy, SourcePz;
  Double_t SourceX, SourceY, SourceZ;
  Double_t ThetaP, PhiP;

  tree_out->Branch("XEdep", &XEdep);
  tree_out->Branch("YEdep", &YEdep);
  tree_out->Branch("ZEdep", &ZEdep);
  tree_out->Branch("XSmeared", &XSmeared);
  tree_out->Branch("YSmeared", &YSmeared);
  tree_out->Branch("CosThetaEdep", &CosThetaEdep);
  tree_out->Branch("CosThetaSmeared", &CosThetaSmeared);
  tree_out->Branch("Ang", &Ang);
  tree_out->Branch("AngSmeared", &AngSmeared);
  tree_out->Branch("E", &E);
  //tree_out->Branch("SourcePid", &SourcePid);
  //tree_out->Branch("SourcePx", &SourcePx);
  //tree_out->Branch("SourcePy", &SourcePy);
  //tree_out->Branch("SourcePz", &SourcePz);
  //tree_out->Branch("ThetaP", &ThetaP);
  //tree_out->Branch("PhiP", &PhiP);
  //tree_out->Branch("SourceX", &SourceX);
  //tree_out->Branch("SourceY", &SourceY);
  //tree_out->Branch("SourceZ", &SourceZ);


  // Temporaries.
  vector<Double_t> X2(nlayer), Y2(nlayer), Z2(nlayer), E2(nlayer);
  Long64_t nvalid = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  double dz = 425 - params->LayerZ[0];

  Long64_t nentry = tree_in->GetEntries();

  for(Long64_t ientry = 0; ientry < nentry; ientry++) {
    if(ientry % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2) << (ientry / (double)nentry) * 100 << "%" << endl;
    }
    tree_in->GetEntry(ientry);

    bool valid = true;
    
    // Simulate detector response.
    E2.assign(E2.size(), 0);
    X2.assign(X2.size(), 0);
    Y2.assign(Y2.size(), 0);
    Z2.assign(Z2.size(), 0);
    Int_t nedep = Edeps->GetEntries();
    //Int_t ntrack = Tracks->GetEntries();

    for(Int_t iedep = 0; iedep < nedep; ++iedep) {
      auto edep = (Edep *)Edeps->UncheckedAt(iedep);
      assert((size_t)edep->Id < E2.size());
      string process = edep->Process >= 0 ? ((Process *)Processes->UncheckedAt(edep->Process))->Name : "";
      //cout << "Processing Edep: id=" << edep->Id << " pid=" << edep->Pid << " trackid=" << edep->trackID << " process=" << process << endl;
      //if (edep->Value < 0.5) continue; 
      E2[edep->Id] += edep->Value;
      X2[edep->Id] += edep->Value * edep->X;
      Y2[edep->Id] += edep->Value * edep->Y;
      Z2[edep->Id] += edep->Value * params->LayerZ[edep->Id];
    }

    for (size_t l = 0; l < nlayer; ++l) {
      
      bool trigger = (E2[0] > 0 && E2[1] > 0) || (E2[1] > 0 && E2[2] > 0) || (E2[2] > 0 && E2[3] > 0);
      if(!trigger) {
        valid = false;
        break;
      }
      
      XEdep[l] = X2[l] / E2[l];
      YEdep[l] = Y2[l] / E2[l];     
      ZEdep[l] = Z2[l] / E2[l] + dz;

      // Simulate detector resolution.
      Double_t radius, phi, deltaphi, newphi;
      Double_t sigma = 0.749;  // mm
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

    CosThetaEdep = GetCosTheta(XEdep, YEdep, ZEdep);
    CosThetaSmeared = GetCosTheta(XSmeared, YSmeared, ZEdep);
    Ang = TMath::ACos(CosThetaEdep);
    AngSmeared = TMath::ACos(CosThetaSmeared);
    //cout << "ang: " << Ang << endl;

    auto event = (::Event *)Event->UncheckedAt(0);
    E = event->E, SourcePid = event->Pid;
    //SourcePx = event->Px, SourcePy = event->Py, SourcePz = event->Pz;
    //ThetaP = TMath::ATan(TMath::Sqrt(SourcePx*SourcePx + SourcePy*SourcePy) / SourcePz);
    //PhiP = TMath::ATan(SourcePy / SourcePx);
    //SourceX = event->X, SourceY = event->Y, SourceZ = event->Z;

    if (!valid) continue;
    nvalid++;

    tree_out->Fill();
  }

  gettimeofday(&end, NULL);
  time_t elapsedTime = 1000000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec);
  printf("time = %lf s\n", elapsedTime / 1e6);
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
