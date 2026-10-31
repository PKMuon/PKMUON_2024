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
#include "../include/garfield_lut.h"

using namespace std;

static const double SIGMA_MM = 0.749;       // fixed experimental resolution (mm)
static const double SIGMA_COORD = SIGMA_MM; // Geant4 positions are mm.
static const double ABOVE_MAX_MARKER = -999.0; // marker for energy above LUT max

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
    const char *outfile = "../../build/root_file/CryMuAna.root",
    const char *lutfile = "rpc_cluster_params.root")
{
  TRandom *rand = new TRandom();

  // ── Load Garfield LUT ──────────────────────────────────────────
  GarfieldLUT lut;
  bool useLUT = lut.Load(lutfile);
  if (!useLUT) {
    cout << "WARNING: Garfield LUT not found, using fixed sigma=" << SIGMA_MM << " mm\n";
  } else {
    cout << "Using Garfield LUT: gasGap=" << lut.gasZ << " cm, E="
         << lut.voltage / lut.gasZ << " V/cm\n";
  }

  TFile *file_in = TFile::Open(infile);
  if (!file_in || !file_in->IsOpen()) {
    cerr << "Failed to open input file: " << infile << endl;
    return;
  }

  TTree *tree_in = (TTree *)file_in->Get("tree");
  TClonesArray *Edeps = NULL;
  TClonesArray *Event = NULL;
  tree_in->SetBranchAddress("Edeps", &Edeps);
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
  vector<Double_t> ClusterE(nlayer), ClusterR(nlayer);
  Double_t CosThetaEdep, CosThetaSmeared;
  Double_t Ang, AngSmeared;
  Double_t E;

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
  tree_out->Branch("ClusterE", &ClusterE);
  tree_out->Branch("ClusterR", &ClusterR);

  struct RawHit { double x, y, e; int pid; };
  struct Cluster { double x, y, e; int pid; double rlut; };

  Long64_t nvalid = 0;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  double dz = 325 - params->LayerZ[0];
  Long64_t nentry = tree_in->GetEntries();

  for(Long64_t ientry = 0; ientry < nentry; ientry++) {
    if(ientry % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2)
           << (ientry / (double)nentry) * 100 << "%" << endl;
    }
    tree_in->GetEntry(ientry);

    bool valid = true;
    // Collect raw hits per layer
    vector<vector<RawHit>> layerHits(nlayer);
    Int_t nedep = Edeps->GetEntries();
    for(Int_t iedep = 0; iedep < nedep; ++iedep) {
      auto edep = (Edep *)Edeps->UncheckedAt(iedep);
      if ((size_t)edep->Id >= nlayer) continue;
      layerHits[edep->Id].push_back({edep->X, edep->Y, edep->Value, edep->Pid});
    }

    for (size_t l = 0; l < nlayer; ++l) {
      auto &hits = layerHits[l];

      // Energy-weighted barycenter (unchanged)
      double sumE = 0, sumX = 0, sumY = 0;
      for (auto &h : hits) { sumE += h.e; sumX += h.e * h.x; sumY += h.e * h.y; }
      if (sumE <= 2.72e-5) { valid = false; break; }
      XEdep[l] = sumX / sumE;
      YEdep[l] = sumY / sumE;
      ZEdep[l] = params->LayerZ[l] + dz;

      // ── Cluster merging on RAW hits (cross-pid supported) ──
      if (useLUT && hits.size() >= 2) {
        // Sort by energy descending (seed = largest energy hit)
        sort(hits.begin(), hits.end(),
             [](const RawHit &a, const RawHit &b) { return a.e > b.e; });

        vector<Cluster> clusters;
        for (auto &h : hits) {
          // Use fixed experimental sigma as the merging criterion
          // (the user wants merge based on sigma=0.749 mm)
          bool merged = false;
          for (auto &c : clusters) {
            double d = GetDistance(h.x, h.y, c.x, c.y);
            // Merge if distance < 2 * sigma  (fixed, not from LUT)
            if (d < 2.0 * SIGMA_COORD) {
              double E_total = c.e + h.e;
              // Cross-pid energy-weighted LUT radius:
              // R = (c.e/E_total)*LUT[c.pid](E_total) + (h.e/E_total)*LUT[h.pid](E_total)
              double R_LUT = lut.InterpolateRMS(c.pid, h.pid, c.e, h.e);
              // Flag if above max
              if (lut.IsAboveMax(E_total)) { c.rlut = ABOVE_MAX_MARKER; }
              else { c.rlut = R_LUT; }
              // Barycenter merge
              c.x = (c.x * c.e + h.x * h.e) / E_total;
              c.y = (c.y * c.e + h.y * h.e) / E_total;
              c.e = E_total;
              merged = true;
              break;
            }
          }
          if (!merged) {
            double r = lut.InterpolateRMS(h.pid, h.e);
            if (lut.IsAboveMax(h.e)) r = ABOVE_MAX_MARKER;
            clusters.push_back({h.x, h.y, h.e, h.pid, r});
          }
        }

        // Dominant cluster = largest energy after merging
        size_t best = 0;
        for (size_t i = 1; i < clusters.size(); ++i)
          if (clusters[i].e > clusters[best].e) best = i;

        auto &dom = clusters[best];
        ClusterE[l] = dom.e;
        ClusterR[l] = dom.rlut == ABOVE_MAX_MARKER ? ABOVE_MAX_MARKER : 10.0 * dom.rlut;

        // ── Smear with FIXED experimental sigma ────────────
        XSmeared[l] = dom.x + rand->Gaus(0, SIGMA_COORD);
        YSmeared[l] = dom.y + rand->Gaus(0, SIGMA_COORD);
      } else {
        // Single hit or no LUT
        ClusterE[l] = sumE;
        if (useLUT) {
          int pid0 = hits.empty() ? 13 : hits[0].pid;
          double r = lut.InterpolateRMS(pid0, sumE);
          ClusterR[l] = lut.IsAboveMax(sumE) ? ABOVE_MAX_MARKER : 10.0 * r;
        } else {
          ClusterR[l] = SIGMA_COORD;
        }
        // Smear with fixed sigma
        XSmeared[l] = XEdep[l] + rand->Gaus(0, SIGMA_COORD);
        YSmeared[l] = YEdep[l] + rand->Gaus(0, SIGMA_COORD);
      }
    }

    if (!valid) continue;

    CosThetaEdep   = GetCosTheta(XEdep,    YEdep,    ZEdep);
    CosThetaSmeared = GetCosTheta(XSmeared, YSmeared, ZEdep);
    Ang      = TMath::ACos(CosThetaEdep);
    AngSmeared = TMath::ACos(CosThetaSmeared);

    auto event = (::Event *)Event->UncheckedAt(0);
    E = event->E;

    nvalid++;
    tree_out->Fill();
  }

  file_out->cd();
  file_out->Write(NULL, TObject::kOverwrite);
  file_out->Close();
  file_in->Close();

  delete rand;
}
