// ============================================================================
// garfield_lut.h — Garfield++ cluster radius LUT loader (analysis-side)
//
// Usage in analysis.C:
//   #include "../include/garfield_lut.h"
//   GarfieldLUT lut;
//   lut.Load("rpc_cluster_params.root");
//   double sigma = lut.InterpolateRMS(pid, energyMeV);           // cm
//   double sigma = lut.InterpolateRMS(pid1, pid2, e1, e2);      // merged
//   if (lut.IsAboveMax(energyMeV)) { /* special handling (>60MeV) */ }
// ============================================================================

#ifndef garfield_lut_h
#define garfield_lut_h

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

struct GarfieldLUT {
  double gasZ = 0, voltage = 0;

  /// Load the ROOT file produced by rpc_cluster_sim.
  bool Load(const std::string &path) {
    auto *f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) { std::cerr << "[LUT] cannot open " << path << "\n"; delete f; return false; }

    auto *t = dynamic_cast<TTree *>(f->Get("clusterParams"));
    if (!t) { std::cerr << "[LUT] no clusterParams tree\n"; f->Close(); delete f; return false; }

    double energy, rms;
    int pid;
    t->SetBranchAddress("energy", &energy);
    t->SetBranchAddress("rms",    &rms);
    t->SetBranchAddress("pid",    &pid);

    data_.clear(); byPid_.clear(); eMax_ = 0;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      data_[{pid, energy}] = rms;
      byPid_[pid].push_back(energy);
      if (energy > eMax_) eMax_ = energy;
    }
    for (auto &[p, v] : byPid_) std::sort(v.begin(), v.end());

    if (auto *info = dynamic_cast<TTree *>(f->Get("info"))) {
      info->SetBranchAddress("gasZ",    &gasZ);
      info->SetBranchAddress("voltage", &voltage);
      info->GetEntry(0);
    }
    f->Close(); delete f;
    ready_ = true;
    return true;
  }

  /// Log-log interpolated RMS cluster radius (cm).
  double InterpolateRMS(int pid, double eMeV) const { return InterpInternal(pid, eMeV); }

  /// Merged-cluster LUT radius: R = (e1/E)·LUT[pid1](E) + (e2/E)·LUT[pid2](E)
  double InterpolateRMS(int pid1, int pid2, double e1, double e2) const {
    double E = e1 + e2;
    if (E <= 0) return 0;
    return (e1 / E) * InterpInternal(pid1, E) + (e2 / E) * InterpInternal(pid2, E);
  }

  bool IsAboveMax(double eMeV) const { return eMeV > eMax_; }
  bool IsReady() const { return ready_; }

private:
  bool ready_ = false;
  double eMax_ = 0;
  std::map<std::tuple<int, double>, double> data_;
  std::map<int, std::vector<double>> byPid_;

  double InterpInternal(int pid, double eMeV) const {
    auto it = byPid_.find(pid);
    if (it == byPid_.end() || it->second.empty()) return 0;

    const auto &v = it->second;
    auto hi = std::upper_bound(v.begin(), v.end(), eMeV);
    if (hi == v.begin()) return data_.at({pid, v.front()});
    if (hi == v.end())   return data_.at({pid, v.back()});

    double e2 = *hi, e1 = *(hi - 1);
    if (e1 == e2) return data_.at({pid, e1});
    double r1 = data_.at({pid, e1}), r2 = data_.at({pid, e2});
    double w = (std::log(eMeV) - std::log(e1)) / (std::log(e2) - std::log(e1));
    w = std::max(0.0, std::min(1.0, w));
    return std::exp((1 - w) * std::log(std::max(r1, 1e-12)) +
                    w * std::log(std::max(r2, 1e-12)));
  }
};

#endif
