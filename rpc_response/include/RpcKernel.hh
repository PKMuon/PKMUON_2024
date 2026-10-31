#ifndef PKMUON_RPC_KERNEL_HH
#define PKMUON_RPC_KERNEL_HH
#include "RpcResponse.hh"
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>

namespace rpc {
struct KernelSample { double charge_e=0, duration_ns=0; bool valid=false; };
class Kernel {
 public:
  void Load(const std::string& path) {
    TFile f(path.c_str(),"READ"); Require(!f.IsZombie(),"Cannot open single-electron kernel");
    auto schema=dynamic_cast<TNamed*>(f.Get("schema"));
    Require(schema && std::string(schema->GetTitle())=="rpc_single_electron_kernel_v1; mm ns elementary_charge",
             "Unknown kernel schema");
    auto tree=dynamic_cast<TTree*>(f.Get("singleElectronKernel")); Require(tree,"Kernel tree missing");
    double depth=0,charge=0,duration=0; Bool_t valid=false;
    Require(tree->SetBranchAddress("depth_mm",&depth)>=0 &&
             tree->SetBranchAddress("electron_induced_e",&charge)>=0 &&
             tree->SetBranchAddress("duration_ns",&duration)>=0 &&
             tree->SetBranchAddress("valid",&valid)>=0,"Invalid kernel branch schema");
    std::map<double,std::vector<KernelSample>> temporary;
    for(Long64_t i=0;i<tree->GetEntries();++i) {
      Require(tree->GetEntry(i)>0,"Kernel read failed");
      Require(std::isfinite(depth) && depth>=0 && std::isfinite(charge) && charge>=0 &&
               std::isfinite(duration) && duration>=0,"Nonphysical kernel row");
      temporary[depth].push_back({charge,duration,valid});
    }
    Require(temporary.size()>=2 && temporary.begin()->first==0,"Insufficient depth coverage");
    data_=std::move(temporary);
  }
  // Nearest depth bin; valid only within coverage. The depth error is <= half a bin.
  // Invalid draws remain visible: never discard them and bias toward small avalanches.
  KernelSample Sample(double distance_mm,std::mt19937_64& rng) const {
    Require(!data_.empty() && std::isfinite(distance_mm) && distance_mm>=data_.begin()->first &&
             distance_mm<=data_.rbegin()->first,"Drift distance outside kernel coverage");
    auto hi=data_.lower_bound(distance_mm);
    if(hi!=data_.begin()) {
      auto lo=std::prev(hi);
      if(distance_mm-lo->first<=hi->first-distance_mm) hi=lo;
    }
    const auto& samples=hi->second;
    return samples[std::uniform_int_distribution<std::size_t>(0,samples.size()-1)(rng)];
  }
 private:
  std::map<double,std::vector<KernelSample>> data_;
};
}
#endif
