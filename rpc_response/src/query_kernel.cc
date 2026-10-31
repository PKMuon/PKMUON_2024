#include "RpcKernel.hh"
#include <chrono>
#include <iostream>
int main(int argc,char** argv) try {
  if(argc!=5) { std::cerr<<"Usage: rpc_query_kernel KERNEL.root DISTANCE_MM DRAWS SEED\n"; return 2; }
  rpc::Kernel kernel; kernel.Load(argv[1]);
  const auto count=std::stoull(argv[3]); rpc::Require(count>0,"Draws must be positive");
  const double distance=std::stod(argv[2]); std::mt19937_64 rng(std::stoull(argv[4]));
  double sum=0; std::uint64_t invalid=0;
  const auto start=std::chrono::steady_clock::now();
  for(std::uint64_t i=0;i<count;++i) {
    const auto s=kernel.Sample(distance,rng); if(!s.valid) { ++invalid; continue; } sum+=s.charge_e;
  }
  std::cout<<"Draws="<<count<<" invalid="<<invalid<<" mean_valid_charge_e="
           <<(count>invalid?sum/(count-invalid):0)<<" elapsed_s="
           <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
  std::cout<<"Conditional valid mean is diagnostic only when invalid>0; no efficiency prediction.\n";
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
