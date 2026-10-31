#include "../include/RpcResponse.hh"
#include <iomanip>
#include <iostream>
int main(int argc,char**argv){
  const auto t=rpc::Transport::Load(argv[1]);
  const double distance=1.,step=.01;
  const uint64_t seed=20261021,limit=10000000;
  std::mt19937_64 rng(seed),check(seed);
  const auto reference=rpc::Grow(1,distance,t,step,limit,check);
  std::ofstream f(argv[2]); f<<std::setprecision(17)<<"t0_ns,t1_ns,n0,n1,current_nA\n";
  uint64_t n=1;double charge=0;
  for(double s=0;s<distance&&n;){
    double ds=std::min(step,distance-s);
    auto next=rpc::Advance(n,t.alpha_per_mm,t.eta_per_mm,ds,rng);
    rpc::Require(next<=limit,"Outside linear validity guard");
    const double current=.1602176634*t.velocity_mm_ns*.5*(n+next)/t.gap_mm;
    f<<s/t.velocity_mm_ns<<','<<(s+ds)/t.velocity_mm_ns<<','<<n<<','<<next<<','<<current<<'\n';
    charge+=.5*(n+next)*ds/t.gap_mm;s+=ds;n=next;
  }
  rpc::Require(std::abs(charge-reference.electron_induced_e)<1.e-9*std::max(1.,charge),"Ramo integral mismatch");
  std::cout<<std::setprecision(17)<<"{\"seed\":"<<seed<<",\"distance_mm\":1,\"step_mm\":0.01,\"electrons_final\":"<<n<<",\"ramo_e\":"<<charge<<",\"ramo_fC\":"<<charge*.0001602176634<<",\"grow_integral_matches\":true}";
}
