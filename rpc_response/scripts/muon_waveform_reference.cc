// Existing 1GeV Heed electron records -> independently sampled linear avalanches.
// R001/R002/R006. Continue above physical-reference guard ONLY with explicit flag.
#include "../include/RpcResponse.hh"
#include <sstream>
#include <iomanip>
#include <iostream>
int main(int argc,char**argv)try{
 if(argc!=5)throw std::runtime_error("TRANSPORT.rpc ELECTRONS.csv OUTPUT_PREFIX EVENTS");
 auto tr=rpc::Transport::Load(argv[1]);int ne=std::stoi(argv[4]);const double dt=.1;const int bins=600;
 std::vector<double> currents(ne*bins,0.);std::vector<int> counts(ne),invalid(ne);
 std::ifstream input(argv[2]);rpc::Require(bool(input),"Cannot open electron source");std::string line;std::getline(input,line);
 std::mt19937_64 rng(20261031);
 while(std::getline(input,line)){
  std::stringstream ss(line);std::vector<std::string>a;for(std::string v;std::getline(ss,v,',');)a.push_back(v);
  int ev=std::stoi(a.at(0));rpc::Require(ev>=0&&ev<ne,"Event index outside range");double depth=std::stod(a.at(1)),birth=std::stod(a.at(2));++counts[ev];
  uint64_t n=1;bool bad=false;
  for(double s=0;s<depth&&n;){double ds=std::min(.01,depth-s);auto next=rpc::Advance(n,tr.alpha_per_mm,tr.eta_per_mm,ds,rng);
   if(n>10000000||next>10000000)bad=true;
   rpc::Require(next<100000000000000ULL,"Numerical population limit, not physical saturation");
   double t0=birth+s/tr.velocity_mm_ns,t1=birth+(s+ds)/tr.velocity_mm_ns;
   rpc::Require(t0>=0&&t1<bins*dt,"Waveform window truncated");
   double i=.1602176634*tr.velocity_mm_ns*.5*(double(n)+double(next))/tr.gap_mm;
   for(int k=int(t0/dt);k<=int(t1/dt);++k){double overlap=std::max(0.,std::min(t1,(k+1)*dt)-std::max(t0,k*dt));currents[ev*bins+k]+=i*overlap/dt;}
   n=next;s+=ds;
  }if(bad)++invalid[ev];
 }
 std::ofstream bin(std::string(argv[3])+".f64",std::ios::binary);bin.write(reinterpret_cast<const char*>(currents.data()),currents.size()*sizeof(double));
 std::ofstream csv(std::string(argv[3])+"_status.csv");csv<<"event,n_electrons,n_outside_linear_guard,charge_fC\n"<<std::setprecision(12);
 for(int e=0;e<ne;++e){double q=0;for(int k=0;k<bins;++k)q+=currents[e*bins+k]*dt*.001;csv<<e<<','<<counts[e]<<','<<invalid[e]<<','<<q<<'\n';}
 std::cout<<"events="<<ne<<" bins="<<bins<<" dt_ns="<<dt<<" seed=20261031; flagged over-guard references retained\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
