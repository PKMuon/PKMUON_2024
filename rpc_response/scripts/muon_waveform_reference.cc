// Existing 1GeV Heed electron records -> independently sampled linear avalanches.
// R001/R002/R006. Population markers are numerical settings, not physical limits.
#include "../include/RpcResponse.hh"
#include <sstream>
#include <iomanip>
#include <iostream>
int main(int argc,char**argv)try{
 if(argc!=5&&argc!=6)throw std::runtime_error("TRANSPORT.rpc ELECTRONS.csv OUTPUT_PREFIX EVENTS [POPULATION_MARKER=1e14]");
 const uint64_t marker=argc==6?std::stoull(argv[5]):100000000000000ULL;
 rpc::Require(marker>0&&marker<1000000000000000ULL,"Population marker must be positive and below numerical safety bound 1e15");
 auto tr=rpc::Transport::Load(argv[1]);int ne=std::stoi(argv[4]);const double dt=.1;const int bins=600;
 rpc::Require(ne>0,"Event count must be positive");
 std::vector<double> currents(ne*bins,0.);std::vector<int> counts(ne),invalid(ne);std::vector<uint64_t> max_population(ne,0);
 std::ifstream input(argv[2]);rpc::Require(bool(input),"Cannot open electron source");std::string line;std::getline(input,line);
 std::mt19937_64 rng(20261031);
 while(std::getline(input,line)){
  std::stringstream ss(line);std::vector<std::string>a;for(std::string v;std::getline(ss,v,',');)a.push_back(v);
  int ev=std::stoi(a.at(0));rpc::Require(ev>=0&&ev<ne,"Event index outside range");double depth=std::stod(a.at(1)),birth=std::stod(a.at(2));++counts[ev];
  uint64_t n=1;bool bad=false;
  for(double s=0;s<depth&&n;){double ds=std::min(.01,depth-s);auto next=rpc::Advance(n,tr.alpha_per_mm,tr.eta_per_mm,ds,rng);
   max_population[ev]=std::max(max_population[ev],std::max(n,next));
   if(n>marker||next>marker)bad=true;
   rpc::Require(next<1000000000000000ULL,"Numerical safety bound, not physical saturation");
   double t0=birth+s/tr.velocity_mm_ns,t1=birth+(s+ds)/tr.velocity_mm_ns;
   rpc::Require(t0>=0&&t1<bins*dt,"Waveform window truncated");
   double i=.1602176634*tr.velocity_mm_ns*.5*(double(n)+double(next))/tr.gap_mm;
   for(int k=int(t0/dt);k<=int(t1/dt);++k){double overlap=std::max(0.,std::min(t1,(k+1)*dt)-std::max(t0,k*dt));currents[ev*bins+k]+=i*overlap/dt;}
   n=next;s+=ds;
  }if(bad)++invalid[ev];
 }
 std::ofstream bin(std::string(argv[3])+".f64",std::ios::binary);bin.write(reinterpret_cast<const char*>(currents.data()),currents.size()*sizeof(double));
 std::ofstream csv(std::string(argv[3])+"_status.csv");csv<<"event,n_electrons,n_above_population_marker,charge_fC,max_family_population\n"<<std::setprecision(12);
 for(int e=0;e<ne;++e){double q=0;for(int k=0;k<bins;++k)q+=currents[e*bins+k]*dt*.001;csv<<e<<','<<counts[e]<<','<<invalid[e]<<','<<q<<','<<max_population[e]<<'\n';}
 std::ofstream metadata(std::string(argv[3])+"_generation.json");metadata<<"{\"events\":"<<ne<<",\"population_marker\":"<<marker<<",\"numerical_bound\":1000000000000000,\"seed\":20261031,\"dt_ns\":0.1,\"bins\":600,\"meaning\":\"user-selected marker only; no physical validity claim; no saturation or space charge\"}\n";
 std::cout<<"events="<<ne<<" bins="<<bins<<" dt_ns="<<dt<<" seed=20261031 population_marker="<<marker<<"; all references retained\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
