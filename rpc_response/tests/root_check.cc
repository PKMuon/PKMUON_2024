#include "RpcResponse.hh"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <set>

int main(int argc,char** argv) try {
  rpc::Require(argc==3,"Usage: rpc_root_check fixture|hits|segments|baseline|kernel PATH");
  const std::string mode=argv[1];
  if(mode=="fixture") {
    TFile f(argv[2],"RECREATE"); TTree t("rpcSegments","Test input, mm ns eV");
    Long64_t event=0; int layer=0,track=1,parent=0,pid=13;
    double x0=0,y0=0,z0=0,x1=0,y1=0,z1=0,t0=0,t1=0,e=0,nonion=0;
    t.Branch("event",&event); t.Branch("layer",&layer); t.Branch("track",&track);
    t.Branch("parent",&parent); t.Branch("pid",&pid);
    t.Branch("x0_mm",&x0); t.Branch("y0_mm",&y0); t.Branch("z0_mm",&z0);
    t.Branch("x1_mm",&x1); t.Branch("y1_mm",&y1); t.Branch("z1_mm",&z1);
    t.Branch("t0_ns",&t0); t.Branch("t1_ns",&t1); t.Branch("edep_ev",&e); t.Branch("nonion_ev",&nonion);
    for(event=0;event<20000;++event) {
      x0=x1=1; y0=y1=2; e=100; t0=t1=0; t.Fill();
      x0=x1=3; y0=y1=4; e=300; t0=t1=.01; t.Fill();
    }
    t.Write(); return 0;
  }
  TFile f(argv[2],"READ"); rpc::Require(!f.IsZombie(),"Cannot open ROOT validation input");
  if(mode=="baseline") {
    auto t=dynamic_cast<TTree*>(f.Get("rpcHits")); rpc::Require(t && t->GetEntries()==800,"Wrong beam readout count");
    Bool_t ambiguous=false,valid=false; double x=0,y=0;
    t->SetBranchAddress("ambiguous",&ambiguous); t->SetBranchAddress("valid_coordinate",&valid);
    t->SetBranchAddress("x_mm",&x); t->SetBranchAddress("y_mm",&y);
    Long64_t accepted=0,flagged=0;
    for(Long64_t i=0;i<t->GetEntries();++i) {
      t->GetEntry(i);
      rpc::Require(valid!=ambiguous,"Ambiguous baseline must not silently emit a valid coordinate");
      if(valid) { ++accepted; rpc::Require(std::isfinite(x+y),"Nonfinite valid coordinate"); }
      else { ++flagged; rpc::Require(std::isnan(x)&&std::isnan(y),"Invalid coordinate must be explicit NaN"); }
    }
    std::cout<<"Baseline PASS: total="<<t->GetEntries()<<" valid="<<accepted<<" ambiguous="<<flagged<<'\n';
    return 0;
  }
  if(mode=="kernel") {
    auto t=dynamic_cast<TTree*>(f.Get("singleElectronKernel")); rpc::Require(t && t->GetEntries()>0,"Missing kernel");
    Bool_t valid=false; double depth=0,charge=0;
    t->SetBranchAddress("valid",&valid); t->SetBranchAddress("depth_mm",&depth);
    t->SetBranchAddress("electron_induced_e",&charge);
    Long64_t invalid=0; std::map<double,std::pair<int,int>> counts;
    for(Long64_t i=0;i<t->GetEntries();++i) {
      t->GetEntry(i); ++counts[depth].first;
      if(!valid) { ++invalid; ++counts[depth].second; }
      rpc::Require(std::isfinite(charge) && charge>=0,"Invalid stored charge");
    }
    rpc::Require(counts.begin()->first==0 && counts.rbegin()->first==2 && counts.size()==41,"Wrong depth grid");
    rpc::Require(f.Get("transport_snapshot")!=nullptr,"Missing portable transport snapshot");
    std::cout<<"Kernel PASS: samples="<<t->GetEntries()<<" invalid="<<invalid<<" depth_bins="<<counts.size()<<'\n';
    for(auto& p:counts) if(p.first==0 || p.first==1 || p.first==1.5 || p.first==2)
      std::cout<<"depth_mm="<<p.first<<" invalid="<<p.second.second<<"/"<<p.second.first<<'\n';
    return 0;
  }
  if(mode=="hits") {
    auto t=dynamic_cast<TTree*>(f.Get("rpcHits")); rpc::Require(t && t->GetEntries()==20000,"Wrong digitized count");
    double x=0,y=0,w=0; t->SetBranchAddress("x_mm",&x); t->SetBranchAddress("y_mm",&y);
    t->SetBranchAddress("ionizing_edep_ev",&w);
    double sx=0,sy=0,sxx=0,syy=0;
    for(Long64_t i=0;i<t->GetEntries();++i) {
      t->GetEntry(i); rpc::Require(w==400,"Energy unit/centroid error");
      x-=2.5; y-=3.5; sx+=x; sy+=y; sxx+=x*x; syy+=y*y;
    }
    const double n=t->GetEntries(),mx=sx/n,my=sy/n;
    const double sigx=std::sqrt(sxx/n-mx*mx),sigy=std::sqrt(syy/n-my*my);
    rpc::Require(std::abs(mx)<.025 && std::abs(my)<.025 && std::abs(sigx-.749)<.025 && std::abs(sigy-.749)<.025,
                 "End-to-end measured sigma or centroid bias check failed");
    std::cout<<"ROOT chain PASS: 20000 hits, sigma_x_mm="<<sigx<<", sigma_y_mm="<<sigy<<'\n'; return 0;
  }
  rpc::Require(mode=="segments","Unknown test mode");
  auto t=dynamic_cast<TTree*>(f.Get("rpcSegments")); rpc::Require(t && t->GetEntries()>0,"No gas segments");
  Long64_t event=0; int layer=0,parent=0; double z0=0,z1=0,center=0,gap=0,e=0;
  t->SetBranchAddress("event",&event); t->SetBranchAddress("layer",&layer); t->SetBranchAddress("parent",&parent);
  t->SetBranchAddress("z0_mm",&z0); t->SetBranchAddress("z1_mm",&z1);
  t->SetBranchAddress("gap_center_z_mm",&center); t->SetBranchAddress("gap_mm",&gap); t->SetBranchAddress("edep_ev",&e);
  std::set<Long64_t> events; std::set<int> layers; Long64_t secondaries=0;
  for(Long64_t i=0;i<t->GetEntries();++i) {
    t->GetEntry(i); events.insert(event); layers.insert(layer); secondaries+=parent!=0;
    rpc::Require(std::abs(gap-2)<1e-8 && std::abs((z0+z1)/2-center)<=gap/2+1e-6 && e>=0,
                 "Gas selection/length/energy unit validation failed");
  }
  rpc::Require(events.size()==200 && layers.size()==4,"Missing beam events or gas layers");
  std::cout<<"Geant4 segments PASS: "<<t->GetEntries()<<" steps; "<<events.size()
           <<" events; "<<layers.size()<<" layers; secondary steps="<<secondaries<<'\n';
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
