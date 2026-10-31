#include "RpcResponse.hh"
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <chrono>
#include <iostream>

int main(int argc,char** argv) try {
  if(argc!=9) {
    std::cerr<<"Usage: rpc_make_kernel TRANSPORT.rpc OUTPUT.root DEPTH_BINS SAMPLES_PER_BIN SEED MAX_ELECTRONS STEP_MM --linear-reference\n"
             <<"Unsaturated single-electron planar Ramo kernel; max electrons is a validity guard, NOT saturation physics.\n";
    return 2;
  }
  auto transport=rpc::Transport::Load(argv[1]);
  const int bins=std::stoi(argv[3]),samples=std::stoi(argv[4]);
  const auto seed=std::stoull(argv[5]),limit=std::stoull(argv[6]); const double step=std::stod(argv[7]);
  // Last argument is an explicit model acknowledgement, preventing accidental production use.
  rpc::Require(std::string(argv[8])=="--linear-reference", "Must acknowledge --linear-reference");
  rpc::Require(bins>0 && samples>0 && step>0 && limit>0,"Invalid kernel grid");
  TFile file(argv[2],"CREATE"); rpc::Require(!file.IsZombie(),"Cannot create kernel; output already exists?");
  TTree tree("singleElectronKernel","Linear birth-death; planar electron-only induced charge");
  double depth=0,charge=0,duration=0; ULong64_t finalElectrons=0; Bool_t valid=false;
  tree.Branch("depth_mm",&depth); tree.Branch("electron_induced_e",&charge);
  tree.Branch("duration_ns",&duration); tree.Branch("final_electrons",&finalElectrons);
  tree.Branch("valid",&valid);
  std::mt19937_64 rng(seed); Long64_t rejected=0;
  const auto begin=std::chrono::steady_clock::now();
  for(int i=0;i<=bins;++i) for(int j=0;j<samples;++j) {
    depth=transport.gap_mm*i/bins; charge=duration=0; finalElectrons=0; valid=true;
    try {
      const auto a=rpc::Grow(1,depth,transport,step,limit,rng);
      charge=a.electron_induced_e; duration=a.duration_ns; finalElectrons=a.final_electrons;
    } catch(const std::exception&) { valid=false; ++rejected; }
    tree.Fill();
  }
  tree.Write();
  TNamed("transport_path",argv[1]).Write();
  std::ifstream source(argv[1]);
  const std::string snapshot((std::istreambuf_iterator<char>(source)),std::istreambuf_iterator<char>());
  TNamed("transport_snapshot",snapshot.c_str()).Write();
  TNamed("model_status","linear_reference_only; no space charge; no segmented weighting field; invalid draws must not be discarded and resampled").Write();
  TNamed("schema","rpc_single_electron_kernel_v1; mm ns elementary_charge").Write();
  TNamed("seed",argv[5]).Write(); TNamed("max_electrons_guard",argv[6]).Write();
  TNamed("step_mm",argv[7]).Write();
  TTree readout("empiricalReadout","Measured isolated-hit coordinate sigma and nominal delay line");
  double sigma=.749,totalDelay=320,pitchX=3.457,pitchY=3.463;
  readout.Branch("sigma_x_mm",&sigma); readout.Branch("sigma_y_mm",&sigma);
  readout.Branch("nominal_delay_ns",&totalDelay); readout.Branch("nominal_pitch_x_mm",&pitchX);
  readout.Branch("nominal_pitch_y_mm",&pitchY); readout.Fill(); readout.Write();
  file.Close();
  std::cout<<"Kernel samples="<<(bins+1)*samples<<" invalid="<<rejected<<" elapsed_s="
           <<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<'\n';
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
