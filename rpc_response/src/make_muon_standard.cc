#include "RpcKernel.hh"
#include <Garfield/ComponentConstant.hh>
#include <Garfield/MediumMagboltz.hh>
#include <Garfield/Random.hh>
#include <Garfield/RandomEngineRoot.hh>
#include <Garfield/Sensor.hh>
#include <Garfield/TrackHeed.hh>
#include <chrono>
#include <filesystem>
#include <iostream>

int main(int argc,char** argv) try {
  if(argc!=8) {
    std::cerr<<"Usage: rpc_make_muon_standard TRANSPORT.rpc KERNEL.root GAS.gas PREFIX EVENTS KINETIC_EV SEED\n";
    return 2;
  }
  const auto transport=rpc::Transport::Load(argv[1]);
  rpc::Kernel kernel; kernel.Load(argv[2]);
  Garfield::MediumMagboltz gas;
  rpc::Require(gas.LoadGasFile(argv[3]),"Cannot load gas file");
  rpc::Require(std::abs(gas.GetTemperature()-transport.temperature_k)<.01 &&
               std::abs(gas.GetPressure()-transport.pressure_torr)<.01,"Gas condition mismatch");
  const std::string prefix=argv[4]; const int events=std::stoi(argv[5]);
  const double kinetic=std::stod(argv[6]); const auto seed=std::stoull(argv[7]);
  rpc::Require(events>0 && kinetic>0,"Invalid standard sample size/energy");
  for(const auto& ext:{".root","_events.csv","_electrons.csv","_clusters.csv","_kernel.csv","_transport.csv","_metadata.json"})
    rpc::Require(!std::filesystem::exists(prefix+ext),"Output exists: "+prefix+ext);
  Garfield::RandomEngineRoot engine(seed); Garfield::Random::SetEngine(engine);
  Garfield::ComponentConstant component;
  const double half=transport.gap_mm/20.; // mm -> cm, half gap
  component.SetArea(-14,-14,-half,14,14,half); component.SetMedium(&gas);
  component.SetElectricField(0,0,10*transport.voltage_v/transport.gap_mm);
  Garfield::Sensor sensor(&component);
  Garfield::TrackHeed heed(&sensor); heed.SetParticle("mu-"); heed.SetKineticEnergy(kinetic);
  rpc::Require(heed.Initialise(&gas),"Heed initialization failed");
  TFile output((prefix+".root").c_str(),"CREATE"); rpc::Require(!output.IsZombie(),"Cannot create muon standard");
  TTree tree("muonStandard","Isolated straight mu-; Heed ionization + linear reference lookup + empirical readout");
  int event=0,nClusters=0,nElectrons=0,nOutside=0,nInvalid=0;
  double transferred=0,xTrue=0,yTrue=0,xRec=0,yRec=0,charge=0,maxDrift=0;
  double xl=0,xr=0,yd=0,yu=0; Bool_t chargeValid=false;
  tree.Branch("event",&event); tree.Branch("n_clusters",&nClusters); tree.Branch("n_electrons",&nElectrons);
  tree.Branch("n_outside_electrons",&nOutside); tree.Branch("n_invalid_avalanches",&nInvalid);
  tree.Branch("energy_transfer_ev",&transferred); tree.Branch("x_true_mm",&xTrue); tree.Branch("y_true_mm",&yTrue);
  tree.Branch("x_rec_mm",&xRec); tree.Branch("y_rec_mm",&yRec);
  tree.Branch("electron_induced_fc",&charge); tree.Branch("charge_valid",&chargeValid);
  tree.Branch("max_nominal_drift_ns",&maxDrift);
  tree.Branch("nominal_txl_ns",&xl); tree.Branch("nominal_txr_ns",&xr);
  tree.Branch("nominal_tyd_ns",&yd); tree.Branch("nominal_tyu_ns",&yu);
  std::ofstream ev(prefix+"_events.csv"),el(prefix+"_electrons.csv"),cl(prefix+"_clusters.csv");
  rpc::Require(bool(ev)&&bool(el)&&bool(cl),"Cannot create CSV outputs");
  ev.precision(12); el.precision(12); cl.precision(12);
  ev<<"event,n_clusters,n_electrons,n_outside,n_invalid,energy_transfer_ev,x_true_mm,y_true_mm,x_rec_mm,y_rec_mm,charge_fc,charge_valid,max_drift_ns,txl_ns,txr_ns,tyd_ns,tyu_ns\n";
  el<<"event,depth_mm,creation_ns,nominal_arrival_ns,charge_fc,kernel_duration_ns,valid\n";
  cl<<"event,z_mm,t_ns,energy_transfer_ev,n_electrons\n";
  // Separate streams: Heed ionization; kernel draws; empirical readout.
  std::mt19937_64 kernelRng(seed+1),readoutRng(seed+2),positionRng(seed+3);
  rpc::NominalDelayLine delay;
  const auto start=std::chrono::steady_clock::now(); Long64_t invalidEvents=0,totalElectrons=0;
  for(event=0;event<events;++event) {
    xTrue=std::uniform_real_distribution<double>(-100,100)(positionRng);
    yTrue=std::uniform_real_distribution<double>(-100,100)(positionRng);
    rpc::Require(heed.NewTrack(xTrue/10,yTrue/10,-half+1.e-8,0,0,0,1),"Heed track generation failed");
    nClusters=nElectrons=nOutside=nInvalid=0; transferred=charge=maxDrift=0;
    for(const auto& cluster:heed.GetClusters()) {
      ++nClusters; transferred+=cluster.energy;
      cl<<event<<','<<10*cluster.z<<','<<cluster.t<<','<<cluster.energy<<','<<cluster.electrons.size()<<'\n';
      for(const auto& electron:cluster.electrons) {
        const double z=10*electron.z;
        if(z< -transport.gap_mm/2 || z>transport.gap_mm/2 || std::abs(electron.x)>14 || std::abs(electron.y)>14) { ++nOutside; continue; }
        ++nElectrons;
        // E points +z: electrons drift toward -z. Positive depth is distance to that plane.
        const double depth=std::clamp(z+transport.gap_mm/2,0.,transport.gap_mm);
        const double arrival=electron.t+depth/transport.velocity_mm_ns;
        maxDrift=std::max(maxDrift,arrival);
        const auto sample=kernel.Sample(depth,kernelRng);
        const double fc=sample.valid?sample.charge_e*1.602176634e-4:std::numeric_limits<double>::quiet_NaN();
        if(!sample.valid) ++nInvalid; else charge+=fc;
        el<<event<<','<<depth<<','<<electron.t<<','<<arrival<<','<<fc<<','<<sample.duration_ns<<','<<sample.valid<<'\n';
      }
    }
    chargeValid=nInvalid==0;
    if(!chargeValid) { ++invalidEvents; charge=std::numeric_limits<double>::quiet_NaN(); }
    const auto ends=delay.Encode(xTrue,yTrue,0,.749,readoutRng);
    xl=ends.xl; xr=ends.xr; yd=ends.yd; yu=ends.yu;
    xRec=delay.calibration.x(xl,xr); yRec=delay.calibration.y(yd,yu);
    tree.Fill(); totalElectrons+=nElectrons;
    ev<<event<<','<<nClusters<<','<<nElectrons<<','<<nOutside<<','<<nInvalid<<','<<transferred<<','
      <<xTrue<<','<<yTrue<<','<<xRec<<','<<yRec<<','<<charge<<','<<chargeValid<<','<<maxDrift<<','
      <<xl<<','<<xr<<','<<yd<<','<<yu<<'\n';
    if((event+1)%1000==0) std::cout<<event+1<<'/'<<events<<'\n';
  }
  tree.Write();
  TNamed("model_status","Heed standalone mu- (delta cascade owned by Heed); independent unsaturated kernel; invalid totals NaN; empirical total sigma 0.749mm; no threshold, space charge or timing validation").Write();
  TNamed("schema","rpc_muon_standard_v1").Write(); TNamed("kinetic_ev",argv[6]).Write(); TNamed("seed",argv[7]).Write();
  TNamed("source_kernel",argv[2]).Write(); TNamed("gas_file",argv[3]).Write();
  TNamed("position_distribution","uniform x,y in [-100,100] mm; perpendicular straight track through 2mm gap; NOT cosmic flux distribution").Write();
  output.Close();
  std::ofstream metadata(prefix+"_metadata.json"); metadata.precision(17);
  metadata<<"{\n\"schema\":\"rpc_muon_standard_v1\",\n\"events\":"<<events
          <<",\n\"muon_kinetic_ev\":"<<kinetic<<",\n\"angle_deg\":0,\n\"seed\":"<<seed
          <<",\n\"sigma_input_mm\":0.749,\n\"gas_gap_mm\":"<<transport.gap_mm
          <<",\n\"temperature_k\":"<<transport.temperature_k<<",\n\"pressure_torr\":"<<transport.pressure_torr
          <<",\n\"gas_voltage_v\":"<<transport.voltage_v<<",\n\"velocity_mm_ns\":"<<transport.velocity_mm_ns
          <<",\n\"alpha_per_mm\":"<<transport.alpha_per_mm<<",\n\"eta_per_mm\":"<<transport.eta_per_mm
          <<",\n\"heed_w_ev\":"<<heed.GetW()<<"\n}\n";
  // Export the source kernel with all invalid rows intact for plotting.
  TFile source(argv[2],"READ"); auto k=dynamic_cast<TTree*>(source.Get("singleElectronKernel"));
  rpc::Require(k,"Missing source kernel tree"); double depth=0,q=0,duration=0; Bool_t valid=false;
  k->SetBranchAddress("depth_mm",&depth); k->SetBranchAddress("electron_induced_e",&q);
  k->SetBranchAddress("duration_ns",&duration); k->SetBranchAddress("valid",&valid);
  std::ofstream kc(prefix+"_kernel.csv"); kc.precision(12); kc<<"depth_mm,charge_fc,duration_ns,valid\n";
  for(Long64_t i=0;i<k->GetEntries();++i) { k->GetEntry(i); kc<<depth<<','<<q*1.602176634e-4<<','<<duration<<','<<valid<<'\n'; }
  std::ofstream tc(prefix+"_transport.csv"); tc.precision(12);
  tc<<"field_kv_cm,velocity_mm_ns,alpha_per_mm,eta_per_mm,dl_sqrt_mm,dt_sqrt_mm\n";
  std::vector<double> fields,bfields,angles; gas.GetFieldGrid(fields,bfields,angles);
  fields.push_back(10*transport.voltage_v/transport.gap_mm); std::sort(fields.begin(),fields.end());
  for(auto f:fields) {
    double vx=0,vy=0,vz=0,a=0,b=0,dl=0,dt=0;
    rpc::Require(gas.ElectronVelocity(0,0,f,0,0,0,vx,vy,vz)&&gas.ElectronTownsend(0,0,f,0,0,0,a)&&
                 gas.ElectronAttachment(0,0,f,0,0,0,b)&&gas.ElectronDiffusion(0,0,f,0,0,0,dl,dt),"Transport export failed");
    tc<<f/1000<<','<<10*std::abs(vz)<<','<<a/10<<','<<b/10<<','<<dl*std::sqrt(10.)<<','<<dt*std::sqrt(10.)<<'\n';
  }
  std::cout<<"Muon standard events="<<events<<" electrons="<<totalElectrons<<" invalid_charge_events="<<invalidEvents
           <<" Heed_W_eV="<<heed.GetW()<<" elapsed_s="
           <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
