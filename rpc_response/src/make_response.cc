#include "RpcResponse.hh"
#include <Garfield/MediumMagboltz.hh>
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <iostream>

int main(int argc, char** argv) try {
  if (argc != 6) {
    std::cerr << "Usage: rpc_make_response GAS_FILE GAP_MM EFFECTIVE_GAS_VOLTAGE OUTPUT.rpc OUTPUT.root\n"
              << "Voltage is across GAS, not automatically the external electrode voltage.\n";
    return 2;
  }
  Garfield::MediumMagboltz gas;
  rpc::Require(gas.LoadGasFile(argv[1]), "Cannot load Magboltz gas file");
  rpc::Transport t;
  t.gas_file=argv[1]; t.gap_mm=std::stod(argv[2]); t.voltage_v=std::stod(argv[3]);
  rpc::Require(t.gap_mm>0 && t.voltage_v>0,"Invalid gap/voltage");
  t.temperature_k=gas.GetTemperature(); t.pressure_torr=gas.GetPressure();
  const double field=10*t.voltage_v/t.gap_mm; // V/mm -> V/cm
  std::vector<double> ef,bf,angles; gas.GetFieldGrid(ef,bf,angles);
  rpc::Require(!ef.empty() && field>=*std::min_element(ef.begin(),ef.end()) &&
               field<=*std::max_element(ef.begin(),ef.end()),"Field outside gas table coverage; no extrapolation allowed");
  double vx=0,vy=0,vz=0,alpha=0,eta=0,dl=0,dt=0;
  rpc::Require(gas.ElectronVelocity(0,0,field,0,0,0,vx,vy,vz),"Missing electron velocity");
  rpc::Require(gas.ElectronTownsend(0,0,field,0,0,0,alpha),"Missing Townsend coefficient");
  rpc::Require(gas.ElectronAttachment(0,0,field,0,0,0,eta),"Missing attachment coefficient");
  rpc::Require(gas.ElectronDiffusion(0,0,field,0,0,0,dl,dt),"Missing diffusion coefficients");
  t.velocity_mm_ns=10*std::sqrt(vx*vx+vy*vy+vz*vz);
  t.alpha_per_mm=alpha/10; t.eta_per_mm=eta/10;
  t.diffusion_l_sqrt_mm=dl*std::sqrt(10.); t.diffusion_t_sqrt_mm=dt*std::sqrt(10.);
  t.Validate();
  TFile out(argv[5],"CREATE"); rpc::Require(!out.IsZombie(),"Cannot create ROOT response (existing files are protected)");
  TTree tree("transport","Uniform-field transport; not a calibrated full RPC response");
  tree.Branch("gap_mm",&t.gap_mm); tree.Branch("voltage_v",&t.voltage_v);
  tree.Branch("temperature_k",&t.temperature_k); tree.Branch("pressure_torr",&t.pressure_torr);
  tree.Branch("velocity_mm_ns",&t.velocity_mm_ns); tree.Branch("alpha_per_mm",&t.alpha_per_mm);
  tree.Branch("eta_per_mm",&t.eta_per_mm);
  tree.Branch("diffusion_l_sqrt_mm",&t.diffusion_l_sqrt_mm);
  tree.Branch("diffusion_t_sqrt_mm",&t.diffusion_t_sqrt_mm);
  tree.Fill(); tree.Write();
  std::string composition;
  for(unsigned i=0;i<gas.GetNumberOfComponents();++i) {
    std::string name; double fraction=0; gas.GetComponent(i,name,fraction);
    composition+=name+":"+std::to_string(fraction)+";";
  }
  TNamed("gas_composition",composition.c_str()).Write();
  TNamed("gas_file",argv[1]).Write();
  TNamed("model_status","transport_only; no space-charge, electrode weighting map or electronics calibration").Write();
  t.Save(argv[4]); out.Close();
  std::cout<<"Transport exported: E="<<field<<" V/cm; v="<<t.velocity_mm_ns
           <<" mm/ns; alpha="<<t.alpha_per_mm<<" /mm; eta="<<t.eta_per_mm<<" /mm\n";
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
