#include "RpcResponse.hh"
#include <Garfield/MediumMagboltz.hh>
#include <iostream>

int main(int argc,char** argv) try {
  if(argc!=6) {
    std::cerr<<"Usage: rpc_make_gas OUTPUT.gas E_MIN_V_CM E_MAX_V_CM N_POINTS N_COLLISIONS\n"
             <<"Confirmed gas: pure R134a, 298.15 K, 760 Torr; gas data only, no assumed operating field.\n";
    return 2;
  }
  const double emin=std::stod(argv[2]),emax=std::stod(argv[3]);
  const int points=std::stoi(argv[4]),collisions=std::stoi(argv[5]);
  rpc::Require(emin>0 && emax>=emin && points>0 && collisions>0,"Invalid Magboltz grid");
  rpc::Require(!std::ifstream(argv[1]).good(),"Output already exists");
  Garfield::MediumMagboltz gas("C2H2F4",100.);
  gas.SetTemperature(298.15); gas.SetPressure(760.);
  gas.SetFieldGrid(emin,emax,points,false);
  gas.GenerateGasTable(collisions);
  rpc::Require(gas.WriteGasFile(argv[1]),"Gas table export failed");
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
