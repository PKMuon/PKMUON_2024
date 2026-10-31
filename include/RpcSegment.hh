#ifndef PKMUON_RPC_SEGMENT_HH
#define PKMUON_RPC_SEGMENT_HH
#include <Rtypes.h>

// Flat ROOT branches: no extra ROOT dictionary required. All units explicit.
struct RpcSegment {
  Long64_t event = 0;
  Int_t layer = 0, track = 0, parent = 0, pid = 0, root_track = 0;
  Double_t x0_mm=0,y0_mm=0,z0_mm=0,x1_mm=0,y1_mm=0,z1_mm=0;
  Double_t t0_ns=0,t1_ns=0,edep_ev=0,nonion_ev=0,kinetic_ev=0,length_mm=0;
  Double_t gap_center_z_mm=0,gap_mm=0;
};
struct RpcPrimary {
  Long64_t event=0;
  Int_t pid=0,source_mode=0,source_proposals=1;
  Double_t kinetic_ev=0,x_mm=0,y_mm=0,z_mm=0,dx=0,dy=0,dz=0,t_ns=0;
};
struct RpcLegacyEdep {
  Long64_t event=0;
  Int_t layer=0,pid=0,process=0;
  Double_t energy_ev=0,x_mm=0,y_mm=0;
  Bool_t legacy_trigger=false;
};
#endif
