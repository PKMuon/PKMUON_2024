#include "RpcResponse.hh"
#include "RpcSpatial.hh"
#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <chrono>
#include <iostream>
#include <set>
#include <filesystem>

// Explicit baseline: deposited-energy centroid + measured TOTAL spatial sigma.
// No claim of calibrated avalanche efficiency, waveform pileup, or two-end pairing.
int main(int argc,char** argv) try {
  if(argc<5 || argc>7) {
    std::cerr<<"Usage: rpc_digitize INPUT.root OUTPUT.root WINDOW_NS SEED [--allow-multi-track | --spatial-multi SIGMA_FACTOR]\n"
             <<"Window must be supplied explicitly. Multiple-track windows get invalid coordinates by default.\n";
    return 2;
  }
  const double window=std::stod(argv[3]); const auto seed=std::stoull(argv[4]);
  const bool allowMulti=argc==6 && std::string(argv[5])=="--allow-multi-track";
  const bool spatial=argc==7 && std::string(argv[5])=="--spatial-multi";
  rpc::Require(argc==5 || allowMulti || spatial,"Unknown option");
  const double spatialFactor=spatial?std::stod(argv[6]):0;
  rpc::Require(!spatial||(std::isfinite(spatialFactor)&&spatialFactor>0),"Invalid spatial sigma factor");
  rpc::Require(window>=0 && std::isfinite(window),"Invalid coincidence window");
  rpc::Require(std::string(argv[1])!=std::string(argv[2]),"Input and output must differ");
  TFile input(argv[1],"READ"); rpc::Require(!input.IsZombie(),"Cannot open input ROOT file");
  rpc::Require(input.Get("rpcSegments")!=nullptr,"Missing rpcSegments: rerun Geant4 with segment exporter");
  // Read and validate all branch types before creating output.
  TTreeReader reader("rpcSegments",&input);
  TTreeReaderValue<Long64_t> event(reader,"event");
  TTreeReaderValue<Int_t> layer(reader,"layer"),track(reader,"track"),parent(reader,"parent"),pid(reader,"pid");
  TTreeReaderValue<Double_t> x0(reader,"x0_mm"),y0(reader,"y0_mm"),z0(reader,"z0_mm");
  TTreeReaderValue<Double_t> x1(reader,"x1_mm"),y1(reader,"y1_mm"),z1(reader,"z1_mm");
  TTreeReaderValue<Double_t> t0(reader,"t0_ns"),t1(reader,"t1_ns"),energy(reader,"edep_ev"),nonion(reader,"nonion_ev");
  const std::string finalPath=argv[2], partialPath=finalPath+".partial";
  rpc::Require(!std::filesystem::exists(finalPath),"Output already exists");
  TFile output(partialPath.c_str(),"CREATE"); rpc::Require(!output.IsZombie(),"Cannot create partial output; existing files are protected");
  TTree hits("rpcHits","Empirical gas-deposit centroid; isolated sigma=0.749 mm; mm/ns/eV");
  Long64_t outEvent=0; Int_t outLayer=0,nSegments=0,nTracks=0;
  double x=0,y=0,time=0,weight=0; Bool_t ambiguous=false,validCoordinate=false;
  double xl=0,xr=0,yd=0,yu=0;
  Int_t windowIndex=0,positionIndex=0,nPositions=1,windowTracks=0;
  Bool_t idealSpatial=spatial,unresolved=false,pairingValidated=false;
  double centerX=0,centerY=0,windowTime=0;
  hits.Branch("event",&outEvent); hits.Branch("layer",&outLayer);
  hits.Branch("x_mm",&x); hits.Branch("y_mm",&y); hits.Branch("t_ns",&time);
  hits.Branch("ionizing_edep_ev",&weight); hits.Branch("n_segments",&nSegments);
  hits.Branch("n_tracks",&nTracks); hits.Branch("ambiguous",&ambiguous);
  hits.Branch("valid_coordinate",&validCoordinate);
  hits.Branch("nominal_txl_ns",&xl); hits.Branch("nominal_txr_ns",&xr);
  hits.Branch("nominal_tyd_ns",&yd); hits.Branch("nominal_tyu_ns",&yu);
  hits.Branch("window_index",&windowIndex); hits.Branch("position_index",&positionIndex);
  hits.Branch("n_positions",&nPositions); hits.Branch("window_n_tracks",&windowTracks);
  hits.Branch("window_t_ns",&windowTime);
  hits.Branch("centroid_x_mm",&centerX); hits.Branch("centroid_y_mm",&centerY);
  hits.Branch("ideal_spatial_model",&idealSpatial);
  hits.Branch("unresolved_component",&unresolved);
  hits.Branch("readout_pairing_validated",&pairingValidated);
  const rpc::NominalDelayLine delay;
  std::mt19937_64 rng(seed); std::vector<rpc::Deposit> deposits;
  Long64_t current=-1,rows=0;
  auto flush=[&]() {
    if(deposits.empty()) return;
    if(spatial) {
      windowIndex=0;
      for(const auto& group:rpc::SpatialReadout(deposits,window,.749*spatialFactor)) {
        nPositions=static_cast<Int_t>(group.positions.size()); windowTracks=group.tracks;
        windowTime=group.t_ns; positionIndex=0;
        for(const auto& position:group.positions) {
          const auto& h=position.hit;
          outEvent=current; outLayer=h.layer; time=h.t_ns; weight=h.weight;
          centerX=h.x_mm; centerY=h.y_mm; unresolved=position.tracks>1;
          ambiguous=false; validCoordinate=true;
          const auto endpoints=delay.Encode(h.x_mm,h.y_mm,h.t_ns,.749,rng);
          xl=endpoints.xl; xr=endpoints.xr; yd=endpoints.yd; yu=endpoints.yu;
          x=delay.calibration.x(xl,xr); y=delay.calibration.y(yd,yu);
          nSegments=static_cast<Int_t>(h.contributors); nTracks=position.tracks;
          hits.Fill(); ++rows; ++positionIndex;
        }
        ++windowIndex;
      }
      deposits.clear(); return;
    }
    windowIndex=0;
    auto response=rpc::Empirical(deposits,window,0,rng);
    for(const auto& h:response) {
      std::set<int> tracks;
      for(const auto& d:deposits)
        if(d.layer==h.layer && d.t_ns>=h.t_ns && d.t_ns<=h.t_ns+window) tracks.insert(d.track);
      ambiguous=tracks.size()>1;
      outEvent=current; outLayer=h.layer; time=h.t_ns; weight=h.weight;
      windowTime=time; windowTracks=static_cast<Int_t>(tracks.size());
      centerX=h.x_mm; centerY=h.y_mm; unresolved=tracks.size()>1;
      validCoordinate=!ambiguous || allowMulti;
      if(validCoordinate) {
        const auto endpoints=delay.Encode(h.x_mm,h.y_mm,h.t_ns,.749,rng);
        xl=endpoints.xl; xr=endpoints.xr; yd=endpoints.yd; yu=endpoints.yu;
        x=delay.calibration.x(xl,xr); y=delay.calibration.y(yd,yu);
      } else {
        x=y=xl=xr=yd=yu=std::numeric_limits<double>::quiet_NaN();
      }
      nSegments=static_cast<Int_t>(h.contributors); nTracks=static_cast<Int_t>(tracks.size()); hits.Fill(); ++rows;
      ++windowIndex;
    }
    deposits.clear();
  };
  const auto start=std::chrono::steady_clock::now();
  while(reader.Next()) {
    rpc::Require(current<0 || *event>=current,"Input event IDs must be ordered");
    if(current>=0 && *event!=current) flush();
    current=*event;
    const double e=*energy-*nonion; if(e<=0) continue;
    rpc::Deposit d; d.layer=*layer; d.track=*track; d.parent=*parent; d.pid=*pid;
    d.x_mm=.5*(*x0+*x1); d.y_mm=.5*(*y0+*y1); d.z_mm=.5*(*z0+*z1);
    d.t_ns=.5*(*t0+*t1); d.energy_ev=e; deposits.push_back(d);
  }
  rpc::Require(reader.GetEntryStatus()==TTreeReader::kEntryBeyondEnd,"ROOT segment read failed or branch schema mismatch");
  flush(); output.cd(); hits.Write();
  TNamed("response_model","empirical_v1; centroid=ionizing_local_edep; sigma_x_mm=0.749; sigma_y_mm=0.749; no efficiency or delay pairing").Write();
  TNamed("nominal_delay_model","8x1507-40A per axis; 320ns; pitch_x=3.457mm; pitch_y=3.463mm; bx=1.062mm; by=0; equal independent endpoint jitter; timing not calibrated").Write();
  TNamed("coincidence_window_ns",argv[3]).Write(); TNamed("random_seed",argv[4]).Write();
  TNamed("multi_track_policy",spatial?"ideal spatial readout: per-track centroids; close connected components merge; separated components each emit a coordinate; simulation-assisted, no waveform pairing":
    (allowMulti?"explicit centroid approximation; ambiguous flag":"flag and output NaN coordinate; includes secondary-track ambiguity")).Write();
  TNamed("spatial_separation_mm",std::to_string(.749*spatialFactor).c_str()).Write();
  output.Close();
  std::filesystem::rename(partialPath,finalPath);
  std::cout<<"Output hits="<<rows<<" elapsed_s="
           <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<'\n';
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
