#include "RpcResponse.hh"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <filesystem>
#include <iostream>

void Export(TFile& file,const std::string& name,const std::string& output,
            const std::vector<std::string>& columns) {
  auto tree=dynamic_cast<TTree*>(file.Get(name.c_str()));
  rpc::Require(tree,"Missing tree: "+name);
  rpc::Require(!std::filesystem::exists(output),"Output exists: "+output);
  std::vector<TLeaf*> leaves;
  for(const auto& column:columns) {
    auto leaf=tree->GetLeaf(column.c_str()); rpc::Require(leaf,"Missing column: "+column); leaves.push_back(leaf);
  }
  std::ofstream csv(output); csv.precision(17);
  for(size_t i=0;i<columns.size();++i) csv<<(i?",":"")<<columns[i]; csv<<'\n';
  for(Long64_t row=0;row<tree->GetEntries();++row) {
    rpc::Require(tree->GetEntry(row)>0,"ROOT read failed");
    for(size_t i=0;i<leaves.size();++i) csv<<(i?",":"")<<leaves[i]->GetValue(); csv<<'\n';
  }
  std::cout<<name<<" rows="<<tree->GetEntries()<<'\n';
}
int main(int argc,char** argv) try {
  if(argc==4 && std::string(argv[1])=="--hits") {
    TFile file(argv[2],"READ"); rpc::Require(!file.IsZombie(),"Cannot open hits file");
    Export(file,"rpcHits",argv[3],{"event","layer","window_index","position_index","n_positions",
      "window_n_tracks","window_t_ns","x_mm","y_mm","t_ns","centroid_x_mm","centroid_y_mm",
      "ionizing_edep_ev","n_segments","n_tracks","ambiguous","valid_coordinate",
      "ideal_spatial_model","unresolved_component","readout_pairing_validated"});
    return 0;
  }
  rpc::Require(argc==4 || argc==5,"Usage: rpc_export_comparison SOURCE.root BASELINE.root OUTPUT_DIRECTORY [LEGACY_ANALYSIS.root]");
  const std::string dir=argv[3]; std::filesystem::create_directories(dir);
  TFile source(argv[1],"READ"),baseline(argv[2],"READ");
  rpc::Require(!source.IsZombie()&&!baseline.IsZombie(),"Cannot open source/baseline");
  Export(source,"rpcSegments",dir+"/segments.csv",{"event","layer","track","parent","pid","root_track",
    "x0_mm","y0_mm","z0_mm","x1_mm","y1_mm","z1_mm","t0_ns","t1_ns","edep_ev","nonion_ev",
    "kinetic_ev","length_mm","gap_center_z_mm","gap_mm"});
  std::vector<std::string> primaryColumns={"event","pid","kinetic_ev","x_mm","y_mm","z_mm","dx","dy","dz","t_ns"};
  auto primary=dynamic_cast<TTree*>(source.Get("rpcPrimaries"));
  if(primary && primary->GetLeaf("source_proposals")) {
    primaryColumns.push_back("source_mode"); primaryColumns.push_back("source_proposals");
  }
  Export(source,"rpcPrimaries",dir+"/primaries.csv",primaryColumns);
  Export(source,"rpcLegacyEdeps",dir+"/legacy_edep.csv",{"event","layer","pid","process","energy_ev","x_mm","y_mm","legacy_trigger"});
  Export(baseline,"rpcHits",dir+"/baseline.csv",{"event","layer","x_mm","y_mm","t_ns","ionizing_edep_ev",
    "n_segments","n_tracks","ambiguous","valid_coordinate","nominal_txl_ns","nominal_txr_ns","nominal_tyd_ns","nominal_tyu_ns"});
  if(argc==5) {
    TFile legacy(argv[4],"READ"); rpc::Require(!legacy.IsZombie(),"Cannot read legacy analysis output");
    Export(legacy,"tree",dir+"/legacy_analysis_summary.csv",{"E","Ang","AngSmeared","CosThetaEdep","CosThetaSmeared"});
  }
  return 0;
} catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
