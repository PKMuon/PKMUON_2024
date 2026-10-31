#include <filesystem>
#include <iostream>
#include "../../analysis/src/analysis.C"
int main(int argc,char** argv) {
  if(argc!=3) { std::cerr<<"Usage: rpc_run_legacy_analysis INPUT.root OUTPUT.root\n"; return 2; }
  if(std::filesystem::exists(argv[2])) { std::cerr<<"Output already exists\n"; return 1; }
  // Force the existing no-LUT Edep fallback. The new kernel is a different schema.
  analysis(argv[1],argv[2],"__comparison_no_legacy_radius_LUT__.root");
  return std::filesystem::exists(argv[2])?0:1;
}
