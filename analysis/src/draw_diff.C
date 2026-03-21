#include <TChain.h>
#include <TCanvas.h>
#include <filesystem>
#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>

using namespace std;
using namespace filesystem;

void draw_diff()
{
  vector<string> rootfiles;
  for(const auto &entry : directory_iterator("../build/root_file")) {
    if(entry.is_directory()) continue;
    string path = entry.path();
    if(regex_match(path, regex("^.*CryMu_\\d+\\.root$"))) {
      rootfiles.emplace_back(std::move(path));
    }
  }
  sort(rootfiles.begin(), rootfiles.end());

  TChain *chain = new TChain("T1");
  for(const string &rootfile : rootfiles) {
    cout << rootfile << endl;
    chain->Add(rootfile.c_str());
  }

  TCanvas *canvas = new TCanvas("diff_x", "diff_x", 4000, 3000);
  canvas->Divide(2, 2);
  for(int i = 0; i < 4; ++i) {
    cout << canvas->GetName() << ": " << i << endl;
    canvas->cd(i + 1);
    chain->Draw(("("
          "(RpcAllX[" + to_string(4 * i + 2) + "] + RpcAllX[" + to_string(4 * i + 3) + "])"
          " - "
          "(RpcTrkX[" + to_string(4 * i + 2) + "] + RpcTrkX[" + to_string(4 * i + 3) + "])"
          ") / 2").c_str(), "RpcTrkComplete");
  }
  canvas->SaveAs((canvas->GetName() + ".pdf"s).c_str());
  delete canvas;

  canvas = new TCanvas("diff_y", "diff_y", 4000, 3000);
  canvas->Divide(2, 2);
  for(int i = 0; i < 4; ++i) {
    cout << canvas->GetName() << ": " << i << endl;
    canvas->cd(i + 1);
    chain->Draw(("("
          "(RpcAllY[" + to_string(4 * i) + "] + RpcAllY[" + to_string(4 * i + 1) + "])"
          " - "
          "(RpcTrkY[" + to_string(4 * i) + "] + RpcTrkY[" + to_string(4 * i + 1) + "])"
          ") / 2").c_str(), "RpcTrkComplete");
  }
  canvas->SaveAs((canvas->GetName() + ".pdf"s).c_str());
  delete canvas;

  delete chain;
}
