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

void draw_momentum()
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

  TCanvas *canvas = new TCanvas("px", "px", 4000, 3000);
  canvas->Divide(2, 2);
  for(int i = 0; i < 4; ++i) {
    cout << canvas->GetName() << ": " << i << endl;
    canvas->cd(i + 1);
    chain->Draw(("("
          "(RpcTrkPx[" + to_string(4 * i + 2) + "] + RpcTrkPx[" + to_string(4 * i + 3) + "])"
          ") / 2").c_str(), "RpcTrkComplete");
  }
  canvas->SaveAs((canvas->GetName() + ".pdf"s).c_str());
  delete canvas;

  canvas = new TCanvas("py", "py", 4000, 3000);
  canvas->Divide(2, 2);
  for(int i = 0; i < 4; ++i) {
    cout << canvas->GetName() << ": " << i << endl;
    canvas->cd(i + 1);
    chain->Draw(("("
          "(RpcTrkPy[" + to_string(4 * i) + "] + RpcTrkPy[" + to_string(4 * i + 1) + "])"
          ") / 2").c_str(), "RpcTrkComplete");
  }
  canvas->SaveAs((canvas->GetName() + ".pdf"s).c_str());
  delete canvas;

  canvas = new TCanvas("pz", "pz", 4000, 3000);
  canvas->Divide(2, 2);
  for(int i = 0; i < 4; ++i) {
    cout << canvas->GetName() << ": " << i << endl;
    canvas->cd(i + 1);
    chain->Draw(("("
          "(RpcTrkPz[" + to_string(4 * i) + "] + RpcTrkPz[" + to_string(4 * i + 1) + "])"
          " + "
          "(RpcTrkPz[" + to_string(4 * i + 2) + "] + RpcTrkPz[" + to_string(4 * i + 3) + "])"
          ") / 4").c_str(), "RpcTrkComplete");
  }
  canvas->SaveAs((canvas->GetName() + ".pdf"s).c_str());
  delete canvas;

  canvas = new TCanvas("E", "E", 4000, 3000);
  canvas->Divide(2, 2);
  for(int i = 0; i < 4; ++i) {
    cout << canvas->GetName() << ": " << i << endl;
    canvas->cd(i + 1);
    chain->Draw(("("
          "(RpcTrkE[" + to_string(4 * i) + "] + RpcTrkE[" + to_string(4 * i + 1) + "])"
          " + "
          "(RpcTrkE[" + to_string(4 * i + 2) + "] + RpcTrkE[" + to_string(4 * i + 3) + "])"
          ") / 4").c_str(), "RpcTrkComplete");
  }
  canvas->SaveAs((canvas->GetName() + ".pdf"s).c_str());
  delete canvas;

  delete chain;
}
