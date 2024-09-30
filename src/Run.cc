// 2020.5.8 by siguang wang (siguang@pku.edu.cn)

#include "Run.hh"
#include "RunMessenger.hh"
#include <TFile.h>
#include <TTree.h>
#include <filesystem>
//#include <syscall.h>
#include <unistd.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

Run::Run()
{
  fRunMessenger = new RunMessenger(this);
  rootName = "CryMu.root";
  _tree = NULL;
  _file = NULL;
}

Run::~Run()
{
  SaveTree();
  delete fRunMessenger;
}

Run *Run::GetInstance()
{
  static Run run;
  return &run;
}

void Run::InitTree()
{
  using namespace std::filesystem;
  auto dirpath = path(rootName.c_str()).parent_path();
  if(!dirpath.empty()) create_directories(dirpath);

  _file = new TFile(rootName, "RECREATE");
  _tree = new TTree("T1", "Simple Out Tree");

  _tree->Branch("PDG",            &PDG,            "PDG/I"             );
  _tree->Branch("RpcTrkPx",       &RpcTrkPx,       "RpcTrkPx[16]/D"    );
  _tree->Branch("RpcTrkPy",       &RpcTrkPy,       "RpcTrkPy[16]/D"    );
  _tree->Branch("RpcTrkPz",       &RpcTrkPz,       "RpcTrkPz[16]/D"    );
  _tree->Branch("RpcTrkE",        &RpcTrkE,        "RpcTrkE[16]/D"     );
  _tree->Branch("RpcTrkEdep",     &RpcTrkEdep,     "RpcTrkEdep[16]/D"  );
  _tree->Branch("RpcTrkX",        &RpcTrkX,        "RpcTrkX[16]/D"     );
  _tree->Branch("RpcTrkY",        &RpcTrkY,        "RpcTrkY[16]/D"     );
  _tree->Branch("RpcTrkZ",        &RpcTrkZ,        "RpcTrkZ[16]/D"     );
  _tree->Branch("RpcTrkStatus",   &RpcTrkStatus,   "RpcTrkStatus[16]/O");
  _tree->Branch("RpcTrkComplete", &RpcTrkComplete, "RpcTrkComplete/O"  );
  _tree->Branch("RpcAllEdep",     &RpcAllEdep,     "RpcAllEdep[16]/D"  );
  _tree->Branch("RpcAllX",        &RpcAllX,        "RpcAllX[16]/D"     );
  _tree->Branch("RpcAllY",        &RpcAllY,        "RpcAllY[16]/D"     );
  _tree->Branch("RpcAllZ",        &RpcAllZ,        "RpcAllZ[16]/D"     );
  _tree->Branch("RpcAllN",        &RpcAllN,        "RpcAllN[16]/i"     );
  _tree->Branch("RpcAllComplete", &RpcAllComplete, "RpcAllComplete/O"  );

  Clear();
}

void Run::SaveTree()
{
  if(!_file) return;
  _file->cd();
  _tree->Write("T1", TObject::kOverwrite);
  _file->Close();
  _tree = NULL;
  _file = NULL;
}

void Run::Fill()
{
  RpcTrkComplete = true;
  for(int i = 0; i < 16; ++i) {
    double Edep = RpcTrkEdep[i];
    if(Edep) {
      RpcTrkPx  [i] /= Edep;
      RpcTrkPy  [i] /= Edep;
      RpcTrkPz  [i] /= Edep;
      RpcTrkE   [i] /= Edep;
      RpcTrkX   [i] /= Edep;
      RpcTrkY   [i] /= Edep;
      RpcTrkZ   [i] /= Edep;
    }
    if(!RpcTrkStatus[i]) RpcTrkComplete = false;
  }
  RpcAllComplete = true;
  for(int i = 0; i < 16; ++i) {
    double Edep = RpcAllEdep[i];
    if(Edep) {
      RpcAllX   [i] /= Edep;
      RpcAllY   [i] /= Edep;
      RpcAllZ   [i] /= Edep;
    }
    RpcAllN[i] = RpcAllIds[i].size();
    if(!RpcAllN[i]) RpcAllComplete = false;
  }
  _tree->Fill();
  Clear();
}

void Run::AutoSave()
{
  _tree->AutoSave("SaveSelf Overwrite");
}

void Run::AddRpcTrkInfo(int i, double Px, double Py, double Pz,
    double E, double Edep, double X, double Y, double Z, int pdgCode)
{
  if(i < 0 || i >= 16) return;
  RpcTrkPx    [i] += Px * Edep;
  RpcTrkPy    [i] += Py * Edep;
  RpcTrkPz    [i] += Pz * Edep;
  RpcTrkE     [i] += E  * Edep;
  RpcTrkEdep  [i] +=      Edep;
  RpcTrkX     [i] += X  * Edep;
  RpcTrkY     [i] += Y  * Edep;
  RpcTrkZ     [i] += Z  * Edep;
  RpcTrkStatus[i]  = true;
  PDG              = pdgCode;
}

void Run::AddRpcAllInfo(int i, int id, double Edep, double X, double Y, double Z)
{
  if(i < 0 || i >= 16) return;
  if(id != 1) {  // not primary
    auto [it, _] = RpcAllLayer.emplace(id, i);
    if(it->second != i) return;  // not belong to this layer
  }
  RpcAllIds   [i].insert(id);
  RpcAllEdep  [i] +=      Edep;
  RpcAllX     [i] += X  * Edep;
  RpcAllY     [i] += Y  * Edep;
  RpcAllZ     [i] += Z  * Edep;
}

void Run::Clear()
{
  for(int i = 0; i < 16; i++) {
    RpcTrkPx    [i] = 0;
    RpcTrkPy    [i] = 0;
    RpcTrkPz    [i] = 0;
    RpcTrkE     [i] = 0;
    RpcTrkEdep  [i] = 0;
    RpcTrkX     [i] = 0;
    RpcTrkY     [i] = 0;
    RpcTrkZ     [i] = 0;
    RpcTrkStatus[i] = false;
  }
  RpcTrkComplete = false;
  RpcAllLayer.clear();
  for(int i = 0; i < 16; i++) {
    RpcAllIds   [i].clear();
    RpcAllEdep  [i] = 0;
    RpcAllX     [i] = 0;
    RpcAllY     [i] = 0;
    RpcAllZ     [i] = 0;
    RpcAllN     [i] = 0;
  }
  RpcAllComplete = false;
  PDG = 0;
}

uint64_t Run::GetThreadId()
{
  return gettid();
}

uint64_t Run::GetSeed()
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::system_clock::now().time_since_epoch()
  ).count() + GetThreadId();
}
