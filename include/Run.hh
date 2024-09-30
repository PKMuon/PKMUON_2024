// Origin: 2020.5.8 by Siguang WANG (siguang@pku.edu.cn)

#ifndef GEANT4_INTRODUCTION_RUN_HH
#define GEANT4_INTRODUCTION_RUN_HH 1

#include "globals.hh"
#include <Rtypes.h>
#include <stdint.h>
#include <unordered_set>
#include <unordered_map>

class TFile;
class TTree;

class RunMessenger;

class Run {
public:
  static Run *GetInstance();
  static uint64_t GetThreadId();
  static uint64_t GetSeed();

  void SetRootName(G4String name) { rootName = name; }

  void InitTree();
  void SaveTree();
  void Fill();
  void AutoSave();

  void AddRpcTrkInfo(int i, double Px, double Py, double Pz, double E,
      double Edep, double X, double Y, double Z, int pdgCode);
  void AddRpcAllInfo(int i, int id, double Edep, double X, double Y, double Z);

private:
  Run();
  ~Run();

  RunMessenger *fRunMessenger;
  G4String rootName;
  TTree *_tree;
  TFile *_file;

  // Altered by other routines.
  Double_t RpcTrkPx[16];
  Double_t RpcTrkPy[16];
  Double_t RpcTrkPz[16];
  Double_t RpcTrkE[16];
  Double_t RpcTrkEdep[16];
  Double_t RpcTrkX[16];
  Double_t RpcTrkY[16];
  Double_t RpcTrkZ[16];
  bool RpcTrkStatus[16];
  std::unordered_map<int, int> RpcAllLayer;
  std::unordered_set<int> RpcAllIds[16];
  Double_t RpcAllEdep[16];
  Double_t RpcAllX[16];
  Double_t RpcAllY[16];
  Double_t RpcAllZ[16];
  Int_t PDG;

  // Mantained by us.
  Bool_t RpcTrkComplete;
  UInt_t RpcAllN[16];
  Bool_t RpcAllComplete;

  void Clear();
};

#endif  // GEANT4_INTRODUCTION_RUN_H
