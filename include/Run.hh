// Origin: 2020.5.8 by Siguang WANG (siguang@pku.edu.cn)

#ifndef GEANT4_INTRODUCTION_RUN_HH
#define GEANT4_INTRODUCTION_RUN_HH 1

#include <Rtypes.h>

#include <map>
#include <vector>

#include "Object.hh"
#include "globals.hh"

class TFile;
class TTree;

class RunMessenger;
class PrimaryGeneratorAction;
class DetectorConstruction;
class G4Step;
class G4Track;

class Run {
public:
  static Run *GetInstance();
  static uint64_t GetThreadId();
  static uint64_t GetSeed();

  void SetRootName(G4String name) { fRootName = name; }

  void InitGeom();
  void InitTree();
  void SaveTree();
  void FillAndReset();
  void AutoSave();

  void AddStep(const G4Step *step);
  void AddTrack(const G4Track *step);

private:
  Run();
  ~Run();

  RunMessenger *fRunMessenger;
  PrimaryGeneratorAction *fPrimaryGeneratorAction;
  DetectorConstruction *fDetectorConstruction;
  G4String fRootName;
  TTree *fTree;
  TFile *fFile;
  G4double fScoringHalfX, fScoringHalfY, fScoringZ;
  std::vector<G4double> fScoringMaxZs;
  std::map<Long64_t, EdepData> fEdepData;
  std::vector<bool> fStatus;
};

#endif  // GEANT4_INTRODUCTION_RUN_H
