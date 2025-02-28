// Origin: 2020.5.8 by Siguang WANG (siguang@pku.edu.cn)

#ifndef GEANT4_INTRODUCTION_RUN_HH
#define GEANT4_INTRODUCTION_RUN_HH 1

#include <Rtypes.h>

#include <map>
#include <string>
#include <vector>

#include "EdepData.hh"
#include "globals.hh"

class TFile;
class TTree;

class RunMessenger;
class PrimaryGeneratorAction;
class DetectorConstruction;
class G4Step;
class G4Track;
class Event;

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
  void AddTrack(const G4Track *);
  void AddStep(const G4Step *);
  Event *GetEvent();

private:
  Run();
  ~Run();

  RunMessenger *fRunMessenger;
  PrimaryGeneratorAction *fPrimaryGeneratorAction;
  DetectorConstruction *fDetectorConstruction;
  G4String fRootName;
  TTree *fTree, *fParams;
  TFile *fFile;
  G4double fScoringHalfX, fScoringHalfY, fScoringZ;
  G4double fScoringOffsetX, fScoringOffsetY;
  std::vector<G4double> fScoringMaxZs;
  std::map<std::string, int> fProcessMap;
  std::map<EdepKey, EdepValue> fEdep;
  std::vector<bool> fStatus;
  G4long fIEvent;

  void BuildProcessMap();
};

#endif  // GEANT4_INTRODUCTION_RUN_H
