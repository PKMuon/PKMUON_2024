//******************************************************************************
// PrimaryGeneratorAction.hh
//
// This class is a class derived from G4VUserPrimaryGeneratorAction for
// constructing the process used to generate incident particles.
//
// 1.00 JMV, LLNL, JAN-2007:  First version.
//******************************************************************************
//
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4DataVector.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "RNGWrapper.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "vector"

#include "GaisserTangGenerator.hh"

class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction();
  void Initialize(const DetectorConstruction *);

public:
  void GeneratePrimaries(G4Event *anEvent);
  G4bool IsPrimary(G4int trackID) const { return trackID > 0 && trackID <= fNPrimary; }

  void SetMinEnergy(G4double e) { if (fGaisserTangGen) fGaisserTangGen->SetMinEnergy(e); }
  void SetMaxEnergy(G4double e) { if (fGaisserTangGen) fGaisserTangGen->SetMaxEnergy(e); }
  
  G4double GetMinEnergy() const { return fGaisserTangGen ? fGaisserTangGen->GetMinEnergy() : 0.0; }
  G4double GetMaxEnergy() const { return fGaisserTangGen ? fGaisserTangGen->GetMaxEnergy() : 0.0; }

private:
  G4ParticleTable *particleTable;
  G4ParticleGun *particleGun;
  PrimaryGeneratorMessenger *gunMessenger;
  GaisserTangGenerator* fGaisserTangGen;
  G4int InputState;
  G4int fNPrimary;
  G4double fDetectorMaxZ, fDetectorHalfX, fDetectorHalfY;
};

#endif
