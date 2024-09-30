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

#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYSetup.h"
#include "CRYUtils.h"
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

class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction(const char *filename);
  ~PrimaryGeneratorAction();
  void Initialize(const DetectorConstruction *);

public:
  void GeneratePrimaries(G4Event *anEvent);
  G4bool IsPrimary(G4int trackID) const { return trackID > 0 && trackID <= fNPrimary; }
  G4double GetHalfTimeWindow() const { return fHalfTimeWindow; }
  void InputCRY();
  void UpdateCRY(std::string *MessInput);
  void CRYFromFile(G4String newValue);

private:
  std::vector<CRYParticle *> *vect;  // vector of generated particles
  G4ParticleTable *particleTable;
  G4ParticleGun *particleGun;
  CRYGenerator *gen;
  PrimaryGeneratorMessenger *gunMessenger;
  G4int InputState;
  G4int fNPrimary;
  G4double fDetectorMinZ, fDetectorHalfX, fDetectorHalfY;
  G4double fHalfTimeWindow;
};

#endif
