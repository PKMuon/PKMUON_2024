#pragma once
#include "G4VPhysicsConstructor.hh"

class DMProcess;

class DMPhysics : public G4VPhysicsConstructor {
public:
  void ConstructParticle() override;
  void ConstructProcess() override;

  void Configure(G4double dmMass, G4double dmXS);
  static DMPhysics *GetInstance() { return fInstance; }

private:
  DMPhysics();
  ~DMPhysics() override;

  static DMPhysics *fInstance;
  static thread_local DMProcess *fProcess;
};
