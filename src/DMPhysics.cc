#include "DMPhysics.hh"

#include <G4MuonMinus.hh>
#include <G4MuonPlus.hh>
#include <G4ProcessManager.hh>

#include "DMProcess.hh"

DMPhysics *DMPhysics::fInstance = new DMPhysics;
thread_local DMProcess *DMPhysics::fProcess;

DMPhysics::DMPhysics() { }

DMPhysics::~DMPhysics() { delete fProcess; }

void DMPhysics::ConstructParticle()
{
  G4MuonPlus::Definition();
  G4MuonMinus::Definition();
}

void DMPhysics::ConstructProcess()
{
  G4ProcessManager *processManager = G4MuonPlus::Definition()->GetProcessManager();
  fProcess = new DMProcess;
  processManager->AddDiscreteProcess(fProcess);
}

void DMPhysics::Configure(G4double dmMass, G4double dmXS) { fProcess->Configure(dmMass, dmXS); }
