#pragma once
#include <G4VDiscreteProcess.hh>
#include <G4ThreeVector.hh>

class Run;

class DMProcess : public G4VDiscreteProcess {
public:
  DMProcess();
  ~DMProcess() override;

  G4double PostStepGetPhysicalInteractionLength(
      const G4Track &track, G4double previousStepSize, G4ForceCondition *condition) override;
  G4VParticleChange *PostStepDoIt(const G4Track &, const G4Step &) override;
  G4double GetCrossSection(const G4double, const G4MaterialCutsCouple *) override;
  G4double MinPrimaryEnergy(const G4ParticleDefinition *, const G4Material *) override;

  // [NOTE] dmMass: GeV/MeV/... dmXS: cm2/barn/...
  void Configure(G4double dmMass, G4double dmXS);

protected:
  G4double GetMeanFreePath(const G4Track &aTrack, G4double previousStepSize, G4ForceCondition *condition) override;

private:
  Run *fRun;
  G4double fDMMass;  // DM particle mass
  G4double fDMXS;  // DM scattering cross section
  G4double fMuonMass;

  G4ThreeVector RandDirection() const;
  G4ThreeVector RandDMMomentum() const;
  void DoScatter(G4ThreeVector &momentum) const;
};
