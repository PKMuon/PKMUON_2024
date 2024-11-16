#include "DMProcess.hh"

#include <TLorentzVector.h>
#include <math.h>

#include <G4MuonMinus.hh>
#include <G4ParticleChange.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include "Run.hh"

DMProcess::DMProcess()
    : G4VDiscreteProcess("DMProcess", fUserDefined), fRun(Run::GetInstance()), fDMMass(0.0), fDMXS(0.0), fMuonMass(0.0)
{
  // empty
}

DMProcess::~DMProcess() { }

G4ThreeVector DMProcess::RandDirection() const
{
  G4double cos_theta = 2 * G4UniformRand() - 1;
  G4double sin_theta = sqrt(1 - cos_theta * cos_theta);
  G4double phi = 2 * M_PI * G4UniformRand();
  G4double sin_phi, cos_phi;
  sincos(phi, &sin_phi, &cos_phi);
  return G4ThreeVector(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
}

G4ThreeVector DMProcess::RandDMMomentum() const
{
  G4double v = 220 * km / s, beta = v / CLHEP::c_light;
  G4double E = fDMMass / sqrt(1 - pow(beta, 2)), p = E * beta;
  return p * RandDirection();
}

void DMProcess::DoScatter(G4ThreeVector &momentum) const
{
  TLorentzVector pMu, pDM;
  pMu.SetPxPyPzE(momentum.x(), momentum.y(), momentum.z(), hypot(momentum.getR(), fMuonMass));
  G4ThreeVector dmMomentum = RandDMMomentum();
  pDM.SetPxPyPzE(dmMomentum.x(), dmMomentum.y(), dmMomentum.z(), hypot(dmMomentum.getR(), fDMMass));

  // Boost to the COM frame.
  TLorentzVector pCOM = pMu + pDM;
  TVector3 boost = pCOM.BoostVector();
  pMu.Boost(-boost);

  // Scatter.
  momentum = pMu.P() * RandDirection();
  pMu.SetPxPyPzE(momentum.x(), momentum.y(), momentum.z(), pMu.E());

  // Boost back to the lab frame.
  pMu.Boost(boost);

  momentum = G4ThreeVector(pMu.Px(), pMu.Py(), pMu.Pz());
}

G4double DMProcess::PostStepGetPhysicalInteractionLength(
    const G4Track &track, G4double previousStepSize, G4ForceCondition *condition)
{
  //if(fDMMass) G4cout << __FUNCTION__ << "(" << fDMMass << ", " << track.GetTrackID() << ")" << G4endl;
  G4double stepLength = GetMeanFreePath(track, previousStepSize, condition) * 0.001;
  //if(fDMMass) {
  //  G4cout << __FUNCTION__ << "(" << fDMMass << ", " << track.GetTrackID() << ") -> " << stepLength << G4endl;
  //}
  return stepLength;
}

G4VParticleChange *DMProcess::PostStepDoIt(const G4Track &track, const G4Step &step)
{
  //if(fDMMass) G4cout << __FUNCTION__ << "(" << fDMMass << ", " << track.GetTrackID() << ")" << G4endl;
  //bool changed = false;
  thread_local G4ParticleChange change;
  change.Initialize(track);
  do {
    if(!fDMMass) break;
    if(fabs(track.GetParticleDefinition()->GetPDGEncoding()) != 13) break;
    G4double xs = GetCrossSection(track.GetKineticEnergy(), track.GetMaterialCutsCouple());
    if(!(xs > 0)) break;

    G4double stepLength = step.GetStepLength();
    G4double logProbKeep = -xs * stepLength;
    G4double logRandom = log(G4UniformRand());
    //G4cout << __FUNCTION__ << "(" << fDMMass << ", " << track.GetTrackID() << "): stepLength=" << stepLength
    //       << " logProbKeep=" << std::fixed << std::setprecision(8) << logProbKeep << " logRandom=" << logRandom
    //       << std::defaultfloat << G4endl;
    if(logRandom < logProbKeep) break;

    //changed = true;
    G4ThreeVector momentum = track.GetMomentum();
    DoScatter(momentum);
    fRun->AddScatter(&track, momentum);

    G4double p = momentum.getR(), E = hypot(p, fMuonMass);
    change.ProposeMomentumDirection(momentum / p);
    change.ProposeEnergy(E - fMuonMass);
    change.ProposeVelocity(p / E * CLHEP::c_light);
  } while(0);
  //if(fDMMass) {
  //  G4cout << __FUNCTION__ << "(" << fDMMass << ", " << track.GetTrackID() << ") -> " << std::boolalpha << changed
  //         << std::noboolalpha << G4endl;
  //}
  return &change;
}

G4double DMProcess::GetCrossSection(
    const G4double energy [[maybe_unused]], const G4MaterialCutsCouple *couple [[maybe_unused]])
{
  //if(fDMMass) G4cout << __FUNCTION__ << "(" << fDMMass << ", " << energy << ")" << G4endl;
  G4double xs = fDMXS;
  //if(fDMMass) G4cout << __FUNCTION__ << "(" << fDMMass << ", " << energy << ") -> " << xs << G4endl;
  return xs;
}

G4double DMProcess::MinPrimaryEnergy(const G4ParticleDefinition *definition, const G4Material *)
{
  //if(fDMMass) G4cout << __FUNCTION__ << "(" << fDMMass << ", " << definition->GetPDGEncoding() << ")" << G4endl;
  G4double minPrimaryEnergy = INFINITY;
  if(fDMMass && fabs(definition->GetPDGEncoding()) == 13) { minPrimaryEnergy = 0.0; }
  //if(fDMMass) {
  //  G4cout << __FUNCTION__ << "(" << fDMMass << ", " << definition->GetPDGEncoding() << ") -> " << minPrimaryEnergy
  //         << G4endl;
  //}
  return minPrimaryEnergy;
}

void DMProcess::Configure(G4double dmMass, G4double dmXS)
{
  fDMMass = dmMass;
  fDMXS = dmXS * (0.3 * GeV / cm3 / dmMass);
  fMuonMass = G4MuonMinus::Definition()->GetPDGMass();
}

G4double DMProcess::GetMeanFreePath(
    const G4Track &track, [[maybe_unused]] G4double previousStepSize, G4ForceCondition *condition)
{
  //if(fDMMass) G4cout << __FUNCTION__ << "(" << fDMMass << ", " << track.GetTrackID() << ")" << G4endl;
  G4double mfp = INFINITY;
  if(!fDMMass) {
    *condition = InActivated;
  } else {
    if(fabs(track.GetParticleDefinition()->GetPDGEncoding()) == 13) {
      G4double xs = GetCrossSection(track.GetKineticEnergy(), track.GetMaterialCutsCouple());
      mfp = 1 / xs;
    }
    *condition = mfp == INFINITY ? NotForced : Forced;
  }
  //if(fDMMass) G4cout << __FUNCTION__ << "(" << fDMMass << ", " << track.GetTrackID() << ") -> " << mfp << G4endl;
  return mfp;
}
