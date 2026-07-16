//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "Object.hh"

#include <math.h>

#include "DetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCuts.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

Track &Track::operator=(const G4Track &track)
{
  auto position = track.GetPosition();
  auto momentum = track.GetMomentum();

  Id = track.GetTrackID();
  Mother = track.GetParentID();
  Pid = track.GetParticleDefinition()->GetPDGEncoding();
  Px = momentum.getX();
  Py = momentum.getY();
  Pz = momentum.getZ();
  E = track.GetTotalEnergy();
  X = position.getX();
  Y = position.getY();
  Z = position.getZ();
  T = track.GetGlobalTime();

  return *this;
}

BoxVolume &BoxVolume::operator=(const std::tuple<const void *, const void *, const void *> &t)
{
  auto &hr = *(G4ThreeVector *)std::get<0>(t);
  auto &r = *(G4ThreeVector *)std::get<1>(t);
  auto rot = (G4RotationMatrix *)std::get<2>(t);
  HalfX = hr.x(), HalfY = hr.y(), HalfZ = hr.z();
  CenterX = r.x(), CenterY = r.y(), CenterZ = r.z();
  if(rot) {
    G4ThreeVector axis = rot->axis();
    Theta = axis.theta(), Phi = axis.phi(), Alpha = rot->delta();
  } else {
    Theta = Phi = Alpha = 0.0;
  }
  return *this;
}

Bool_t BoxVolume::Test(const G4Track *track) const
{
  G4ThreeVector axis;
  axis.setRThetaPhi(1.0, Theta, Phi);
  G4RotationMatrix rotation;
  rotation.rotate(Alpha, axis);
  auto r = rotation.inverse() * track->GetPosition();
  if(r.x() < CenterX - HalfX || r.x() > CenterX + HalfX) return kFALSE;
  if(r.y() < CenterY - HalfY || r.y() > CenterY + HalfY) return kFALSE;
  if(r.z() < CenterZ - HalfZ || r.z() > CenterZ + HalfZ) return kFALSE;
  return kTRUE;
}

Params &Params::operator=(const DetectorConstruction &detectorConstruction)
{
  //const G4MaterialCutsCouple *couple = detectorConstruction.GetScoringVolume()->GetMaterialCutsCouple();
  //const G4Material *material = couple->GetMaterial();
  //G4ProductionCuts *cuts = couple->GetProductionCuts();

  //GammaCut = cuts->GetProductionCut("gamma");
  //ElectronCut = cuts->GetProductionCut("e-");
  //PositronCut = cuts->GetProductionCut("e+");
  //ProtonCut = cuts->GetProductionCut("proton");
  GammaCut = 0.0;
  ElectronCut = 0.0;
  PositronCut = 0.0;
  ProtonCut = 0.0;

  //GammaThreshold = G4RToEConvForGamma().Convert(GammaCut, material);
  //ElectronThreshold = G4RToEConvForElectron().Convert(ElectronCut, material);
  //PositronThreshold = G4RToEConvForPositron().Convert(PositronCut, material);
  //ProtonThreshold = G4RToEConvForProton().Convert(ProtonCut, material);
  GammaThreshold = 0.0;
  ElectronThreshold = 0.0;
  PositronThreshold = 0.0;
  ProtonThreshold = 0.0;

  ScoringVolumes = detectorConstruction.GetScoringVolumes();
  return *this;
}

Scatter &Scatter::operator=(const std::tuple<const G4Track *, const G4DynamicParticle *, const G4DynamicParticle *> &t)
{
  auto [muon, lp, ln] = t;
  Id = muon->GetTrackID();
  const G4DynamicParticle *particles[3] = { muon->GetDynamicParticle(), lp, ln };
  for(size_t i = 0; i < 3; ++i) {
    Pid[i] = particles[i]->GetParticleDefinition()->GetPDGEncoding();
    auto momentum = particles[i]->GetMomentum();
    Px[i] = momentum.getX();
    Py[i] = momentum.getY();
    Pz[i] = momentum.getZ();
    E[i] = particles[i]->GetTotalEnergy();
  }
  auto position = muon->GetPosition();
  X = position.getX();
  Y = position.getY();
  Z = position.getZ();
  T = muon->GetGlobalTime();
  return *this;
}
