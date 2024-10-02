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

#include "G4LogicalVolume.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCuts.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
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

Cuts &Cuts::operator=(const G4LogicalVolume &volume)
{
  G4Material *material = volume.GetMaterial();
  G4ProductionCuts *cuts = volume.GetMaterialCutsCouple()->GetProductionCuts();

  GammaCut = cuts->GetProductionCut("gamma");
  ElectronCut = cuts->GetProductionCut("e-");
  PositronCut = cuts->GetProductionCut("e+");
  ProtonCut = cuts->GetProductionCut("proton");

  GammaThreshold = G4RToEConvForGamma().Convert(GammaCut, material);
  ElectronThreshold = G4RToEConvForElectron().Convert(ElectronCut, material);
  PositronThreshold = G4RToEConvForPositron().Convert(PositronCut, material);
  ProtonThreshold = G4RToEConvForProton().Convert(ProtonCut, material);

  return *this;
}

std::vector<Double_t> EdepData::fScoringZs;

Scatter &Scatter::operator=(const std::tuple<Double_t, Double_t, const G4Track *, const G4Track *, const G4Track *> &t)
{
  auto [prob, xs, muon, lp, ln] = t;
  Prob = prob;
  XS = xs;
  const G4Track *tracks[3] = { muon, lp, ln };
  for(size_t i = 0; i < 3; ++i) {
    Id[i] = tracks[i]->GetTrackID();
    Pid[i] = tracks[i]->GetParticleDefinition()->GetPDGEncoding();
    auto momentum = tracks[i]->GetMomentum();
    Px[i] = momentum.getX();
    Py[i] = momentum.getY();
    Pz[i] = momentum.getZ();
    E[i] = tracks[i]->GetTotalEnergy();
  }
  auto position = tracks[0]->GetPosition();
  X = position.getX();
  Y = position.getY();
  Z = position.getZ();
  T = tracks[0]->GetGlobalTime();
  return *this;
}
