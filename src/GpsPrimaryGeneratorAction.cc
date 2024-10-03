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

#include "GpsPrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SPSPosDistribution.hh"

GpsPrimaryGeneratorAction::GpsPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fGeneralParticleSource(nullptr)
{
  fGeneralParticleSource = new G4GeneralParticleSource();
}

GpsPrimaryGeneratorAction::~GpsPrimaryGeneratorAction() { delete fGeneralParticleSource; }

void GpsPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  fGeneralParticleSource->GeneratePrimaryVertex(anEvent);
}

void GpsPrimaryGeneratorAction::Initialize(DetectorConstruction *detectorConstruction)
{
  auto posDist = fGeneralParticleSource->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType("Plane");
  posDist->SetPosDisShape("Square");
  posDist->SetCentreCoords({ 0, 0, detectorConstruction->GetDetectorMinZ() });
  posDist->SetHalfX(detectorConstruction->GetDetectorHalfX());
  posDist->SetHalfY(detectorConstruction->GetDetectorHalfY());
}

void GpsPrimaryGeneratorAction::SetTotalEnergy(G4double energy)
{
  auto eneDist = fGeneralParticleSource->GetCurrentSource()->GetEneDist();
  energy -= fGeneralParticleSource->GetCurrentSource()->GetParticleDefinition()->GetPDGMass();
  if(energy < 0) throw std::invalid_argument("energy less than mass");
  eneDist->SetMonoEnergy(energy);
}
