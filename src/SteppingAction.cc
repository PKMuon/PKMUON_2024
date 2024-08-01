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
//
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia

#include "SteppingAction.hh"
#include "Run.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

SteppingAction::SteppingAction() : fScoringVolume(nullptr) { }
SteppingAction::~SteppingAction() { }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // Reject steps in non-scoring volumes.
  G4LogicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if(!fScoringVolume) {  // memoization (this term is not a typo)
    fScoringVolume = (
        (const DetectorConstruction *)G4RunManager::GetRunManager()->GetUserDetectorConstruction()
    )->GetScoringVolume();
    fScoringZRanges = (
        (const DetectorConstruction *)G4RunManager::GetRunManager()->GetUserDetectorConstruction()
    )->GetScoringZRanges();
    G4cout << "Scoring volume: " << fScoringVolume->GetName() << G4endl;
    G4cout << "Scoring ZRanges:" << G4endl;
    for(auto [zm, dz] : fScoringZRanges) {
      G4cout << " * " << zm << " +/- " << dz << " mm" << G4endl;
    }
  }
  if(volume != fScoringVolume) return;

  // Reject steps without energy deposition.
  G4double energy = aStep->GetTotalEnergyDeposit();
  G4double totalenergy = aStep->GetTrack()->GetTotalEnergy();
  if(!(energy > 0)) return;

  // Reject secondary tracks.
  G4int iTrkID = aStep->GetTrack()->GetTrackID();
  G4int iTrkparentID = aStep->GetTrack()->GetParentID();
  if(!(iTrkID==1 && iTrkparentID==0)) return;

  // Compute the hit point of this step.
  G4StepPoint* prePoint  = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  G4double x = (prePoint->GetPosition().x() + postPoint->GetPosition().x()) / 2.;
  G4double y = (prePoint->GetPosition().y() + postPoint->GetPosition().y()) / 2.;
  G4double z = (prePoint->GetPosition().z() + postPoint->GetPosition().z()) / 2.;

  // Compute the detector layer (if any) hit by this step.
  int igem = -1;
  for(int i = 0; i < (int)fScoringZRanges.size(); ++i) {
    if(std::fabs(z/mm - fScoringZRanges[i].first) <= fScoringZRanges[i].second) {
      igem = i; break;
    }
  }
  if(igem == -1) return;

  // Avoid duplicate hits of the same layer.
  if(!Run::GetInstance()->TestAndSetRpcTrkStatus(igem)) return;

  // Get momentum of the track.
  G4ThreeVector curDirection = aStep->GetPreStepPoint()->GetMomentumDirection();
  G4double px = curDirection.x();
  G4double py = curDirection.y();
  G4double pz = curDirection.z();

  // Record the hit info.
  Run::GetInstance()->SetRpcTrkInfo(igem, px, py, pz, totalenergy/MeV, energy/MeV, x/mm, y/mm, z/mm);
}
