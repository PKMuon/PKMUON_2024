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

#include "TrackingAction.hh"

//#include "G4VProcess.hh"
#include "Run.hh"

TrackingAction::TrackingAction() { }

TrackingAction::~TrackingAction() { }

void TrackingAction::PreUserTrackingAction([[maybe_unused]] const G4Track *track)
{
  Run::GetInstance()->AddTrack(track);
  //if(const G4VProcess *process = track->GetCreatorProcess()) {
  //  G4cout << __PRETTY_FUNCTION__ << ": " << track->GetTrackID() << ": "
  //         << track->GetParticleDefinition()->GetParticleName() << ", " << process->GetProcessName() << G4endl;
  //}
}

void TrackingAction::PostUserTrackingAction([[maybe_unused]] const G4Track *track) { }
