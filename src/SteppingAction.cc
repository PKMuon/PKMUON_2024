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

#include "SteppingAction.hh"

#include "EventAction.hh"
#include "G4Step.hh"
#include "Randomize.hh"
#include "Run.hh"
#include "ScatterProcess.hh"

SteppingAction::SteppingAction(EventAction *eventAction)
    : fEventAction(eventAction), fScoringHalfX(0), fScoringHalfY(0), fScoringZ(0)
{
  memset(fScatterProbability, 0, sizeof fScatterProbability);
  memset(fScatterProcess, 0, sizeof fScatterProcess);
}

SteppingAction::~SteppingAction() { }

void SteppingAction::UserSteppingAction([[maybe_unused]] const G4Step *step)
{
  Run::GetInstance()->AddStep(step);

  G4Track *lp_track = step->GetTrack(), *ln_track;
  //G4cout << "Step of track " << lp_track->GetTrackID() << ": " << lp_track->GetPosition().getZ() << G4endl;
  if(lp_track->GetTrackID() == 1 && lp_track->GetPosition().getZ() >= fEventAction->GetScatterZ()) {
    double r = G4UniformRand(), s = 0;
    for(size_t i = 0; i < 2; ++i) {
      s += fScatterProbability[i];
      if(s <= r) continue;
      G4Track *muon_track = lp_track;
      double xs = fScatterProcess[i]->Scatter(lp_track, ln_track);
      Run::GetInstance()->AddScatter(fScatterProbability[i], xs, muon_track, lp_track, ln_track);
      break;
    }
  }
}

void SteppingAction::SetMupTargetEnToLL(size_t index, double probability, const char *points_file)
{
  if(fScatterProbability[1 - index] + probability > 1) throw std::logic_error("probability sum exceeds 1");
  fScatterProbability[index] = probability;
  fScatterProcess[index] = new MupTargetEnToLL(11 + index * 2, points_file);
}
