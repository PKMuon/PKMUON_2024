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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli,University of Wollongong, Australia
// Contributions by F. Ambroglini INFN Perugia, Italy
//

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include <map>
class G4LogicalVolume;
class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
   ~SteppingAction();
    void UserSteppingAction(const G4Step*);

   private:
    G4LogicalVolume*  fScoringVolume;
    G4LogicalVolume*  fScoringVolume2;
    G4LogicalVolume*  fScoringVolume3;

   private:
    void printDaughters(const G4LogicalVolume* mother);
    double Z1 = 0;
    double Z2 = 13.34;
    double Z3 = 13.44;
    double deltaZ = 0.1;
    double deltaZ3 = 0.01;

};
#endif
