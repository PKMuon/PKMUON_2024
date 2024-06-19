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
    G4LogicalVolume*  fScoringVolume4;

   private:
    void printDaughters(const G4LogicalVolume* mother);
    /*
    double Z1 = -513.34;
    double Z2 = -500.;
    double Z3 = 513.24;
    double Z4 = 526.58;
    double deltaZ = 0.1;
    */
     G4double Z1 = -459.5;
     G4double Z2 = -451.9; 
     G4double Z3 = -445.8;
     G4double Z4 = -438.2;
     G4double Z5 = -259.5;
     G4double Z6 = -251.9;
     G4double Z7 = -245.8;
     G4double Z8 = -238.2;
     G4double Z9 = 240.5;
     G4double Z10 = 248.1;
     G4double Z11 = 254.2;
     G4double Z12 = 261.8;
     G4double Z13 = 440.5;
     G4double Z14 = 448.1;
     G4double Z15 = 454.2;
     G4double Z16 = 461.8;

     G4double deltaZ = 0.1;

};
#endif
