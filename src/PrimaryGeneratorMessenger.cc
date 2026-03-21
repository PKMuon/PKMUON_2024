//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: PrimaryGeneratorMessenger.cc,v 1.8 2002/12/16 16:37:27 maire Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                          PrimaryGeneratorAction* Gun)
  :G4UImessenger(),
  fAction(Gun),
  fGunDir(nullptr),
  fMinEnergyCmd(nullptr),
  fMaxEnergyCmd(nullptr)
{
  fGunDir = new G4UIdirectory("/gun/");
  fGunDir->SetGuidance("Primary generator control");
   
  fMinEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/setMinEnergy", this);
  fMinEnergyCmd->SetGuidance("Set minimum energy for Gaisser-Tang generator");
  fMinEnergyCmd->SetParameterName("MinEnergy", false);
  fMinEnergyCmd->SetRange("MinEnergy>0");
  fMinEnergyCmd->SetUnitCategory("Energy");
  fMinEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fMaxEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/setMaxEnergy", this);
  fMaxEnergyCmd->SetGuidance("Set maximum energy for Gaisser-Tang generator");
  fMaxEnergyCmd->SetParameterName("MaxEnergy", false);
  fMaxEnergyCmd->SetRange("MaxEnergy>0");
  fMaxEnergyCmd->SetUnitCategory("Energy");
  fMaxEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fMinEnergyCmd;
  delete fMaxEnergyCmd;
  delete fGunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
    if (command == fMinEnergyCmd) {
        G4double value = fMinEnergyCmd->GetNewDoubleValue(newValue);
        G4cout << "========================================" << G4endl;
        G4cout << "Setting MinEnergy to: " << value/GeV << " GeV" << G4endl;
        fAction->SetMinEnergy(value / GeV);
        G4cout << "MinEnergy set to: " << fAction->GetMinEnergy() << " GeV" << G4endl;
        G4cout << "========================================" << G4endl;
    }
    if (command == fMaxEnergyCmd) {
        G4double value = fMaxEnergyCmd->GetNewDoubleValue(newValue);
        G4cout << "========================================" << G4endl;
        G4cout << "Setting MaxEnergy to: " << value/GeV << " GeV" << G4endl;
        fAction->SetMaxEnergy(value / GeV);
        G4cout << "MaxEnergy set to: " << fAction->GetMaxEnergy() << " GeV" << G4endl;
        G4cout << "========================================" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

