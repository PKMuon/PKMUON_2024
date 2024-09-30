
//
// 2020.5.8 by siguang wang (siguang@pku.edu.cn  PKU)
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionMessenger.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* primary)
  : G4UImessenger(),  fPrimaryGeneratorAction(primary)
{
  fFileNameDir = new G4UIdirectory("/primary/");
  fFileNameDir->SetGuidance("provide the input particle file full name");
          
  fSetInputFileNameCmd = new G4UIcmdWithAString("/primary/SetFileName",this);
  fSetInputFileNameCmd->SetGuidance("input the primary particle info");
  fSetInputFileNameCmd->SetParameterName("fileName",true);
  fSetInputFileNameCmd->SetDefaultValue ("xxx.txt");
  fSetInputFileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
  delete fSetInputFileNameCmd; 
  delete fFileNameDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
                 
  if (command == fSetInputFileNameCmd) {
    G4cout << "\n---> input primary particle file name: " << newValues << G4endl;
    fPrimaryGeneratorAction->SetInputFullName(newValues);
  }   

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

   
