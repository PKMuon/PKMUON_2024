// Origin: 2020.5.8 by siguang wang (siguang@pku.edu.cn PKU)

#include "RunMessenger.hh"

#include "G4RunManager.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"
#include "SteppingAction.hh"

class RunMessenger::Driver {
public:
  Driver(RunMessenger *messenger);
  ~Driver();
  void SetNewValue(G4UIcommand *, G4String);

private:
  PrimaryGeneratorAction *fPrimaryGeneratorAction;
  SteppingAction *fSteppingAction;
  G4UIdirectory *fScatterDir;
  G4UIcommand *fSetMupTargetEnToEECmd;
  G4UIcommand *fSetMupTargetEnToMuMuCmd;
  G4UIcmdWithADoubleAndUnit *fSetTotalEnergyCmd;
};

RunMessenger::RunMessenger(Run *run) : G4UImessenger(), fRun(run)
{
  fFileNameDir = new G4UIdirectory("/rlt/");
  fFileNameDir->SetGuidance("Interact with ROOT library.");

  fSetFileNameCmd = new G4UIcmdWithAString("/rlt/SetFileName", this);
  fSetFileNameCmd->SetGuidance("Set output pathname.");
  fSetFileNameCmd->SetParameterName("fileName", true);
  fSetFileNameCmd->SetDefaultValue("rlt.root");
  fSetFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDriver = new Driver(this);
}

RunMessenger::~RunMessenger()
{
  delete fSetFileNameCmd;
  delete fFileNameDir;
  delete fDriver;
}

void RunMessenger::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if(cmd == fSetFileNameCmd) {
    G4cout << "\n---> root name from file: " << val << G4endl;
    fRun->SetRootName(val);
  } else {
    fDriver->SetNewValue(cmd, val);
  }
}

RunMessenger::Driver::Driver(RunMessenger *messenger)
{
  fPrimaryGeneratorAction = (PrimaryGeneratorAction *)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  fSteppingAction = (SteppingAction *)G4RunManager::GetRunManager()->GetUserSteppingAction();

  fScatterDir = new G4UIdirectory("/scatter");
  fScatterDir->SetGuidance("Control scattering processes.");

  fSetMupTargetEnToEECmd = new G4UIcommand("/scatter/mupTargetEnToEE", messenger);
  fSetMupTargetEnToEECmd->SetParameter(new G4UIparameter("probability", 'd', false));
  fSetMupTargetEnToEECmd->GetParameter(0)->SetParameterRange("probability >= 0.0 && probability <= 1.0");
  fSetMupTargetEnToEECmd->SetParameter(new G4UIparameter("points_file", 's', false));
  fSetMupTargetEnToEECmd->SetGuidance("Configure MupTargetEnToEE process.");
  fSetMupTargetEnToEECmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSetMupTargetEnToMuMuCmd = new G4UIcommand("/scatter/mupTargetEnToMuMu", messenger);
  fSetMupTargetEnToMuMuCmd->SetParameter(new G4UIparameter("probability", 'd', false));
  fSetMupTargetEnToMuMuCmd->GetParameter(0)->SetParameterRange("probability >= 0.0 && probability <= 1.0");
  fSetMupTargetEnToMuMuCmd->SetParameter(new G4UIparameter("points_file", 's', false));
  fSetMupTargetEnToMuMuCmd->SetGuidance("Configure MupTargetEnToMuMu process.");
  fSetMupTargetEnToMuMuCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSetTotalEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/totalEnergy", messenger);
  fSetTotalEnergyCmd->SetGuidance("Set total energy.");
  fSetTotalEnergyCmd->SetParameterName("TotalEnergy", false);
  fSetTotalEnergyCmd->SetUnitCategory("Energy");
  fSetTotalEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunMessenger::Driver::~Driver()
{
  delete fSetTotalEnergyCmd;
  delete fSetMupTargetEnToMuMuCmd;
  delete fSetMupTargetEnToEECmd;
  delete fScatterDir;
}

void RunMessenger::Driver::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if(cmd == fSetTotalEnergyCmd) {
    //fPrimaryGeneratorAction->SetTotalEnergy(fSetTotalEnergyCmd->GetNewDoubleValue(val));
    throw std::runtime_error("setting total energy not compatible with CRY");
  } else if(cmd == fSetMupTargetEnToEECmd || cmd == fSetMupTargetEnToMuMuCmd) {
    G4Tokenizer next(val);
    G4String probability_s = next();
    G4String points_file = next();
    G4String trailing = next();
    if(probability_s.empty() || points_file.empty() || !trailing.empty()) {
      throw std::runtime_error("expect 2 arguments");
    }
    G4double probability = std::stod(probability_s);
    if(cmd == fSetMupTargetEnToEECmd) {
      fSteppingAction->SetMupTargetEnToEE(probability, points_file);
    } else {
      fSteppingAction->SetMupTargetEnToMuMu(probability, points_file);
    }
  }
}
