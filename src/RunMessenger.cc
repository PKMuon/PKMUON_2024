// Origin: 2020.5.8 by siguang wang (siguang@pku.edu.cn PKU)

#include "RunMessenger.hh"

#include "G4RunManager.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "GpsPrimaryGeneratorAction.hh"
#include "MupTargetEnToLLPhysics.hh"
#include "Run.hh"

class RunMessenger::Driver {
public:
  Driver(RunMessenger *messenger);
  ~Driver();
  void SetNewValue(G4UIcommand *, G4String);

private:
  GpsPrimaryGeneratorAction *fGpsPrimaryGeneratorAction;

  G4UIdirectory *fScatterDir;
  G4UIcommand *fSetMupTargetEnToLLCmd;
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
  fGpsPrimaryGeneratorAction =
      (GpsPrimaryGeneratorAction *)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

  fScatterDir = new G4UIdirectory("/scatter/");
  fScatterDir->SetGuidance("Control scattering processes.");

  fSetMupTargetEnToLLCmd = new G4UIcommand("/scatter/mupTargetEnToLL", messenger);
  fSetMupTargetEnToLLCmd->SetParameter(new G4UIparameter("pid", 'i', false));
  fSetMupTargetEnToLLCmd->SetParameter(new G4UIparameter("points_file", 's', false));
  fSetMupTargetEnToLLCmd->SetParameter(new G4UIparameter("xssf", 'd', false));
  fSetMupTargetEnToLLCmd->SetGuidance("Configure MupTargetEnToLL process.");
  fSetMupTargetEnToLLCmd->AvailableForStates(G4State_Idle);

  fSetTotalEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gps/totalEnergy", messenger);
  fSetTotalEnergyCmd->SetGuidance("Set total energy.");
  fSetTotalEnergyCmd->SetParameterName("TotalEnergy", false);
  fSetTotalEnergyCmd->SetUnitCategory("Energy");
  fSetTotalEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunMessenger::Driver::~Driver()
{
  delete fSetTotalEnergyCmd;
  delete fSetMupTargetEnToLLCmd;
  delete fScatterDir;
}

void RunMessenger::Driver::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if(cmd == fSetTotalEnergyCmd) {
    fGpsPrimaryGeneratorAction->SetTotalEnergy(fSetTotalEnergyCmd->GetNewDoubleValue(val));
  } else if(cmd == fSetMupTargetEnToLLCmd) {
    G4Tokenizer next(val);
    G4String pid_s = next();
    G4String points_file = next();
    G4String xssf_s = next();
    G4String trailing = next();
    if(pid_s.empty() || points_file.empty() || xssf_s.empty() || !trailing.empty()) {
      throw std::runtime_error("expect 3 arguments");
    }
    G4int pid = stoi(pid_s);
    G4double xssf = stod(xssf_s);
    MupTargetEnToLLPhysics::GetInstance()->Configure(pid, points_file, xssf);
  }
}
