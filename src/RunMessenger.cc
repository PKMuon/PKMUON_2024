// Origin: 2020.5.8 by siguang wang (siguang@pku.edu.cn PKU)

#include "RunMessenger.hh"

#include "DMPhysics.hh"
#include "G4RunManager.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UnitsTable.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"

class RunMessenger::Driver {
public:
  Driver(RunMessenger *messenger);
  ~Driver();
  void SetNewValue(G4UIcommand *, G4String);

private:
  PrimaryGeneratorAction *fPrimaryGeneratorAction;

  G4UIdirectory *fScatterDir;
  G4UIcommand *fSetDMCmd;
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

  fScatterDir = new G4UIdirectory("/scatter/");
  fScatterDir->SetGuidance("Control scattering processes.");

  fSetDMCmd = new G4UIcommand("/scatter/DM", messenger);
  fSetDMCmd->SetParameter(new G4UIparameter("mass", 'd', false));
  fSetDMCmd->SetParameter(new G4UIparameter("mass_unit", 's', false));
  fSetDMCmd->SetParameter(new G4UIparameter("xs", 'd', false));
  fSetDMCmd->SetParameter(new G4UIparameter("xs_unit", 's', false));
  fSetDMCmd->SetGuidance("Configure DM process.");
  fSetDMCmd->AvailableForStates(G4State_Idle);

  fSetTotalEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gun/totalEnergy", messenger);
  fSetTotalEnergyCmd->SetGuidance("Set total energy.");
  fSetTotalEnergyCmd->SetParameterName("TotalEnergy", false);
  fSetTotalEnergyCmd->SetUnitCategory("Energy");
  fSetTotalEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunMessenger::Driver::~Driver()
{
  delete fSetTotalEnergyCmd;
  delete fSetDMCmd;
  delete fScatterDir;
}

void RunMessenger::Driver::SetNewValue(G4UIcommand *cmd, G4String val)
{
  if(cmd == fSetTotalEnergyCmd) {
    throw std::runtime_error("setting total energy not compatible with CRY");
  } else if(cmd == fSetDMCmd) {
    G4Tokenizer next(val);
    G4String mass_s = next();
    G4String mass_unit_s = next();
    G4String xs_s = next();
    G4String xs_unit_s = next();
    G4String trailing = next();
    if(mass_s.empty() || mass_unit_s.empty() || xs_s.empty() || xs_unit_s.empty() || !trailing.empty()) {
      throw std::runtime_error("expect 4 arguments");
    }
    G4double mass = stod(mass_s) * G4UnitDefinition::GetValueOf(mass_unit_s);
    G4double xs = stod(xs_s) * G4UnitDefinition::GetValueOf(xs_unit_s);
    DMPhysics::GetInstance()->Configure(mass, xs);
  }
}
