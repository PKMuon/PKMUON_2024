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

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G4AutoDelete.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4RunManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "MupTargetEnToLLPhysics.hh"
#include "Run.hh"

int main(int argc, char **argv)
{
  // Detect interactive mode (if no arguments) and define UI session.
  G4UIExecutive *ui = NULL;
  if(argc == 1) ui = new G4UIExecutive(argc, argv);

  // Choose the random engine.
  auto engine = new CLHEP::RanecuEngine;
  G4Random::setTheEngine(engine);  // ownership kept
  G4AutoDelete::Register(engine);
  uint64_t seed = Run::GetSeed();
  G4cout << "seed=" << seed << G4endl;
  G4Random::setTheSeed(seed);

  G4RunManager *runManager = new G4RunManager;

  // Set mandatory initialization classes.
  runManager->SetUserInitialization(new DetectorConstruction);
  G4VModularPhysicsList *physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  physicsList->RegisterPhysics(MupTargetEnToLLPhysics::GetInstance());
  runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new ActionInitialization);

  runManager->Initialize();

  // Initialize visualization.
  G4VisManager *visManager = nullptr;

  // Get the pointer to the User Interface manager.
  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  if(ui) {  // interactive mode
    visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  } else {  // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  delete visManager;
  delete runManager;
}
