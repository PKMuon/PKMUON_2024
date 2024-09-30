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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicalConstants.hh"  //twopi
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4GenericTrap.hh"
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),fScoringVolume(nullptr),fScoringVolume2(nullptr)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
     delete vacum;
     delete air;
     delete pb;
     delete cu;
     delete cuLess;
     delete pbMore;
     delete kapton;
     delete kaptonLess;
     delete gasMixture;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4cout << "Defining the materials" << G4endl;
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  // Build materials
  vacum   = nistManager->FindOrBuildMaterial("G4_Galactic");
  air   = nistManager->FindOrBuildMaterial("G4_AIR");
  pb   = nistManager->FindOrBuildMaterial("G4_Pb");
  fe   = nistManager->FindOrBuildMaterial("G4_Fe");
  w   = nistManager->FindOrBuildMaterial("G4_W");
  cu   = nistManager->FindOrBuildMaterial("G4_Cu");

  G4cout<<"define high density Pb"<<G4endl;
  // define high density Pb
  G4int nComponent;
  G4double fracmass=1.0;
  pbMore   = new G4Material("pbMore",11.34*1*g/cm3,nComponent=1);
  pbMore -> AddMaterial(pb,fracmass);

  G4cout<<"define low density Cu"<<G4endl;
  // define low density Cu
  G4double ratio_outter = 0.83342659; // 1-(2700*pi)/(29400*sqrt(3)) 
  G4double density = ratio_outter * 8.94 * g/cm3;
  cuLess   = new G4Material("cuLess",density,nComponent=1);
  cuLess -> AddMaterial(cu,fracmass);

  G4cout<<"define Kapton"<<G4endl;
  // Define Kapton material
  density = 1.42 * g/cm3; // Density of Kapton
  kapton = new G4Material("Kapton", density, 3);
  kapton->AddElement(nistManager->FindOrBuildElement("H"), 0.0273);
  kapton->AddElement(nistManager->FindOrBuildElement("C"), 0.7213);
  kapton->AddElement(nistManager->FindOrBuildElement("O"), 0.2514);

  G4double ratio_inner = 0.85887530; // 1-(2700*pi+1875*pi)/2/(29400*sqrt(3))
  density = ratio_inner * 1.42 * g/cm3; // Density of Kapton
  kaptonLess = new G4Material("KaptonLess", density, 3);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("H"), 0.0273);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("C"), 0.7213);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("O"), 0.2514);

  G4cout<<"define mixgass"<<G4endl;
  // Define 70%Ar and 30%CO2 gas mixture
  density = 1.822 * mg/cm3; // Density of the gas mixture
  gasMixture = new G4Material("ArCO2GasMixture", density, 3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("Ar"), 0.7);
  gasMixture->AddElement(nistManager->FindOrBuildElement("C"), 0.3/3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("O"), 0.3/3*2);

  G4cout<<"define detector material"<<G4endl;
  // Define detector material
  Drift_cathode_Mat = cu;
  Gem_inner_Mat = kaptonLess;
  Gem_outter_Mat = cuLess;
  Readout_plate_Mat = cu;
  Shell_Mat = kapton;
  Readout_bar_Mat = cu;
  Gem_Mat = gasMixture;
  world_Mat = vacum;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineConstants()
{
  Gem_outter_x=1*m;
  Gem_outter_y=1*m;
  Gem_outter_z=5*um;

  Gem_inner_x=1*m;
  Gem_inner_y=1*m;
  Gem_inner_z=50*um;

  drift_cathode_x=1*m;
  drift_cathode_y=1*m;
  drift_cathode_z=0.1*mm;

  readout_bar_gap_x=210*um;
  readout_bar_x=150*um;
  readout_bar_y=1*m;
  readout_bar_z=0.1*mm;

  readout_plate_x=1*m;
  readout_plate_y=1*m;
  readout_plate_z=0.1*mm;

  shell_x=Gem_outter_x;
  shell_y=Gem_outter_y;
  shell_z=100*um;

  gap1 = 4.8*mm;
  gap2 = 2*mm;
  Gem_x=1*m;
  Gem_y=1*m;
  Gem_z=
        gap1
        +num_Gem_inner*gap2
        +(num_Gem_outter*Gem_outter_z+num_Gem_inner*Gem_inner_z)
        +drift_cathode_z
        +readout_bar_z
        +readout_plate_z;

  Box_x=1*m;
  Box_y=1*m;
  Box_z=1*m;

  num_Gem=2;

  experimentalHall_x=1.1*Box_x;
  experimentalHall_y=1.1*Box_y;
  experimentalHall_z=1.1*(Gem_z*num_Gem);
 
  Pbplate_x=100*mm;
  Pbplate_y=100*mm;
  Pbplate_z=0.1*mm;

  Pbbox_x=0.5*m;
  Pbbox_y=0.5*m;
  //Pbbox_z=0.5*m;
  //Pbbox_z=20*mm;
  Pbbox_z=1*mm;

/*
  Pbbox_x=20*mm;
  Pbbox_z=0.5*m;
  Pbbox_y=0.5*m;
*/

/*
  Pbbox_x=0.3*m;
  Pbbox_y=0.3*m;
  Pbbox_z=0.3*m;
*/

  Febox_x=0.3*m;
  Febox_y=0.3*m;
  Febox_z=0.3*m;

  Wbox_x=0.3*m;
  Wbox_y=0.3*m;
  Wbox_z=0.3*m;

  Vacbox_x=0.5*m-40*mm;
  Vacbox_y=0.5*m-40*mm;
  Vacbox_z=0.5*m-40*mm;

  pku_box_x = 30*cm;
  pku_box_y = 40*cm;
  pku_box_z = 3*cm;

  pku_bar_x = 3*cm;
  pku_bar_y1 = 20*cm;
  pku_bar_y2 = 40*cm;
  pku_bar_z = 3*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{

    /*******************************
   * Gem outter plate       *
   *******************************/
  G4VSolid* Gem_outter_box
    = new G4Box("Gem_outter_box",             // World Volume
                Gem_outter_x/2,        // x size
                Gem_outter_y/2,        // y size
                Gem_outter_z/2);       // z size

  G4LogicalVolume* Gem_outter_Log
    = new G4LogicalVolume(Gem_outter_box,
			  Gem_outter_Mat,
                          "Gem_outter_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Gem_outter_Log->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Gem inner plate       *
   *******************************/
  G4VSolid* Gem_inner_box
    = new G4Box("Gem_inner_box",             // World Volume
                Gem_inner_x/2,        // x size
                Gem_inner_y/2,        // y size
                Gem_inner_z/2);       // z size

  G4LogicalVolume* Gem_inner_Log
    = new G4LogicalVolume(Gem_inner_box,
			  Gem_inner_Mat,
                          "Gem_inner_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Gem_inner_Log->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Drift cathode       *
   *******************************/
  G4VSolid* drift_cathode_box 
    = new G4Box("dricath_box",             // World Volume
                drift_cathode_x/2,        // x size
                drift_cathode_y/2,        // y size
                drift_cathode_z/2);       // z size
  
  G4LogicalVolume* driftcathodeLog 
    = new G4LogicalVolume(drift_cathode_box,
			  Drift_cathode_Mat,
                          "dricathLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  driftcathodeLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Readout bar       *
   *******************************/
  G4VSolid* readout_bar_box 
    = new G4Box("readoutbar_box",             // World Volume
                readout_bar_x/2,        // x size
                readout_bar_y/2,        // y size
                readout_bar_z/2);       // z size

  G4LogicalVolume* readoutbarLog
    = new G4LogicalVolume(readout_bar_box,
			  Readout_bar_Mat,
                          "readoutbarLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  readoutbarLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Readout plate       *
   *******************************/
  G4VSolid* readout_plate_box 
    = new G4Box("readoutplate_box",             // World Volume
                readout_plate_x/2,        // x size
                readout_plate_y/2,        // y size
                readout_plate_z/2);       // z size
  
  G4LogicalVolume* readoutplateLog 
    = new G4LogicalVolume(readout_plate_box,
			  Readout_plate_Mat,
                          "readoutplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  readoutplateLog->SetVisAttributes(G4VisAttributes::GetInvisible());


    /*******************************
   * The Gem       *
   *******************************/

  G4VSolid* Gem_box 
    = new G4Box("Gem_box",             // World Volume
                Gem_x/2,        // x size
                Gem_y/2,        // y size
                Gem_z/2);       // z size
  
  G4LogicalVolume* GemLog 
    = new G4LogicalVolume(Gem_box,
			  Gem_Mat,
                          "GemLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  //GemLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    // put Gem outter plate in Gem
    //

  Zcenter = -0.5*Gem_z + drift_cathode_z + gap1 - gap2 - 0.5*Gem_outter_z;
  for (int i=0; i<num_Gem_outter; i++){
    if(i%2==0) Zcenter = Zcenter + gap2 + Gem_outter_z;
    if(i%2==1) Zcenter = Zcenter + Gem_inner_z + Gem_outter_z;
    G4cout<<i<<" Gem outter plate with Zcenter = "<<Zcenter<<" mm"<<G4endl;
    G4String name = "Gem_outter_" + G4String(std::to_string(i).c_str());
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    Gem_outter_Log,             //its logical volume
                    name,                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    1);          //overlaps checking
  }

    // put Gem inner plate in Gem
    //

  Zcenter = -0.5*Gem_z + drift_cathode_z + gap1 - gap2 - Gem_outter_z - 0.5*Gem_inner_z;
  for (int i=0; i<num_Gem_inner; i++){
    Zcenter = Zcenter + gap2 + 2*Gem_outter_z + Gem_inner_z;
    G4cout<<i<<" Gem inner plate with Zcenter = "<<Zcenter<<" mm"<<G4endl;
    G4String name = "Gem_inner_" + G4String(std::to_string(i).c_str());
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    Gem_inner_Log,             //its logical volume
                    name,                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    1);          //overlaps checking
  }

  // put drift cathode in Gem
  //
  Zcenter = -0.5*(Gem_z-drift_cathode_z);
  G4VPhysicalVolume* driftcathodePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    driftcathodeLog,             //its logical volume
                    "driftcathode",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

  // put readout bar in Gem
  //
  const G4int num_readout_bar=Gem_x/readout_bar_gap_x;

  Zcenter = 0.5*Gem_z-readout_plate_z-0.5*readout_bar_z;
  Xcenter = -0.5*(Gem_x-readout_bar_x)-readout_bar_gap_x;
  G4cout<<" Gem readout bar with Zcenter = "<<Zcenter<<" mm"<<G4endl;
  for (int i=0; i<num_readout_bar; i++){
    Xcenter = Xcenter + readout_bar_gap_x;
    //G4cout<<i<<" Gem readout bar with Xcenter = "<<Xcenter<<" mm"<<G4endl;
    G4String name = "Gem_readoutbar_" + G4String(std::to_string(i).c_str());
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,Zcenter),                    //at position
                    readoutbarLog,             //its logical volume
                    name,                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    0);          //overlaps checking
  }

  // put readout plate in Gem
  //
  Zcenter = 0.5*(Gem_z-readout_plate_z);
  G4VPhysicalVolume* readoutplatePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    readoutplateLog,             //its logical volume
                    "readoutplate",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

    /*******************************
   * The Experimental Hall       *
   *******************************/

  G4VSolid* experimentalHall_box 
    = new G4Box("expHall_box",             // World Volume
                experimentalHall_x/2,        // x size
                experimentalHall_y/2,        // y size
                experimentalHall_z/2);       // z size
  
  G4LogicalVolume* experimentalHallLog 
    = new G4LogicalVolume(experimentalHall_box,
                          //air,
			  world_Mat,
                          "expHallLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  G4VPhysicalVolume* experimentalHallPhys 
    = new G4PVPlacement(0,
                        G4ThreeVector(),   //at (0,0,0)
                        "expHall",
                        experimentalHallLog,
                        0,
                        false, 
                        0);
    
  experimentalHallLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  // put Gem in world
  //
  Zcenter = -0.5*Gem_z;
  G4cout<<"1st Gem with Zcenter = "<<Zcenter<<G4endl;
  new G4PVPlacement(0,                       //no rotation
              G4ThreeVector(0,0,Zcenter),                    //at position
              GemLog,             //its logical volume
              "Gem",                //its name
              experimentalHallLog,                //its mother  volume
              false,                   //no boolean operation
              0);               //copy number

  Zcenter = 0.5*Gem_z;
  G4cout<<"2nd Gem with Zcenter = "<<Zcenter<<G4endl;
  new G4PVPlacement(0,                       //no rotation
              G4ThreeVector(0,0,Zcenter),                    //at position
              GemLog,             //its logical volume
              "Gem",                //its name
              experimentalHallLog,                //its mother  volume
              false,                   //no boolean operation
              1);               //copy number
  G4VisAttributes* GemLogAtt = new G4VisAttributes(G4Colour(233/256.,206/256.,238./256,0.1));
  GemLog->SetVisAttributes(GemLogAtt);

  G4VSolid* shell_box = new G4Box("shell_box",shell_x/2,  shell_y/2,  shell_z/2); 
  G4LogicalVolume* shellLog = new G4LogicalVolume(shell_box, Shell_Mat,"shellLog",0, 0, 0);

  Zcenter = -Gem_z-0.5*shell_z;
  G4cout<<"Shell with Zcenter = "<<Zcenter<<G4endl;
  new G4PVPlacement(0,                       //no rotation
              G4ThreeVector(0,0,Zcenter),                    //at position
              shellLog,             //its logical volume
              "Shell",                //its name
              experimentalHallLog,                //its mother  volume
              false,                   //no boolean operation
              0);               //copy number

  fScoringVolume = GemLog;
  fScoringVolume2 = readoutplateLog;

  // visualization attributes ------------------------------------------------

  auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  experimentalHallLog->SetVisAttributes(visAttributes);

  return experimentalHallPhys;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4cout << "Construt the DetectorGeometry" <<G4endl;
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // define a material
  DefineMaterials();
    
  // define some constant
  DefineConstants();

  // Define volumes
  return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
