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

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "GeometryConfig.hh"
#include "GpsPrimaryGeneratorAction.hh"

// geometry
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"

// visualization
#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include <unordered_set>

DetectorConstruction::DetectorConstruction(int o)
    : fOptions(o),
      fWorld(NULL)
{
  if(fOptions) { throw std::invalid_argument("options unimplemented"); }
  fLogicalVolumeStore = G4LogicalVolumeStore::GetInstance();
  fPhysicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
}

DetectorConstruction::~DetectorConstruction()
{
  // empty
}

static std::vector<std::string> split(const std::string &str, char c)
{
  std::vector<std::string> v;
  size_t p = 0, q = 0;
  while((q = str.find(c, p)) != str.npos) {
    v.push_back(str.substr(p, q - p));
    p = q + 1;
  }
  v.push_back(str.substr(p));
  return v;
}

void DetectorConstruction::DefineMaterials()
{
  std::vector<std::string> paths = {
    "../config/hiaf_material.yaml",
  };
  char *p = getenv("MUPOS_MATERIAL_CONFIG");
  if(p) { paths = split(p, ':'); }
  for(const std::string &path : paths) { GeometryConfig::LoadMaterials(path.c_str()); }
}

void DetectorConstruction::DefineVolumes()
{
  std::vector<std::string> paths = {
    "../config/mwdc_up_detector.yaml",
    "../config/mwdc_down_detector.yaml",
    "../config/hiaf_layout.yaml",
  };
  char *p = getenv("MUPOS_VOLUME_CONFIG");
  if(p) { paths = split(p, ':'); }
  for(const std::string &path : paths) { GeometryConfig::LoadVolumes(path.c_str()); }

  fWorld = new G4PVPlacement(0, { 0, 0, 0 }, fLogicalVolumeStore->GetVolume("world"), "world", 0, false, 0, true);

  std::unordered_set<std::string> scoringNames {
    "mwdc_up_gas_wires_0",
    "mwdc_up_gas_wires_1",
    "mwdc_up_gas_wires_2",
    "mwdc_up_gas_wires_3",
    "mwdc_up_gas_wires_5",
    "mwdc_up_gas_wires_6",
    "mwdc_up_gas_wires_7",
    "mwdc_up_gas_wires_8",
    "mwdc_down_gas_wires_0",
    "mwdc_down_gas_wires_1",
    "mwdc_down_gas_wires_2",
    "mwdc_down_gas_wires_3",
    "mwdc_down_gas_wires_5",
    "mwdc_down_gas_wires_6",
    "mwdc_down_gas_wires_7",
    "mwdc_down_gas_wires_8",
  };
  WalkVolume(fWorld, [this, &scoringNames](G4VPhysicalVolume *volume,
        const G4ThreeVector &r, const G4RotationMatrix &rm) {
    if(scoringNames.count(volume->GetLogicalVolume()->GetName()) == 0) return;
    auto box = dynamic_cast<G4Box *>(volume->GetLogicalVolume()->GetSolid());
    G4ThreeVector hr { box->GetXHalfLength(), box->GetYHalfLength(), box->GetZHalfLength() };
    fScoringVolumes.emplace_back();
    fScoringVolumes.back() = { &hr, &r, &rm };
  });
  sort(fScoringVolumes.begin(), fScoringVolumes.end(), [](const BoxVolume &x, const BoxVolume &y) {
    return x.CenterZ < y.CenterZ;  // [XXX] Rotation not considered; OK for HIAF.
  });
  G4cout << "Scoring volumes:" << G4endl;
  for(size_t i = 0; i < fScoringVolumes.size(); ++i) {
    G4cout << "  * " << fScoringVolumes[i].CenterZ << ": " << fScoringVolumes[i].Alpha / deg << G4endl;
  }
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4SolidStore::GetInstance()->Clean();
  fLogicalVolumeStore->Clean();
  fPhysicalVolumeStore->Clean();

  DefineMaterials();
  DefineVolumes();
  PrintVolumes(NULL);

  ((GpsPrimaryGeneratorAction *)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())->Initialize(this);
  return fWorld;
}

void DetectorConstruction::PrintVolumes(G4VPhysicalVolume *volume) const
{
  size_t depth = 0;
  WalkVolume(
      volume,
      [&depth](G4VPhysicalVolume *v) {
        G4cout << std::string(2 * depth++, ' ');
        G4cout << v->GetName() << " - " << v->GetLogicalVolume()->GetName() << " - "
               << v->GetLogicalVolume()->GetSolid()->GetName() << G4endl;
      },
      [&depth](G4VPhysicalVolume *) { --depth; });
}

static void WalkVolume(G4LogicalVolume *volume, const std::function<void(G4LogicalVolume *)> &enter,
    const std::function<void(G4LogicalVolume *)> &leave)
{
  if(enter) { enter(volume); }
  for(size_t i = 0; i < volume->GetNoDaughters(); ++i) {
    WalkVolume(volume->GetDaughter(i)->GetLogicalVolume(), enter, leave);
  }
  if(leave) { leave(volume); }
}

void DetectorConstruction::WalkVolume(G4LogicalVolume *volume, const std::function<void(G4LogicalVolume *)> &enter,
    const std::function<void(G4LogicalVolume *)> &leave) const
{
  if(volume == NULL) { volume = fWorld->GetLogicalVolume(); }
  if(volume == NULL) { return; }
  ::WalkVolume(volume, enter, leave);
}

static void WalkVolume(G4VPhysicalVolume *volume, const std::function<void(G4VPhysicalVolume *)> &enter,
    const std::function<void(G4VPhysicalVolume *)> &leave)
{
  if(enter) { enter(volume); }
  G4LogicalVolume *logical = volume->GetLogicalVolume();
  for(size_t i = 0; i < logical->GetNoDaughters(); ++i) { WalkVolume(logical->GetDaughter(i), enter, leave); }
  if(leave) { leave(volume); }
}

void DetectorConstruction::WalkVolume(G4VPhysicalVolume *volume, const std::function<void(G4VPhysicalVolume *)> &enter,
    const std::function<void(G4VPhysicalVolume *)> &leave) const
{
  if(volume == NULL) { volume = fWorld; }
  if(volume == NULL) { return; }
  ::WalkVolume(volume, enter, leave);
}

void DetectorConstruction::WalkVolume(G4VPhysicalVolume *volume,
    const std::function<void(G4VPhysicalVolume *, const G4ThreeVector &, const G4RotationMatrix &)> &enter,
    const std::function<void(G4VPhysicalVolume *, const G4ThreeVector &, const G4RotationMatrix &)> &leave) const
{
  G4ThreeVector r = { 0, 0, 0 };
  G4RotationMatrix rm = { 0, 0, 0 };
  WalkVolume(
      volume,
      [&r, &rm, &enter](G4VPhysicalVolume *v) {
        r += rm * v->GetObjectTranslation();
        if(G4RotationMatrix *rotation = v->GetObjectRotation()) { rm = rm * *rotation; }
        if(enter) { enter(v, r, rm); }
      },
      [&r, &rm, &leave](G4VPhysicalVolume *v) {
        if(leave) { leave(v, r, rm); }
        if(G4RotationMatrix *rotation = v->GetObjectRotation()) { rm = rm * rotation->inverse(); }
        r -= rm * v->GetObjectTranslation();
      });
}
