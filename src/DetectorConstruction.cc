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
#include "PrimaryGeneratorAction.hh"

// geometry
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction(int o)
    : fOptions(o),
      fWorld(NULL),
      fScoringHalfX(DETECTOR_TYPE_COUNT, 0.0),
      fScoringHalfY(DETECTOR_TYPE_COUNT, 0.0),
      fScoringHalfZ(DETECTOR_TYPE_COUNT, 0.0)
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
    "../config/silicon_material.yaml",
    "../config/newrpc_material.yaml",
  };
  char *p = getenv("MUPOS_MATERIAL_CONFIG");
  if(p) { paths = split(p, ':'); }
  for(const std::string &path : paths) { GeometryConfig::LoadMaterials(path.c_str()); }
}

void DetectorConstruction::DefineVolumes()
{
  std::vector<std::string> paths = {
    "../config/newrpc_readout.yaml",
    "../config/newrpc.yaml",
    "../config/silicon_detector.yaml",
    "../config/newrpc_silicon_layout.yaml",
  };
  if(char *p = getenv("MUPOS_VOLUME_CONFIG")) paths = split(p, ':');
  for(const auto &path : paths) GeometryConfig::LoadVolumes(path.c_str());

  fWorld = new G4PVPlacement(0, { 0, 0, 0 }, fLogicalVolumeStore->GetVolume("world"), "world", 0, false, 0, true);
  fScoringTypes.clear();
  fScoringZs.clear();
  fScoringRotations.clear();

  {
    fScoringVolume[DETECTOR_TYPE_NEWRPC] = fLogicalVolumeStore->GetVolume("newrpc_gas");
    G4LogicalVolume *newrpc_electrode = fLogicalVolumeStore->GetVolume("newrpc_electrode");
    fScoringHalfX[DETECTOR_TYPE_NEWRPC] = dynamic_cast<G4Box *>(newrpc_electrode->GetSolid())->GetXHalfLength();
    fScoringHalfY[DETECTOR_TYPE_NEWRPC] = dynamic_cast<G4Box *>(newrpc_electrode->GetSolid())->GetYHalfLength();

    G4double electrodeHalfZ = dynamic_cast<G4Box *>(newrpc_electrode->GetSolid())->GetZHalfLength();
    std::vector<G4double> electrodeZs(0, 0.0);
    WalkVolume(fWorld,
        [newrpc_electrode, &electrodeZs](G4VPhysicalVolume *volume, const G4ThreeVector &r, const G4RotationMatrix &) {
          if(volume->GetLogicalVolume() != newrpc_electrode) { return; }
          electrodeZs.push_back(r.z());
        });
    sort(electrodeZs.begin(), electrodeZs.end());
    fScoringHalfZ[DETECTOR_TYPE_NEWRPC] = (electrodeZs.at(1) - electrodeZs.at(0)) * 0.5 - electrodeHalfZ;

    for(size_t i = 0; 2 * i < electrodeZs.size(); ++i) {
      fScoringTypes.push_back(DETECTOR_TYPE_NEWRPC);
      fScoringZs.push_back((electrodeZs[2 * i] + electrodeZs[2 * i + 1]) * 0.5);
      fScoringRotations.push_back(G4RotationMatrix());
    }
  }

  {
    fScoringVolume[DETECTOR_TYPE_SILICON] = fLogicalVolumeStore->GetVolume("silicon_strip");
    G4LogicalVolume *silicon_plane = fLogicalVolumeStore->GetVolume("silicon_plane");
    fScoringHalfX[DETECTOR_TYPE_SILICON] = dynamic_cast<G4Box *>(silicon_plane->GetSolid())->GetXHalfLength();
    fScoringHalfY[DETECTOR_TYPE_SILICON] = dynamic_cast<G4Box *>(silicon_plane->GetSolid())->GetYHalfLength();
    fScoringHalfZ[DETECTOR_TYPE_SILICON] = dynamic_cast<G4Box *>(silicon_plane->GetSolid())->GetZHalfLength();

    WalkVolume(fWorld, [silicon_plane, this](G4VPhysicalVolume *v, const G4ThreeVector &r, const G4RotationMatrix &rm) {
      if(v->GetLogicalVolume() != silicon_plane) return;
      fScoringTypes.push_back(DETECTOR_TYPE_SILICON);
      fScoringZs.push_back(r.z());
      fScoringRotations.push_back(rm);
    });
  }

  {
    std::vector<size_t> indexes(fScoringZs.size());
    for(size_t i = 0; i < fScoringZs.size(); ++i) { indexes[i] = i; }
    sort(indexes.begin(), indexes.end(), [this](size_t i, size_t j) { return fScoringZs[i] < fScoringZs[j]; });
    std::vector<G4int> types(fScoringTypes.size());
    std::vector<G4double> zs(fScoringZs.size());
    std::vector<G4RotationMatrix> rotations(fScoringRotations.size());
    for(size_t i = 0; i < fScoringZs.size(); ++i) {
      types[i] = fScoringTypes[indexes[i]];
      zs[i] = fScoringZs[indexes[i]];
      rotations[i] = fScoringRotations[indexes[i]];
    }
    fScoringTypes = types;
    fScoringZs = zs;
    fScoringRotations = rotations;
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

  ((PrimaryGeneratorAction *)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction())->Initialize(this);
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

G4double DetectorConstruction::GetDetectorMinZ() const
{
  G4double z = 1.0 / 0.0;
  WalkVolume(NULL, [&z](G4VPhysicalVolume *volume, const G4ThreeVector &r, const G4RotationMatrix &) {
    const G4String &name = volume->GetLogicalVolume()->GetName();
    std::vector<G4String> names = { "silicon_module_xy", "newrpc" };
    if(find(names.begin(), names.end(), name) == names.end()) return;
    auto box = dynamic_cast<G4Box *>(volume->GetLogicalVolume()->GetSolid());
    z = std::min(z, r.z() - box->GetZHalfLength());
  });
  return z;
}

G4double DetectorConstruction::GetDetectorHalfX() const
{
  return *std::max_element(fScoringHalfX.begin(), fScoringHalfX.end());
}

G4double DetectorConstruction::GetDetectorHalfY() const
{
  return *std::max_element(fScoringHalfY.begin(), fScoringHalfY.end());
}
