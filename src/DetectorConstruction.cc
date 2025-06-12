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

DetectorConstruction::DetectorConstruction(int o)
    : fOptions(o),
      fWorld(NULL),
      fSiliconStrip(NULL),
      fScoringHalfX(0.0),
      fScoringHalfY(0.0),
      fScoringHalfZ(0.0),
      fStripInterval(0.0),
      fNSiliconStrips(0)
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
  };
  char *p = getenv("MUPOS_MATERIAL_CONFIG");
  if(p) { paths = split(p, ':'); }
  for(const std::string &path : paths) { GeometryConfig::LoadMaterials(path.c_str()); }
}

void DetectorConstruction::DefineVolumes()
{
  std::vector<std::string> paths = {
    "../config/silicon_detector.yaml",
    "../config/silicon_layout.yaml",
  };
  char *p = getenv("MUPOS_VOLUME_CONFIG");
  if(p) { paths = split(p, ':'); }
  for(const std::string &path : paths) { GeometryConfig::LoadVolumes(path.c_str()); }

  fWorld = new G4PVPlacement(0, { 0, 0, 0 }, fLogicalVolumeStore->GetVolume("world"), "world", 0, false, 0, true);
  fSiliconStrip = fLogicalVolumeStore->GetVolume("silicon_strip");
  G4LogicalVolume *silicon_plane = fLogicalVolumeStore->GetVolume("silicon_plane");
  fScoringHalfX = dynamic_cast<G4Box *>(silicon_plane->GetSolid())->GetXHalfLength();
  fScoringHalfY = dynamic_cast<G4Box *>(silicon_plane->GetSolid())->GetYHalfLength();
  fScoringHalfZ = dynamic_cast<G4Box *>(silicon_plane->GetSolid())->GetZHalfLength();

  fScoringZs.assign(0, 0.0);
  fScoringRotations.assign(0, G4RotationMatrix());
  G4VPhysicalVolume *silicon_plane_physical = NULL;
  WalkVolume(fWorld,
      [silicon_plane, this, &silicon_plane_physical](
          G4VPhysicalVolume *volume, const G4ThreeVector &r, const G4RotationMatrix &rm) {
        if(volume->GetLogicalVolume() != silicon_plane) { return; }
        fScoringZs.push_back(r.z());
        fScoringRotations.push_back(rm);
        silicon_plane_physical = volume;
      });
  {
    std::vector<size_t> indexes(fScoringZs.size());
    for(size_t i = 0; i < fScoringZs.size(); ++i) { indexes[i] = i; }
    sort(indexes.begin(), indexes.end(), [this](size_t i, size_t j) { return fScoringZs[i] < fScoringZs[j]; });
    std::vector<G4double> zs(fScoringZs.size());
    std::vector<G4RotationMatrix> rotations(fScoringRotations.size());
    for(size_t i = 0; i < fScoringZs.size(); ++i) {
      zs[i] = fScoringZs[indexes[i]];
      rotations[i] = fScoringRotations[indexes[i]];
    }
    fScoringZs = zs;
    fScoringRotations = rotations;
  }
  G4cout << "Scoring volumes:" << G4endl;
  for(size_t i = 0; i < fScoringZs.size(); ++i) {
    G4cout << "  * " << fScoringZs[i] << ": " << fScoringRotations[i].delta() / deg << G4endl;
  }

  fNSiliconStrips = 0;
  G4double x0 = -1, x1 = -1;
  WalkVolume(silicon_plane_physical,
      [this, &x0, &x1](G4VPhysicalVolume *volume, const G4ThreeVector &r, const G4RotationMatrix &) {
        if(volume->GetLogicalVolume() != fSiliconStrip) { return; }
        ++fNSiliconStrips;
        if(x0 < 0) {
          x0 = r.x();
        } else if(x1 < 0) {
          x1 = r.x();
        }
      });
  assert(x0 >= 0 && x1 >= 0);
  fStripInterval = x1 - x0;
  assert(fStripInterval > 0);
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4SolidStore::GetInstance()->Clean();
  fLogicalVolumeStore->Clean();
  fPhysicalVolumeStore->Clean();

  DefineMaterials();
  DefineVolumes();
  //PrintVolumes(NULL);

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

G4double DetectorConstruction::GetDetectorMinZ() const
{
  G4double z = 1.0 / 0.0;
  WalkVolume(NULL, [&z](G4VPhysicalVolume *volume, const G4ThreeVector &r, const G4RotationMatrix &) {
    if(volume->GetLogicalVolume()->GetName() != "silicon_module_xy") { return; }
    auto box = dynamic_cast<G4Box *>(volume->GetLogicalVolume()->GetSolid());
    z = std::min(z, r.z() - box->GetZHalfLength());
  });
  return z;
}

G4double DetectorConstruction::GetDetectorHalfX() const { return GetScoringHalfX(); }
G4double DetectorConstruction::GetDetectorHalfY() const { return GetScoringHalfY(); }
