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

#include "G4AutoDelete.hh"
#include "G4RunManager.hh"
#include "GeometryConfig.hh"
#include "PrimaryGeneratorAction.hh"

// geometry
#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"

// visualization
#include "G4Color.hh"
#include "G4VisAttributes.hh"

// EM field
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4EqMagElectricField.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"

DetectorConstruction::DetectorConstruction(int o)
    : fOptions(o),
      fWorld(NULL),
      fElectrodeHalfX(0.0),
      fElectrodeHalfY(0.0),
      fElectrodeHalfZ(0.0),
      fScoringHalfZ(0.0),
      fScoringGasVolume(NULL)
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
    "../config/rpc_material.yaml",
  };
  char *p = getenv("MUPOS_MATERIAL_CONFIG");
  if(p) { paths = split(p, ':'); }
  for(const std::string &path : paths) { GeometryConfig::LoadMaterials(path.c_str()); }
}

void DetectorConstruction::DefineVolumes()
{
  std::vector<std::string> paths = {
    "../config/rpc_readout.yaml",
    "../config/rpc.yaml",
    "../config/layout.yaml",
  };
  char *p = getenv("MUPOS_VOLUME_CONFIG");
  if(p) { paths = split(p, ':'); }
  for(const std::string &path : paths) { GeometryConfig::LoadVolumes(path.c_str()); }

  fWorld = new G4PVPlacement(0, { 0, 0, 0 }, fLogicalVolumeStore->GetVolume("world"), "world", 0, false, 0, true);
  G4LogicalVolume *rpc_electrode = fLogicalVolumeStore->GetVolume("rpc_electrode");
  fElectrodeHalfX = dynamic_cast<G4Box *>(rpc_electrode->GetSolid())->GetXHalfLength();
  fElectrodeHalfY = dynamic_cast<G4Box *>(rpc_electrode->GetSolid())->GetYHalfLength();
  fElectrodeHalfZ = dynamic_cast<G4Box *>(rpc_electrode->GetSolid())->GetZHalfLength();
  fElectrodeZs.assign(0, 0.0);
  WalkVolume(
      fWorld, [rpc_electrode, this](G4VPhysicalVolume *volume, const G4ThreeVector &r, const G4RotationMatrix &) {
        if(volume->GetLogicalVolume() != rpc_electrode) { return; }
        fElectrodeZs.push_back(r.z());
      });
  sort(fElectrodeZs.begin(), fElectrodeZs.end());
  fScoringHalfZ = (fElectrodeZs.at(1) - fElectrodeZs.at(0)) * 0.5 - fElectrodeHalfZ;
  fScoringZs.resize(fElectrodeZs.size() / 2);
  for(size_t i = 0; i < fScoringZs.size(); ++i) {
    fScoringZs[i] = (fElectrodeZs[2 * i] + fElectrodeZs[2 * i + 1]) * 0.5;
  }
  fScoringGasVolume = fLogicalVolumeStore->GetVolume("rpc_gas");
}

void DetectorConstruction::DefineFields()
{
  // Find all unique physical occurrences of rpc_electric.
  G4String name = "rpc_electric";
  std::vector<G4VPhysicalVolume *> rpc_electrics;
  WalkVolume(NULL, [&name, &rpc_electrics](G4VPhysicalVolume *v) {
    if(v->GetLogicalVolume()->GetName() == name) { rpc_electrics.push_back(v); }
  });
  std::sort(rpc_electrics.begin(), rpc_electrics.end());
  rpc_electrics.erase(std::unique(rpc_electrics.begin(), rpc_electrics.end()), rpc_electrics.end());

  // Determine electric field volume.
  G4double z = fScoringHalfZ * 2;
  auto electric = new G4Box("electric", fElectrodeHalfX, fElectrodeHalfY, fScoringHalfZ);

  // Turn on electric field.
  G4double step = z * 0.01;
  G4double U = 10100 * volt, E = U / z;
  auto field = new G4UniformElectricField(G4ThreeVector{ 0, 0, E });
  //auto equation_of_motion = new G4EqMagElectricField(field);
  //auto stepper = new G4ClassicalRK4(equation_of_motion);
  //auto int_driver = new G4MagInt_Driver(step, stepper, stepper->GetNumberOfVariables());
  //auto chord_finder = new G4ChordFinder(int_driver);
  //auto manager = new G4FieldManager(field, chord_finder);
  auto magField = new G4UniformMagField(G4ThreeVector{ 0, 0, 0 });
  auto manager = new G4FieldManager(field);  // no ownership transfer
  manager->CreateChordFinder(magField);      // no ownership transfer
  G4AutoDelete::Register(field);
  G4AutoDelete::Register(magField);

  for(G4VPhysicalVolume *rpc_electric : rpc_electrics) {
    // Split rpc_electric into two parts, w/ and w/o electric field.
    PrintVolumes(rpc_electric);
    rpc_electric = PartitionVolume(
        rpc_electric, [&electric, &name](G4VSolid *solid, const G4ThreeVector &r, const G4RotationMatrix &rm) {
          // r_s = r + rm * r'_s  <=>  r'_s = - (rm^-1 * r) + rm^-1 * r_s
          static size_t g_id;
          size_t id = g_id++;
          std::vector<G4VSolid *> parts;
          parts.reserve(2);
          auto rotation = new G4RotationMatrix(rm.inverse());
          auto rps = -(*rotation * r);
          name = "part_" + std::to_string(id) + "_0_" + solid->GetName();
          parts.push_back(new G4IntersectionSolid(name, solid, electric, rotation, rps));  // rotation owned still
          name = "part_" + std::to_string(id) + "_1_" + solid->GetName();
          parts.push_back(new G4SubtractionSolid(name, solid, electric, rotation, rps));  // rotation owned still
          G4AutoDelete::Register(rotation);
          return parts;
        });
    PrintVolumes(rpc_electric);
    G4VPhysicalVolume *electric_volume = rpc_electric->GetLogicalVolume()->GetDaughter(0);
    electric_volume->GetLogicalVolume()->SetFieldManager(manager, true);
    {
      G4VisAttributes attr;
      G4Color color;
      G4Color::GetColour("green", color), color.SetAlpha(0.2);
      attr.SetColor(color), attr.SetForceSolid();
      electric_volume->GetLogicalVolume()->SetVisAttributes(attr);
    }

    // Force step limit in field areas.
    auto limits = new G4UserLimits(step);
    WalkVolume(electric_volume->GetLogicalVolume(), [limits](G4LogicalVolume *volume) {
      G4cout << "Setting step limit for " << volume->GetName() << G4endl;
      volume->SetUserLimits(limits);  // no ownership transfer
    });
    G4AutoDelete::Register(limits);
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
  DefineFields();
  PrintVolumes(NULL);

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

G4VPhysicalVolume *DetectorConstruction::PartitionVolume(G4VPhysicalVolume *volume,
    const std::function<std::vector<G4VSolid *>(G4VSolid *, const G4ThreeVector &, const G4RotationMatrix &)>
        &partition) const
{
  std::vector<std::vector<std::vector<G4LogicalVolume *>>> stack(1);
  WalkVolume(
      volume, [&stack](G4VPhysicalVolume *, const G4ThreeVector &, const G4RotationMatrix &) { stack.emplace_back(); },
      [&stack, &partition](G4VPhysicalVolume *v, const G4ThreeVector &r, const G4RotationMatrix &rm) {
        static size_t g_id;
        size_t id = g_id++;

        // Pop children from stack.
        auto children = std::move(stack.back());  // [nchild, npart]
        stack.pop_back();
        size_t nchild = children.size();

        // Add self as a child of parent.
        stack.back().emplace_back();
        auto &result = stack.back().back();  // self, [npart]

        // Partition solid.
        auto solids = partition(v->GetLogicalVolume()->GetSolid(), r, rm);
        size_t npart = solids.size();
        for(auto &child : children) {
          if(child.size() != npart) { throw std::logic_error("inconsistent number of partitions"); }
        }

        // Create npart volumes as result.
        result.reserve(npart);
        for(size_t ipart = 0; ipart < npart; ++ipart) {
          if(solids[ipart] == NULL) {
            for(size_t ichild = 0; ichild < nchild; ++ichild) {
              if(children[ichild][ipart]) { throw std::logic_error("null volume contains non-null child"); }
            }
            result.push_back(NULL);
            continue;
          }
          // Create logical volume.
          G4String name =
              "part_" + std::to_string(id) + "_" + std::to_string(ipart) + "_" + v->GetLogicalVolume()->GetName();
          result.push_back(new G4LogicalVolume(solids[ipart], v->GetLogicalVolume()->GetMaterial(), name));
          result.back()->SetVisAttributes(v->GetLogicalVolume()->GetVisAttributes());
          // Add children.
          for(size_t ichild = 0; ichild < nchild; ++ichild) {
            if(children[ichild][ipart] == NULL) { continue; }
            G4VPhysicalVolume *dau = v->GetLogicalVolume()->GetDaughter(ichild);
            name = "part_" + std::to_string(id) + "_" + std::to_string(ipart) + "_" + dau->GetName();
            new G4PVPlacement(dau->GetRotation(), dau->GetTranslation(), children[ichild][ipart], name, result.back(),
                false, 0, true);
          }
        }
      });

  // Create group.
  static size_t g_id;
  size_t id = g_id++;
  G4String name = "group_" + std::to_string(id) + "_" + volume->GetLogicalVolume()->GetName();
  auto solid = volume->GetLogicalVolume()->GetSolid();
  auto material = volume->GetLogicalVolume()->GetMaterial();
  auto logical = new G4LogicalVolume(solid, material, name);
  logical->SetVisAttributes(volume->GetLogicalVolume()->GetVisAttributes());
  for(size_t ichild = 0; ichild < stack[0][0].size(); ++ichild) {
    if(stack[0][0][ichild] == NULL) { continue; }
    name = "group_" + std::to_string(id) + "_" + std::to_string(ichild) + "_" + volume->GetName();
    new G4PVPlacement(0, { 0, 0, 0 }, stack[0][0][ichild], name, logical, false, 0, true);
  }

  // Replace physical volume.
  name = "group_" + std::to_string(id) + "_" + volume->GetName();
  auto rotation = volume->GetObjectRotation();
  auto translation = volume->GetObjectTranslation();
  auto mother = volume->GetMotherLogical();
  mother->RemoveDaughter(volume);
  delete volume;
  return new G4PVPlacement(rotation, translation, logical, name, mother, false, 0, true);
}

G4double DetectorConstruction::GetDetectorMinZ() const
{
  G4double z = 1.0 / 0.0;
  WalkVolume(NULL, [&z](G4VPhysicalVolume *volume, const G4ThreeVector &r, const G4RotationMatrix &) {
    if(volume->GetLogicalVolume()->GetName() != "rpc") { return; }
    auto box = dynamic_cast<G4Box *>(volume->GetLogicalVolume()->GetSolid());
    z = std::min(z, r.z() - box->GetZHalfLength());
  });
  return z;
}

G4double DetectorConstruction::GetDetectorHalfX() const
{
  //return dynamic_cast<G4Box *>(fLogicalVolumeStore->GetVolume("rpc")->GetSolid())->GetXHalfLength();  // more precise
  return GetScoringHalfX();  // faster
}

G4double DetectorConstruction::GetDetectorHalfY() const
{
  //return dynamic_cast<G4Box *>(fLogicalVolumeStore->GetVolume("rpc")->GetSolid())->GetYHalfLength();  // more precise
  return GetScoringHalfY();  // faster
}
