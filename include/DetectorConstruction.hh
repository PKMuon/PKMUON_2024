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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include <vector>
#include <functional>

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4LogicalVolumeStore;
class G4PhysicalVolumeStore;

#define DETECTOR_OPTION_SCORING_ONLY  0b00000001
#define DETECTOR_OPTION_VACUUM_ENV    0b00000010

class DetectorConstruction : public G4VUserDetectorConstruction {
public:

  DetectorConstruction(int options = 0);
  ~DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();

  // Call these methods after Construct().
  G4double GetScoringHalfX() const { return fScoringHalfX; }
  G4double GetScoringHalfY() const { return fScoringHalfY; }
  G4double GetScoringHalfZ() const { return fScoringHalfZ; }
  std::vector<G4double> GetScoringZs() const { return fScoringZs; }
  G4double GetDetectorMinZ() const;
  G4double GetDetectorHalfX() const;
  G4double GetDetectorHalfY() const;

  // Hierarchic options.
  void PrintVolumes(G4VPhysicalVolume *) const;
  void WalkVolume(G4LogicalVolume *volume,
      const std::function<void(G4LogicalVolume *)> &enter,
      const std::function<void(G4LogicalVolume *)> &leave = nullptr) const;
  void WalkVolume(G4VPhysicalVolume *volume,
      const std::function<void(G4VPhysicalVolume *)> &enter,
      const std::function<void(G4VPhysicalVolume *)> &leave = nullptr) const;
  void WalkVolume(G4VPhysicalVolume *volume,
      const std::function<void(G4VPhysicalVolume *, const G4ThreeVector &, const G4RotationMatrix &)> &enter,
      const std::function<void(G4VPhysicalVolume *, const G4ThreeVector &, const G4RotationMatrix &)> &leave = nullptr) const;
  G4VPhysicalVolume *PartitionVolume(G4VPhysicalVolume *volume,
      const std::function<std::vector<G4VSolid *>(G4VSolid *, const G4ThreeVector &, const G4RotationMatrix &)> &partition) const;

private:
  void DefineMaterials();
  void DefineVolumes();
  void DefineFields();

  const int fOptions;
  G4LogicalVolumeStore *fLogicalVolumeStore;
  G4PhysicalVolumeStore *fPhysicalVolumeStore;
  G4VPhysicalVolume *fWorld;
  G4LogicalVolume *fElectrodeVolume;
  G4double fScoringHalfX, fScoringHalfY, fScoringHalfZ;
  std::vector<G4double> fScoringZs;
};

#endif
