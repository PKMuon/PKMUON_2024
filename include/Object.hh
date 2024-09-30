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

#ifndef Object_h
#define Object_h 1

#include <TObject.h>

#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "EdepData.hh"

class DetectorConstruction;
class G4Track;
class G4DynamicParticle;

class Track : public TObject {
public:
  Track &operator=(const G4Track &);

  Int_t Id;
  Int_t Mother;
  Int_t Pid;
  Double_t Px;
  Double_t Py;
  Double_t Pz;
  Double_t E;
  Double_t X;
  Double_t Y;
  Double_t Z;
  Double_t T;

  ClassDef(Track, 1);
};

class Params : public TObject {
public:
  Params &operator=(const DetectorConstruction &);

  Double_t GammaCut;
  Double_t GammaThreshold;
  Double_t ElectronCut;
  Double_t ElectronThreshold;
  Double_t PositronCut;
  Double_t PositronThreshold;
  Double_t ProtonCut;
  Double_t ProtonThreshold;

  std::vector<double> LayerZ;
  Double_t CellX;
  Double_t CellY;
  Int_t HalfNCellX;
  Int_t HalfNCellY;

  ClassDef(Params, 1);
};

class Edep : public TObject {
public:
  Edep &operator=(std::pair<const EdepKey, Double_t> &p)
  {
    auto &[key, value] = p;
    std::tie(Id, Pid, Process) = key.Tuple();
    Value = value;
    return *this;
  }

  Int_t Id;
  Int_t Pid;
  Int_t Process;
  Double_t Value;

  ClassDef(Edep, 1);
};

class Scatter : public TObject {
public:
  Scatter &operator=(const std::tuple<const G4Track *, const G4DynamicParticle *, const G4DynamicParticle *> &t);

  Int_t Id;
  Int_t Pid[3];
  Double_t Px[3];
  Double_t Py[3];
  Double_t Pz[3];
  Double_t E[3];
  Double_t X;
  Double_t Y;
  Double_t Z;
  Double_t T;

  ClassDef(Scatter, 1);
};

class Process : public TObject {
public:
  Process &operator=(const std::pair<Int_t, const std::string &> &t)
  {
    std::tie(Id, Name) = t;
    return *this;
  }

  Int_t Id;
  std::string Name;

  ClassDef(Process, 1);
};

#endif
