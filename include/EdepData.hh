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

#include <Rtypes.h>

#include <tuple>

#ifndef EdepData_h
#define EdepData_h 1

struct EdepKey {
public:
  Int_t Id;
  Int_t Pid;
  Int_t Process;

  auto Tuple() { return std::tie(Id, Pid, Process); }
  auto Tuple() const { return std::tie(Id, Pid, Process); }

  // Typically ROOT is not built with C++20, unfortunately.
  //auto operator<=>(const EdepKey &) const = default;
};

inline bool operator<(const EdepKey &lhs, const EdepKey &rhs) { return lhs.Tuple() < rhs.Tuple(); }
inline bool operator==(const EdepKey &lhs, const EdepKey &rhs) { return lhs.Tuple() == rhs.Tuple(); }

namespace std {

template<>
struct hash<EdepKey> {
  size_t operator()(const EdepKey &key) const
  {
    return hash<int>()(key.Id) ^ hash<int>()(key.Pid) ^ hash<int>()(key.Process);
  }
};

}  // namespace std

struct EdepValue {
public:
  Double_t Value;
  Double_t X;
  Double_t Y;

  auto Tuple() && { return std::tie(Value, X, Y); }
  auto Tuple() const && { return std::tie(Value, X, Y); }

  EdepValue &Add(Double_t v, Double_t x, Double_t y)
  {
    Value += v, X += v * x, Y += v * y;
    return *this;
  }

  EdepValue &&Finish()
  {
    if(Value > 0) X /= Value, Y /= Value;
    return std::move(*this);
  }
};

#endif
