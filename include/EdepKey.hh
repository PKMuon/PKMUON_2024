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

#include <TString.h>

#ifndef EdepKey_h
#define EdepKey_h 1

struct EdepKey {
public:
  Int_t Id;
  Int_t Pid;
  TString Process;

  // Typically ROOT is not built with C++20, unfortunately.
  //auto operator<=>(const EdepKey &) const = default;
};

inline bool operator<(const EdepKey &lhs, const EdepKey &rhs)
{
  return lhs.Id < rhs.Id || (lhs.Id == rhs.Id && lhs.Pid < rhs.Pid)
      || (lhs.Id == rhs.Id && lhs.Pid == rhs.Pid && lhs.Process < rhs.Process);
}

inline bool operator==(const EdepKey &lhs, const EdepKey &rhs)
{
  return lhs.Id == rhs.Id && lhs.Pid == rhs.Pid && lhs.Process == rhs.Process;
}

namespace std {

template<>
struct hash<EdepKey> {
  size_t operator()(const EdepKey &key) const
  {
    return hash<size_t>{}(((size_t)key.Id << 32) | key.Pid) ^ hash<const char *>{}(key.Process);
  }
};

}  // namespace std

#endif
