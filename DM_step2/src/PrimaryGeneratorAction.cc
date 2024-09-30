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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fPrimaryGeneratorActionMessenger(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  G4String particleName;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="mu-");
  fParticleGun->SetParticleDefinition(particle);

  input_filename = "xxx.txt";
  fPrimaryGeneratorActionMessenger = new PrimaryGeneratorActionMessenger(this);

//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0.,1.));
//  fParticleGun->SetParticleEnergy(3*GeV);  //change this
//  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,-8*mm)); //change this

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  input_stream.close();
  delete fParticleGun;
  delete fPrimaryGeneratorActionMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    G4int i = anEvent->GetEventID();
    if(i==0){
      G4cout<<"input_filename = "<<input_filename<<G4endl;
      input_stream.open(input_filename, std::ios::in);
      if (input_stream.fail()) {
          G4cout << "Sorry, couldn’t open file" << G4endl;
      }
    }

    G4int eventid;
    G4double x, y, px, py, pz, energy;
    input_stream >> eventid >> x >> y >> px >> py >> pz >> energy;
//G4cout<<"evtid = "<<eventid<<" , x = "<<x<<" , y = "<<y<<", px = "<<px<<" , py = "<<py<<" , pz = "<<pz<<" , energy = "<<energy<<G4endl;

    Run::GetInstance()->SetEvtID(eventid);
    fParticleGun->SetParticleEnergy(energy * MeV);
    //----------------------------------------------------------
    fParticleGun->SetParticlePosition(G4ThreeVector(x * mm, y * mm, -13.45 * mm));
    //----------------------------------------------------------
    G4ThreeVector v(px,py,pz);
    G4ThreeVector unitv = v.unit();
    fParticleGun->SetParticleMomentumDirection(unitv);
    //-------------------------------------------------------

    fParticleGun->GeneratePrimaryVertex(anEvent);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
