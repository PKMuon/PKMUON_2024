#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Object.hh"
#include "Run.hh"
#include "PrimaryGeneratorMessenger.hh"

#include <iomanip>
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"

using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    particleGun = new G4ParticleGun();  
    fGaisserTangGen = new GaisserTangGenerator(); 

    gunMessenger = new PrimaryGeneratorMessenger(this);
    particleTable = G4ParticleTable::GetParticleTable();
    
    fDetectorMaxZ = 0.0;
    fDetectorHalfX = 0.0;
    fDetectorHalfY = 0.0;

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete particleGun;
    delete gunMessenger;
    delete fGaisserTangGen;
}

void PrimaryGeneratorAction::Initialize(const DetectorConstruction* detectorConstruction)
{
    fDetectorMaxZ = detectorConstruction->GetDetectorMaxZ();
    fDetectorHalfX = detectorConstruction->GetDetectorHalfX();
    fDetectorHalfY = detectorConstruction->GetDetectorHalfY();
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    G4double energy, theta, phi;
    fGaisserTangGen->GenerateEvent(energy, theta, phi);
    
    particleGun->SetParticleDefinition(particleTable->FindParticle(13));
    
    particleGun->SetParticleEnergy(energy * GeV);

    G4double r = std::sqrt(fDetectorHalfX*fDetectorHalfX + fDetectorHalfY*fDetectorHalfY + fDetectorMaxZ*fDetectorMaxZ);
    G4double P_theta = std::acos(2.0 * G4UniformRand() - 1.0);
    G4double P_phi = 2.0 * CLHEP::pi * G4UniformRand();
    G4double x = r * std::sin(P_theta) * std::cos(P_phi);
    G4double y = r * std::sin(P_theta) * std::sin(P_phi);
    G4double z = r * std::cos(P_theta);
    
    particleGun->SetParticlePosition(G4ThreeVector(-x, z, y));
    
    G4double sinTheta = std::sin(theta);
    G4double dirX = sinTheta * std::cos(phi);
    G4double dirY = sinTheta * std::sin(phi);
    G4double dirZ = -std::cos(theta);  
    
    particleGun->SetParticleMomentumDirection(G4ThreeVector(-dirX, dirZ, dirY));
    
    particleGun->SetParticleTime(0.0);
    
    particleGun->GeneratePrimaryVertex(anEvent);
    
    Event* event = Run::GetInstance()->GetEvent();
    event->Reset();
    event->Pid = particleGun->GetParticleDefinition()->GetPDGEncoding();
    
    G4double mass = particleGun->GetParticleDefinition()->GetPDGMass();
    G4double e = particleGun->GetParticleEnergy() + mass;
    G4ThreeVector v = std::sqrt(e*e - mass*mass) * particleGun->GetParticleMomentumDirection();
    
    event->Px = -v.x();
    event->Py = v.z();
    event->Pz = v.y();
    event->E = e - mass;
    event->X = -x;
    event->Y = z;
    event->Z = y;
    event->T = particleGun->GetParticleTime();
}