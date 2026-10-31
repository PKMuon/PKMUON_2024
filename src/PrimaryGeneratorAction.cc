#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "Object.hh"
#include "Run.hh"
#include "PrimaryGeneratorMessenger.hh"

#include <iomanip>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
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
    fPlaneSourceZ=std::numeric_limits<double>::infinity();
    detectorConstruction->WalkVolume(nullptr,[this](G4VPhysicalVolume *volume,const G4ThreeVector& r,const G4RotationMatrix&) {
        if(volume->GetLogicalVolume()->GetName()!="newrpc") return;
        auto box=dynamic_cast<G4Box*>(volume->GetLogicalVolume()->GetSolid());
        if(!box) throw std::runtime_error("Plane source requires box RPC envelopes");
        if(r.z()-box->GetZHalfLength()<fPlaneSourceZ) {
            fPlaneSourceZ=r.z()-box->GetZHalfLength();
            fPlaneHalfX=box->GetXHalfLength(); fPlaneHalfY=box->GetYHalfLength();
        }
    });
    fPlaneSourceZ-=1*mm;
    auto store=G4LogicalVolumeStore::GetInstance();
    auto half=[&](const char* name) { return dynamic_cast<G4Box*>(store->GetVolume(name)->GetSolid())->GetZHalfLength(); };
    const auto offset=half("newrpc_xy_readout_board")+2*half("newrpc_insulating_film")+
                      2*half("newrpc_electrode")+2*half("newrpc_glass")+half("newrpc_gas");
    fPlaneGasCenters.clear();
    for(auto z:detectorConstruction->GetScoringZs()) fPlaneGasCenters.push_back(z+offset);
    if(std::getenv("MUPOS_GAISSER_PLANE")) G4cout<<"RPC plane source: z="<<fPlaneSourceZ/mm
        <<" mm; size="<<2*fPlaneHalfX/mm<<" x "<<2*fPlaneHalfY/mm<<" mm; straight full-gap intersections required"<<G4endl;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Opt-in reproducible vertical beam for RPC validation; production source unchanged.
    if (std::getenv("MUPOS_RPC_PENCIL_BEAM")) {
        particleGun->SetParticleDefinition(particleTable->FindParticle(13));
        particleGun->SetParticleEnergy(1*GeV);
        particleGun->SetParticlePosition({0,0,-1*m});
        particleGun->SetParticleMomentumDirection({0,0,1});
        particleGun->SetParticleTime(0);
        particleGun->GeneratePrimaryVertex(anEvent);
        Run::GetInstance()->RecordPrimary(particleGun,1);
        auto event = Run::GetInstance()->GetEvent();
        event->Reset(); event->Pid=13; event->E=1*GeV; event->Z=-1*m;
        const auto mass=particleGun->GetParticleDefinition()->GetPDGMass();
        event->Pz=std::sqrt(1*GeV*(1*GeV+2*mass));
        return;
    }
    if(std::getenv("MUPOS_GAISSER_PLANE")) {
        G4double energy=0,theta=0,phi=0,x=0,y=0; G4ThreeVector direction;
        int proposals=0;
        for(;;) {
            if(++proposals>1000000) throw std::runtime_error("Plane source failed geometric acceptance");
            fGaisserTangGen->GenerateEvent(energy,theta,phi);
            x=fPlaneHalfX*(2*G4UniformRand()-1); y=fPlaneHalfY*(2*G4UniformRand()-1);
            direction={std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)};
            bool accepted=direction.z()>0;
            for(auto center:fPlaneGasCenters) for(double face:{center-1*mm,center+1*mm}) {
                const auto distance=(face-fPlaneSourceZ)/direction.z();
                accepted &= distance>0 && std::abs(x+distance*direction.x())<fDetectorHalfX-1.e-6*mm &&
                            std::abs(y+distance*direction.y())<fDetectorHalfY-1.e-6*mm;
            }
            if(accepted) break;
        }
        particleGun->SetParticleDefinition(particleTable->FindParticle(13));
        particleGun->SetParticleEnergy(energy*GeV); particleGun->SetParticlePosition({x,y,fPlaneSourceZ});
        particleGun->SetParticleMomentumDirection(direction); particleGun->SetParticleTime(0);
        particleGun->GeneratePrimaryVertex(anEvent); Run::GetInstance()->RecordPrimary(particleGun,2,proposals);
        auto event=Run::GetInstance()->GetEvent(); event->Reset(); event->Pid=13; event->E=energy*GeV;
        const auto mass=particleGun->GetParticleDefinition()->GetPDGMass();
        const auto momentum=std::sqrt(event->E*(event->E+2*mass))*direction;
        event->Px=momentum.x(); event->Py=momentum.y(); event->Pz=momentum.z();
        event->X=x; event->Y=y; event->Z=fPlaneSourceZ; event->T=0;
        return;
    }
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
    const bool zenithZ = std::getenv("MUPOS_GAISSER_Z_AXIS") != nullptr;
    if (zenithZ) {
        // Explicit opt-in coordinate convention; keep Gaisser-Tang energy/zenith law.
        particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
        particleGun->SetParticleMomentumDirection(G4ThreeVector(dirX,dirY,dirZ));
    }
    
    particleGun->SetParticleTime(0.0);
    
    particleGun->GeneratePrimaryVertex(anEvent);
    Run::GetInstance()->RecordPrimary(particleGun,zenithZ?3:0);
    
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
    if (zenithZ) {
        const auto actualPosition=particleGun->GetParticlePosition();
        const auto actualMomentum=std::sqrt(e*e-mass*mass)*particleGun->GetParticleMomentumDirection();
        event->X=actualPosition.x(); event->Y=actualPosition.y(); event->Z=actualPosition.z();
        event->Px=actualMomentum.x(); event->Py=actualMomentum.y(); event->Pz=actualMomentum.z();
    }
}
