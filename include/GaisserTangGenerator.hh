//******************************************************************************
// Gaisser-Tang model
// Tang et al. (2006)
//******************************************************************************

#ifndef GaisserTangGenerator_h
#define GaisserTangGenerator_h 1

#include "globals.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include <cmath>
#include <vector>
#include <TH2F.h>
#include <TRandom.h>

class GaisserTangGenerator {
public:
    GaisserTangGenerator();
    ~GaisserTangGenerator();
    
    void GenerateEvent(G4double &energy, G4double &theta, G4double &phi);
    
    void SetMinEnergy(G4double e) { 
        fMinEnergy = e;
        fNeedUpdate = true;
    }
    void SetMaxEnergy(G4double e) { 
        fMaxEnergy = e;
        fNeedUpdate = true;
    }
    
    G4double GetMinEnergy() const { return fMinEnergy; }
    G4double GetMaxEnergy() const { return fMaxEnergy; }
    
    G4double GetDifferentialFlux(G4double E, G4double cosTheta);

    
private:

    G4double ComputeCosThetaStar(G4double cosTheta);
    G4double ComputeDelta(G4double cosThetaStar);
    G4double ComputeAT(G4double E, G4double cosTheta, G4double cosThetaStar, G4double Eeff);

    void InitializeSegmentData();
    
    static const G4double fPI;
    static const G4double fEpsilonPi;     
    static const G4double fEpsilonK;      
    static const G4double fGammaIndex;   
    
    G4double fMuMass;         
    G4double fMinEnergy;       
    G4double fMaxEnergy;
    bool fNeedUpdate;     
           
    TH2F* fFluxHistogram;
    TRandom* fRandom;

};

#endif