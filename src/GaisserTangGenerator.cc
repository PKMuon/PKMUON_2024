//******************************************************************************
// Gaisser-Tang model
// Tang et al. (2006)
//******************************************************************************

#include "GaisserTangGenerator.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>
#include <TH2F.h>
#include <TRandom.h>
#include <TCanvas.h>

const G4double GaisserTangGenerator::fPI = 3.14159265358979323846;
const G4double GaisserTangGenerator::fEpsilonPi = 115.0;     // GeV
const G4double GaisserTangGenerator::fEpsilonK = 850.0;      // GeV
const G4double GaisserTangGenerator::fGammaIndex = 2.7;

GaisserTangGenerator::GaisserTangGenerator()
    : fMuMass(0.10566),
      fMinEnergy(0.1), 
      fMaxEnergy(1000.0),
      fNeedUpdate(true),
      fFluxHistogram(nullptr) {
    
    
    G4ParticleDefinition* muon = G4ParticleTable::GetParticleTable()->FindParticle(13);
    if (muon) {
        fMuMass = muon->GetPDGMass() / GeV;
    } else {
        G4cout << "GaisserTangGenerator: Warning - Cannot get muon mass from Geant4, using default = " 
               << fMuMass << " GeV" << G4endl;
    }
    
    InitializeSegmentData();
}

GaisserTangGenerator::~GaisserTangGenerator() {    
    fFluxHistogram = nullptr;
}

G4double GaisserTangGenerator::ComputeCosThetaStar(G4double cosTheta) {
    G4double x = cosTheta;
    G4double p1 = 0.102573;
    G4double p2 = -0.068287;
    G4double p3 = 0.958633;
    G4double p4 = 0.0407253;
    G4double p5 = 0.817285;
    
    G4double numerator = x*x + p1*p1 + p2 * std::pow(x, p3) + p4 * std::pow(x, p5);
    G4double denominator = 1 + p1*p1 + p2 + p4;
    
    return std::sqrt(numerator / denominator);
}

G4double GaisserTangGenerator::ComputeDelta(G4double cosThetaStar) {
    G4double hF = 950.0;           // g/cm^2
    G4double lambdaN = 90.0;        // g/cm^2
    G4double stoppingPower = 2.06e-3; // GeV per g/cm^2
    
    return stoppingPower * (hF / cosThetaStar - lambdaN);
}

G4double GaisserTangGenerator::ComputeAT(G4double E, G4double cosTheta, 
                                          G4double cosThetaStar, G4double Eeff) {
    if (E > 100.0 / cosThetaStar) {
        return 1.0;
    } else {
        return 1.1 * std::pow((90.0 * std::sqrt(cosTheta + 0.001) / 1030.0), 
                               4.5 / (Eeff * cosThetaStar));
    }
}

G4double GaisserTangGenerator::GetDifferentialFlux(G4double E, G4double cosTheta) {
    if (cosTheta <= 0.0 || E <= 0.0) return 0.0;
    
    G4double cosThetaStar = ComputeCosThetaStar(cosTheta);
    G4double Eeff = E;
    G4double EforPower = E;
    G4double Delta = ComputeDelta(cosThetaStar);
    
    if (E <= 100.0 / cosThetaStar) {
        Eeff = E + Delta;
    }
    
    if (E <= 1.0 / cosThetaStar) {
        EforPower = (3.0 * E + 7.0 / cosThetaStar) / 10.0;
        Eeff = EforPower + Delta;
    }
    
    G4double AT = ComputeAT(E, cosTheta, cosThetaStar, Eeff);
    G4double rc = (E > 100.0 / cosThetaStar) ? 0.0 : 1.0e-4;
    
    G4double factor1 = 1.0 / (1.0 + 1.1 * Eeff * cosThetaStar / fEpsilonPi);
    G4double factor2 = 0.054 / (1.0 + 1.1 * Eeff * cosThetaStar / fEpsilonK);
    
    G4double flux = AT * 0.14 * std::pow(EforPower, -fGammaIndex) * (factor1 + factor2 + rc);
    
    return flux;
}

void GaisserTangGenerator::InitializeSegmentData() {

    if (fFluxHistogram) {
        fFluxHistogram->SetDirectory(0);
        delete fFluxHistogram;
        fFluxHistogram = nullptr;
    }

    G4int nBinsE = 200;      
    G4int nBinsTheta = 100;  
    
    std::vector<G4double> eBins(nBinsE + 1);
    G4double logEMin = std::log10(fMinEnergy);
    G4double logEMax = std::log10(fMaxEnergy);
    for (G4int i = 0; i <= nBinsE; ++i) {
        G4double logE = logEMin + (logEMax - logEMin) * i / nBinsE;
        eBins[i] = std::pow(10.0, logE);
    }
    
    std::vector<G4double> thetaBins(nBinsTheta + 1);
    for (G4int i = 0; i <= nBinsTheta; ++i) {
        thetaBins[i] = static_cast<G4double>(i) / nBinsTheta;
    }
    
    TString histName = "hFluxDistribution";
    TString histTitle = "Gaisser-Tang Flux Distribution;log_{10}(E/GeV);cos#theta;Flux";
    fFluxHistogram = new TH2F(histName, histTitle, 
                              nBinsE, &eBins[0], 
                              nBinsTheta, &thetaBins[0]);
    
    G4double totalFlux = 0.0;
    for (G4int iE = 1; iE <= nBinsE; ++iE) {
        G4double E = fFluxHistogram->GetXaxis()->GetBinCenter(iE);
        
        for (G4int iTheta = 1; iTheta <= nBinsTheta; ++iTheta) {
            G4double cosTheta = fFluxHistogram->GetYaxis()->GetBinCenter(iTheta);
            
            G4double flux = GetDifferentialFlux(E, cosTheta);

            G4double dE = fFluxHistogram->GetXaxis()->GetBinWidth(iE);
            G4double dCosTheta = fFluxHistogram->GetYaxis()->GetBinWidth(iTheta);
            G4double weight = flux * dE * dCosTheta;
           
            fFluxHistogram->SetBinContent(iE, iTheta, weight);
            totalFlux += weight;
        }
    }
    
    fFluxHistogram->Scale(1.0 / totalFlux);
/**/
    if (!fNeedUpdate) {
        TCanvas* canvas = new TCanvas("canvas", "Gaisser-Tang Flux Distribution", 800, 600);
        canvas->SetLogx();  
        canvas->SetLogz();  
        
        fFluxHistogram->GetXaxis()->SetTitle("Energy E_{#mu} (GeV)");
        fFluxHistogram->GetYaxis()->SetTitle("cos #theta");
        fFluxHistogram->GetZaxis()->SetTitle("Probability Density");
        
        fFluxHistogram->Draw("COLZ");
        
        canvas->Update();
        
        TString filename = TString::Format("GaisserTang_Flux2D.png");
        
        canvas->SaveAs(filename);
        
        delete canvas;
    }    
}

void GaisserTangGenerator::GenerateEvent(G4double &energy, G4double &theta, G4double &phi) {

    if (fNeedUpdate) {
        InitializeSegmentData();
        fNeedUpdate = false;
    }

    if (!fFluxHistogram) {
        G4cerr << "Error: Flux histogram not initialized!" << G4endl;
        return;
    }
    
    G4double x, y;
    fFluxHistogram->GetRandom2(x, y);
    
    energy = x;
    theta = std::acos(y);
    phi = 2.0 * CLHEP::pi * G4UniformRand();

}