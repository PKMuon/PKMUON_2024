//******************************************************************************
// PrimaryGeneratorAction.cc
//
// 1.00 JMV, LLNL, Jan-2007:  First version.
//******************************************************************************
//

#include "PrimaryGeneratorAction.hh"

#include <iomanip>

#include "DetectorConstruction.hh"
using namespace std;

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#ifndef CRY_DATA
#define CRY_DATA "../data"
#endif /* CRY_DATA */

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::PrimaryGeneratorAction(const char *inputfile)
{
  // define a particle gun
  particleGun = new G4ParticleGun();

  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(inputfile, std::ios::in);
  char buffer[1000];

  if(inputFile.fail()) {
    if(*inputfile != 0) {  //....only complain if a filename was given
      G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << G4endl;
    }
    InputState = -1;
  } else {
    std::string setupString("");
    while(!inputFile.getline(buffer, 1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup = new CRYSetup(setupString, CRY_DATA);

    gen = new CRYGenerator(setup);

    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState = 0;
  }
  // create a vector to store the CRY particle properties
  vect = new std::vector<CRYParticle *>;

  // Create the table containing all particle names
  particleTable = G4ParticleTable::GetParticleTable();

  // Create the messenger file
  gunMessenger = new PrimaryGeneratorMessenger(this);

  fNPrimary = -1;
  fDetectorMinZ = NAN;
  fDetectorHalfX = NAN;
  fDetectorHalfY = NAN;
  fHalfTimeWindow = 25 * ns;
}

//----------------------------------------------------------------------------//
PrimaryGeneratorAction::~PrimaryGeneratorAction() { }

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::Initialize(const DetectorConstruction *detectorConstruction)
{
  fDetectorMinZ = detectorConstruction->GetDetectorMinZ();
  fDetectorHalfX = detectorConstruction->GetDetectorHalfX();
  fDetectorHalfY = detectorConstruction->GetDetectorHalfY();
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::InputCRY() { InputState = 1; }

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::UpdateCRY(std::string *MessInput)
{
  CRYSetup *setup = new CRYSetup(*MessInput, CRY_DATA);

  gen = new CRYGenerator(setup);

  // set random number generator
  RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
  setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
  InputState = 0;
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::CRYFromFile(G4String newValue)
{
  // Read the cry input file
  std::ifstream inputFile;
  inputFile.open(newValue, std::ios::in);
  char buffer[1000];

  if(inputFile.fail()) {
    G4cout << "Failed to open input file " << newValue << G4endl;
    G4cout << "Make sure to define the cry library on the command line" << G4endl;
    InputState = -1;
  } else {
    std::string setupString("");
    while(!inputFile.getline(buffer, 1000).eof()) {
      setupString.append(buffer);
      setupString.append(" ");
    }

    CRYSetup *setup = new CRYSetup(setupString, CRY_DATA);

    gen = new CRYGenerator(setup);

    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState = 0;
  }
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  if(InputState != 0) {
    G4String *str = new G4String("CRY library was not successfully initialized");
    //G4Exception(*str);
    G4Exception("PrimaryGeneratorAction", "1", RunMustBeAborted, *str);
  }
  G4String particleName;
  vect->clear();
  gen->genEvent(vect);

  ////....debug output
  //G4cout << "\nEvent=" << anEvent->GetEventID() << " "
  //  << "CRY generated nparticles=" << vect->size()
  //  << G4endl;

  if(__builtin_expect(vect->empty(), false)) return;
  G4double t0 = vect->at(G4UniformRand() * vect->size())->t() * s;
  fNPrimary = 0;
  for(unsigned j = 0; j < vect->size(); j++) {
    particleName = CRYUtils::partName((*vect)[j]->id());

    ////....debug output
    //cout << scientific << setprecision(2) << "  " << setw(12) << left << particleName
    //  << "\tC=" << (*vect)[j]->charge()
    //  << "\tE[MeV]=" << (*vect)[j]->ke()
    //  << "\tX[m]=" << showpos << G4ThreeVector((*vect)[j]->x(), (*vect)[j]->y(), (*vect)[j]->z())
    //  << "\tA=" << G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w()) << noshowpos
    //  << "\tT[s]=" << (*vect)[j]->t()
    //  << endl;

    G4double t = (*vect)[j]->t() * s;
    if(fabs(t - t0) <= fHalfTimeWindow) {
      particleGun->SetParticleDefinition(particleTable->FindParticle((*vect)[j]->PDGid()));
      particleGun->SetParticleEnergy((*vect)[j]->ke() * MeV);
      //particleGun->SetParticlePosition(G4ThreeVector((*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->z()*m +
      //fDetectorMinZ));
      particleGun->SetParticlePosition({
          fDetectorHalfX * (2 * G4UniformRand() - 1),
          fDetectorHalfY * (2 * G4UniformRand() - 1),
          fDetectorMinZ,
      });
      particleGun->SetParticleMomentumDirection(G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), -(*vect)[j]->w()));
      particleGun->SetParticleTime((*vect)[j]->t() * s);
      particleGun->GeneratePrimaryVertex(anEvent);
      ++fNPrimary;
    }
    delete(*vect)[j];
  }
}
