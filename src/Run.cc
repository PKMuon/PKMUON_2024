// 2020.5.8 by siguang wang (siguang@pku.edu.cn)

#include "Run.hh"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <syscall.h>
#include <unistd.h>

#include <filesystem>

#include "DetectorConstruction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "Object.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunMessenger.hh"

Run::Run()
{
  fRunMessenger = new RunMessenger(this);
  fPrimaryGeneratorAction = (PrimaryGeneratorAction *)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  fDetectorConstruction = (DetectorConstruction *)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  fRootName = "CryMu.root";
  fTree = NULL;
  fFile = NULL;
}

Run::~Run()
{
  SaveTree();
  delete fRunMessenger;
}

Run *Run::GetInstance()
{
  static Run run;
  return &run;
}

void Run::InitGeom()
{
  fScoringHalfX = fDetectorConstruction->GetScoringHalfX();
  fScoringHalfY = fDetectorConstruction->GetScoringHalfY();
  G4double scoringHalfZ = fDetectorConstruction->GetScoringHalfZ();
  const std::vector<G4double> &scoringZs = fDetectorConstruction->GetScoringZs();
  fScoringZ = scoringHalfZ * 2;
  EdepData::fScoringZs.assign(scoringZs.begin(), scoringZs.end());
  fScoringMaxZs = scoringZs;
  for(G4double &z : fScoringMaxZs) z += scoringHalfZ;
  fStatus.resize(fScoringMaxZs.size());
}

void Run::InitTree()
{
  using namespace std::filesystem;
  auto dirpath = path(fRootName.c_str()).parent_path();
  if(!dirpath.empty()) { create_directories(dirpath); }

  fFile = TFile::Open(fRootName, "RECREATE");
  fTree = new TTree("tree", "tree");
  fTree->Branch("Tracks", new TClonesArray("Track"));
  fTree->Branch("Edeps", new TClonesArray("Edep"));

  TClonesArray Cuts("Cuts");
  TTree *cuts = new TTree("cuts", "cuts");
  cuts->Branch("Cuts", &Cuts);
  *((::Cuts *)Cuts.ConstructedAt(0)) = *G4LogicalVolumeStore::GetInstance()->GetVolume("world");
  cuts->Fill();
  cuts->Write(NULL, cuts->kOverwrite);
  cuts->SetBranchAddress("Cuts", NULL);
}

void Run::SaveTree()
{
  if(!fFile) { return; }
  fFile->cd();
  fTree->Write(NULL, TObject::kOverwrite);
  delete *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();
  fFile->Close();
  fTree = NULL;
  fFile = NULL;
}

void Run::FillAndReset()
{
  auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  auto Edeps = *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();

  if(all_of(fStatus.begin(), fStatus.end(), [](bool b) { return b; })) {
    for(auto &[id, data] : fEdepData) { *(::Edep *)Edeps->ConstructedAt(Edeps->GetEntries()) = { id, data.EndSum() }; }
    fTree->Fill();
    Edeps->Clear();
  }
  fStatus.assign(fStatus.size(), false);

  Tracks->Clear();
  fEdepData.clear();
}

void Run::AutoSave() { fTree->AutoSave("SaveSelf Overwrite"); }

void Run::AddStep(const G4Step *step)
{
  const G4ThreeVector &r = step->GetTrack()->GetPosition();
  G4double x = r.x(), y = r.y(), z = r.z();
  if(fabs(x) > fScoringHalfX || fabs(y) > fScoringHalfY) return;
  auto ub = std::upper_bound(fScoringMaxZs.begin(), fScoringMaxZs.end(), z);
  if(ub == fScoringMaxZs.end()) return;
  if(z < *ub - fScoringZ) return;

  G4double edep = step->GetTotalEnergyDeposit();
  if(edep == 0) return;

  Long64_t layer = ub - fScoringMaxZs.begin();
  Long64_t pid = (uint32_t)step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  Long64_t id = (layer << 32) | pid;
  fStatus[layer] = true;
  fEdepData[id].Add(edep, x, y);
}

void Run::AddTrack(const G4Track *track)
{
  //G4cout << __PRETTY_FUNCTION__ << ": " << track->GetTrackID()
  //  << "(" << track->GetParentID() << ")"
  //  << ": primary=" << fPrimaryGeneratorAction->IsPrimary(track->GetTrackID())
  //  << G4endl;
  auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  *(Track *)Tracks->ConstructedAt(Tracks->GetEntries()) = *track;
}

uint64_t Run::GetThreadId()
{
#ifdef __APPLE__
  uint64_t tid;
  pthread_threadid_np(NULL, &tid);
  return tid;
#else  /* __APPLE__ */
  int64_t tid = syscall(SYS_gettid);
  if(tid < 0) {  // probably ENOSYS
    perror("gettid");
    exit(EXIT_FAILURE);
  }
  return tid;
#endif /* __APPLE__ */
}

uint64_t Run::GetSeed()
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch())
             .count()
      + GetThreadId();
}
