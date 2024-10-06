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
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "Object.hh"
#include "RunMessenger.hh"

Run::Run()
{
  fRunMessenger = new RunMessenger(this);
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
  G4double scoringHalfZ = fDetectorConstruction->GetScoringHalfZ();
  const std::vector<G4double> &scoringZs = fDetectorConstruction->GetScoringZs();

  fScoringHalfX = fDetectorConstruction->GetScoringHalfX();
  fScoringHalfY = fDetectorConstruction->GetScoringHalfY();
  fScoringZ = scoringHalfZ * 2;
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
  fTree->Branch("Scatters", new TClonesArray("Scatter"));

  // The params tree is only accessed here.
  TClonesArray Params("Params");
  TTree *params = new TTree("params", "params");
  params->Branch("Params", &Params);
  *((::Params *)Params.ConstructedAt(0)) = *fDetectorConstruction;
  fCellX = ((::Params *)Params.UncheckedAt(0))->CellX;
  fCellY = ((::Params *)Params.UncheckedAt(0))->CellY;
  fNCellX = ((::Params *)Params.UncheckedAt(0))->HalfNCellX * 2;
  fNCellY = ((::Params *)Params.UncheckedAt(0))->HalfNCellY * 2;
  fScoringOffsetX = -((::Params *)Params.UncheckedAt(0))->HalfNCellX * fCellX;
  fScoringOffsetY = -((::Params *)Params.UncheckedAt(0))->HalfNCellY * fCellY;
  params->Fill();
  params->Write(NULL, params->kOverwrite);
  params->SetBranchAddress("Params", NULL);
}

void Run::SaveTree()
{
  if(!fFile) { return; }
  fFile->cd();
  fTree->Write(NULL, TObject::kOverwrite);
  delete *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Scatters")->GetAddress();
  fFile->Close();
  fTree = NULL;
  fFile = NULL;
}

void Run::FillAndReset()
{
  auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  auto Edeps = *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();
  auto Scatters = *(TClonesArray **)fTree->GetBranch("Scatters")->GetAddress();

  // Sort the tracks by ID.
  std::vector<Track *> tracks;
  tracks.resize(Tracks->GetEntries());
  for(size_t i = 0; i < tracks.size(); ++i) tracks[i] = (Track *)(*Tracks)[i];
  sort(tracks.begin(), tracks.end(), [](Track *a, Track *b) { return a->Id < b->Id; });
  for(size_t i = 0; i < tracks.size(); ++i) (*Tracks)[i] = tracks[i];

  // Export Edeps.
  if(all_of(fStatus.begin(), fStatus.end(), [](bool b) { return b; })) {
    for(auto &[id, value] : fEdep) { *(::Edep *)Edeps->ConstructedAt(Edeps->GetEntries()) = { id, value }; }
    fTree->Fill();
    Edeps->Clear();
  }
  fStatus.assign(fStatus.size(), false);

  Tracks->Clear();
  fEdep.clear();
  Scatters->Clear();
}

void Run::AutoSave() { fTree->AutoSave("SaveSelf Overwrite"); }

void Run::AddStep(const G4Step *step)
{
  const G4ThreeVector &r = step->GetTrack()->GetPosition();
  G4double x = r.x(), y = r.y(), z = r.z();
  if(fabs(x) >= fScoringHalfX || fabs(y) >= fScoringHalfY) return;
  auto ub = std::upper_bound(fScoringMaxZs.begin(), fScoringMaxZs.end(), z);
  if(ub == fScoringMaxZs.end()) return;
  if(z < *ub - fScoringZ) return;

  G4double edep = step->GetTotalEnergyDeposit();
  if(edep == 0) return;

  Long64_t zid = ub - fScoringMaxZs.begin();
  Long64_t xid = (x - fScoringOffsetX) / fCellX;
  Long64_t yid = (y - fScoringOffsetY) / fCellY;
  Long64_t rid = zid * (fNCellX * fNCellY) + xid * fNCellY + yid;
  Long64_t pid = (uint32_t)step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  Long64_t id = (rid << 32) | pid;
  fStatus[zid] = true;
  fEdep[id] += edep;
}

void Run::AddTrack(const G4Track *track)
{
  //G4cout << __PRETTY_FUNCTION__ << ": " << track->GetTrackID()
  //  << "(" << track->GetParentID() << ")"
  //  << G4endl;
  auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  *(Track *)Tracks->ConstructedAt(Tracks->GetEntries()) = *track;
}

void Run::AddScatter(const G4Track *muon, const G4DynamicParticle *lp, const G4DynamicParticle *ln)
{
  auto Scatters = *(TClonesArray **)fTree->GetBranch("Scatters")->GetAddress();
  *(Scatter *)Scatters->ConstructedAt(Scatters->GetEntries()) = { muon, lp, ln };
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
