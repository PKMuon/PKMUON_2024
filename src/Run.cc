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
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
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
  fIEvent = 0;
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
  //fTree->Branch("Tracks", new TClonesArray("Track"));
  fTree->Branch("Edeps", new TClonesArray("Edep"));
  fTree->Branch("Event", new TClonesArray("Event"));
  fTree->Branch("Scatters", new TClonesArray("Scatter"));
  (*(TClonesArray **)fTree->GetBranch("Event")->GetAddress())->ConstructedAt(0);

  fParams = new TTree("params", "params");
  fParams->Branch("Params", new TClonesArray("Params"));
  (*(TClonesArray **)fParams->GetBranch("Params")->GetAddress())->ConstructedAt(0);
  fParams->Branch("Processes", new TClonesArray("Process"));

  BuildProcessMap();
}

void Run::SaveTree()
{
  if(!fFile) { return; }
  fFile->cd();

  fTree->Write(NULL, TObject::kOverwrite);
  //delete *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Event")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Scatters")->GetAddress();
  fTree = NULL;

  Params *params = (Params *)(*(TClonesArray **)fParams->GetBranch("Params")->GetAddress())->At(0);
  params->NEvent = fIEvent;
  *params = *fDetectorConstruction;
  TClonesArray *Processes = *(TClonesArray **)fParams->GetBranch("Processes")->GetAddress();
  for(auto &[name, id] : fProcessMap) *(::Process *)Processes->ConstructedAt(Processes->GetEntries()) = { id, name };
  fParams->Fill();
  fParams->Write(NULL, TObject::kOverwrite);
  delete *(TClonesArray **)fParams->GetBranch("Params")->GetAddress();
  delete *(TClonesArray **)fParams->GetBranch("Processes")->GetAddress();
  fParams = NULL;

  fFile->Close();
  fFile = NULL;
}

void Run::FillAndReset()
{
  //auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  auto Edeps = *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();
  auto Scatters = *(TClonesArray **)fTree->GetBranch("Scatters")->GetAddress();

  //// Sort the tracks by ID.
  //std::vector<Track *> tracks;
  //tracks.resize(Tracks->GetEntries());
  //for(size_t i = 0; i < tracks.size(); ++i) tracks[i] = (Track *)(*Tracks)[i];
  //sort(tracks.begin(), tracks.end(), [](Track *a, Track *b) { return a->Id < b->Id; });
  //for(size_t i = 0; i < tracks.size(); ++i) (*Tracks)[i] = tracks[i];

  // Export Edeps.
  if(all_of(fStatus.begin(), fStatus.end(), [](bool b) { return b; })) {
    for(auto &edep : fEdep) { *(::Edep *)Edeps->ConstructedAt(Edeps->GetEntries()) = edep; }
    fTree->Fill();
    Edeps->Clear();
  }
  fStatus.assign(fStatus.size(), false);

  //Tracks->Clear();
  fEdep.clear();
  Scatters->Clear();
  ++fIEvent;
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

  Int_t zid = ub - fScoringMaxZs.begin();
  Int_t pid = (uint32_t)step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  Int_t process = -1;
  if(const G4VProcess *p = step->GetTrack()->GetCreatorProcess()) {
    auto it = fProcessMap.find(p->GetProcessName());
    if(it != fProcessMap.end()) process = it->second;
  }
  fStatus[zid] = true;
  fEdep[{ zid, pid, process }].Add(edep, x, y);
}

void Run::AddTrack([[maybe_unused]] const G4Track *track)
{
  //G4cout << __PRETTY_FUNCTION__ << ": " << track->GetTrackID()
  //  << "(" << track->GetParentID() << ")"
  //  << ": primary=" << fPrimaryGeneratorAction->IsPrimary(track->GetTrackID())
  //  << G4endl;
  //auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  //*(Track *)Tracks->ConstructedAt(Tracks->GetEntries()) = *track;
}

void Run::AddScatter(const G4Track *track, const G4ThreeVector &momentum)
{
  auto Scatters = *(TClonesArray **)fTree->GetBranch("Scatters")->GetAddress();
  *(Scatter *)Scatters->ConstructedAt(Scatters->GetEntries()) = {track, momentum.x(), momentum.y(), momentum.z()};
}

void Run::BuildProcessMap()
{
  auto &iter = *G4ParticleTable::GetParticleTable()->GetIterator();
  iter.reset();
  while(iter()) {
    G4ParticleDefinition *particle = iter.value();
    G4ProcessManager *processManager = particle->GetProcessManager();
    if(!processManager) continue;
    G4ProcessVector *processList = processManager->GetProcessList();
    if(!processList) continue;
    for(size_t i = 0; i < processList->size(); ++i) {
      G4VProcess *process = (*processList)[i];
      //G4cout << __PRETTY_FUNCTION__ << ": " << particle->GetParticleName() << ", " << process->GetProcessName() <<
      //G4endl;
      fProcessMap[process->GetProcessName()] = 0;  // Delay numbering to the end.
    }
  }
  size_t i = 0;
  for(auto &[name, id] : fProcessMap) {
    id = i++;
    G4cout << __PRETTY_FUNCTION__ << ": " << std::setw(3) << id << " " << name << G4endl;
  }
}

Event *Run::GetEvent() { return (Event *)(*(TClonesArray **)fTree->GetBranch("Event")->GetAddress())->At(0); }

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
