// 2020.5.8 by siguang wang (siguang@pku.edu.cn)

#include "Run.hh"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <syscall.h>
#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "DetectorConstruction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ParticleGun.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "Object.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunMessenger.hh"

// ── LUT check & Garfield++ auto-invoke ────────────────────────────────────
static bool CheckAndGenerateLUT()
{
  const char *lutPath = "rpc_cluster_params.root";
  if (std::ifstream(lutPath).good()) {
    G4cout << "[GarfieldLUT] Found existing LUT: " << lutPath << G4endl;
    return true;
  }
  G4cout << "[GarfieldLUT] LUT not found. Attempting Garfield++ auto-generation...\n";
  // Try to run garfield_rpc/rpc_cluster_sim from build dir
  const char *cmd =
#ifdef _WIN32
    "cd garfield_rpc\\build 2>nul && rpc_cluster_sim.exe "
    "--material ../../config/newrpc_material.yaml "
    "--volume ../../config/newrpc.yaml "
    "--output rpc_cluster_params.root";
#else
    "(cd garfield_rpc/build 2>/dev/null && ./rpc_cluster_sim "
    "--material ../../config/newrpc_material.yaml "
    "--volume ../../config/newrpc.yaml "
    "--output rpc_cluster_params.root)";
#endif
  int ret = std::system(cmd);
  if (ret != 0) {
    G4cerr << "[GarfieldLUT] Garfield++ invocation failed (ret=" << ret
           << "). Continuing without LUT.\n";
    return false;
  }
  if (!std::ifstream(lutPath).good()) {
    G4cerr << "[GarfieldLUT] Still no LUT after Garfield++ run.\n";
    return false;
  }
  G4cout << "[GarfieldLUT] LUT generated successfully.\n";
  return true;
}

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

  // ── Garfield++ LUT check ────────────────────────────────────────
  // The old radius LUT is not a calibrated signal response. Never launch costly
  // Garfield jobs implicitly in each Geant4 process.
  if (std::getenv("RPC_LEGACY_AUTOGENERATE")) fLutReady = CheckAndGenerateLUT();
  auto store = G4LogicalVolumeStore::GetInstance();
  auto gas = dynamic_cast<G4Box *>(store->GetVolume("newrpc_gas")->GetSolid());
  auto glass = dynamic_cast<G4Box *>(store->GetVolume("newrpc_glass")->GetSolid());
  auto film = dynamic_cast<G4Box *>(store->GetVolume("newrpc_insulating_film")->GetSolid());
  auto electrode = dynamic_cast<G4Box *>(store->GetVolume("newrpc_electrode")->GetSolid());
  if (!gas || !glass || !film || !electrode) throw std::runtime_error("RPC segment export requires configured box layers");
  fRpcGasHalfZ = gas->GetZHalfLength();
  // Current stack: readout board | film | electrode | glass | active gas.
  // GetScoringZs are board centers, NOT gas centers. Keep legacy branches intact.
  auto board = dynamic_cast<G4Box *>(store->GetVolume("newrpc_xy_readout_board")->GetSolid());
  if (!board) throw std::runtime_error("RPC segment export requires a box readout board");
  const double offset = board->GetZHalfLength() + 2*film->GetZHalfLength() +
                        2*electrode->GetZHalfLength() + 2*glass->GetZHalfLength() + fRpcGasHalfZ;
  fRpcGapCenters.clear();
  for (auto boardZ : scoringZs) fRpcGapCenters.push_back(boardZ + offset);
  std::sort(fRpcGapCenters.begin(), fRpcGapCenters.end(), std::greater<G4double>());
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
  fTree->Branch("Event", new TClonesArray("Event"));
  fRpcSegments = new TTree("rpcSegments", "Actual gas steps before legacy event trigger; mm ns eV");
  fRpcSegments->Branch("event", &fRpcSegment.event);
  fRpcSegments->Branch("layer", &fRpcSegment.layer);
  fRpcSegments->Branch("track", &fRpcSegment.track);
  fRpcSegments->Branch("parent", &fRpcSegment.parent);
  fRpcSegments->Branch("pid", &fRpcSegment.pid);
  fRpcSegments->Branch("root_track", &fRpcSegment.root_track);
#define RPC_BRANCH(name) fRpcSegments->Branch(#name, &fRpcSegment.name)
  RPC_BRANCH(x0_mm); RPC_BRANCH(y0_mm); RPC_BRANCH(z0_mm);
  RPC_BRANCH(x1_mm); RPC_BRANCH(y1_mm); RPC_BRANCH(z1_mm);
  RPC_BRANCH(t0_ns); RPC_BRANCH(t1_ns); RPC_BRANCH(edep_ev); RPC_BRANCH(nonion_ev);
  RPC_BRANCH(kinetic_ev); RPC_BRANCH(length_mm); RPC_BRANCH(gap_center_z_mm); RPC_BRANCH(gap_mm);
#undef RPC_BRANCH
  fRpcPrimaries = new TTree("rpcPrimaries", "Actual generated primary coordinates and directions; all events");
#define RPC_PRIMARY(name) fRpcPrimaries->Branch(#name, &fRpcPrimary.name)
  RPC_PRIMARY(event); RPC_PRIMARY(pid); RPC_PRIMARY(kinetic_ev); RPC_PRIMARY(x_mm); RPC_PRIMARY(y_mm);
  RPC_PRIMARY(z_mm); RPC_PRIMARY(dx); RPC_PRIMARY(dy); RPC_PRIMARY(dz); RPC_PRIMARY(t_ns);
  RPC_PRIMARY(source_mode); RPC_PRIMARY(source_proposals);
#undef RPC_PRIMARY
  fRpcLegacyEdeps = new TTree("rpcLegacyEdeps", "Original Edep map before four-layer trigger; mm eV");
#define RPC_LEGACY(name) fRpcLegacyEdeps->Branch(#name, &fRpcLegacyEdep.name)
  RPC_LEGACY(event); RPC_LEGACY(layer); RPC_LEGACY(pid); RPC_LEGACY(process);
  RPC_LEGACY(energy_ev); RPC_LEGACY(x_mm); RPC_LEGACY(y_mm); RPC_LEGACY(legacy_trigger);
#undef RPC_LEGACY
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
  fRpcSegments->Write(NULL, TObject::kOverwrite);
  fRpcPrimaries->Write(NULL, TObject::kOverwrite);
  fRpcLegacyEdeps->Write(NULL, TObject::kOverwrite);
  delete *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();
  delete *(TClonesArray **)fTree->GetBranch("Event")->GetAddress();
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
  fRpcSegments = nullptr;
  fRpcPrimaries = fRpcLegacyEdeps = nullptr;
}

void Run::FillAndReset()
{
  auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  auto Edeps = *(TClonesArray **)fTree->GetBranch("Edeps")->GetAddress();

  //// Sort the tracks by ID.
  std::vector<Track *> tracks;
  tracks.resize(Tracks->GetEntries());
  for(size_t i = 0; i < tracks.size(); ++i) tracks[i] = (Track *)(*Tracks)[i];
  sort(tracks.begin(), tracks.end(), [](Track *a, Track *b) { return a->Z > b->Z; });
  for(size_t i = 0; i < tracks.size(); ++i) (*Tracks)[i] = tracks[i];

  // Export Edeps.
  //bool trigger = (fStatus[0] && fStatus[1]) || (fStatus[2] && fStatus[3]) || (fStatus[1] && fStatus[2]);
  bool trigger = fStatus[0] && fStatus[1] && fStatus[2] && fStatus[3];
  for(const auto &p:fEdep) {
    const auto &key=p.first; const auto &value=p.second;
    if(value.Value<=0) continue;
    fRpcLegacyEdep.event=fIEvent; fRpcLegacyEdep.layer=key.Id;
    fRpcLegacyEdep.pid=key.Pid; fRpcLegacyEdep.process=key.Process;
    fRpcLegacyEdep.energy_ev=value.Value/eV;
    fRpcLegacyEdep.x_mm=value.X/value.Value/mm; fRpcLegacyEdep.y_mm=value.Y/value.Value/mm;
    fRpcLegacyEdep.legacy_trigger=trigger; fRpcLegacyEdeps->Fill();
  }
  if(trigger) {
    for(auto &edep : fEdep) { *(::Edep *)Edeps->ConstructedAt(Edeps->GetEntries()) = edep; }
    fTree->Fill();
    Edeps->Clear();
  }
  fStatus.assign(fStatus.size(), false);

  Tracks->Clear();
  fEdep.clear();
  fRpcParents.clear();
  ++fIEvent;
}

void Run::AutoSave() { fTree->AutoSave("SaveSelf Overwrite"); }

void Run::AddStep(const G4Step *step)
{
  fRpcParents[step->GetTrack()->GetTrackID()]=step->GetTrack()->GetParentID();
  RecordRpcSegment(step);
  const G4ThreeVector &r = step->GetTrack()->GetPosition();
  G4double x = r.x(), y = r.y(), z = r.z();
  G4double edep = step->GetTotalEnergyDeposit();
  auto ub = std::upper_bound(fScoringMaxZs.begin(), fScoringMaxZs.end(), z);

  bool xyRefuseCondition = fabs(x) >= fScoringHalfX || fabs(y) >= fScoringHalfY;
  bool zRefuseCondition = ub == fScoringMaxZs.end() || z < *ub - fScoringZ;
  bool energyRefuseCondition = edep <= 0;
  //bool pidRefuseCondition = step->GetTrack()->GetParticleDefinition()->GetPDGCharge() == 0;
  bool refuseCondition = xyRefuseCondition || zRefuseCondition || energyRefuseCondition;

  if (refuseCondition) return;

  Int_t zid = fScoringMaxZs.size() - 1 - (ub - fScoringMaxZs.begin());
  Int_t pid = (uint32_t)step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  Int_t process = -1;
  if(const G4VProcess *p = step->GetTrack()->GetCreatorProcess()) {
    auto it = fProcessMap.find(p->GetProcessName());
    if(it != fProcessMap.end()) process = it->second;
  }
  fStatus[zid] = true;
  fEdep[{ zid, pid, process}].Add(edep, x, y);

  auto Tracks = *(TClonesArray **)fTree->GetBranch("Tracks")->GetAddress();
  Track *existingTrack = nullptr;
  Int_t trackid = step->GetTrack()->GetTrackID();

  for (Int_t i = 0; i < Tracks->GetEntries(); ++i) {
    Track *track = (Track *)(*Tracks)[i];
    auto ubTrack = std::upper_bound(fScoringMaxZs.begin(), fScoringMaxZs.end(), track->Z);
    Int_t trackZid = fScoringMaxZs.size() - 1 - (ubTrack - fScoringMaxZs.begin()); 

    if (track->Id == trackid && trackZid == zid) {
      existingTrack = track;
      break;
    }
  }

  if (existingTrack) {
    G4double totalEdep = existingTrack->E + edep;
    existingTrack->X = (existingTrack->X * existingTrack->E + x * edep) / totalEdep;
    existingTrack->Y = (existingTrack->Y * existingTrack->E + y * edep) / totalEdep;
    existingTrack->Z = (existingTrack->Z * existingTrack->E + z * edep) / totalEdep;
    existingTrack->E = totalEdep;
  } else *(Track *)Tracks->ConstructedAt(Tracks->GetEntries()) = *(step->GetTrack());
  
}

void Run::RecordRpcSegment(const G4Step *step)
{
  if (!fRpcSegments) return;
  const auto pre = step->GetPreStepPoint();
  const auto post = step->GetPostStepPoint();
  if (pre->GetMaterial() != fDetectorConstruction->GetScoringGasVolume()->GetMaterial()) return;
  const auto mid = (pre->GetPosition() + post->GetPosition()) * .5;
  if (std::abs(mid.x()) >= fScoringHalfX || std::abs(mid.y()) >= fScoringHalfY) return;
  for (size_t i=0; i<fRpcGapCenters.size(); ++i) {
    if (std::abs(mid.z()-fRpcGapCenters[i]) > fRpcGasHalfZ) continue;
    const auto track = step->GetTrack();
    fRpcSegment.event=fIEvent; fRpcSegment.layer=static_cast<Int_t>(i);
    fRpcSegment.track=track->GetTrackID(); fRpcSegment.parent=track->GetParentID();
    fRpcSegment.pid=track->GetParticleDefinition()->GetPDGEncoding();
    int root=track->GetTrackID();
    while(fRpcParents.count(root) && fRpcParents.at(root)>0) root=fRpcParents.at(root);
    fRpcSegment.root_track=root;
    fRpcSegment.x0_mm=pre->GetPosition().x()/mm; fRpcSegment.y0_mm=pre->GetPosition().y()/mm;
    fRpcSegment.z0_mm=pre->GetPosition().z()/mm; fRpcSegment.x1_mm=post->GetPosition().x()/mm;
    fRpcSegment.y1_mm=post->GetPosition().y()/mm; fRpcSegment.z1_mm=post->GetPosition().z()/mm;
    fRpcSegment.t0_ns=pre->GetGlobalTime()/ns; fRpcSegment.t1_ns=post->GetGlobalTime()/ns;
    fRpcSegment.edep_ev=step->GetTotalEnergyDeposit()/eV;
    fRpcSegment.nonion_ev=step->GetNonIonizingEnergyDeposit()/eV;
    fRpcSegment.kinetic_ev=pre->GetKineticEnergy()/eV; fRpcSegment.length_mm=step->GetStepLength()/mm;
    fRpcSegment.gap_center_z_mm=fRpcGapCenters[i]/mm; fRpcSegment.gap_mm=2*fRpcGasHalfZ/mm;
    fRpcSegments->Fill(); return;
  }
}

void Run::RecordPrimary(G4ParticleGun *gun, G4int mode, G4int proposals)
{
  if(!fRpcPrimaries) return;
  const auto p=gun->GetParticlePosition(), d=gun->GetParticleMomentumDirection();
  fRpcPrimary.event=fIEvent; fRpcPrimary.pid=gun->GetParticleDefinition()->GetPDGEncoding();
  fRpcPrimary.kinetic_ev=gun->GetParticleEnergy()/eV;
  fRpcPrimary.x_mm=p.x()/mm; fRpcPrimary.y_mm=p.y()/mm; fRpcPrimary.z_mm=p.z()/mm;
  fRpcPrimary.dx=d.x(); fRpcPrimary.dy=d.y(); fRpcPrimary.dz=d.z();
  fRpcPrimary.t_ns=gun->GetParticleTime()/ns;
  fRpcPrimary.source_mode=mode; fRpcPrimary.source_proposals=proposals; fRpcPrimaries->Fill();
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
      //G4cout << __PRETTY_FUNCTION__ << ": " << particle->GetParticleName() << ", " << process->GetProcessName() << G4endl;
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
