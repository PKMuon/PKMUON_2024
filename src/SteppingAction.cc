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

#include "SteppingAction.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "Run.hh"

SteppingAction::SteppingAction() { }

SteppingAction::~SteppingAction() { }

void SteppingAction::UserSteppingAction(const G4Step *step)
{
  SteppingScatter(step);
  Run::GetInstance()->AddStep(step);
}

void SteppingAction::SteppingScatter(const G4Step *step)
{
  const G4Track *trk = step->GetTrack();
  if(trk->GetParentID() != 0) return;
  // const std::string& name = trk->GetDefinition()->GetParticleName();
  // if (!(name == "mu+" || name == "mu-")) return;
  const int pdg = trk->GetDefinition()->GetPDGEncoding();
  if(pdg != 13 && pdg != -13) return;

  //获取次级粒子列表（返回类型为G4Track*（s））
  const auto *secs = step->GetSecondaryInCurrentStep();
  if(!secs || secs->empty()) return;

  const G4Track *e_trk = nullptr;
  // int e_cnt = 0,  other_cnt = 0;
  int e_num = 0;

  for(const auto *s : *secs) {
    if(!s) continue;

    // 只看当前 μ 子直接产生的次级
    if(s->GetParentID() != trk->GetTrackID()) {
      // other_cnt++; continue;
      continue;
    }

    //只看次级粒子为e的step
    const int spdg = s->GetDefinition()->GetPDGEncoding();
    if(spdg == 11) {
      e_trk = s;
      e_num++;
    } else {
      continue;
    }
  }
  if(!e_trk) return;
  // if (spdg == 11||spdg==-11) {  e_trk = s; e_cnt++; }      // 电子
  // else {other_cnt++;}}

  // 要求：恰好 1 个 e
  // if (!(e_cnt == 1 && other_cnt == 0)) return;

  // —— 出射 μ 用 post-step（不是次级）
  const G4StepPoint *post = step->GetPostStepPoint();
  const G4ParticleDefinition *muon = trk->GetDefinition();
  auto *mu_out = new G4DynamicParticle(muon, post->GetMomentum(), post->GetKineticEnergy());

  // —— 出射 e 用次级 track 自带的 DynamicParticle
  const G4DynamicParticle *e_out = e_trk->GetDynamicParticle();

  // 记录
  if(e_num >= 1) { Run::GetInstance()->AddScatter(trk, mu_out, e_out); }

  delete mu_out;
}

//测试接口
// void SteppingAction::SteppingScatter(const G4Step *step) {

// const G4Track* trk = step->GetTrack();
// const auto* secs = step->GetSecondaryInCurrentStep();
// if (!secs || secs->empty()) return;
// const G4Track* sec = (*secs)[0];

// const G4StepPoint* post = step->GetPostStepPoint();
// const G4ParticleDefinition* muon = trk->GetDefinition();
// auto* mu_out = new G4DynamicParticle(muon, post->GetMomentum(), post->GetKineticEnergy());

// // —— 出射 e 用次级 track 自带的 DynamicParticle
// const G4DynamicParticle* e_out= sec->GetDynamicParticle();

// // 记录
// Run::GetInstance()->AddScatter(trk, mu_out, e_out);

// delete mu_out;
//   }
