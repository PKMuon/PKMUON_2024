//2020.5.8 by Siguang WANG ( siguang@pku.edu.cn)
//

#ifndef GEANT4_INTRODUCTION_RUN_HH
#define GEANT4_INTRODUCTION_RUN_HH 1



#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "globals.hh"


class RunMessenger;

class Run {
public:
  Run();
  static Run* GetInstance();
  virtual ~Run();
  void initTree();
  void saveTree();
  void Fill();
  void AddEnergy1(double Eng1);
 
  void RecPrimPartEng(double fGunEnergy){_GunEng =fGunEnergy; };
  void SetROOTFullName(G4String NewRootFileFullName){ rootFileName = NewRootFileFullName;};
 
  void  AddX(G4double inX);
  void  AddY(G4double inY);
  void  AddZ(G4double inZ);
  void  AddR(G4double inR);

  void  AddReadoutEdepX(G4double X);
  void  AddReadoutEdepY(G4double Y);
  void  AddReadoutEdepZ(G4double Z);
  void  AddReadoutEdep(G4double E);
  void  AddReadoutE(G4double E);
  void  AddReadoutTrkid(G4int i);
  void  AddReadoutTrkparentid(G4int i);
  void AddPbEdepX(G4double X);
  void AddPbEdepY(G4double Y);
  void AddPbEdepZ(G4double Z);
  void  AddPbTrkid(G4int i);

  void  AddPx(G4double X);
  void  AddPy(G4double Y);
  void  AddPz(G4double Z);

  void  SetPos(G4double inX, G4double inY,G4double inZ);
  void  SetPxyz(G4double inX, G4double inY,G4double inZ);
  void  SetPosF(G4double fiX, G4double fiY,G4double fiZ);
  void  SetPxyzF(G4double fiX, G4double fiY,G4double fiZ);
  void  AddEdep(G4double inE); 
  void AddTrkID(G4int iTrkID);
  void SetCosTheta(G4double inCosTh){cosTh=inCosTh;};
  Double_t GetCosTheta(){return cosTh;};
  void SetStatus(bool instatus){status=instatus;};
  bool GetStatus(){return status;};
  TTree *GetTree(){return _tree;};
  void AutoSave(){_tree->AutoSave("SaveSelf Overwrite");};
  void AutoSave(TString opt){_tree->AutoSave(opt.Data());};
  void FillEdep(G4double Edep, G4double iR);
  void SetEdepZero();
  void ClearAll();
  G4double GetPx(){return px0;};
  G4double GetPy(){return py0;};
  G4double GetPz(){return pz0;};
  void SetGemTrkID(int i, int ID);
  void SetGemTrkMID(int i, int MID);
  void SetGemTrkPx(int i, double Px);
  void SetGemTrkPy(int i, double Py);
  void SetGemTrkPz(int i, double Pz);
  void SetGemTrkE(int i, double E);
  void SetGemTrkEdep(int i, double Edep);
  void SetGemTrkX(int i, double X);
  void SetGemTrkY(int i, double Y);
  void SetGemTrkZ(int i, double Z);
  void SetGemTrkStatus(int i, bool status);
  bool GetGemTrkStatus(int i){return GemTrkStatus[i];};
  void SetGemTruthZ(int i, double Z);
  double GetGemTruthZ(int i){return GemTruthZ[i];};
  
private:
  G4String rootFileName;
  TTree *_tree;
  TFile *_file;
  Double_t _energyDisFirst;
  Double_t _GunEng,totEdep;
  Double_t X0,Y0,Z0,px0,py0,pz0;
  Double_t X1,Y1,Z1,px1,py1,pz1;
  Double_t cosTh;
  bool status;
  std::vector<Double_t> vX; //Edep X
  std::vector<Double_t> vY; //Edep Y
  std::vector<Double_t> vZ; //Edep Z
  std::vector<Double_t> vEdep; //Edep value
  std::vector<Double_t> vReadoutPosX; //Edep X in readout bar
  std::vector<Double_t> vReadoutPosY; //Edep Y in readout bar
  std::vector<Double_t> vReadoutPosZ; //Edep Z in readout bar
  std::vector<Double_t> vReadoutE; //Total Energy in readout bar
  std::vector<Double_t> vReadoutEdep; //Edep Energy in readout bar
  std::vector<Int_t> vReadoutTrkid; //Trk id in readout bar
  std::vector<Int_t> vReadoutTrkparentid; //Trk id in readout bar
  std::vector<Double_t> vPbPosX; //Edep X in Pb Box
  std::vector<Double_t> vPbPosY; //Edep Y in Pb Box
  std::vector<Double_t> vPbPosZ; //Edep Z in Pb Box
  std::vector<Int_t> vPbTrkid; //Edep Energy in readout bar
  std::vector<Double_t> vPx;
  std::vector<Double_t> vPy;
  std::vector<Double_t> vPz;
  std::vector<Int_t> vTrkID; //TrackID

  Int_t GemTrkID[3]={0};
  Int_t GemTrkMID[3]={0};
  Double_t GemTrkPx[3]={0};
  Double_t GemTrkPy[3]={0};
  Double_t GemTrkPz[3]={0};
  Double_t GemTrkE[3]={0};
  Double_t GemTrkEdep[3]={0};
  Double_t GemTrkX[3]={0};
  Double_t GemTrkY[3]={0};
  Double_t GemTrkZ[3]={0};
  bool GemTrkStatus[3]={false};
  
  Double_t GemTruthZ[3]={0};

  RunMessenger* fRunMessenger;
    
};

#endif //GEANT4_INTRODUCTION_RUN_H
