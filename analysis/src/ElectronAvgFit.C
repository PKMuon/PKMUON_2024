// root -l ElectronAvgFit.C
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <cmath>
#include <array>
#include <algorithm>
#include <iostream>

static const std::array<double,6> layerZ = { -607.5, -426.5, -244.5, 237.5, 418.5, 600.5 };

struct P3 { double x,y,z; };

// 用 z 做自变量的线性回归，返回单位方向向量
static std::array<double,3> FitDir(const std::array<P3,3>& pts){
  double zbar = (pts[0].z + pts[1].z + pts[2].z)/3.0;
  double Szz=0, Sxz=0, Syz=0;
  for (auto& p: pts){
    double dz = p.z - zbar;
    Szz += dz*dz;
    Sxz += p.x*dz;
    Syz += p.y*dz;
  }
  if (Szz==0) return {0,0,1};
  double bx = Sxz/Szz, by = Syz/Szz; // x = ax + bx z, y = ay + by z
  double vx=bx, vy=by, vz=1.0;
  double m = std::sqrt(vx*vx+vy*vy+vz*vz); if (m==0) return {0,0,1};
  return {vx/m, vy/m, vz/m};
}

void ElectronAvgFit(const char* inFile  = "../../build/root_file/CryMu.root",
                    const char* treeName= "tree",
                    const char* outFile = "../../build/root_file/e_e_avg_fit.root")
{
  // 输入
  TFile fin(inFile);
  if (fin.IsZombie()) { std::cerr<<"Cannot open "<<inFile<<"\n"; return; }
  TTree* tin = (TTree*)fin.Get(treeName);
  if (!tin) { std::cerr<<"No tree '"<<treeName<<"'\n"; return; }

  // 分支
  TTreeReader r(tin);
  TTreeReaderArray<int>    pid(r, "Edeps.Pid");
  TTreeReaderArray<int>    lid(r, "Edeps.Id");
  TTreeReaderArray<double> ex (r, "Edeps.X");
  TTreeReaderArray<double> ey (r, "Edeps.Y");

  // 输出
  TFile fout(outFile, "RECREATE");
  TTree tout("tree","angles from per-layer averaged electrons");
  Long64_t EventID=0;
  double thetaEE = std::nan(""), thetaEE_deg = std::nan("");
  bool has_e_top=false, has_e_bot=false;
  tout.Branch("EventID", &EventID);
  tout.Branch("thetaEE", &thetaEE);
  tout.Branch("thetaEE_deg", &thetaEE_deg);
  tout.Branch("has_e_top", &has_e_top);
  tout.Branch("has_e_bot", &has_e_bot);

  TH1D h("h_thetaEE","thetaEE (e_in vs e_out);thetaEE [rad];Entries",120,0,0.6);

  for (EventID=0; r.Next(); ++EventID){
    // 每层电子加和/计数
    double sumX_e[6]={0}, sumY_e[6]={0}; int cnt_e[6]={0};

    size_t n = std::min({(size_t)pid.GetSize(), (size_t)lid.GetSize(),
                         (size_t)ex.GetSize(),  (size_t)ey.GetSize()});
    for (size_t i=0;i<n;++i){
      int L = lid[(int)i]; if (L<0 || L>=6) continue;
      int P = pid[(int)i];
      if (P==11 || P==-11){
        sumX_e[L] += ex[(int)i];
        sumY_e[L] += ey[(int)i];
        cnt_e[L]  += 1;
      }
    }

    has_e_top = (cnt_e[0]>0 && cnt_e[1]>0 && cnt_e[2]>0);
    has_e_bot = (cnt_e[3]>0 && cnt_e[4]>0 && cnt_e[5]>0);

    thetaEE = std::nan(""); thetaEE_deg = std::nan("");

    if (has_e_top && has_e_bot){
      std::array<P3,3> topPts{
        P3{sumX_e[0]/cnt_e[0], sumY_e[0]/cnt_e[0], layerZ[0]},
        P3{sumX_e[1]/cnt_e[1], sumY_e[1]/cnt_e[1], layerZ[1]},
        P3{sumX_e[2]/cnt_e[2], sumY_e[2]/cnt_e[2], layerZ[2]}
      };
      std::array<P3,3> botPts{
        P3{sumX_e[3]/cnt_e[3], sumY_e[3]/cnt_e[3], layerZ[3]},
        P3{sumX_e[4]/cnt_e[4], sumY_e[4]/cnt_e[4], layerZ[4]},
        P3{sumX_e[5]/cnt_e[5], sumY_e[5]/cnt_e[5], layerZ[5]}
      };
      auto vin  = FitDir(topPts);
      auto vout = FitDir(botPts);

      // 和之前一致：把 z 方向统一向下游（z>0）
      auto fixz=[&](std::array<double,3>& v){ if(v[2]<0){ v[0]=-v[0]; v[1]=-v[1]; v[2]=-v[2]; } };
      fixz(vin); fixz(vout);

      double d = vin[0]*vout[0]+vin[1]*vout[1]+vin[2]*vout[2];
      if (d> 1) d= 1; if (d<-1) d=-1;
      thetaEE = std::acos(d);
      thetaEE_deg = thetaEE * 180.0/M_PI;

      h.Fill(thetaEE);
    }

    tout.Fill();
  }

  h.Write();
  tout.Write();
  fout.Close();
  std::cout<<"Done. Wrote tree 'tree' and hist 'h_thetaEE' to "
           << outFile << "\n";
}

