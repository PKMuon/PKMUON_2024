// root -l MuonElectronAvgFit.C
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <cmath>
#include <array>
#include <iostream>
#include <algorithm>

static const std::array<double,6> layerZ = {
  -607.5, -426.5, -244.5, 237.5, 418.5, 600.5
};

struct P3 { double x,y,z; };

// 用 z 做自变量的线性回归：x=ax+bx z, y=ay+by z
// 返回单位方向向量 (bx, by, 1) 归一化
static std::array<double,3> FitDirByLinReg(const std::array<P3,3>& pts) {
  double zbar = (pts[0].z + pts[1].z + pts[2].z) / 3.0;
  double Szz = 0.0, Sxz = 0.0, Syz = 0.0;
  for (auto &p: pts) {
    const double dz = p.z - zbar;
    Szz += dz*dz;
    Sxz += p.x * dz;
    Syz += p.y * dz;
  }
  if (Szz == 0.0) return {0,0,1};

  const double bx = Sxz / Szz;
  const double by = Syz / Szz;

  double vx = bx, vy = by, vz = 1.0;
  const double m = std::sqrt(vx*vx + vy*vy + vz*vz);
  if (m == 0.0) return {0,0,1};
  return {vx/m, vy/m, vz/m};
}

void MuonElectronAvgFit(const char* inFile  = "../../build/root_file/CryMu.root",
                        const char* treeName= "tree",
                        const char* outFile = "../../build/root_file/mu_e_avg_fit.root")
{
  // 输入
  TFile fin(inFile);
  if (fin.IsZombie()) { std::cerr << "Cannot open " << inFile << "\n"; return; }
  TTree* tin = (TTree*)fin.Get(treeName);
  if (!tin) { std::cerr << "No TTree '" << treeName << "'\n"; return; }

  // 读分支
  TTreeReader r(tin);
  TTreeReaderArray<int>    pid(r, "Edeps.Pid");
  TTreeReaderArray<int>    lid(r, "Edeps.Id");
  TTreeReaderArray<double> ex (r, "Edeps.X");
  TTreeReaderArray<double> ey (r, "Edeps.Y");

  // 输出
  TFile fout(outFile, "RECREATE");
  TTree tout("tree", "angles from per-layer averaged mu/e points");
  Long64_t EventID = 0;
  double theta12 = std::nan(""), theta12_deg = std::nan("");
  double theta13 = std::nan(""), theta13_deg = std::nan("");
  double theta23 = std::nan(""), theta23_deg = std::nan("");
  bool has_mu_top=false, has_mu_bot=false, has_e_bot=false;

  tout.Branch("EventID", &EventID);
  tout.Branch("theta12", &theta12);     tout.Branch("theta12_deg", &theta12_deg);
  tout.Branch("theta13", &theta13);     tout.Branch("theta13_deg", &theta13_deg);
  tout.Branch("theta23", &theta23);     tout.Branch("theta23_deg", &theta23_deg);
  tout.Branch("has_mu_top", &has_mu_top);
  tout.Branch("has_mu_bot", &has_mu_bot);
  tout.Branch("has_e_bot",  &has_e_bot);

  TH1D h12("h_theta12", "theta12 (mu_in vs mu_out);theta12 [rad];Entries", 120, 0, 0.3);
  TH1D h13("h_theta13", "theta13 (mu_in vs e_out);theta13 [rad];Entries",   120, 0, 0.6);
  TH1D h23("h_theta23", "theta23 (mu_out vs e_out);theta23 [rad];Entries",  120, 0, 0.6);

  // 事件循环
  for (EventID = 0; r.Next(); ++EventID) {
    // 每层 μ/e 的累加与计数
    double sumX_mu[6] = {0}, sumY_mu[6] = {0}; int cnt_mu[6] = {0};
    double sumX_e [6] = {0}, sumY_e [6] = {0}; int cnt_e [6] = {0};

    const size_t n = std::min({ (size_t)pid.GetSize(),
                                (size_t)lid.GetSize(),
                                (size_t)ex.GetSize(),
                                (size_t)ey.GetSize() });
    for (size_t i=0; i<n; ++i) {
      const int L = lid[(int)i];
      if (L < 0 || L >= 6) continue;
      const int P = pid[(int)i];
      if (P==13 || P==-13) { // mu
        sumX_mu[L] += ex[(int)i];
        sumY_mu[L] += ey[(int)i];
        cnt_mu[L]  += 1;
      } else if (P==11 || P==-11) { // e
        sumX_e[L] += ex[(int)i];
        sumY_e[L] += ey[(int)i];
        cnt_e[L]  += 1;
      }
    }

    // 条件：上三层要有 μ；下三层 μ/e 分别看各自是否有
    has_mu_top = (cnt_mu[0]>0 && cnt_mu[1]>0 && cnt_mu[2]>0);
    has_mu_bot = (cnt_mu[3]>0 && cnt_mu[4]>0 && cnt_mu[5]>0);
    has_e_bot  = (cnt_e [3]>0 && cnt_e [4]>0 && cnt_e [5]>0);

    // 拟合入射 μ 方向
    std::array<double,3> v_in{0,0,1};
    if (has_mu_top) {
      std::array<P3,3> topPts{
        P3{sumX_mu[0]/cnt_mu[0], sumY_mu[0]/cnt_mu[0], layerZ[0]},
        P3{sumX_mu[1]/cnt_mu[1], sumY_mu[1]/cnt_mu[1], layerZ[1]},
        P3{sumX_mu[2]/cnt_mu[2], sumY_mu[2]/cnt_mu[2], layerZ[2]}
      };
      v_in = FitDirByLinReg(topPts);
    }

    // 拟合出射 μ 方向
    std::array<double,3> v_out_mu{0,0,1};
    if (has_mu_bot) {
      std::array<P3,3> botPtsMu{
        P3{sumX_mu[3]/cnt_mu[3], sumY_mu[3]/cnt_mu[3], layerZ[3]},
        P3{sumX_mu[4]/cnt_mu[4], sumY_mu[4]/cnt_mu[4], layerZ[4]},
        P3{sumX_mu[5]/cnt_mu[5], sumY_mu[5]/cnt_mu[5], layerZ[5]}
      };
      v_out_mu = FitDirByLinReg(botPtsMu);
    }

    // 拟合出射 e 方向
    std::array<double,3> v_out_e{0,0,1};
    if (has_e_bot) {
      std::array<P3,3> botPtsE{
        P3{sumX_e[3]/cnt_e[3], sumY_e[3]/cnt_e[3], layerZ[3]},
        P3{sumX_e[4]/cnt_e[4], sumY_e[4]/cnt_e[4], layerZ[4]},
        P3{sumX_e[5]/cnt_e[5], sumY_e[5]/cnt_e[5], layerZ[5]}
      };
      v_out_e = FitDirByLinReg(botPtsE);
    }

    // 为避免 π 翻转，方向统一取 z>0
    auto fixz = [](std::array<double,3>& v){ if (v[2]<0){ v[0]=-v[0]; v[1]=-v[1]; v[2]=-v[2]; } };
    fixz(v_in); fixz(v_out_mu); fixz(v_out_e);

    // 夹角工具
    auto angle = [](const std::array<double,3>& a, const std::array<double,3>& b){
      double d = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
      if (d> 1) d= 1; if (d<-1) d=-1; return std::acos(d);
    };

    // 计算并填充（存在才算）
    theta12 = std::nan(""); theta12_deg = std::nan("");
    theta13 = std::nan(""); theta13_deg = std::nan("");
    theta23 = std::nan(""); theta23_deg = std::nan("");

    if (has_mu_top && has_mu_bot) {
      theta12 = angle(v_in, v_out_mu);
      theta12_deg = theta12 * 180.0 / M_PI;
      h12.Fill(theta12);
    }
    if (has_mu_top && has_e_bot) {
      theta13 = angle(v_in, v_out_e);
      theta13_deg = theta13 * 180.0 / M_PI;
      h13.Fill(theta13);
    }
    if (has_mu_bot && has_e_bot) {
      theta23 = angle(v_out_mu, v_out_e);
      theta23_deg = theta23 * 180.0 / M_PI;
      h23.Fill(theta23);
    }

    tout.Fill();
  }

  // 写出
  h12.Write(); h13.Write(); h23.Write();
  tout.Write();
  fout.Close();
  std::cout << "Done. Wrote tree 'tree' and hists {h_theta12,h_theta13,h_theta23} to "
            << outFile << "\n";
}

