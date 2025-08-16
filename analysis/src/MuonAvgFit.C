// root -l MuonAvgFit_AvgPerLayer.C
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <cmath>
#include <array>
#include <iostream>
#include <algorithm>

// 每层的几何 Z（单位：mm；和你之前保持一致）
static const std::array<double,6> layerZ = {
  -607.5, -426.5, -244.5, 237.5, 418.5, 600.5
};

struct P3 { double x,y,z; };

// 用 z 做自变量的线性回归：x = ax + bx*z, y = ay + by*z
// 返回归一化方向向量 v = (bx, by, 1) / ||...||
static std::array<double,3> FitDirByLinReg(const std::array<P3,3>& pts) {
  // 计算 z 的均值、Szz
  double zbar = (pts[0].z + pts[1].z + pts[2].z) / 3.0;
  double Szz = 0.0;
  for (auto &p: pts) { double dz = p.z - zbar; Szz += dz*dz; }
  if (Szz == 0) return {0,0,1}; // 极端情况，退化成 z 轴

  // bx = Sxz / Szz, by = Syz / Szz
  double Sxz = 0.0, Syz = 0.0;
  for (auto &p: pts) {
    double dz = p.z - zbar;
    Sxz += (p.x) * dz;
    Syz += (p.y) * dz;
  }
  double bx = Sxz / Szz;
  double by = Syz / Szz;

  // 方向向量 (bx, by, 1) 并归一
  double vx = bx, vy = by, vz = 1.0;
  double m = std::sqrt(vx*vx + vy*vy + vz*vz);
  if (m == 0) return {0,0,1};
  return {vx/m, vy/m, vz/m};
}

void MuonAvgFit(const char* inFile  = "../../build/root_file/CryMu.root",
                            const char* treeName= "tree",
                            const char* outFile = "../../build/root_file/mu_avg.root")
{
  // 1) 打开输入
  TFile fin(inFile);
  if (fin.IsZombie()) { std::cerr << "Cannot open " << inFile << "\n"; return; }
  TTree* tin = (TTree*)fin.Get(treeName);
  if (!tin) { std::cerr << "No TTree '" << treeName << "' in file\n"; return; }

  // 2) 分支读取（split 叶子）
  TTreeReader r(tin);
  TTreeReaderArray<int>    pid(r, "Edeps.Pid");
  TTreeReaderArray<int>    lid(r, "Edeps.Id");
  TTreeReaderArray<double> ex (r, "Edeps.X");
  TTreeReaderArray<double> ey (r, "Edeps.Y");

  // 3) 输出文件/树/直方图
  TFile fout(outFile, "RECREATE");
  TTree tout("tree", "theta12 from per-layer muon-averaged hits");
  Long64_t EventID = 0;
  double theta12 = std::nan("");     // 弧度
  double theta12_deg = std::nan(""); // 度
  tout.Branch("EventID",     &EventID);
  tout.Branch("theta12",     &theta12);
  tout.Branch("theta12_deg", &theta12_deg);

  #TH1D h("h_theta12", "theta12;theta12 [rad];Entries", 120, 0, 0.3);
  TH1D h("h_theta", "theta12;theta12 [rad];Entries", 120, 0, 0.3);

  // 4) 事件循环
  for (EventID = 0; r.Next(); ++EventID) {
    // 每层 μ 命中的累加和与计数
    double sumX[6] = {0}, sumY[6] = {0};
    int    cnt [6] = {0};

    const size_t n = std::min({ (size_t)pid.GetSize(),
                                (size_t)lid.GetSize(),
                                (size_t)ex.GetSize(),
                                (size_t)ey.GetSize() });
    for (size_t i=0; i<n; ++i) {
      int P = pid[(int)i];
      if (P!=13 && P!=-13) continue;       // 只要 μ
      int L = lid[(int)i];
      if (L<0 || L>=6) continue;           // 只关心 0..5 六层
      sumX[L] += ex[(int)i];
      sumY[L] += ey[(int)i];
      cnt [L] += 1;
    }

    // 要求上三层和下三层每层至少有 1 个 μ 命中
    bool top_ok = (cnt[0]>0 && cnt[1]>0 && cnt[2]>0);
    bool bot_ok = (cnt[3]>0 && cnt[4]>0 && cnt[5]>0);
    if (!(top_ok && bot_ok)) continue;

    // 形成三层的“平均点”
    std::array<P3,3> topPts, botPts;
    for (int k=0; k<3; ++k) {
      topPts[k] = { sumX[k]/cnt[k], sumY[k]/cnt[k], layerZ[k] };
    }
    for (int k=0; k<3; ++k) {
      int L = 3+k;
      botPts[k] = { sumX[L]/cnt[L], sumY[L]/cnt[L], layerZ[L] };
    }

    // 5) 用线性回归拟合两段方向
    auto vin  = FitDirByLinReg(topPts);
    auto vout = FitDirByLinReg(botPts);

    // 方向统一朝 +z，避免 π 的跳变
    if (vin[2]  < 0) { vin[0]  = -vin[0];  vin[1]  = -vin[1];  vin[2]  = -vin[2];  }
    if (vout[2] < 0) { vout[0] = -vout[0]; vout[1] = -vout[1]; vout[2] = -vout[2]; }

    // 6) 夹角
    double dot = vin[0]*vout[0] + vin[1]*vout[1] + vin[2]*vout[2];
    if (dot >  1) dot = 1;
    if (dot < -1) dot = -1;
    theta12 = std::acos(dot);                // rad
    theta12_deg = theta12 * 180.0 / M_PI;    // deg

    // 填充
    tout.Fill();
    h.Fill(theta12);
  }

  // 7) 写出
  h.Write();
  tout.Write();
  fout.Close();
  std::cout << "Done. Wrote tree 'tree' and hist 'h_theta12' to " << outFile << "\n";
}

