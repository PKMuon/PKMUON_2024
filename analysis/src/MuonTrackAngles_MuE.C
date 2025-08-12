// #include <TFile.h>
// #include <TTree.h>
// #include <TTreeReader.h>
// #include <TTreeReaderArray.h>
// #include <iostream>
// #include <map>
// #include <vector>
// #include <set>
// #include <array>
// #include <cmath>
// #include <algorithm>
// #include <limits>

// struct Hit { double x,y,z; };

// static const std::map<int,double> layerZ = {
//   {0, -607.5},
//   {1, -426.5},
//   {2, -244.5},
//   {3,  237.5},
//   {4,  418.5},
//   {5,  600.5},
// };

// // PCA 主方向拟合（幂迭代）
// static std::array<double,3> fitDir(const std::vector<Hit>& pts) {
//   double cx=0, cy=0, cz=0;
//   for (auto& p: pts) cx+=p.x, cy+=p.y, cz+=p.z;
//   cx/=pts.size(); cy/=pts.size(); cz/=pts.size();

//   auto S = [&](const std::array<double,3>& v){
//     std::array<double,3> r{0,0,0};
//     for (auto& p: pts) {
//       double dx=p.x-cx, dy=p.y-cy, dz=p.z-cz;
//       double s = dx*v[0]+dy*v[1]+dz*v[2];
//       r[0]+=dx*s; r[1]+=dy*s; r[2]+=dz*s;
//     }
//     return r;
//   };

//   std::array<double,3> v{1,0,0};
//   for (int it=0; it<20; ++it) {
//     v = S(v);
//     double m = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
//     if (m==0) break;
//     v[0]/=m; v[1]/=m; v[2]/=m;
//   }
//   return v;
// }

// static inline double clamp(double x, double a, double b){ return x<a? a : (x>b? b : x); }
// static inline double angle_between(const std::array<double,3>& a, const std::array<double,3>& b){
//   double c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
//   return std::acos(clamp(c, -1.0, 1.0)); // 弧度
// }

// void MuonTrackAngles_MuE(const char* inFile="../../build/root_file/CryMu.root",
//                     const char* treeName="tree",
//                     const char* outFile="../../build/root_file/mu_costheta_filtered.root")
// {
//   // 打开输入
//   TFile fin(inFile);
//   if (fin.IsZombie()) { std::cerr<<"Cannot open "<<inFile<<"\n"; return; }
//   TTree* tin = (TTree*)fin.Get(treeName);
//   if (!tin) { std::cerr<<"No tree: "<<treeName<<"\n"; return; }

//   // 读分支
//   TTreeReader r(tin);
//   TTreeReaderArray<int>    pid(r, "Edeps.Pid");
//   TTreeReaderArray<int>    lid(r, "Edeps.Id");
//   TTreeReaderArray<double> ex (r, "Edeps.X");
//   TTreeReaderArray<double> ey (r, "Edeps.Y");

//   // 输出
//   TFile fout(outFile, "RECREATE");
//   TTree tout("tree","angles (rad) with same filter as previous code");
//   Long64_t EventID = 0;
//   double theta12 = std::numeric_limits<double>::quiet_NaN(); // 入 μ vs 出 μ
//   double theta13 = std::numeric_limits<double>::quiet_NaN(); // 入 μ vs 出 e
//   double theta23 = std::numeric_limits<double>::quiet_NaN(); // 出 μ vs 出 e
//   tout.Branch("EventID", &EventID);
//   tout.Branch("theta12", &theta12);
//   tout.Branch("theta13", &theta13);
//   tout.Branch("theta23", &theta23);

//   // 事件循环
//   for (EventID=0; r.Next(); ++EventID) {

//     // 每层收集: 该层出现过的 PID 集合 & 该层的 μ / e 命中点
//     std::map<int, std::set<int>> pidsPerLayer;
//     std::map<int, std::vector<Hit>> muHitsPerLayer;
//     std::map<int, std::vector<Hit>> eHitsPerLayer;

//     const size_t n = std::min({ (size_t)pid.GetSize(),
//                                 (size_t)lid.GetSize(),
//                                 (size_t)ex.GetSize(),
//                                 (size_t)ey.GetSize() });

//     for (size_t i=0; i<n; ++i) {
//       int L = lid[(int)i];
//       auto itZ = layerZ.find(L);
//       if (itZ == layerZ.end()) continue; // 只关心 0..5 层
//       int P = pid[(int)i];
//       pidsPerLayer[L].insert(P);
//       if (P==13 || P==-13) {
//         muHitsPerLayer[L].push_back({ ex[(int)i], ey[(int)i], itZ->second });
//       } else if (P==11 || P==-11) {
//         eHitsPerLayer[L].push_back({ ex[(int)i], ey[(int)i], itZ->second });
//       }
//     }

//     // 过滤条件：完全沿用你“刚才那份”
//     // 顶三层(0,1,2)：三层都有命中，且 PID 只能是 ±13
//     bool top_ok = true;
//     for (int L=0; L<=2; ++L) {
//       auto it = pidsPerLayer.find(L);
//       if (it == pidsPerLayer.end()) { top_ok=false; break; }
//       for (int P : it->second) { if (!(P==13 || P==-13)) { top_ok=false; break; } }
//       if (!top_ok) break;
//     }
//     // 底三层(3,4,5)：三层都有命中，且 PID 只能是 ±13 或 ±11
//     bool bot_ok = true;
//     for (int L=3; L<=5; ++L) {
//       auto it = pidsPerLayer.find(L);
//       if (it == pidsPerLayer.end()) { bot_ok=false; break; }
//       for (int P : it->second) {
//         if (!(P==13 || P==-13 || P==11 || P==-11)) { bot_ok=false; break; }
//       }
//       if (!bot_ok) break;
//     }
//     if (!(top_ok && bot_ok)) continue;

//     // 拟合用点：入 μ（0-2）；出 μ（3-5）；出 e（3-5）
//     std::vector<Hit> inMu, outMu, outE;
//     for (int L=0; L<=2; ++L) {
//       auto it = muHitsPerLayer.find(L);
//       if (it != muHitsPerLayer.end()) {
//         inMu.insert(inMu.end(), it->second.begin(), it->second.end());
//       }
//     }
//     for (int L=3; L<=5; ++L) {
//       auto itM = muHitsPerLayer.find(L);
//       if (itM != muHitsPerLayer.end()) {
//         outMu.insert(outMu.end(), itM->second.begin(), itM->second.end());
//       }
//       auto itE = eHitsPerLayer.find(L);
//       if (itE != eHitsPerLayer.end()) {
//         outE.insert(outE.end(), itE->second.begin(), itE->second.end());
//       }
//     }

//     // 要求与“刚才那份”一致：入/出 μ 拟合都得 ≥3 点，否则直接跳过事件
//     if (inMu.size() < 3 || outMu.size() < 3) continue;

//     // 方向
//     auto vIn   = fitDir(inMu);
//     auto vOutM = fitDir(outMu);
//     // theta12：一定能算
//     theta12 = angle_between(vIn, vOutM);

//     // theta13/theta23：只有当对应段 ≥3 点时才算，否则保留 NaN
//     if (outE.size() >= 3) {
//       auto vOutE = fitDir(outE);
//       theta13 = angle_between(vIn,   vOutE);
//       theta23 = angle_between(vOutM, vOutE);
//     } else {
//       theta13 = std::numeric_limits<double>::quiet_NaN();
//       theta23 = std::numeric_limits<double>::quiet_NaN();
//     }

//     tout.Fill();
//   }

//   tout.Write();
//   fout.Close();
//   std::cout << "Saved to " << outFile
//             << " (branches: theta12/theta13/theta23 in radians; selection identical to previous code)\n";
// }

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>

struct Hit { double x,y,z; };

// 层中心 Z（单位随你的 X/Y 一致）
static const std::map<int,double> layerZ = {
  {0, -607.5},
  {1, -426.5},
  {2, -244.5},
  {3,  237.5},
  {4,  418.5},
  {5,  600.5},
};

// PCA 主方向拟合（幂迭代）
static std::array<double,3> fitDir(const std::vector<Hit>& pts) {
  double cx=0, cy=0, cz=0;
  for (auto& p: pts) cx+=p.x, cy+=p.y, cz+=p.z;
  cx/=pts.size(); cy/=pts.size(); cz/=pts.size();

  auto S = [&](const std::array<double,3>& v){
    std::array<double,3> r{0,0,0};
    for (auto& p: pts) {
      double dx=p.x-cx, dy=p.y-cy, dz=p.z-cz;
      double s = dx*v[0]+dy*v[1]+dz*v[2];
      r[0]+=dx*s; r[1]+=dy*s; r[2]+=dz*s;
    }
    return r;
  };

  std::array<double,3> v{1,0,0};
  for (int it=0; it<20; ++it) {
    v = S(v);
    double m = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    if (m==0) break;
    v[0]/=m; v[1]/=m; v[2]/=m;
  }
  return v;
}

static inline double clamp(double x, double a, double b){ return x<a? a : (x>b? b : x); }
static inline double angle_between(const std::array<double,3>& a, const std::array<double,3>& b){
  double c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return std::acos(clamp(c, -1.0, 1.0)); // 弧度
}

// 统一方向到 +z，避免 PCA 方向符号翻转
static inline void align_to_plus_z(std::array<double,3>& v) {
  if (v[2] < 0) { v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2]; }
}

void MuonTrackAngles_MuE(const char* inFile ="../../build/root_file/CryMu.root",
                         const char* treeName="tree",
                         const char* outFile="../../build/root_file/mu_angles_mue.root")
{
  // 输入
  TFile fin(inFile);
  if (fin.IsZombie()) { std::cerr<<"Cannot open "<<inFile<<"\n"; return; }
  TTree* tin = (TTree*)fin.Get(treeName);
  if (!tin) { std::cerr<<"No tree: "<<treeName<<"\n"; return; }

  // 分支
  TTreeReader r(tin);
  TTreeReaderArray<int>    pid(r, "Edeps.Pid");
  TTreeReaderArray<int>    lid(r, "Edeps.Id");
  TTreeReaderArray<double> ex (r, "Edeps.X");
  TTreeReaderArray<double> ey (r, "Edeps.Y");

  // 输出
  TFile fout(outFile, "RECREATE");
  TTree tout("tree","theta12/13/23 (rad) with mu-only top and mu/e-only bottom");
  Long64_t EventID = 0;
  double theta12 = std::numeric_limits<double>::quiet_NaN(); // 入 μ vs 出 μ
  double theta13 = std::numeric_limits<double>::quiet_NaN(); // 入 μ vs 出 e
  double theta23 = std::numeric_limits<double>::quiet_NaN(); // 出 μ vs 出 e
  tout.Branch("EventID", &EventID);
  tout.Branch("theta12", &theta12);
  tout.Branch("theta13", &theta13);
  tout.Branch("theta23", &theta23);

  // 事件循环
  for (EventID=0; r.Next(); ++EventID) {

    // 每层：出现过的 PID 集合；每层 μ/e 命中点
    std::map<int, std::set<int>>    pidsPerLayer;
    std::map<int, std::vector<Hit>> muHitsPerLayer;
    std::map<int, std::vector<Hit>> eHitsPerLayer;

    const size_t n = std::min({ (size_t)pid.GetSize(),
                                (size_t)lid.GetSize(),
                                (size_t)ex.GetSize(),
                                (size_t)ey.GetSize() });

    for (size_t i=0; i<n; ++i) {
      int L = lid[(int)i];
      auto itZ = layerZ.find(L);
      if (itZ == layerZ.end()) continue; // 只关心 0..5 层
      int P = pid[(int)i];
      pidsPerLayer[L].insert(P);//把当前命中点的粒子 PID（P）加入到 L 层的 PID 集合中，如果之前没出现过，就记录下来；如果已经出现过，就忽略（去重）
      if (P==13 || P==-13) {
        muHitsPerLayer[L].push_back({ ex[(int)i], ey[(int)i], itZ->second });
      } else if (P==11 || P==-11) {
        eHitsPerLayer[L].push_back({ ex[(int)i], ey[(int)i], itZ->second });
      }
    }

    // 过滤条件（与“刚才那份”一致）：
    // 顶三层(0,1,2)：三层都有命中，且 PID 只能是 ±13
    bool top_ok = true;
    for (int L=0; L<=2; ++L) {
      auto it = pidsPerLayer.find(L);
      if (it == pidsPerLayer.end()) { top_ok=false; break; }
      for (int P : it->second) { if (!(P==13 || P==-13)) { top_ok=false; break; } }
      if (!top_ok) break;
    }
    // 底三层(3,4,5)：三层都有命中，且 PID 只能是 ±13 或 ±11
    bool bot_ok = true;
    for (int L=3; L<=5; ++L) {
      auto it = pidsPerLayer.find(L);
      if (it == pidsPerLayer.end()) { bot_ok=false; break; }
      for (int P : it->second) {
        if (!(P==13 || P==-13 || P==11 || P==-11)) { bot_ok=false; break; }
      }
      if (!bot_ok) break;
    }
    if (!(top_ok && bot_ok)) continue;

    // 聚合点：入 μ（0-2）；出 μ（3-5）；出 e（3-5）
    std::vector<Hit> inMu, outMu, outE;
    for (int L=0; L<=2; ++L) {
      auto it = muHitsPerLayer.find(L);
      if (it != muHitsPerLayer.end())
        inMu.insert(inMu.end(), it->second.begin(), it->second.end());
    }
    for (int L=3; L<=5; ++L) {
      auto itM = muHitsPerLayer.find(L);
      if (itM != muHitsPerLayer.end())
        outMu.insert(outMu.end(), itM->second.begin(), itM->second.end());
      auto itE = eHitsPerLayer.find(L);
      if (itE != eHitsPerLayer.end())
        outE.insert(outE.end(), itE->second.begin(), itE->second.end());
    }

    // 与你之前口径一致：入/出 μ 各 ≥3 点才计入事件
    if (inMu.size() < 3 || outMu.size() < 3) continue;

    // 拟合并统一方向到 +z
    auto vIn   = fitDir(inMu);
    auto vOutM = fitDir(outMu);
    align_to_plus_z(vIn);
    align_to_plus_z(vOutM);

    // theta12 一定能算
    theta12 = angle_between(vIn, vOutM);

    // 出射 e 若足够（≥3点）则计算 theta13/theta23，否则置 NaN
    if (outE.size() >= 3) {
      auto vOutE = fitDir(outE);
      align_to_plus_z(vOutE);
      theta13 = angle_between(vIn,   vOutE);
      theta23 = angle_between(vOutM, vOutE);
    } else {
      theta13 = theta23 = std::numeric_limits<double>::quiet_NaN();
    }

    tout.Fill();
  }

  tout.Write();
  fout.Close();
  std::cout << "Saved to " << outFile
            << " (branches: theta12/theta13/theta23 in radians; selection identical to your previous code)\n";
}
