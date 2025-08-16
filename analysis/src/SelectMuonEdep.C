#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <iostream>

// 只保留 muon(±13) 的 Edep 命中坐标，写到新树 "mu_edep"
void SelectMuonEdep(const char* inFile  = "../../build/root_file/CryMu.root",
                    const char* treeName= "tree",
                    const char* outFile = "../../build/root_file/mu_edep.root")
{
  // 1) 打开输入
  TFile fin(inFile);
  if (fin.IsZombie()) { std::cerr << "Cannot open " << inFile << "\n"; return; }
  TTree* tin = (TTree*)fin.Get(treeName);
  if (!tin) { std::cerr << "No TTree '" << treeName << "' in file\n"; return; }

  // 2) 读取 split 叶子
  TTreeReader r(tin);
  TTreeReaderArray<int>    id   (r, "Edeps.Id");
  TTreeReaderArray<int>    pid  (r, "Edeps.Pid");
  TTreeReaderArray<double> x    (r, "Edeps.X");
  TTreeReaderArray<double> y    (r, "Edeps.Y");

  // 可选：能量
  TTreeReaderArray<double> val  (r, "Edeps.Value");

  // 3) 输出文件与树（扁平化逐命中）
  TFile fout(outFile, "RECREATE");
  Long64_t oEvent = 0;   // 事件号（方便回溯）
  int      oLayer = 0;   // 层/单元 Id
  int      oPID   = 0;   // PID（只会是 ±13）
  double   oX=0, oY=0;   // 坐标
  double   oEdep=0;      // 能量（可选）

  TTree tout("mu_edep","muon-only Edep hits (flattened)");
  tout.Branch("Event", &oEvent);
  tout.Branch("Layer", &oLayer);
  tout.Branch("PID",   &oPID);
  tout.Branch("X",     &oX);
  tout.Branch("Y",     &oY);
  tout.Branch("Edep",  &oEdep);  // 如不需要可注释掉

  // 4) 逐事件、逐命中过滤 μ 子
  Long64_t evt = 0;
  while (r.Next()) {
    const int n = std::min({id.GetSize(), pid.GetSize(), x.GetSize(), y.GetSize(), val.GetSize()});
    for (int i = 0; i < n; ++i) {
      const int p = pid[i];
      if (p == 13 || p == -13) {       // μ⁻ or μ⁺
        oEvent = evt;
        oLayer = id[i];
        oPID   = p;
        oX     = x[i];
        oY     = y[i];
        oEdep  = val[i];
        tout.Fill();
      }
    }
    ++evt;
  }

  tout.Write();
  fout.Close();
  std::cout << "Saved muon-only hits to " << outFile
            << " (tree: mu_edep)\n";
}

// #include <TFile.h>
// #include <TTree.h>
// #include <TTreeReader.h>
// #include <TTreeReaderArray.h>
// #include <iostream>
// #include <map>
// #include <vector>
// #include <cmath>

// struct Hit {
//     double x, y, z;
// };

// static const std::map<int,double> layerZ = {
//     {0, -607.5},
//     {1, -426.5},
//     {2, -244.5},
//     {3,  237.5},
//     {4,  418.5},
//     {5,  600.5}
// };

// static std::array<double,3> fitDir(const std::vector<Hit>& pts) {
//     double cx=0, cy=0, cz=0;
//     for (auto& p : pts) {
//         cx += p.x; cy += p.y; cz += p.z;
//     }
//     cx /= pts.size(); cy /= pts.size(); cz /= pts.size();

//     auto S = [&](const std::array<double,3>& v) {
//         std::array<double,3> r{0,0,0};
//         for (auto& p : pts) {
//             double dx = p.x - cx;
//             double dy = p.y - cy;
//             double dz = p.z - cz;
//             double dotp = dx*v[0] + dy*v[1] + dz*v[2];
//             r[0] += dx*dotp;
//             r[1] += dy*dotp;
//             r[2] += dz*dotp;
//         }
//         return r;
//     };

//     std::array<double,3> v{1,0,0};
//     for (int it=0; it<20; ++it) {
//         v = S(v);
//         double m = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
//         if (m==0) break;
//         v[0] /= m; v[1] /= m; v[2] /= m;
//     }
//     return v;
// }

// void SelectMuonEdep(const char* inFile="../../build/root_file/CryMu.root",
//                        const char* treeName="tree",
//                        const char* outFile="../../build/root_file/mu_costheta.root")
// {
//     TFile fin(inFile);
//     if (fin.IsZombie()) { std::cerr << "Cannot open " << inFile << "\n"; return; }
//     TTree* tin = (TTree*)fin.Get(treeName);
//     if (!tin) { std::cerr << "No tree: " << treeName << "\n"; return; }

//     TTreeReader r(tin);
//     TTreeReaderArray<int>    pid(r, "Edeps.Pid");
//     TTreeReaderArray<int>    lid(r, "Edeps.Id");
//     TTreeReaderArray<double> ex (r, "Edeps.X");
//     TTreeReaderArray<double> ey (r, "Edeps.Y");

//     TFile fout(outFile,"RECREATE");
//     TTree tout("tree","with costheta and theta_rad");
//     double costheta=0;
//     double theta_rad=0; // 新增：弧度制角度
//     tout.Branch("costheta",&costheta);
//     tout.Branch("theta_rad",&theta_rad);

//     while (r.Next()) {
//         std::map<int,std::vector<Hit>> hitsPerLayer;

//         const size_t n = std::min({(size_t)pid.GetSize(), (size_t)lid.GetSize(), (size_t)ex.GetSize(), (size_t)ey.GetSize()});
//         for (size_t i=0; i<n; ++i) {
//             if (pid[i] == 13 || pid[i] == -13) {
//                 int layer = lid[i];
//                 auto it = layerZ.find(layer);
//                 if (it != layerZ.end()) {
//                     hitsPerLayer[layer].push_back({ex[i], ey[i], it->second});
//                 }
//             }
//         }

//         bool fullMuon = true;
//         for (int layer=0; layer<=5; ++layer) {
//             if (hitsPerLayer.find(layer) == hitsPerLayer.end()) {
//                 fullMuon = false;
//                 break;
//             }
//         }
//         if (!fullMuon) continue;

//         std::vector<Hit> inHits, outHits;
//         for (int layer=0; layer<=2; ++layer) {
//             auto& v = hitsPerLayer[layer];
//             inHits.insert(inHits.end(), v.begin(), v.end());
//         }
//         for (int layer=3; layer<=5; ++layer) {
//             auto& v = hitsPerLayer[layer];
//             outHits.insert(outHits.end(), v.begin(), v.end());
//         }

//         auto vin  = fitDir(inHits);
//         auto vout = fitDir(outHits);

//         costheta = vin[0]*vout[0] + vin[1]*vout[1] + vin[2]*vout[2];
//         if (costheta > 1) costheta = 1;     // 防止数值误差
//         if (costheta < -1) costheta = -1;

//         theta_rad = std::acos(costheta);    // 弧度

//         tout.Fill();
//     }

//     tout.Write();
//     fout.Close();
//     std::cout << "Saved " << outFile << "\n";
// }

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

// void SelectMuonEdep(const char* inFile="../../build/root_file/CryMu.root",
//                                 const char* treeName="tree",
//                                 const char* outFile="../../build/root_file/mu_costheta_filtered.root")
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
//   TTree tout("tree","costheta with event-level mu/e filter");
//   Long64_t EventID = 0;
//   double costheta = NAN;
//   double theta_rad = NAN;
//   tout.Branch("EventID",   &EventID);
//   tout.Branch("costheta",  &costheta);
//   tout.Branch("theta_rad", &theta_rad);

//   // 事件循环
//   for (EventID=0; r.Next(); ++EventID) {

//     // 每层收集: 该层出现过的 PID 集合 & 该层的 μ 命中点
//     std::map<int, std::set<int>> pidsPerLayer;
//     std::map<int, std::vector<Hit>> muHitsPerLayer;

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
//       if (P==13 || P==-13) { // 只把 μ 命中点用于方向拟合
//         muHitsPerLayer[L].push_back({ ex[(int)i], ey[(int)i], itZ->second });
//       }
//     }

//     // 过滤条件：
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
//     if (!(top_ok && bot_ok)) {  // 不满足条件则跳过事件
//       continue;
//     }

//     // 拟合方向：入射=0/1/2 层的 μ 命中；出射=3/4/5 层的 μ 命中
//     std::vector<Hit> inHits, outHits;
//     for (int L=0; L<=2; ++L) {
//       auto it = muHitsPerLayer.find(L);
//       if (it != muHitsPerLayer.end()) {
//         inHits.insert(inHits.end(), it->second.begin(), it->second.end());
//       }
//     }
//     for (int L=3; L<=5; ++L) {
//       auto it = muHitsPerLayer.find(L);
//       if (it != muHitsPerLayer.end()) {
//         outHits.insert(outHits.end(), it->second.begin(), it->second.end());
//       }
//     }

//     // 要求方向拟合至少各有 3 个点（更稳）
//     if (inHits.size() < 3 || outHits.size() < 3) { continue; }

//     auto vin  = fitDir(inHits);
//     auto vout = fitDir(outHits);

//     // 计算角度
//     costheta = vin[0]*vout[0] + vin[1]*vout[1] + vin[2]*vout[2];
//     if (costheta >  1) costheta = 1;
//     if (costheta < -1) costheta = -1;
//     theta_rad = std::acos(costheta);

//     tout.Fill();
//   }

//   tout.Write();
//   fout.Close();
//   std::cout << "Saved filtered events with costheta/theta_rad to "
//             << outFile << std::endl;
// }
