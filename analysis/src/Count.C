#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

// 统计每层每种粒子的命中数，并写出一个结果树
void Count(const char* inFile = "../../build/root_file/CryMu.root",
                      const char* treeName = "tree",
                      const char* outFile = "../../build/root_file/your_pid_counts.root")
{
  // 1) 打开输入树
  TFile fin(inFile);
  if (fin.IsZombie()) { std::cerr << "Cannot open " << inFile << "\n"; return; }
  TTree* t = (TTree*)fin.Get(treeName);
  if (!t) { std::cerr << "No TTree '" << treeName << "'\n"; return; }

  // 2) 读取 split 的叶子（Edeps.Id / Edeps.Pid）
  TTreeReader r(t);
  TTreeReaderArray<int> lid (r, "Edeps.Id");
  TTreeReaderArray<int> lpid(r, "Edeps.Pid");

  // 3) 统计：layer -> (pid -> count)
  std::map<int, std::map<int, long long>> counts;

  while (r.Next()) {
    const int n = std::min(lid.GetSize(), lpid.GetSize());
    for (int i = 0; i < n; ++i) {
      counts[lid[i]][lpid[i]]++;
    }
  }

  // 4) 打印结果（按层号、再按 |count| 降序）
  auto pdg = TDatabasePDG::Instance();
  std::cout << "\n=== Per-layer PID counts ===\n";
  for (auto &kv : counts) {
    int layer = kv.first;
    // 排序：按 count 降序
    std::vector<std::pair<int,long long>> v(kv.second.begin(), kv.second.end());
    std::sort(v.begin(), v.end(), [](auto& a, auto& b){ return a.second > b.second; });

    std::cout << "Layer " << layer << ":\n";
    for (auto &pp : v) {
      int pid = pp.first;
      long long c = pp.second;
      const auto* p = pdg->GetParticle(pid);
      const char* name = p ? p->GetName() : "?";
      std::cout << "  PID " << pid << " (" << name << "): " << c << "\n";
    }
  }

  // 5) 写到新 root 文件（便于后续画图）
  TFile fout(outFile, "RECREATE");
  int    oLayer = 0;
  int    oPID   = 0;
  long long oCount = 0;
  TString oName;

  TTree tout("layer_pid_counts", "per-layer PID counts");
  tout.Branch("Layer", &oLayer);
  tout.Branch("PID",   &oPID);
  tout.Branch("Count", &oCount);
  tout.Branch("Name",  &oName);

  for (auto &kv : counts) {
    oLayer = kv.first;
    for (auto &pp : kv.second) {
      oPID   = pp.first;
      oCount = pp.second;
      const auto* p = pdg->GetParticle(oPID);
      oName = p ? p->GetName() : "?";
      tout.Fill();
    }
  }

  tout.Write();
  fout.Close();
  std::cout << "\nSaved result tree 'layer_pid_counts' to " << outFile << "\n";
}

