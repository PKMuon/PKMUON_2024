#include <Math/Vector3D.h>
#include <TFile.h>
#include <TTree.h>
#include <math.h>

#include <fstream>
#include <tuple>
#include <vector>

#include "../../include/Object.hh"

using namespace std;

Long64_t PoCA_sim(const char *infile = "../../build/root_file/CryMuAna.root",
    const char *outfile = "../../build/root_file/CryMuPoca.root")
{
  // Input file and tree.
  TFile *file_in = new TFile(infile);
  TTree *tree_in = (TTree *)file_in->Get("tree");
  vector<Double_t> *XEdep = NULL, *YEdep = NULL, *ZEdep = NULL, *XSmeared = NULL, *YSmeared = NULL;
  tree_in->SetBranchAddress("XEdep", &XEdep);
  tree_in->SetBranchAddress("YEdep", &YEdep);
  tree_in->SetBranchAddress("ZEdep", &ZEdep);
  tree_in->SetBranchAddress("XSmeared", &XSmeared);
  tree_in->SetBranchAddress("YSmeared", &YSmeared);

  // Output file and tree.
  TFile *file_out = new TFile(outfile, "RECREATE");
  TTree *tree_out = new TTree("tree", "tree");
  tree_out = tree_in->CloneTree(0);  // copy 0 entries
  Double_t DPoCAEdep, XPoCAEdep, YPoCAEdep, ZPoCAEdep;
  Double_t DPoCASmeared, XPoCASmeared, YPoCASmeared, ZPoCASmeared;
  tree_out->Branch("DPoCAEdep", &DPoCAEdep);
  tree_out->Branch("XPoCAEdep", &XPoCAEdep);
  tree_out->Branch("YPoCAEdep", &YPoCAEdep);
  tree_out->Branch("ZPoCAEdep", &ZPoCAEdep);
  tree_out->Branch("DPoCASmeared", &DPoCASmeared);
  tree_out->Branch("XPoCASmeared", &XPoCASmeared);
  tree_out->Branch("YPoCASmeared", &YPoCASmeared);
  tree_out->Branch("ZPoCASmeared", &ZPoCASmeared);

  Long64_t nentry = tree_in->GetEntries();
  for(Long64_t ientry = 0; ientry < nentry; ientry++) {
    if(ientry % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2) << (ientry / (double)nentry) * 100 << "%" << endl;
    }
    tree_in->GetEntry(ientry);

    auto args = {
      tie(*XEdep, *YEdep, *ZEdep, DPoCAEdep, XPoCAEdep, YPoCAEdep, ZPoCAEdep),
      tie(*XSmeared, *YSmeared, *ZEdep, DPoCASmeared, XPoCASmeared, YPoCASmeared, ZPoCASmeared),
    };
    for(auto &[x, y, z, d_poca, x_poca, y_poca, z_poca] : args) {
      ROOT::Math::XYZVector a, b, va, vb, uva, uvb, vab;
      ROOT::Math::XYZVector vn, uvn;
      ROOT::Math::XYZVector vm;

      size_t n = x.size();
      assert(n >= 4 && n % 2 == 0);
      size_t i1 = 0, i2 = n / 2 - 1, i3 = n / 2, i4 = n - 1;

      a.SetCoordinates(x[i2], y[i2], z[i2]);  // 点 a 定义为打在第二层探测器上的位置
      b.SetCoordinates(x[i3], y[i3], z[i3]);  // 点 b 定义为打在第三层探测器上的位置
      va.SetCoordinates(x[i2] - x[i1], y[i2] - y[i1], z[i2] - z[i1]);  // va 是未归一化的过点 a 的径迹方向矢量
      vb.SetCoordinates(x[i4] - x[i3], y[i4] - y[i3], z[i4] - z[i3]);  // vb 是未归一化的过点 b 的径迹方向矢量

      uva = va.Unit(), uvb = vb.Unit();
      vn = va.Cross(vb), uvn = vn.Unit();  // 计算公垂线方向向量
      vab = b - a;
      d_poca = vab.Dot(uvn);  // 两条直线的距离
      a += d_poca * uvn;      // 将 a 平移以 uvn 为法向量，包含 b 的平面上

      /*
       * Move a -> A.
       *
       * a ----A--> [u]va
       *       |
       *       b
       *        \
       *         \
       *          * [u]vb
       */
      vab = b - a;              // 在该平面内重新计算 vab
      a += uva * uva.Dot(vab);  // 将 a 平移到 [b 到过 a 点的径迹上的垂足]

      /*
       * Compute m.
       *
       *   --m-a--> [u]va
       *      \|
       *       b
       *        \
       *         \
       *          * [u]vb
       */
      vab = b - a;                                   // 计算 b 到过 a 点的径迹之间的距离 ab
      vm = b - uvb * (vab.Dot(vab) / uvb.Dot(vab));  // 轨迹交点
      vm -= 0.5 * d_poca * uvn;                      // PoCA 点
      x_poca = vm.x(), y_poca = vm.y(), z_poca = vm.z();
    }

    tree_out->Fill();
  }

  tree_out->Write(NULL, TObject::kOverwrite);
  file_out->Close();
  file_in->Close();
  return nentry;
}
