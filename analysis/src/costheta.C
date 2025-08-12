#include <Math/Vector3D.h>
#include <TFile.h>
#include <TTree.h>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <tuple>
#include <vector>

#include "../../include/Object.hh"

using namespace std;

Long64_t costheta(const char *infile = "../../build/root_file/CryMuAna.root",
                  const char *outfile = "../../build/root_file/CryMuPoca.root")
{
  // Input file and tree.
  TFile *file_in = new TFile(infile);
  TTree *tree_in = (TTree *)file_in->Get("tree");
  vector<Double_t> *XEdep = nullptr, *YEdep = nullptr, *ZEdep = nullptr, *XSmeared = nullptr, *YSmeared = nullptr;
  tree_in->SetBranchAddress("XEdep", &XEdep);
  tree_in->SetBranchAddress("YEdep", &YEdep);
  tree_in->SetBranchAddress("ZEdep", &ZEdep);
  tree_in->SetBranchAddress("XSmeared", &XSmeared);
  tree_in->SetBranchAddress("YSmeared", &YSmeared);

  // Output file and tree.
  TFile *file_out = new TFile(outfile, "RECREATE");
  TTree *tree_out = tree_in->CloneTree(0);  // 复制结构但不复制条目

  // 新增：输出变量
  Double_t DPoCAEdep=0, XPoCAEdep=0, YPoCAEdep=0, ZPoCAEdep=0, CosThetaEdep=0;
  Double_t DPoCASmeared=0, XPoCASmeared=0, YPoCASmeared=0, ZPoCASmeared=0, CosThetaSmeared=0;

  tree_out->Branch("DPoCAEdep",      &DPoCAEdep);
  tree_out->Branch("XPoCAEdep",      &XPoCAEdep);
  tree_out->Branch("YPoCAEdep",      &YPoCAEdep);
  tree_out->Branch("ZPoCAEdep",      &ZPoCAEdep);
  tree_out->Branch("CosThetaEdep",   &CosThetaEdep);      // 新增

  tree_out->Branch("DPoCASmeared",   &DPoCASmeared);
  tree_out->Branch("XPoCASmeared",   &XPoCASmeared);
  tree_out->Branch("YPoCASmeared",   &YPoCASmeared);
  tree_out->Branch("ZPoCASmeared",   &ZPoCASmeared);
  tree_out->Branch("CosThetaSmeared",&CosThetaSmeared);   // 新增

  Long64_t nentry = tree_in->GetEntries();
  for(Long64_t ientry = 0; ientry < nentry; ++ientry) {
    if(ientry % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2)
           << (ientry / (double)nentry) * 100 << "%" << endl;
    }
    tree_in->GetEntry(ientry);

    // 打两套：Edep 与 Smeared
    auto args = {
      tie(*XEdep,    *YEdep,    *ZEdep,    DPoCAEdep,    XPoCAEdep,    YPoCAEdep,    ZPoCAEdep,    CosThetaEdep),
      tie(*XSmeared, *YSmeared, *ZEdep,    DPoCASmeared, XPoCASmeared, YPoCASmeared, ZPoCASmeared, CosThetaSmeared)
    };

    for (auto &[x, y, z, d_poca, x_poca, y_poca, z_poca, cos_theta] : args) {
      ROOT::Math::XYZVector a, b, va, vb, uva, uvb, vab;
      ROOT::Math::XYZVector vn, uvn, vm;

      const size_t n = x.size();
      assert(n >= 4 && n % 2 == 0);

      // 这里你选择了中间四层来构造两段直线
      size_t i1 = n/2 - 1, i2 = n/2 - 2, i3 = n/2, i4 = n/2 + 1;

      a.SetCoordinates(x[i2], y[i2], z[i2]);
      b.SetCoordinates(x[i3], y[i3], z[i3]);
      va.SetCoordinates(x[i2] - x[i1], y[i2] - y[i1], z[i2] - z[i1]);
      vb.SetCoordinates(x[i4] - x[i3], y[i4] - y[i3], z[i4] - z[i3]);

      // 计算 cosθ（注意零长度防护）
      const double magA = va.R();
      const double magB = vb.R();
      if (magA > 0 && magB > 0) {
        cos_theta = (va.Dot(vb)) / (magA * magB);
        // 数值稳定性：夹角余弦理论上在 [-1,1]
        if (cos_theta > 1.0)  cos_theta = 1.0;
        if (cos_theta < -1.0) cos_theta = -1.0;
      } else {
        cos_theta = std::numeric_limits<double>::quiet_NaN();
      }

      // 以下是你原本的 PoCA 计算
      uva = va.Unit(); uvb = vb.Unit();
      vn  = va.Cross(vb);
      const double vnMag = vn.R();
      if (vnMag == 0) {
        // 两线平行或几乎平行：PoCA 距离按 0 处理，并取中点投影
        d_poca = 0.0;
        ROOT::Math::XYZVector mid = 0.5*(a+b);
        x_poca = mid.x(); y_poca = mid.y(); z_poca = mid.z();
      } else {
        uvn = vn.Unit();         // 公垂线方向
        vab = b - a;
        d_poca = vab.Dot(uvn);   // 线间最短距离（有符号）
        a += d_poca * uvn;       // 把 a 平移到包含 b 的那个平面

        // 在该平面内，把 a 投到“过 a 的轨迹”上（求垂足）
        vab = b - a;
        a  += uva * uva.Dot(vab);

        // 从 b 沿 vb 的反方向，求到达这条线的最近点，再回退半个 d_poca 即 PoCA
        vab = b - a;
        vm  = b - uvb * (vab.Dot(vab) / uvb.Dot(vab));
        vm -= 0.5 * d_poca * uvn;

        x_poca = vm.x(); y_poca = vm.y(); z_poca = vm.z();
      }
    }

    tree_out->Fill();
  }

  tree_out->Write(nullptr, TObject::kOverwrite);
  file_out->Close();
  file_in->Close();
  return nentry;
}
