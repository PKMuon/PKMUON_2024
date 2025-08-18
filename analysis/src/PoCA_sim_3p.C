#include <Math/Vector3D.h>
#include <TFile.h>
#include <TTree.h>
#include <math.h>

#include <fstream>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include "../../include/Object.hh"

using namespace std;
// 声明三点拟合直线方向向量的函数原型
std::vector<Double_t> FitLine3Points(const std::vector<Double_t>& x, 
                                     const std::vector<Double_t>& y, 
                                     const std::vector<Double_t>& z);
Long64_t PoCA_sim_3p(const char *infile = "../../build/root_file/CryMuAna.root",
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
      ROOT::Math::XYZVector a, b, uva, uvb, vab;  //va,vb
      ROOT::Math::XYZVector vn, uvn;
      ROOT::Math::XYZVector vm;

      size_t n = x.size();
      assert(n >= 4 && n % 2 == 0);
      size_t i1 = n / 2 - 1, i2 = n / 2 -2, i3 = n / 2 , i4 = n / 2 + 1;//size_t i1 = 0, i2 = n / 2 - 1, i3 = n / 2, i4 = n - 1;

      a.SetCoordinates(x[i2], y[i2], z[i2]);//a.SetCoordinates(x[i2], y[i2], z[i2]);  // 点 a 定义为打在第二层探测器上的位置
      b.SetCoordinates(x[i3], y[i3], z[i3]);//b.SetCoordinates(x[i3], y[i3], z[i3]);  // 点 b 定义为打在第三层探测器上的位置
//      va.SetCoordinates(x[i2] - x[i1], y[i2] - y[i1], z[i2] - z[i1]);//va.SetCoordinates(x[i2] - x[i1], y[i2] - y[i1], z[i2] - z[i1]);  // va 是未归一化的过点 a 的径迹方向矢量
//      vb.SetCoordinates(x[i4] - x[i3], y[i4] - y[i3], z[i4] - z[i3]);//vb.SetCoordinates(x[i4] - x[i3], y[i4] - y[i3], z[i4] - z[i3]);  // vb 是未归一化的过点 b 的径迹方向矢量
// 前半部分：取3个点（索引0,1,2）拟合方向向量va
vector<Double_t> x_va = {x[0], x[1], x[2]};
vector<Double_t> y_va = {y[0], y[1], y[2]};
vector<Double_t> z_va = {z[0], z[1], z[2]};
//vector<Double_t> va = FitLine3Points(x_va, y_va, z_va);
std::vector<Double_t> va_vec = FitLine3Points(x_va, y_va, z_va);
TVector3 va(va_vec[0], va_vec[1], va_vec[2]);  // 转换为 TVector3

// 后半部分：取3个点（索引3,4,5）拟合方向向量vb
vector<Double_t> x_vb = {x[3], x[4], x[5]};
vector<Double_t> y_vb = {y[3], y[4], y[5]};
vector<Double_t> z_vb = {z[3], z[4], z[5]};
//vector<Double_t> vb = FitLine3Points(x_vb, y_vb, z_vb);
std::vector<Double_t> vb_vec = FitLine3Points(x_vb, y_vb, z_vb);
TVector3 vb(vb_vec[0], vb_vec[1], vb_vec[2]);  // 转换为 TVector3


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

// 输入三个点的坐标，返回拟合直线的方向向量
vector<Double_t> FitLine3Points(const vector<Double_t>& x, const vector<Double_t>& y, const vector<Double_t>& z) {
    assert(x.size() == 3 && y.size() == 3 && z.size() == 3); // 确保输入3个点

    // Step 1: 计算质心
    Double_t cx = (x[0] + x[1] + x[2]) / 3;
    Double_t cy = (y[0] + y[1] + y[2]) / 3;
    Double_t cz = (z[0] + z[1] + z[2]) / 3;

    // Step 2: 构建协方差矩阵 (3x3)
    Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; ++i) {
        Double_t dx = x[i] - cx;
        Double_t dy = y[i] - cy;
        Double_t dz = z[i] - cz;
        C(0,0) += dx*dx;  // Cxx
        C(0,1) += dx*dy;  // Cxy
        C(0,2) += dx*dz;  // Cxz
        C(1,1) += dy*dy;  // Cyy
        C(1,2) += dy*dz;  // Cyz
        C(2,2) += dz*dz;  // Czz
    }
    C(1,0) = C(0,1); // 对称矩阵
    C(2,0) = C(0,2);
    C(2,1) = C(1,2);
    C /= 2; // 除以 (n-1)，n=3

    // Step 3: 特征值分解，取最大特征值对应的特征向量
    Eigen::EigenSolver<Eigen::Matrix3d> solver(C);
    Eigen::Vector3d eigenvalues = solver.eigenvalues().real();
    Eigen::Matrix3d eigenvectors = solver.eigenvectors().real();

    // 找到最大特征值的索引
    int max_idx = 0;
    if (eigenvalues[1] > eigenvalues[max_idx]) max_idx = 1;
    if (eigenvalues[2] > eigenvalues[max_idx]) max_idx = 2;

    // 返回方向向量（未归一化，后续计算余弦值时自动归一化）
    return {eigenvectors(0, max_idx), eigenvectors(1, max_idx), eigenvectors(2, max_idx)};
}
