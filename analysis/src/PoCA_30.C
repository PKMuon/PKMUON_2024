#include <Math/Vector3D.h>
#include <TFile.h>
#include <TTree.h>
#include <math.h>
#include <fstream>
#include <tuple>
#include <vector>
#include <iomanip>  // 用于进度条格式化
#include <cassert>

#include "../../include/Object.hh"

using namespace std;

// 函数名已修改为PoCA_30
Long64_t PoCA_30(const char *infile = "../../build/root_file/e_trace_output.root",
    const char *outfile = "../../build/root_file/CryMuPoca.root")
{
  // 输入文件和Tree
  TFile *file_in = new TFile(infile);
  if (!file_in->IsOpen()) {
    cerr << "Error: 无法打开输入文件 " << infile << endl;
    return -1;
  }
  TTree *tree_in = (TTree *)file_in->Get("tree");
  if (!tree_in) {
    cerr << "Error: 无法找到输入Tree" << endl;
    return -1;
  }

  // 输入分支地址
  vector<Double_t> *XEdep = NULL, *YEdep = NULL, *ZEdep = NULL, *XSmeared = NULL, *YSmeared = NULL;
  tree_in->SetBranchAddress("XEdep", &XEdep);
  tree_in->SetBranchAddress("YEdep", &YEdep);
  tree_in->SetBranchAddress("ZEdep", &ZEdep);
  tree_in->SetBranchAddress("XSmeared", &XSmeared);
  tree_in->SetBranchAddress("YSmeared", &YSmeared);

  // 输出文件1: CryMuPoca.root (原始需求)
  TFile *file_out = new TFile(outfile, "RECREATE");
  TTree *tree_out = tree_in->CloneTree(0);  // 克隆结构
  tree_out->SetDirectory(file_out);  // 显式关联到输出文件
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

  // 输出文件2: e_Poca_output.root (筛选后的原始数据)
  TFile *file_poca_out = new TFile("../../build/root_file/e_Poca_output.root", "RECREATE");
  TTree *tree_poca_out = tree_in->CloneTree(0);  // 克隆原始结构
  tree_poca_out->SetDirectory(file_poca_out);  // 显式关联到输出文件

  // 事件循环
  Long64_t nentry = tree_in->GetEntries();
  for (Long64_t ientry = 0; ientry < nentry; ientry++) {
    if (ientry % 1000 == 0) {
      cout << "Processing progress: " << fixed << setprecision(2) 
           << (ientry / (double)nentry) * 100 << "%" << endl;
    }
    tree_in->GetEntry(ientry);

    // 保存Edep组PoCA坐标用于筛选
    double edep_x_poca = 0.0, edep_y_poca = 0.0, edep_z_poca = 0.0;

    auto args = {
      tie(*XEdep, *YEdep, *ZEdep, DPoCAEdep, XPoCAEdep, YPoCAEdep, ZPoCAEdep),
      tie(*XSmeared, *YSmeared, *ZEdep, DPoCASmeared, XPoCASmeared, YPoCASmeared, ZPoCASmeared),
    };
    int count = 0;  // 区分两组数据

    for (auto &[x, y, z, d_poca, x_poca, y_poca, z_poca] : args) {
      ROOT::Math::XYZVector a, b, va, vb, uva, uvb, vab;
      ROOT::Math::XYZVector vn, uvn;
      ROOT::Math::XYZVector vm;

      size_t n = x.size() - 25;
      assert(n >= 4 && n % 2 == 0);  // 确保数据量合法
      size_t i1 = 1, i2 = 2, i3 = 4, i4 = 3;

      a.SetCoordinates(x[i2], y[i2], z[i2]);
      b.SetCoordinates(x[i3], y[i3], z[i3]);
      va.SetCoordinates(x[i2]-x[i1], y[i2]-y[i1], z[i2]-z[i1]);
      vb.SetCoordinates(x[i4]-x[i3], y[i4]-y[i3], z[i4]-z[i3]);

      uva = va.Unit(); uvb = vb.Unit();
      vn = va.Cross(vb); uvn = vn.Unit();
      vab = b - a;
      d_poca = vab.Dot(uvn);
      a += d_poca * uvn;

      vab = b - a;
      a += uva * uva.Dot(vab);

      vab = b - a;
      vm = b - uvb * (vab.Dot(vab)/uvb.Dot(vab));
      vm -= 0.5 * d_poca * uvn;
      x_poca = vm.x(); y_poca = vm.y(); z_poca = vm.z();

      // 保存第一组(Edep)的PoCA坐标
      if (count == 0) {
        edep_x_poca = x_poca;
        edep_y_poca = y_poca;
        edep_z_poca = z_poca;
      }
      count++;
    }

    // 填充原始输出树
    tree_out->Fill();

    // 筛选并填充PoCA输出树
    if (edep_x_poca > -60 && edep_x_poca < 60 &&
        edep_y_poca > -60 && edep_y_poca < 60 &&
        edep_z_poca > -15 && edep_z_poca < 15) {
      tree_poca_out->Fill();
    }
  }

  // 写入并关闭文件（严格按文件关联关系操作）
  file_out->cd();  // 切换到输出文件目录
  tree_out->Write("", TObject::kOverwrite);
  file_out->Close();

  file_poca_out->cd();  // 切换到筛选文件目录
  tree_poca_out->Write("", TObject::kOverwrite);
  file_poca_out->Close();

  file_in->Close();

  // 清理动态内存
  delete file_in;
  delete file_out;
  delete file_poca_out;

  return nentry;
}

// 添加ROOT命令行直接执行入口
#ifndef __CINT__
int main() {
  return PoCA_30();
}
#endif