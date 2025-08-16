// BuildLayerZ_FromPeaks.C
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

void LayerZ(const char* file="../../build/root_file/CryMuAna.root",
                           const char* tree="tree",
                           int nbins=2000, double zmin=-1000, double zmax=1000,
                           double sigma=2.0, double thresh=0.05)
{
  TFile fin(file);
  if(fin.IsZombie()){ std::cerr<<"open "<<file<<" failed\n"; return; }
  TTree* t = (TTree*)fin.Get(tree);
  if(!t){ std::cerr<<"no tree "<<tree<<"\n"; return; }

  TH1D h("hZ","ZEdep;Z;counts", nbins, zmin, zmax);
  t->Draw("ZEdep >> hZ", "", "goff");   // 把 ZEdep 投到直方图

  TSpectrum sp(50);                     // 最多找 50 个峰，按需调
  int npeaks = sp.Search(&h, sigma, "", thresh);
  std::cout << "found peaks: " << npeaks << "\n";
  double* x = sp.GetPositionX();

  std::vector<double> peaks(x, x+npeaks);
  std::sort(peaks.begin(), peaks.end());

  // 有时 TSpectrum 会给出左右相近的重复峰，做一次简单去重
  std::vector<double> uniq;
  const double merge_window = (zmax-zmin)/nbins * 3; // 3个bin 以内算重复
  for(double z : peaks){
    if(uniq.empty() || std::fabs(z-uniq.back()) > merge_window)
      uniq.push_back(z);
  }

  std::cout << "unique peak positions (Z):\n";
  std::cout << std::fixed; 
  for(double z : uniq) std::cout << "  " << z << "\n";

  // 如果你需要生成一个顺序层号 → Z 的映射，可直接输出：
  std::cout << "std::unordered_map<int,double> layerZ = {\n";
  for(size_t i=0;i<uniq.size();++i)
    std::cout << "  {" << (int)i << ", " << uniq[i] << "},\n";
  std::cout << "};\n";

  // 想看一下图：
  // TCanvas c; h.Draw(); c.Update();
}

