#include "TFile.h"
#include "TH1.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using std::operator""s;

void CryPKUGps(const char *histname, const char *parname)
{
  TFile *file = TFile::Open("test.root", "READ");
  if(!file || file->IsZombie()) {
    std::cerr << "无法打开 ROOT 文件!" << std::endl;
    return;
  }

  TH1 *ke = dynamic_cast<TH1 *>(file->Get(("ke_"s + histname).c_str()));
  TH1 *ct = dynamic_cast<TH1 *>(file->Get(("ct_"s + histname).c_str()));
  if(!ke || !ct) {
    std::cerr << histname << ": 无法获取直方图!" << std::endl;
    file->Close();
    return;
  }

  std::string opath = "CryGps_"s + parname + ".mac";
  std::ofstream output_file(opath);
  if(!output_file) {
    std::cerr << histname << ": 无法创建输出文件!" << std::endl;
    file->Close();
    return;
  }
  output_file << std::scientific << std::setprecision(18);

  output_file << "# ----------------------------------------------\n";
  output_file << "# General Particle Source (GPS) Settings\n";
  output_file << "# ----------------------------------------------\n";

  output_file << "\n";
  output_file << "/gps/particle " << parname << "\n";
  output_file << "/gps/ang/type user\n";
  output_file << "/gps/hist/type theta\n";
  for(int i = 1; i <= ct->GetNbinsX(); ++i) {
    float bin_center_ct = ct->GetBinCenter(i);
    float bin_content_ct = ct->GetBinContent(i);
    float angle_in_radians = M_PI - std::acos(bin_center_ct);
    output_file << "/gps/hist/point " << angle_in_radians << " " << bin_content_ct << "\n";
  }

  output_file << "\n";
  output_file << "/gps/ene/type Arb\n";
  output_file << "/gps/hist/type arb\n";
  for(int i = 1; i <= ke->GetNbinsX(); ++i) {
    float bin_center_ke = ke->GetBinCenter(i);
    float bin_content_ke = ke->GetBinContent(i);
    output_file << "/gps/hist/point " << pow(10, bin_center_ke) << " " << bin_content_ke << "\n";
  }
  output_file << "/gps/hist/inter Lin\n";

  output_file.close();
  file->Close();
  std::cout << "数据已保存至 " << opath << std::endl;
}

void CryPKUGps()
{
  // https://agenda.infn.it/event/19376/contributions/97267/attachments/64642/78312/08_-_Pandola_-_Physics_-_1.pdf
  CryPKUGps("-13",    "mu+");
  CryPKUGps("+13",    "mu-");
  CryPKUGps("-11",    "e+");
  CryPKUGps("+11",    "e-");
  CryPKUGps("+2112",  "neutron");
  CryPKUGps("-2112",  "anti_neutron");
  CryPKUGps("+211",   "pi+");
  CryPKUGps("-211",   "pi-");
  CryPKUGps("+22",    "gamma");
  CryPKUGps("+2212",  "proton");
  CryPKUGps("-2212",  "anti_proton");
  exit(0);
}
