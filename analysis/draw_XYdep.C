#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include <iostream>

void draw_XYdep() {
    // 1. 打开ROOT文件
    TFile *file = TFile::Open("../build/root_file/poca_CryMu.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: 无法打开文件 ../build/root_file/poca_CryMu.root" << std::endl;
        return;
    }

    // 2. 获取TTree对象
    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: 无法找到树对象 'tree'" << std::endl;
        file->Close();
        return;
    }

    // 3. 创建画布并设置样式
    TCanvas *c1 = new TCanvas("c1", "XPoCAEdep vs YPoCAEdep 二维分布", 800, 600);
    gStyle->SetOptStat(0);  // 关闭统计信息框

    // 4. 绘制二维密度图（颜色编码）
    //tree->Draw("YPoCAEdep:XPoCAEdep", "", "COLZ");  // COLZ: 颜色+Z轴颜色条
    tree->Draw("YEdep[2]:XEdep[2]", "", "COLZ");

    // 5. 获取自动生成的直方图并美化
    TH2F *h2 = (TH2F*)gPad->GetPrimitive("htemp");
    if (h2) {
        h2->SetTitle("XPoCAEdep vs YPoCAEdep 二维密度分布");
        h2->GetXaxis()->SetTitle("XPoCAEdep");
        h2->GetYaxis()->SetTitle("YPoCAEdep");
        h2->GetZaxis()->SetTitle("事件数");  // 设置颜色条标题
    }

    // 6. 保存图像（支持PNG/PDF/ROOT等格式）
    c1->SaveAs("png/XvsY_PoCAEdep_2D.png");
    std::cout << "图像已保存为 XvsY_PoCAEdep_2D.png" << std::endl;

    // 7. 清理资源
    file->Close();
    delete file;
    delete c1;
}
