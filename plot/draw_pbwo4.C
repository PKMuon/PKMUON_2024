#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"

void draw_pbwo4() {
    // 设置画布和样式
    TCanvas *c1 = new TCanvas("c1", "5x5 Id热图", 800, 800);
    gStyle->SetOptStat(0);  // 关闭统计框
    gStyle->SetPalette(kRainBow);  // 彩虹色标（深色表示高事例数）
    gStyle->SetPadGridX(1);  // 显示网格线
    gStyle->SetPadGridY(1);

    // 打开ROOT文件
    TFile *file = TFile::Open("../build/root_file/CryMu_0.root");
    if (!file || file->IsZombie()) {
        printf("无法打开文件!\n");
        return;
    }

    // 获取树（假设树名为"tree"，若实际名称不同需修改）
    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        printf("无法找到树!\n");
        return;
    }

    // 创建5x5直方图（x: 列索引0-4, y: 行索引0-4）
    TH2D *h2 = new TH2D("h2", "Edeps.Id热图(6-30);列索引;行索引",
                        5, -0.5, 4.5,  // x轴：5个bin，覆盖0-4整数坐标
                        5, -0.5, 4.5); // y轴：5个bin，覆盖0-4整数坐标

    // 填充直方图：Id=6~30映射到5x5网格，统计事例数
    tree->Draw("int((Edeps.Id-6)/5):int((Edeps.Id-6)%5) >> h2",  // y:x坐标（行:列）
               "Edeps.Id>=6 && Edeps.Id<=30",              // 筛选Id范围
               "COLZ");                                    // 彩色填充+色标

    // 美化坐标轴标签（显示实际Id范围）
    h2->GetXaxis()->SetTitle("x");
    h2->GetYaxis()->SetTitle("y");
    h2->SetTitle("5x5 Id");

    // 保存图片（支持png/pdf/eps等格式）
    c1->SaveAs("png/pbwo4_heatmap.png");
    printf("热图已保存为 pbwo4_heatmap.png\n");

    // 清理内存
    delete h2;
    file->Close();
    delete file;
    delete c1;
    
}
