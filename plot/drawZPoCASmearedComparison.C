#include <math.h>
#include <fstream>
using namespace std;

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TColor.h"

int drawZPoCASmearedComparison() {
    // 打开两个ROOT文件
    TFile *f1 = TFile::Open("../build/root_file/CryMu.root", "read");
    TFile *f2 = TFile::Open("../build/root_file_2p/CryMu.root", "read");
    
    // 获取树（根据用户提供的结构，树名为tree）
    TTree *tree1 = (TTree*)f1->Get("tree");
    TTree *tree2 = (TTree*)f2->Get("tree");
    
    // 创建画布
    TCanvas *c = new TCanvas("c", "ZPoCASmeared-Events", 800, 600);
    c->SetGrid();
    
    // 创建直方图（设置合适的范围和分箱数）
    TH1D *h2 = new TH1D("h2", "addPb_EdepX", 100, -150, 150); // 根据示例数据调整范围
    TH1D *h1 = new TH1D("h1", "", 100, -150, 150);
    
    // 设置直方图样式
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(2);
    h2->SetLineColor(kRed);
    h2->SetLineWidth(2);
    h2->GetXaxis()->SetTitle("Edeps.X");
    h2->GetYaxis()->SetTitle("Events");
    h2->SetStats(0); // 关闭统计框
    
    // 填充直方图
    //tree1->Draw("ZPoCASmeared>>h1", "", "HIST");
    //tree2->Draw("ZPoCASmeared>>h2", "", "HIST SAME");
    tree2->Draw("Edeps.X>>h2", "Edeps.Id==4 && Edeps.Pid==11", "HIST");
    tree1->Draw("Edeps.X>>h1", "Edeps.Id==4 && Edeps.Pid==11", "HIST SAME");

    // 添加图例
    TLegend *leg = new TLegend(0.76, 0.76, 0.9, 0.9); // 缩小总面积（x1,y1,x2,y2）
    leg->SetTextSize(0.035); // 放大字体（默认约0.025，根据画布大小调整）
    leg->SetMargin(0.15);    // 减小图例内部边距（可选，让文字更紧凑）
    leg->AddEntry(h1, "add Lead", "l");
    leg->AddEntry(h2, "no Lead", "l");
    leg->Draw();
    
    // 保存图片
    //c->SaveAs("png/addPb_ZPoCASmeared_Comparison.png");
    c->SaveAs("png/addPb_EdepX.png");
    // 清理资源
    f1->Close();
    f2->Close();
    delete leg;
    delete c;
    delete h1;
    delete h2;
    
    return 0;
}
