#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

void e_trace(const char *infile = "../../build/root_file/CryMuAna.root", const char *outfile = "../../build/root_file/e_trace_output.root") {
    // 打开输入文件和树
    TFile *file_in = TFile::Open(infile, "READ");
    if (!file_in || file_in->IsZombie()) {
        cerr << "Error: 无法打开输入文件 " << infile << endl;
        return;
    }

    TTree *tree_in = (TTree*)file_in->Get("tree");
    if (!tree_in) {
        cerr << "Error: 无法找到输入树 'tree'" << endl;
        file_in->Close();
        return;
    }

    // 设置输入分支地址
    vector<Double_t> *TotalEdep = nullptr;
    vector<Double_t> *XEdep = nullptr;
    vector<Double_t> *YEdep = nullptr;
    tree_in->SetBranchAddress("TotalEdep", &TotalEdep);
    tree_in->SetBranchAddress("XEdep", &XEdep);
    tree_in->SetBranchAddress("YEdep", &YEdep);

    // 创建输出文件和树（复制原结构）
    TFile *file_out = TFile::Open(outfile, "RECREATE");
    TTree *tree_out = tree_in->CloneTree(0); // 复制分支结构，不复制数据

    // 添加新分支：总沉积点位置（X_reco, Y_reco）
    Double_t X_reco, Y_reco;
    tree_out->Branch("X_reco", &X_reco, "X_reco/D");
    tree_out->Branch("Y_reco", &Y_reco, "Y_reco/D");

    // 事件循环
    Long64_t nentries = tree_in->GetEntries();
    Long64_t passed = 0;
    for (Long64_t i = 0; i < nentries; ++i) {
        tree_in->GetEntry(i);

        // 1. 查找TotalEdep[6..30]中的最大值及对应ID
        Double_t max_edep = -1;
        Int_t max_id = -1;
        for (Int_t id = 6; id <= 30; ++id) { // 5x5方阵ID范围：6-30
            if (id >= (Int_t)TotalEdep->size()) continue; // 索引越界检查
            Double_t edep = (*TotalEdep)[id];
            if (edep > max_edep) {
                max_edep = edep;
                max_id = id;
            }
        }

        // 2. 筛选最大值 > 150 的事例
        if (max_edep <= 220) continue;
        passed++;

        // 3. 计算中心ID在5x5方阵中的坐标 (row, col)
        Int_t offset = max_id - 6; // 转换为0-24的方阵索引
        Int_t row = offset / 5;    // 行坐标（0-4）
        Int_t col = offset % 5;    // 列坐标（0-4）

        // 4. 3x3网格加权平均计算
        Double_t sum_edep = 0, sum_x = 0, sum_y = 0;
        for (Int_t dr = -1; dr <= 1; ++dr) { // 行方向偏移（-1,0,1）
            for (Int_t dc = -1; dc <= 1; ++dc) { // 列方向偏移（-1,0,1）
                Int_t r = row + dr;
                Int_t c = col + dc;

                // 检查是否在5x5方阵内（防止越界）
                if (r < 0 || r >= 5 || c < 0 || c >= 5) continue;

                // 计算当前网格的ID
                Int_t current_id = 6 + r * 5 + c;

                // 跳过无能量沉积的板
                if (current_id >= (Int_t)TotalEdep->size()) continue;
                Double_t edep = (*TotalEdep)[current_id];
                if (edep <= 0) continue;

                // 累加加权坐标（XEdep/YEdep以TotalEdep为权重）
                sum_edep += edep;
                sum_x += (*XEdep)[current_id] * edep;
                sum_y += (*YEdep)[current_id] * edep;
            }
        }

        // 计算加权平均位置（避免除零）
        X_reco = sum_edep > 0 ? sum_x / sum_edep : 0;
        Y_reco = sum_edep > 0 ? sum_y / sum_edep : 0;

        // 写入筛选后的事件
        tree_out->Fill();
    }

    // 保存输出文件
    file_out->cd();
    tree_out->Write();
    cout << "处理完成：共 " << nentries << " 个事例，筛选出 " << passed << " 个有效事例" << endl;

    // 清理资源
    file_out->Close();
    file_in->Close();
    delete file_out;
    delete file_in;
}
