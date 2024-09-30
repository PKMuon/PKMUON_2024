#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <iostream>
#include <map>
#include <vector>

void Hist(const char* infile = "../../build/root_file/CryMu_new.root") {
    TFile *file_in = new TFile(infile);
    if (!file_in->IsOpen()) {
        std::cerr << "Error: Could not open file " << infile << std::endl;
        return;
    }
    std::cout << "Opening ROOT File." << std::endl;

    TTree *tree_in = (TTree *)file_in->Get("T1");
    if (!tree_in) {
        std::cerr << "Error: Could not find tree T1 in file " << infile << std::endl;
        file_in->Close();
        return;
    }

    double costheta_smeared;
    Int_t pdgCode;

    if (tree_in->SetBranchAddress("costheta_smeared", &costheta_smeared) < 0) {
        std::cerr << "Error: could not set branch address for 'costheta_smeared'" << std::endl;
        file_in->Close();
        return;
    }
    if (tree_in->SetBranchAddress("PDG", &pdgCode) < 0) {
        std::cerr << "Error: could not set branch address for 'PDG'" << std::endl;
        file_in->Close();
        return;
    }

    Long64_t nentry = tree_in->GetEntries();
    if (nentry <= 0) {
        std::cerr << "Error: No entries found in the tree" << std::endl;
        file_in->Close();
        return;
    }

    // 使用 map 存储不同 PDG 类型的角度数据和数量
    std::map<int, std::vector<double>> angles;
    std::map<int, int> particleCounts; // 存储每种粒子的数量

    std::map<int, std::string> pdgNames = {
        {13, "mu-"},
        {-13, "mu+"},
        {11, "e-"},
        {-11, "e+"},
        {211, "pi+"},
        {-211, "pi-"},
        {2212, "p+"},
        {-2212, "p-"},
        {22, "gamma"},
        {2112, "n"}
    };

    // 遍历所有事件
    for (Long64_t ientry = 0; ientry < nentry; ++ientry) {
        tree_in->GetEntry(ientry);
        if (costheta_smeared > 0) { // 只保存正值
            angles[pdgCode].push_back(costheta_smeared);
            particleCounts[pdgCode]++; 
        }
    }

    std::cout << "Particle counts:" << std::endl;
    for (const auto& entry : particleCounts) {
        int pdg = entry.first;
        int count = entry.second;
        std::string particleName = (pdgNames.count(pdg) ? pdgNames[pdg] : Form("PDG %d", pdg));
        std::cout << particleName << " (PDG: " << pdg << "), Count: " << count << std::endl;
    }

    TCanvas *canvas = new TCanvas("canvas", "CRY through RPC", 800, 600);
    TLegend *legend = new TLegend(0.1, 0.5, 0.5, 0.9); 
    gStyle->SetOptStat(0); 

    int colorIndex = 1; 

    TH1D *h_total = new TH1D("h_total", "AllParticle", 100, 0, 1);

    // 遍历 map 中的不同 PDG 类型
    for (const auto& entry : angles) {
        int pdg = entry.first;
        const std::vector<double>& angleData = entry.second;

        // 创建对应 PDG 类型的直方图
        TH1D *h = new TH1D(Form("h_%d", pdg), "CRY through RPC", 100, 0, 1);
        h->SetXTitle("costheta"); 
        for (double angle : angleData) {
            h->Fill(angle);
            h_total->Fill(angle); // 填充到总直方图
        }

        h->SetLineColor(colorIndex % 10 + 1);
        h->SetLineWidth(2);
        h->Draw("HIST SAME"); 
        canvas->Update(); 

        std::string legendEntry = (pdgNames.count(pdg) ? pdgNames[pdg] : Form("PDG %d", pdg));
        legend->AddEntry(h, legendEntry.c_str(), "l");

        colorIndex++;
    }

    // 归一化总直方图
    int totalCount = h_total->GetEntries(); 
    if (totalCount > 0) {
        h_total->Scale(1.0 / totalCount); 
    }

    h_total->SetXTitle("costheta");
    h_total->SetLineColor(kBlack); 
    h_total->SetLineWidth(2);      
    h_total->Draw("HIST SAME");          

    canvas->SetLogy();

    double totalMaxY = h_total->GetMaximum();

    if (totalMaxY > 0) {
        h_total->SetMaximum(totalMaxY * 5); 
        h_total->SetMinimum(5e-7);             
    }

    // 设置每个 PDG 直方图的 y 轴范围并进行归一化
    for (const auto& entry : angles) {
        int pdg = entry.first;
        TH1D *h = (TH1D *)canvas->FindObject(Form("h_%d", pdg));
        if (h) {
            h->Scale(1.0 / totalCount); 
            h->SetMaximum(totalMaxY * 5); 
            h->SetMinimum(5e-7);            
        }
    }

    legend->AddEntry(h_total, "AllParticle", "l");
    legend->Draw();

    canvas->SaveAs("scattering_angle_distributions.pdf");

    file_in->Close();
    delete legend;
    delete canvas;
    delete h_total; 

    std::cout << "Plotting complete." << std::endl;
}
