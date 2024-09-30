#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"
#include "TList.h"

void rebinHistograms() {
    // 打开包含直方图的ROOT文件
    TFile *file = new TFile("../Hist.root", "READ");

    // 创建新的ROOT文件，用于保存重新分bin后的直方图
    TFile *newFile = new TFile("Hist.root", "RECREATE");

    // 自定义不均匀分bin的边界数组
    Double_t customBins[52];
    customBins[0]=-1;
    double binwidth=(1.-0.)/50.;
    for (int i = 0; i <= 50; ++i) {
        customBins[i+1] = 0 + (double)i*binwidth;
    }
    //for(int i=0; i<52; i++) cout<<"bin "<<i<<": "<<customBins[i]<<endl;

    // 循环处理ROOT文件中的所有直方图
    TIter nextkey(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)nextkey())) {
        TObject *obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom(TH1::Class())) {
            TH1F *histogram = (TH1F*)obj;

            // 创建新的直方图，重新分bin
            TH1F *newHistogram = new TH1F(histogram->GetName(), histogram->GetTitle(), 51, customBins);
            newHistogram->SetBinContent(1,histogram->Integral(1,50));
            for (int i = 51; i <= 100; i++) newHistogram->SetBinContent(i-49, histogram->GetBinContent(i));
            histogram->Print();
            newHistogram->Print();

            // 将新的直方图保存到新的ROOT文件中
            newHistogram->Write();
        }
    }

    // 关闭输入文件和新文件
    file->Close();
    newFile->Close();
}

