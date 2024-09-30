#ifndef DRAW_H
#define DRAW_H

#pragma cling add_include_path("/opt/homebrew/opt/boost/include")
#pragma cling add_library_path("/opt/homebrew/opt/boost/lib")
#pragma cling load("libboost_filesystem.dylib")
#pragma cling load("libboost_system.dylib")

// root headers are listed below
#include "TCanvas.h"
#include "TColor.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"

#include <vector>
#include <array>

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <fstream>
#include "vector.h"

// 4个GEM读出层的Z坐标
double Z1 = -513.34;
double Z2 = -500.;
double Z3 = 513.24;
double Z4 = 526.58;
using namespace std;
namespace fs = boost::filesystem;
// 使用默认的z方向间距 Lz=Z3-Z2 和 Nz，注释掉 #define SIM 表示自定义z方向间距 Lz_ 和 Nz_，且 z=0 在中央，但保证 dz 不变
// #define SIM
double Lz_ = 1000.; // 自定义z方向间距
double Z2_ = -Lz_ / 2;
double Z3_ = Lz_ / 2;

// 分格数，默认为50，不要比这个小
#define Nx 500
#define Ny 50
#define Nz 20 // Nz 要能被 Z2 和 Z3 整除
#define SPACE_N (Nx * Ny * Nz)

// 图的全宽（mm），不要比物体区域大太多
#define Lx 1000
#define Ly 1000
#define Lz 1000

#define INDEX(i, j, k) ((i) + Nx * (j) + Nx * Ny * (k))
//#define Lz (Z3 - Z2)

template <typename T>
class Matrix;

void DumpStep(std::vector<double> &lambda,
              const char *title,
              const char *directory,
              int step, int flag = 0);

void DumpPOCAData(const Matrix<double> &L,
                  const Matrix<double> &W00,
                  const Matrix<double> &W01,
                  const Matrix<double> &W11,
                  const std::vector<std::array<double, 4>> D,
                  const std::vector<double> pr2,
                  const char *directory);

void DrawL(std::vector<double> &stigma,
           const char *title,
           int pcount, int iter, int shrink_size);

void DrawS(std::vector<double> &stigma,
           const char *title,
           int pcount, int iter);
void DrawCut(std::vector<double> &lambda, const char *title, int iter);
extern double lambda_none[SPACE_N];

// 绘制一维图像
// 例如 x 方向的投影就是将所有 x 方向索引值相等的所有子区域的 lambda 值累加起来
void DrawProjection1D(vector<double> &lambda,
                      const char *title,
                      const char *filename)
{
#ifndef SIM
    //double dz = (Z3 - Z2) / Nz;
    double dz = (Z3_ - Z2_) / Nz;
    int Nz_ = (Z3_ - Z2_) / dz; // 向零取整（舍去小数部分）
    int N0 = (Z2_ - Z2) / dz; // 向零取整（舍去小数部分）
#endif
#ifdef SIM
    Z2_ = Z2;
    Z3_ = Z3;
    int Nz_ = Nz;
    int N0 = 0;
#endif
    TH3D graph3d("result", title,
                 Nx, -Lx / 2, Lx / 2, // x 方向范围：[-Lx/2, Lx/2]
                 Ny, -Ly / 2, Ly / 2, // y 方向范围：[-Ly/2, Ly/2]
                 Nz_, Z2_, Z3_);      // z 方向范围：[Z2_, Z3_]
    // standard values: Fe=0.00117, Pb=0.00367, W=0.00587
    TH3D std_graph3d("standard", title, // I don't know what this mean
                     Nx, -Lx / 2, Lx / 2,
                     Ny, -Ly / 2, Ly / 2,
                     Nz_, Z2_, Z3_);

    double rms = 0;

    for (int k = 0; k < Nz_; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
            {
                double value = lambda[INDEX(i, j, N0 + k)];
                std_graph3d.SetBinContent(i + 1, j + 1, k, 1e-20); // 原来z方向的坐标是 Nz_ - k，原因未知，这会导致z方向镜像翻转
                graph3d.SetBinContent(i + 1, j + 1, k, value + 1e-20); // 原来z方向的坐标是 Nz_ - k，原因未知，这会导致z方向镜像翻转
                rms += value * value;
            }

    cout << "sqrt(rms) = " << sqrt(rms / SPACE_N) << '\n';

    TCanvas canvas("result-canvas", title, 1600, 800);
    canvas.Divide(2, 1);

    canvas.cd(1);
    auto std_x_panel = std_graph3d.Project3D("x");
    std_x_panel->SetLineColor(kRed);
    std_x_panel->SetStats(false);

    auto x_panel = graph3d.Project3D("x");
    x_panel->SetLineColor(kBlue);
    x_panel->SetStats(false);

    if (x_panel->GetMaximum() > std_x_panel->GetMaximum())
    {
        x_panel->Draw("HIST");
        std_x_panel->Draw("HIST SAME");
        x_panel->GetXaxis()->SetTitle("X/mm");
        x_panel->GetYaxis()->SetTitle("#lambda");
        x_panel->GetXaxis()->CenterTitle();
        x_panel->GetYaxis()->CenterTitle();
    }
    else
    {
        std_x_panel->Draw("HIST");
        x_panel->Draw("HIST SAME");
        std_x_panel->GetXaxis()->SetTitle("X/mm");
        std_x_panel->GetYaxis()->SetTitle("#lambda");
        std_x_panel->GetXaxis()->CenterTitle();
        std_x_panel->GetYaxis()->CenterTitle();
    }

    canvas.cd(2);
    auto std_z_panel = std_graph3d.Project3D("z");
    std_z_panel->SetLineColor(kRed);
    std_z_panel->Draw("HIST");
    std_z_panel->SetStats(false);

    auto z_panel = graph3d.Project3D("z");
    z_panel->Draw("SAME HIST");
    z_panel->SetStats(false);

    if (z_panel->GetMaximum() > std_z_panel->GetMaximum())
    {
        z_panel->Draw("HIST");
        std_z_panel->Draw("HIST SAME");
        z_panel->GetXaxis()->SetTitle("Z/mm");
        z_panel->GetYaxis()->SetTitle("#lambda");
        z_panel->GetXaxis()->CenterTitle();
        z_panel->GetYaxis()->CenterTitle();
    }
    else
    {
        std_z_panel->Draw("HIST");
        z_panel->Draw("HIST SAME");
        std_z_panel->GetXaxis()->SetTitle("Z/mm");
        std_z_panel->GetYaxis()->SetTitle("#lambda");
        std_z_panel->GetXaxis()->CenterTitle();
        std_z_panel->GetYaxis()->CenterTitle();
    }

    for (int i = 1; i <= 2; ++i)
    {
        auto *pad = canvas.GetPad(i);
        pad->SetLeftMargin(0.14);
        pad->SetTopMargin(0.12);
        pad->SetBottomMargin(0.12);
    }

    canvas.Draw();
    canvas.SaveAs(filename);
}

// 绘制二维图像
void DrawProjection2D(vector<double> &lambda,
                      const char *title,
                      const char *filename, int flag)
{
#ifndef SIM
    //double dz = (Z3 - Z2) / Nz;
    double dz = (Z3_ - Z2_) / Nz;
    int Nz_ = (Z3_ - Z2_) / dz; // 向零取整（舍去小数部分）
    int N0 = (Z2_ - Z2) / dz; // 向零取整（舍去小数部分）
#endif
#ifdef SIM
    Z2_ = Z2;
    Z3_ = Z3;
    int Nz_ = Nz;
    int N0 = 0;
#endif
    TH3D graph3d("result", title,
                 Nx, -Lx / 2, Lx / 2,
                 Ny, -Ly / 2, Ly / 2,
                 Nz_, Z2_, Z3_);
    TStyle *default_style = new TStyle("Default", "Default Style");
    default_style->SetPalette(55);

    for (int k = 0; k < Nz_; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
            {
                // remove most points in 3d graph (I don't know why)
                double value = lambda[INDEX(i, j, N0 + k)];
                graph3d.SetBinContent(i + 1, j + 1, k, ((value > 1e-3) ? value : 1e-20)); // 原来z方向的坐标是 Nz_ - k，原因未知，这会导致z方向镜像翻转
            }

    // graph3d.Scale(10000 / graph3d.Integral());

    TCanvas canvas("result-canvas", title, 1600, 1600);
    canvas.Divide(2, 2);
    canvas.cd(1);
    //graph3d.Draw("LEGO2");

    // we wanna draw another rotatable version in a root session
    // but TFile always had a linker error
    // TFile *stigma = new TFile("sig-data.root","recreate");
    // graph3d.Write();
    // stigma->Close();

    // do not use the Projection function of root since it always produces empty results
    // creating a new histogram can easily solve the problem
    canvas.cd(2);                         // X-Y Graph
    delete gROOT->FindObject("xy_panel"); // prevent warning: "Potential memory leak"
    TH2D *xy_panel = new TH2D("xy_panel", "xy_panel", Nx, -Lx / 2, Lx / 2, Ny, -Ly / 2, Ly / 2);
    for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            double value = 0;
            for (int k = 0; k < Nz_; ++k)
            {
                //value += lambda[INDEX(j, i, N0 + k)];
                value += lambda[INDEX(i, j, N0 + k)];
            }
            xy_panel->SetBinContent(i + 1, j + 1, ((value > 1e-20) ? value : 1e-20));
        }
    if (flag)
        xy_panel->SetMaximum(8);
    xy_panel->Draw("colz");
    xy_panel->SetStats(false);
    xy_panel->GetXaxis()->SetTitle("X/mm");
    xy_panel->GetYaxis()->SetTitle("Y/mm");
    xy_panel->GetXaxis()->CenterTitle();
    xy_panel->GetYaxis()->CenterTitle();

    canvas.cd(3);                         // X-Z Graph
    delete gROOT->FindObject("xz_panel"); // prevent warning: "Potential memory leak"
    TH2D *xz_panel = new TH2D("xz_panel", "xz_panel", Nx, -Lx / 2, Lx / 2, Nz_, Z2_, Z3_);
    for (int j = 0; j < Nz_; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            double value = 0;
            for (int k = 0; k < Ny; ++k)
            {
                value += lambda[INDEX(i, k, j + N0)]; // 原来z方向的坐标是 Nz_ + N0 - j - 1，原因未知，这会导致z方向镜像翻转
            }
            xz_panel->SetBinContent(i + 1, j + 1, ((value > 1e-20) ? value : 1e-20));
        }
    if (flag)
        xz_panel->SetMaximum(8);
    xz_panel->Draw("colz");
    xz_panel->SetStats(false);
    xz_panel->GetXaxis()->SetTitle("X/mm");
    xz_panel->GetYaxis()->SetTitle("Z/mm");
    xz_panel->GetXaxis()->CenterTitle();
    xz_panel->GetYaxis()->CenterTitle();

    canvas.cd(4);                         // Y-Z Graph
    delete gROOT->FindObject("yz_panel"); // prevent warning: "Potential memory leak"
    TH2D *yz_panel = new TH2D("yz_panel", "yz_panel", Ny, -Ly / 2, Ly / 2, Nz_, Z2_, Z3_);
    for (int j = 0; j < Nz_; ++j)
        for (int i = 0; i < Ny; ++i)
        {
            double value = 0;
            for (int k = 0; k < Nx; ++k)
            {
                value += lambda[INDEX(k, i, j + N0)]; // 原来z方向的坐标是 Nz_ + N0 - j - 1，原因未知，这会导致z方向镜像翻转
            }
            yz_panel->SetBinContent(i + 1, j + 1, ((value > 1e-20) ? value : 1e-20));
        }
    if (flag)
        yz_panel->SetMaximum(8);
    yz_panel->Draw("colz");
    yz_panel->SetStats(false);
    yz_panel->GetXaxis()->SetTitle("Y/mm");
    yz_panel->GetYaxis()->SetTitle("Z/mm");
    yz_panel->GetXaxis()->CenterTitle();
    yz_panel->GetYaxis()->CenterTitle();

    for (int i = 1; i <= 4; ++i)
    {
        auto *pad = canvas.GetPad(i);
        pad->SetLeftMargin(0.15);
        pad->SetRightMargin(0.17);
        pad->SetTopMargin(0.12);
        pad->SetBottomMargin(0.11);
    }

    // canvas.TAttPad::SetRightMargin(0.11);
    canvas.SetTopMargin(0.12);
    canvas.SetBottomMargin(0.11);

    canvas.Draw();
    canvas.SaveAs(filename);
}

void DumpStep(std::vector<double> &lambda,
              const char *title,
              const char *directory,
              int step, int flag)
{
    try
    {
        if (!fs::exists(directory))
            fs::create_directory(directory);
        fs::path dir(directory);
        string filename = "lambda." + to_string(step);
        ofstream fout((dir / filename).string());
        if (fout)
        {
            fout.precision(12);
            fout << Nx << ' ' << Ny << ' ' << Nz << '\n'               // 第一行是分格数 Nx、Ny、Nz（默认都为30）
                 << Lx << ' ' << Ly << '\n'                            // 第二行是绘图区域 x、y 边长 Lx、Ly（默认都为160）
                 << Z1 << ' ' << Z2 << ' ' << Z3 << ' ' << Z4 << '\n'; // 第三行是四块 RPC 板的 z 方向位置
            for (auto value : lambda)
                fout << value << '\n';
        }
        fout.close();
        // 绘制一维图像
        DrawProjection1D(lambda, title,
                         (dir / ("1d." + filename + ".eps")).string().c_str());
        // 绘制二维图像
        DrawProjection2D(lambda, title,
                         (dir / ("2d." + filename + ".eps")).string().c_str(), flag);
    }
    catch (...)
    {
        cerr << "Can not dump for step " << step << endl;
    }
}

void DumpMatrix(const Matrix<double> &data, ofstream &out)
{
    out.precision(12);
    for (auto &path : data.data())
    {
        for (auto &index : path)
            out << index.first << ' ' << index.second << ' ';
        out << '\n';
    }
}

void DumpPOCAData(const Matrix<double> &L,
                  const Matrix<double> &W00,
                  const Matrix<double> &W01,
                  const Matrix<double> &W11,
                  const vector<array<double, 4>> D,
                  const vector<double> pr2,
                  const char *directory)
{
    try
    {
        if (!fs::exists(directory))
            fs::create_directory(directory);
        fs::path dir(directory);
        ofstream fout((dir / "L.txt").string());
        if (fout)
            DumpMatrix(L, fout);
        fout.close();

        fout.open((dir / "W00.txt").string());
        if (fout)
            DumpMatrix(W00, fout);
        fout.close();

        fout.open((dir / "W01.txt").string());
        if (fout)
            DumpMatrix(W01, fout);
        fout.close();

        fout.open((dir / "W11.txt").string());
        if (fout)
            DumpMatrix(W11, fout);
        fout.close();

        fout.open((dir / "pr2.txt").string());
        if (fout)
        {
            for (auto &v : pr2)
                fout << v << '\n';
        }
        fout.close();
        fout.open((dir / "D.txt").string());
        if (fout)
        {
            for (auto &v : D)
                fout << v[0] << ' ' << v[1] << ' '
                     << v[2] << ' ' << v[3] << '\n';
        }
        fout.close();
    }
    catch (...)
    {
        cerr << "Can not dump length " << endl;
    }
}

// a new draw function, this kinda useful
void DrawL(std::vector<double> &stigma,
           const char *title,
           int pcount, int iter, int shrink_size)
{
    // calculate the biggest value
    double max_value = -1;
    for (int i = 0; i < pcount; i++)
    {
        if (fabs(stigma[i]) > max_value)
            max_value = fabs(stigma[i]);
    }

    // you should always create the canvas first
    TCanvas canvas("result-canvas", title, 1600, 1600);
    TH1D graph1d("result", title,
                 1000, 0, max_value / shrink_size);

    for (int k = 0; k < pcount; k++)
    {
        // graph1d.SetBinContent(k,fabs(stigma[k]));
        graph1d.Fill(fabs(stigma[k]));
    }
    graph1d.Draw();
    canvas.Draw();

    string titled(title);
    if (!fs::exists(title))
        fs::create_directory(title);
    string savename = titled + "/" + titled + "." + to_string(iter) + ".eps";
    canvas.SaveAs(savename.c_str());
}

// Draw a stigma[k]-k graph, not a filled histogram
// use iter=0 as default
void DrawS(std::vector<double> &stigma,
           const char *title,
           int pcount, int iter)
{
    // you should always create the canvas first
    TCanvas canvas("result-canvas", title, 1600, 1600);
    TH1D graph1d("result", title,
                 pcount, 0, pcount);
    graph1d.SetStats(kFALSE);

    for (int k = 1; k <= pcount; k++)
    {
        graph1d.SetBinContent(k, fabs(stigma[k - 1]));
    }
    graph1d.Draw();
    canvas.Draw();

    string titled(title);
    if (!fs::exists(title))
        fs::create_directory(title);
    string savename = titled + "/" + titled + "." + to_string(iter) + ".eps";
    canvas.SaveAs(savename.c_str());
}

// draw the cut surface of y=0, which is middle, compare with projectiony
void DrawCut(std::vector<double> &lambda, const char *title, int iter)
{
    TCanvas canvas("result-canvas", title, 1600, 1600);
    delete gROOT->FindObject("cut_panel"); // prevent warning: "Potential memory leak"
    TH2D *cut_panel = new TH2D("cut_panel", "cut_panel", Nx, -Lx / 2, Lx / 2, Nz, Z2, Z3);
    for (int j = 0; j < Nz; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            cut_panel->SetBinContent(i + 1, j + 1, lambda[INDEX(i, Ny / 2, j)] + 1e-20);
        }
    TStyle *default_style = new TStyle("Default", "Default Style");
    default_style->SetPalette(55);

    // cut_panel->Scale(100 / cut_panel->GetMaximum());
    cut_panel->Draw("colz");
    cut_panel->SetStats(false);
    cut_panel->GetXaxis()->SetTitle("Y/mm");
    cut_panel->GetYaxis()->SetTitle("X/mm");
    cut_panel->GetXaxis()->CenterTitle();
    cut_panel->GetYaxis()->CenterTitle();
    canvas.SetLeftMargin(0.15);
    canvas.SetRightMargin(0.17);
    canvas.SetTopMargin(0.12);
    canvas.SetBottomMargin(0.11);
    canvas.Draw();

    string titled(title);
    if (!fs::exists(title))
        fs::create_directory(title);
    string savename = titled + "/" + titled + "." + to_string(iter) + ".eps";
    canvas.SaveAs(savename.c_str());
}

#endif /* DRAW_H */
