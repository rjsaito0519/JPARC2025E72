// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"

Config& conf = Config::getInstance();

void analyze(TString path, TString particle){    
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0/256, 160.0/256, 44.0/256);
    
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x"); // x軸のタイトルサイズ
    gStyle->SetTitleSize(0.06, "y"); // y軸のタイトルサイズ
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);


    // +-----------+
    // | load file |
    // +-----------+
    auto *f = new TFile(path.Data());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << path << std::endl;
        return;
    }

    const char* treeName = nullptr;
    TString in_or_out = "";
    if (f->Get("bcout")) {
        treeName = "bcout";
        in_or_out = "2";
    }
    else if (f->Get("bcin")) {
        treeName = "bcin";
        in_or_out = "1";
    }
    else {
        std::cerr << "Error: neither 'bcout' nor 'bcin' tree exists.\n";
        return;
    }
    TTreeReader reader(treeName, f);
    TTreeReaderValue<unsigned int> run_number(reader, "run_number");
    reader.SetEntry(0);
    Int_t run_num = *run_number;

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = Form("%s/root/run%05d_BLC%s_residual_%s.root", OUTPUT_DIR.Data(), run_num, in_or_out.Data(), particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- tdc ----------
    // hname = f'{name}_Track_Residual{key}_plane{i}{mh.beamflag_for_param}' '_vs_DriftLength
    TH1D *h_blca_residual[conf.num_of_ch.at("blc")];
    TH1D *h_blcb_residual[conf.num_of_ch.at("blc")];
    TH2D *h_blca_XT[conf.num_of_ch.at("blc")];
    TH2D *h_blcb_XT[conf.num_of_ch.at("blc")];
    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++ ) {
        h_blca_residual[i] = (TH1D*)f->Get(Form("BLC%sa_Track_Residual_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        h_blcb_residual[i] = (TH1D*)f->Get(Form("BLC%sb_Track_Residual_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        h_blca_XT[i] = (TH2D*)f->Get(Form("BLC%sa_Track_Residual_vs_DriftLength_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        h_blcb_XT[i] = (TH2D*)f->Get(Form("BLC%sb_Track_Residual_vs_DriftLength_plane%d_%s", in_or_out.Data(), i, particle.Data()));
    }
 
    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_BLC%s_residual_%s.pdf", OUTPUT_DIR.Data(), run_num, in_or_out.Data(), particle.Data());

    // -- container -----
    std::vector<FitResult> blca_residual;
    std::vector<FitResult> blcb_residual;

    auto c_blc = new TCanvas("blc", "", 1500, 1200);
    c_blc->Divide(cols, rows);
    c_blc->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++) {
        if (nth_pad > max_pads) {
            c_blc->Print(pdf_path);
            c_blc->Clear();
            c_blc->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result = ana_helper::residual_fit(h_blca_residual[i], c_blc, nth_pad);
        blca_residual.push_back(result);
        c_blc->cd(++nth_pad);
        h_blca_XT[i]->Draw("colz");

        result = ana_helper::residual_fit(h_blcb_residual[i], c_blc, ++nth_pad);
        blcb_residual.push_back(result);
        c_blc->cd(++nth_pad);
        h_blcb_XT[i]->Draw("colz");
    }
    c_blc->Print(pdf_path);
    c_blc->Print(pdf_path + "]"); // end
    delete c_blc;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> residual_p0_val;
    std::vector<Double_t> residual_p0_err;
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("residual_p0_val", &residual_p0_val);
    tree->Branch("residual_p0_err", &residual_p0_err);
    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++) {
        ch = i;
        residual_p0_val.clear();
        residual_p0_err.clear();
    
        // -- residual -----
        residual_p0_val.push_back(blca_residual[i].par[1]);
        residual_p0_err.push_back(blca_residual[i].err[1]);

        residual_p0_val.push_back(blcb_residual[i].par[1]);
        residual_p0_err.push_back(blcb_residual[i].err[1]);

        tree->Fill();
    }    
    
    fout->cd();
    tree->Write();
    fout->Close(); 
}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <root file path> <particle>" << std::endl;
        return 1;
    }
    TString path = argv[1];
    TString particle = argv[2];
    if (particle != "Pi" && particle != "K") {
        std::cerr << "Error: Unexpected particle name: " << particle << std::endl;
        return 1;
    }

    analyze(path, particle);
    return 0;
}