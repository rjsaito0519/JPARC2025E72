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
    TTreeReader reader("hodo", f);
    TTreeReaderValue<unsigned int> run_number(reader, "run_number");
    reader.SetEntry(0);
    Int_t run_num = *run_number;

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = Form("%s/root/run%05d_BH2_PHC_%s.root", OUTPUT_DIR.Data(), run_num, particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- BTOF vs dE ----------
    TH2D *h_bh2_btof_vs_de[2][conf.num_of_ch.at("bh2")];
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++ ) h_bh2_btof_vs_de[0][i] = (TH2D*)f->Get(Form("BH2_seg%dU_BTOF_vs_DeltaE_%s", i, particle.Data()));
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++ ) h_bh2_btof_vs_de[1][i] = (TH2D*)f->Get(Form("BH2_seg%dD_BTOF_vs_DeltaE_%s", i, particle.Data()));

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_BH2_PHC_%s.pdf", OUTPUT_DIR.Data(), run_num, particle.Data());

    // -- container -----
    std::vector<FitResult> phc_up;
    std::vector<FitResult> phc_down;

    auto c_bh2 = new TCanvas("bh2", "", 1500, 1200);
    c_bh2->Divide(cols, rows);
    c_bh2->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++) {
        if (nth_pad > max_pads) {
            c_bh2->Print(pdf_path);
            c_bh2->Clear();
            c_bh2->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result;
        // -- UP -----
        result = ana_helper::phc_fit(h_bh2_btof_vs_de[0][i], c_bh2, nth_pad);
        phc_up.push_back(result);
        nth_pad++;

        // -- DOWN -----
        result = ana_helper::phc_fit(h_bh2_btof_vs_de[1][i], c_bh2, nth_pad);
        phc_down.push_back(result);
        nth_pad++;
    }
    c_bh2->Print(pdf_path);
    c_bh2->Print(pdf_path + "]"); // end
    delete c_bh2;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> p0_val, p1_val, p2_val; 
    std::vector<Double_t> p0_err, p1_err, p2_err; 
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("p0_val", &p0_val);
    tree->Branch("p1_val", &p1_val);
    tree->Branch("p2_val", &p2_val);
    tree->Branch("p0_err", &p0_err);
    tree->Branch("p1_err", &p1_err);
    tree->Branch("p2_err", &p2_err);
    
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++) {
        ch = i;
        p0_val.clear(); p1_val.clear(); p2_val.clear();
        p0_err.clear(); p1_err.clear(); p2_err.clear();

        // -- p0 -----
        p0_val.push_back(phc_up[i].par[0]);
        p0_val.push_back(phc_down[i].par[0]);
        p0_err.push_back(phc_up[i].err[0]);
        p0_err.push_back(phc_down[i].err[0]);
    
        // -- p1 -----
        p1_val.push_back(phc_up[i].par[1]);
        p1_val.push_back(phc_down[i].par[1]);
        p1_err.push_back(phc_up[i].err[1]);
        p1_err.push_back(phc_down[i].err[1]);
    
        // -- p2 -----
        p2_val.push_back(phc_up[i].par[2]);
        p2_val.push_back(phc_down[i].par[2]);
        p2_err.push_back(phc_up[i].err[2]);
        p2_err.push_back(phc_down[i].err[2]);
    
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

    conf.detector = "bh2";
    conf.phc_de_range_min = 0.5;
    analyze(path, particle);
    return 0;
}