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
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_path = Form("%s/run%05d_T1_PHC_%s.root", root_base_dir.Data(), run_num, particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- TOF vs dE ----------
    TH2D *h_t1_tof_vs_de[1][conf.num_of_ch.at("t1")];
    TH2D *h_t1_ctof_vs_de[1][conf.num_of_ch.at("t1")];
    for (Int_t i = 0; i < conf.num_of_ch.at("t1"); i++ ) {
        h_t1_tof_vs_de[0][i]  = (TH2D*)f->Get(Form("T1_seg%dU_TOF_vs_DeltaE_%s", i, particle.Data()));
        h_t1_ctof_vs_de[0][i] = (TH2D*)f->Get(Form("T1_seg%dU_CTOF_vs_DeltaE_%s", i, particle.Data()));
        // If single ended without U suffix, one might need to adjust here. But assuming U is used for single side.
    }
    
    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2; // Can adjust layout if needed
    Int_t max_pads = rows * cols;
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_T1_PHC_%s.pdf", img_base_dir.Data(), run_num, particle.Data());

    // -- container -----
    std::vector<FitResult> phc_up;
    
    auto c_t1 = new TCanvas("t1", "", 1500, 1200);
    c_t1->Divide(cols, rows);
    c_t1->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("t1"); i++) {
        if (nth_pad > max_pads) {
            c_t1->Print(pdf_path);
            c_t1->Clear();
            c_t1->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result;
        // -- UP -----
        if (h_t1_tof_vs_de[0][i]) {
            result = ana_helper::phc_fit(h_t1_tof_vs_de[0][i], c_t1, nth_pad);
            phc_up.push_back(result);
            nth_pad++;
            c_t1->cd(nth_pad);
            if(h_t1_ctof_vs_de[0][i]) h_t1_ctof_vs_de[0][i]->Draw("colz");
            nth_pad++;
        } else {
            // Handle missing histogram
            FitResult dummy; for(int k=0;k<3;k++){dummy.par.push_back(0); dummy.err.push_back(0);}
            phc_up.push_back(dummy);
            nth_pad+=2;
        }
    }
    c_t1->Print(pdf_path);
    c_t1->Print(pdf_path + "]"); // end
    delete c_t1;

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
    
    for (Int_t i = 0; i < conf.num_of_ch.at("t1"); i++) {
        ch = i;
        p0_val.clear(); p1_val.clear(); p2_val.clear();
        p0_err.clear(); p1_err.clear(); p2_err.clear();

        // -- p0 -----
        p0_val.push_back(phc_up[i].par[0]);
        p0_err.push_back(phc_up[i].err[0]);
    
        // -- p1 -----
        p1_val.push_back(phc_up[i].par[1]);
        p1_err.push_back(phc_up[i].err[1]);
    
        // -- p2 -----
        p2_val.push_back(phc_up[i].par[2]);
        p2_err.push_back(phc_up[i].err[2]);
    
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

    conf.detector = "t1";
    conf.phc_de_range_min = 0.15;
    conf.phc_de_range_max = 2.5;
    analyze(path, particle);
    return 0;
}