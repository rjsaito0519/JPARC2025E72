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
    TString output_path = Form("%s/root/run%05d_T0_Offset_%s.root", OUTPUT_DIR.Data(), run_num, particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- tdc ----------
    TH1D *h_t0_offset[conf.num_of_ch.at("t0")];
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_offset[i] = (TH1D*)f->Get(Form("T0_seg%d_TimeOffset_%s", i, particle.Data()));

    // -- set tdc range ----------
    TH1D *h_sum_tdc = (TH1D*)h_t0_offset[0]->Clone("h_sum_tdc");
    h_sum_tdc->Reset(); 
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++) h_sum_tdc->Add(h_t0_offset[i]);
    ana_helper::set_tdc_search_range(h_sum_tdc);

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_T0_Offset_%s.pdf", OUTPUT_DIR.Data(), run_num, particle.Data());

    // -- container -----
    std::vector<FitResult> t0_offset;

    auto c_t0 = new TCanvas("t0", "", 1500, 1200);
    c_t0->Divide(cols, rows);
    c_t0->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++) {
        if (nth_pad > max_pads) {
            c_t0->Print(pdf_path);
            c_t0->Clear();
            c_t0->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result = ana_helper::t0_offset_fit(h_t0_offset[i], c_t0, nth_pad);
        t0_offset.push_back(result);
        nth_pad++;
    }
    c_t0->Print(pdf_path);
    c_t0->Print(pdf_path + "]"); // end
    delete c_t0;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> offset_p0_val; 
    std::vector<Double_t> offset_p0_err; 
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("offset_p0_val", &offset_p0_val);
    tree->Branch("offset_p0_err", &offset_p0_err);
    
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++) {
        ch = i;
        offset_p0_val.clear();
        offset_p0_err.clear();
        
        offset_p0_val.push_back(t0_offset[i].par[1]);
        offset_p0_err.push_back(t0_offset[i].err[1]);
    
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

    conf.detector = "t0";
    analyze(path, particle);
    return 0;
}