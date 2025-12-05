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
#include "params.h"

Config& conf = Config::getInstance();

void analyze(TString path, TString counter, TString particle){    
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

    TString counter_upper = counter;
    counter_upper.ToUpper();

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
    TString output_path = Form("%s/root/run%05d_%s_HDPRM_%s.root", OUTPUT_DIR.Data(), run_num, counter_upper.Data(), particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- tdc ----------
    TH1D *h_tdc[2][conf.num_of_ch.at(counter.Data())];
    for (Int_t i = 0; i < conf.num_of_ch.at(counter.Data()); i++ ) {
        if (counter == "bac") {
            h_tdc[0][i] = (TH1D*)f->Get(Form("%s_TDC_seg4U_%s", counter_upper.Data(), particle.Data()));
        } else if (counter == "kvc") {
            h_tdc[0][i] = (TH1D*)f->Get(Form("%s_TDC_seg%dS_%s", counter_upper.Data(), i, particle.Data()));
        } else {
            h_tdc[0][i] = (TH1D*)f->Get(Form("%s_TDC_seg%dU_%s", counter_upper.Data(), i, particle.Data()));
            h_tdc[1][i] = (TH1D*)f->Get(Form("%s_TDC_seg%dD_%s", counter_upper.Data(), i, particle.Data()));
        }
    }

    // -- set tdc range ----------
    TH1D *h_sum_tdc = (TH1D*)h_tdc[0][0]->Clone("h_sum_tdc");
    h_sum_tdc->Reset(); 
    for (Int_t i = 0; i < conf.num_of_ch.at(counter.Data()); i++) {
        h_sum_tdc->Add(h_tdc[0][i]);
        if (counter != "bac" && counter != "kvc") h_sum_tdc->Add(h_tdc[1][i]);
    }
    ana_helper::set_tdc_search_range(h_sum_tdc);

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_%s_HDPRM_%s.pdf", OUTPUT_DIR.Data(), run_num, counter_upper.Data(), particle.Data());

    // -- container -----
    std::vector<FitResult> tdc_up;
    std::vector<FitResult> tdc_down;

    auto c_tdc = new TCanvas(counter.Data(), "", 1500, 1200);
    c_tdc->Divide(cols, rows);
    c_tdc->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at(counter.Data()); i++) {
        if (nth_pad > max_pads) {
            c_tdc->Print(pdf_path);
            c_tdc->Clear();
            c_tdc->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result;
        TString key;
        // -- UP -----
        result = ana_helper::tdc_fit(h_tdc[0][i], c_tdc, nth_pad);
        tdc_up.push_back(result);
        nth_pad++;

        if (counter != "bac" && counter != "kvc" && counter != "t1" && counter != "sac3" && counter != "sfv") {
            // -- DOWN -----
            result = ana_helper::tdc_fit(h_tdc[1][i], c_tdc, nth_pad);
            tdc_down.push_back(result);
            nth_pad++;
        }
    }
    c_tdc->Print(pdf_path);
    c_tdc->Print(pdf_path + "]"); // end
    delete c_tdc;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> tdc_p0_val; 
    std::vector<Double_t> tdc_p0_err; 
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("tdc_p0_val", &tdc_p0_val);
    tree->Branch("tdc_p0_err", &tdc_p0_err);
    
    for (Int_t i = 0; i < conf.num_of_ch.at(counter.Data()); i++) {
        ch = i;
        tdc_p0_val.clear();
        tdc_p0_err.clear();
        
        // -- tdc -----
        tdc_p0_val.push_back(tdc_up[i].par[1]);
        tdc_p0_err.push_back(tdc_up[i].err[1]);
        if (counter != "bac" && counter != "kvc" && counter != "t1" && counter != "sac3" && counter != "sfv") {
            tdc_p0_val.push_back(tdc_down[i].par[1]);
            tdc_p0_err.push_back(tdc_down[i].err[1]);
        }
    
        tree->Fill();
    }
    
    fout->cd();
    tree->Write();
    fout->Close(); 
}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <root file path> <counter> <particle>" << std::endl;
        return 1;
    }
    TString path = argv[1];
    TString counter = argv[2];
    TString particle = argv[3];
    if (particle != "Pi" && particle != "K") {
        std::cerr << "Error: Unexpected particle name: " << particle << std::endl;
        return 1;
    }

    conf.detector = counter;
    analyze(path, counter, particle);
    return 0;
}