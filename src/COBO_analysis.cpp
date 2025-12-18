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
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TGraph.h>
#include <TPolyLine.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"

Config& conf = Config::getInstance();

void analyze(Int_t run_num){
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
    // -- hodo -----
    TString root_file_path_hodo = Form("%s/hodo_run%05d_Pi_hdphc.root", DATA_DIR.Data(), run_num);
    auto *f_hodo = new TFile( root_file_path_hodo.Data() );
    if (!f_hodo || f_hodo->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path_hodo << std::endl;
        return;
    }
    TTreeReader reader_hodo("hodo", f_hodo);
    Int_t total_entry = reader_hodo.GetEntries();
    TTreeReaderValue<unsigned int> run_number(reader_hodo, "run_number");
    TTreeReaderValue<unsigned int> evnum_hodo(reader_hodo, "event_number");

    TTreeReaderValue<std::vector<Double_t>> cobo_raw_seg(reader_hodo, "cobo_raw_seg");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> cobo_tdc(reader_hodo, "cobo_tdc_u");
    
    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    // -- COBO -----    
    TH1D *h_cobo[conf.num_of_ch.at("cobo")];
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cobo"); ch++) {
        TString name  = Form("COBO_tdc_diff_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d COBO (TDC-diff) ch%d;TDC;", run_num, ch + 1);
        h_cobo[ch] = new TH1D(name, title, 3000, 80000, 83000);
    }


    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    std::vector<std::vector<Double_t>> cobo_diff(
        conf.num_of_ch.at("cobo")
    );


    Int_t evnum = 0;
    reader_hodo.Restart();
    while (reader_hodo.Next()){ displayProgressBar(++evnum, total_entry);

        for (Int_t i = 0, n_i = (*cobo_raw_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*cobo_raw_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("cobo")) {
                for (Int_t j = 1, n_j = (*cobo_tdc)[i].size(); j < n_j; j++) {
                    Double_t diff = (*cobo_tdc)[i][j-1] - (*cobo_tdc)[i][j];
                    cobo_diff[index].push_back( diff );
                    h_cobo[index]->Fill( diff ); 
                }
            }
        }

    }

    // +------+
    // | plot |
    // +------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_COBO.pdf", OUTPUT_DIR.Data(), run_num);

    auto c_cobo = new TCanvas("", "", 1500, 1200);
    c_cobo->Divide(cols, rows);
    c_cobo->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("cobo"); i++) {
        if (nth_pad > max_pads) {
            c_cobo->Print(pdf_path);
            c_cobo->Clear();
            c_cobo->Divide(cols, rows);
            nth_pad = 1;
        }

        c_cobo->cd(nth_pad);
        h_cobo[i]->Draw();
        nth_pad++;
    }
    c_cobo->Print(pdf_path);
    c_cobo->Print(pdf_path + "]"); // end
    delete c_cobo;

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString output_path = Form("%s/root/run%05d_COBO_diff.root", OUTPUT_DIR.Data(), run_num);
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> diff;
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("diff", &diff);
    
    for (Int_t i = 0; i < conf.num_of_ch.at("cobo"); i++) {
        ch = i;
        diff = cobo_diff[i];

        tree->Fill();
    }
    
    fout->cd();
    tree->Write();
    fout->Close(); 

}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
        return 1;
    }
    Int_t run_num = std::atoi(argv[1]);

    analyze(run_num);
    return 0;
}