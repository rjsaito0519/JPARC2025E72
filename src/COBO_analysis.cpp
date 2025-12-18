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
    TH1D *h_cobo_time[conf.num_of_ch.at("cobo")];
    TH1D *h_cobo_tdc[conf.num_of_ch.at("cobo")];
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cobo"); ch++) {
        {
            TString name  = Form("COBO_time_diff_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d COBO (diff) ch%d;Time [ns];", run_num, ch + 1);
            h_cobo_time[ch] = new TH1D(name, title, 100, 79.5, 80.5);
        }
        {
            TString name  = Form("COBO_tdc_diff_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d COBO (diff) ch%d;TDC [ch];", run_num, ch + 1);
            h_cobo_tdc[ch] = new TH1D(name, title, 1000, 81500.0, 82500.0);
        }
    }
    auto *h_cobo_1st_hit = new TH1D(
        Form("COBO_1st_hit_%d", run_num),
        Form("run%05d COBO 1st hit;TDC [ch];", run_num),
        1000, 1.E6, 2.E6
    );


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
                if ((*cobo_tdc)[i][0]) h_cobo_1st_hit->Fill((*cobo_tdc)[i][0]);
                for (Int_t j = 1, n_j = (*cobo_tdc)[i].size(); j < n_j; j++) {
                    Double_t diff = (*cobo_tdc)[i][j-1] - (*cobo_tdc)[i][j];
                    cobo_diff[index].push_back( diff );
                    h_cobo_tdc[index]->Fill( diff );
                    h_cobo_time[index]->Fill( diff*9.765625e-04 );
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

        c_cobo->cd(nth_pad)->SetLogy(1);
        h_cobo_tdc[i]->Draw();
        nth_pad++;

        c_cobo->cd(nth_pad)->SetLogy(1);
        h_cobo_time[i]->Draw();
        nth_pad++;
    }
    if (nth_pad > max_pads) {
        c_cobo->Print(pdf_path);
        c_cobo->Clear();
        c_cobo->Divide(cols, rows);
        nth_pad = 1;
    }
    c_cobo->cd(nth_pad)->SetLogy(1);
    h_cobo_1st_hit->Draw();

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

        h_cobo_tdc[i]->Write();
        h_cobo_time[i]->Write();        
    }
    h_cobo_1st_hit->Write();
    
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