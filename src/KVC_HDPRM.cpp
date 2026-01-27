// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <unordered_map>

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

    TString counter = "kvc";
    TString counter_upper = "KVC";

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
    Int_t num_ch = 8;
    // -- adc ----------
    TH1D *h_adc[5][num_ch]; // a,b,c,d,S
    // -- tdc ----------
    TH1D *h_tdc[num_ch];   // S

    for (Int_t i = 0; i < num_ch; i++ ) {
        h_adc[0][i] = (TH1D*)f->Get(Form("%s_ADC_seg%da_%s", counter_upper.Data(), i, particle.Data()));
        h_adc[1][i] = (TH1D*)f->Get(Form("%s_ADC_seg%db_%s", counter_upper.Data(), i, particle.Data()));
        h_adc[2][i] = (TH1D*)f->Get(Form("%s_ADC_seg%dc_%s", counter_upper.Data(), i, particle.Data()));
        h_adc[3][i] = (TH1D*)f->Get(Form("%s_ADC_seg%dd_%s", counter_upper.Data(), i, particle.Data()));
        h_adc[4][i] = (TH1D*)f->Get(Form("%s_ADC_seg%dS_%s", counter_upper.Data(), i, particle.Data()));
        
        h_tdc[i]    = (TH1D*)f->Get(Form("%s_TDC_seg%dS_%s", counter_upper.Data(), i, particle.Data()));
    }

    // -- set tdc range ----------
    TH1D *h_sum_tdc = (TH1D*)h_tdc[0]->Clone("h_sum_tdc");
    h_sum_tdc->Reset(); 
    for (Int_t i = 0; i < num_ch; i++) {
        h_sum_tdc->Add(h_tdc[i]);
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

    auto c_kvc = new TCanvas("kvc", "", 1500, 1200);
    c_kvc->Divide(cols, rows);
    c_kvc->Print(pdf_path + "["); // start
    
    std::vector<FitResult> tdc_res;
    std::vector<std::vector<FitResult>> adc_res(num_ch, std::vector<FitResult>());
    
    // suffix map for keys
    const std::vector<std::string> suffix = {"a", "b", "c", "d", "S"};

    nth_pad = 1;
    for (Int_t i = 0; i < num_ch; i++) {
        if (nth_pad > max_pads) {
            c_kvc->Print(pdf_path);
            c_kvc->Clear();
            c_kvc->Divide(cols, rows);
            nth_pad = 1;
        }

        // -- TDC Fit --
        FitResult tRes = ana_helper::tdc_fit(h_tdc[i], c_kvc, nth_pad);
        tdc_res.push_back(tRes);
        nth_pad++;

        // -- ADC Fit (4 channels) --
        for (Int_t j = 0; j < 4; j++) {
            if (nth_pad > max_pads) {
                c_kvc->Print(pdf_path);
                c_kvc->Clear();
                c_kvc->Divide(cols, rows);
                nth_pad = 1;
            }
            
            TString key = Form("kvc-%d-%s", i, suffix[j].c_str());
            conf.hdprm_pedestal_range_right = param::hdprm_params.count(key.Data()) ? param::hdprm_params.at(key.Data())[0] : 2048.0;

            FitResult aRes = ana_helper::pedestal_fit(h_adc[j][i], c_kvc, nth_pad);
            adc_res[i].push_back(aRes);
            nth_pad++;
        }
    }
    c_kvc->Print(pdf_path);
    c_kvc->Print(pdf_path + "]"); // end
    delete c_kvc;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> adc_p0_val;
    std::vector<Double_t> adc_p0_err; 
    std::vector<Double_t> tdc_p0_val; 
    std::vector<Double_t> tdc_p0_err; 

    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("adc_p0_val", &adc_p0_val);
    tree->Branch("adc_p0_err", &adc_p0_err);
    tree->Branch("tdc_p0_val", &tdc_p0_val);
    tree->Branch("tdc_p0_err", &tdc_p0_err);
    
    for (Int_t i = 0; i < num_ch; i++) {
        ch = i;
        adc_p0_val.clear();
        adc_p0_err.clear();
        tdc_p0_val.clear();
        tdc_p0_err.clear();

        // -- pedestal (4 values) -----
        for(int j=0; j<4; ++j){
            adc_p0_val.push_back(adc_res[i][j].par[1]);
            adc_p0_err.push_back(adc_res[i][j].err[1]);
        }

        // -- tdc (1 value) -----
        tdc_p0_val.push_back(tdc_res[i].par[1]);
        tdc_p0_err.push_back(tdc_res[i].err[1]);

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

    conf.detector = "kvc";
    analyze(path, particle);
    return 0;
}
