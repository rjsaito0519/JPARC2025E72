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
    // -- adc ----------
    TH1D *h_adc[conf.num_of_UorD.at(counter.Data())][conf.num_of_ch.at(counter.Data())];
    for (Int_t i = 0; i < conf.num_of_ch.at(counter.Data()); i++ ) {
        if (counter == "bac") {
            h_adc[0][i] = (TH1D*)f->Get(Form("%s_ADC_seg%dU_%s", counter_upper.Data(), i, particle.Data()));
        } else if (counter == "kvc") {
            h_adc[0][i] = (TH1D*)f->Get(Form("%s_ADC_seg%da_%s", counter_upper.Data(), i, particle.Data()));
            h_adc[1][i] = (TH1D*)f->Get(Form("%s_ADC_seg%db_%s", counter_upper.Data(), i, particle.Data()));
            h_adc[2][i] = (TH1D*)f->Get(Form("%s_ADC_seg%dc_%s", counter_upper.Data(), i, particle.Data()));
            h_adc[3][i] = (TH1D*)f->Get(Form("%s_ADC_seg%dd_%s", counter_upper.Data(), i, particle.Data()));
        }
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_%s_HDPRM_%s.pdf", OUTPUT_DIR.Data(), run_num, counter_upper.Data(), particle.Data());

    // -- container -----
    std::vector<std::vector<FitResult>> pedestal_cont(conf.num_of_ch.at(counter.Data()), std::vector<FitResult>());

    const std::unordered_map<std::string, std::string> suffix{
        { "bac-0", "u" },
        { "kvc-0", "a" },
        { "kvc-1", "b" },
        { "kvc-2", "c" },
        { "kvc-3", "d" },        
    };
    
    auto c_pedestal = new TCanvas("cvc", "", 1500, 1200);
    c_pedestal->Divide(cols, rows);
    c_pedestal->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at(counter.Data()); i++) {
        for (Int_t UorD = 0; UorD < conf.num_of_UorD.at(counter.Data()); UorD++) {
            if (nth_pad > max_pads) {
                c_pedestal->Print(pdf_path);
                c_pedestal->Clear();
                c_pedestal->Divide(cols, rows);
                nth_pad = 1;
            }

            FitResult result;
            TString key;

            key = Form("%s-%d-%s", counter.Data(), i, suffix.at(Form("%s-%d", counter.Data(), i)).c_str());
            conf.hdprm_pedestal_range_right = param::hdprm_params.count(key.Data()) ? param::hdprm_params.at(key.Data())[0] : 2048.0;
            result = ana_helper::pedestal_fit(h_adc[UorD][i], c_pedestal, nth_pad);
            pedestal_cont[i].push_back(result);
            nth_pad++;
        }
    }
    c_pedestal->Print(pdf_path);
    c_pedestal->Print(pdf_path + "]"); // end
    delete c_pedestal;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch, UorD;
    std::vector<Double_t> adc_p0_val;
    std::vector<Double_t> adc_p0_err; 
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("UorD", &UorD, "UorD/I");
    tree->Branch("adc_p0_val", &adc_p0_val);
    tree->Branch("adc_p0_err", &adc_p0_err);
    
    for (Int_t i = 0; i < conf.num_of_ch.at(counter.Data()); i++) {
        for (Int_t j = 0; j < conf.num_of_UorD.at(counter.Data()); j++) {
            ch = i;
            UorD = j;
            adc_p0_val.clear();
            adc_p0_err.clear();
            
            // -- pedestal -----
            adc_p0_val.push_back(pedestal_cont[i][j].par[1]);
            adc_p0_err.push_back(pedestal_cont[i][j].err[1]);
        
            tree->Fill();
        }
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