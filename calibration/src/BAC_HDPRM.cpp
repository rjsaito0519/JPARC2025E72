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

    TString counter = "bac";
    TString counter_upper = "BAC";

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
    conf.run_num = run_num;

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_suffix = (particle == "") ? "" : Form("_%s", particle.Data());
    TString output_path = Form("%s/run%05d_%s_HDPRM%s.root", root_base_dir.Data(), run_num, counter_upper.Data(), output_suffix.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // For Cherenkov, we prefer all-event histograms (no suffix) for pedestals
    auto get_h = [&](TFile* f, TString base_name, TString part) -> TH1D* {
        TH1D* h = (TH1D*)f->Get(base_name); // Try no-suffix first
        if (!h && part != "") h = (TH1D*)f->Get(base_name + "_" + part);
        return h;
    };

    // -- adc ----------
    TH1D *h_adc[5]; // BAC seg0-3 U + Sum (seg4)
    for (Int_t i = 0; i < 5; i++ ) {
        h_adc[i] = get_h(f, Form("%s_ADC_seg%dU", counter_upper.Data(), i), particle);
    }
    // -- tdc ----------
    // BAC has only seg4U TDC (Sum)
    TH1D *h_tdc_sum = get_h(f, Form("%s_TDC_seg4U", counter_upper.Data()), particle);

    // -- set tdc range ----------
    TH1D *h_sum_tdc_calc = (TH1D*)h_tdc_sum->Clone("h_sum_tdc_calc");
    ana_helper::set_tdc_search_range(h_sum_tdc_calc);

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 3;
    Int_t max_pads = rows * cols;
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_%s_HDPRM%s.pdf", img_base_dir.Data(), run_num, counter_upper.Data(), output_suffix.Data());

    auto c_bac = new TCanvas("bac", "", 1500, 1200);
    c_bac->Divide(cols, rows);
    c_bac->Print(pdf_path + "["); // start
    
    // -- TDC fit (once) --
    FitResult tdc_res;
    c_bac->Clear();
    c_bac->Divide(cols, rows);
    tdc_res = ana_helper::tdc_fit(h_tdc_sum, c_bac, 1);
    c_bac->Print(pdf_path); // Print TDC page

    // -- ADC fit (loop) --
    std::vector<FitResult> adc_res;
    nth_pad = 1;
    c_bac->Clear();
    c_bac->Divide(cols, rows);

    for (Int_t i = 0; i < 5; i++) {
        if (nth_pad > max_pads) {
            c_bac->Print(pdf_path);
            c_bac->Clear();
            c_bac->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result;
        TString key = (i < 4) ? Form("bac-%d-u", i) : "bac-sum";
        conf.hdprm_pedestal_range_right = param::hdprm_params.count(key.Data()) ? param::hdprm_params.at(key.Data())[0] : 120.0;
        
        result = ana_helper::pedestal_fit(h_adc[i], c_bac, nth_pad);
        adc_res.push_back(result);
        nth_pad++;
    }
    c_bac->Print(pdf_path);
    c_bac->Print(pdf_path + "]"); // end
    delete c_bac;

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
    
    for (Int_t i = 0; i < 5; i++) {
        ch = i;
        adc_p0_val.clear();
        adc_p0_err.clear();
        tdc_p0_val.clear();
        tdc_p0_err.clear();

        // -- pedestal -----
        adc_p0_val.push_back(adc_res[i].par[1]);
        adc_p0_err.push_back(adc_res[i].err[1]);

        // -- tdc (common) -----
        // BAC TDC is essentially the sum (seg 4). 
        // We fill it for all segments including segment 4.
        tdc_p0_val.push_back(tdc_res.par[1]);
        tdc_p0_err.push_back(tdc_res.err[1]);

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

    conf.detector = "bac";
    analyze(path, particle);
    return 0;
}
