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

// HTOF HDPRM: TDC-only using All-event histograms (no beam-particle suffix).
void analyze(TString path){
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0/256, 160.0/256, 44.0/256);
    
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x");
    gStyle->SetTitleSize(0.06, "y");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    const TString tag = "All";

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
    TString output_path = Form("%s/run%05d_HTOF_HDPRM_%s.root", root_base_dir.Data(), run_num, tag.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    // All-event histos: HTOF_TDC_seg%dU (no _Pi/_K suffix)
    TH1D *h_htof_tdc[3][conf.num_of_ch.at("htof")];
    for (Int_t i = 0; i < conf.num_of_ch.at("htof"); i++ ) {
        h_htof_tdc[0][i] = (TH1D*)f->Get(Form("HTOF_TDC_seg%dU", i));
        h_htof_tdc[1][i] = (TH1D*)f->Get(Form("HTOF_TDC_seg%dD", i));
        h_htof_tdc[2][i] = (TH1D*)f->Get(Form("HTOF_TDC_seg%dS", i));
        if (!h_htof_tdc[0][i] || !h_htof_tdc[1][i] || !h_htof_tdc[2][i]) {
            std::cerr << "Error: Missing HTOF_TDC All-event hist for seg " << i << std::endl;
            fout->Close();
            return;
        }
    }

    // -- set tdc range ----------
    TH1D *h_sum_tdc[2];
    h_sum_tdc[0] = (TH1D*)h_htof_tdc[0][0]->Clone("h_sum_tdc_indiv");
    h_sum_tdc[0]->Reset();
    h_sum_tdc[1] = (TH1D*)h_htof_tdc[2][0]->Clone("h_sum_tdc_sum");
    h_sum_tdc[1]->Reset();
    for (Int_t i = 0; i < conf.num_of_ch.at("htof"); i++) {
        h_sum_tdc[0]->Add(h_htof_tdc[0][i]);
        h_sum_tdc[0]->Add(h_htof_tdc[1][i]);
        h_sum_tdc[1]->Add(h_htof_tdc[2][i]);
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    Int_t nth_pad = 1;
    Int_t rows = 3, cols = 3;
    Int_t max_pads = rows * cols;
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_HTOF_HDPRM_%s.pdf", img_base_dir.Data(), run_num, tag.Data());

    std::vector<FitResult> tdc_up;
    std::vector<FitResult> tdc_down;
    std::vector<FitResult> tdc_sum;

    auto c_htof = new TCanvas("htof", "", 1500, 1200);
    c_htof->Divide(cols, rows);
    c_htof->Print(pdf_path + "[");
    nth_pad = 1;
    for (Int_t i = conf.htof_adc_exist_seg.first; i < conf.htof_adc_exist_seg.second+1; i++) {
        if (nth_pad > max_pads) {
            c_htof->Print(pdf_path);
            c_htof->Clear();
            c_htof->Divide(cols, rows);
            nth_pad = 1;
        }

        FitResult result;
        ana_helper::set_tdc_search_range(h_sum_tdc[0]);
        // -- UP -----
        result = ana_helper::tdc_fit(h_htof_tdc[0][i], c_htof, nth_pad);
        tdc_up.push_back(result);
        nth_pad++;

        // -- DOWN -----
        result = ana_helper::tdc_fit(h_htof_tdc[1][i], c_htof, nth_pad);
        tdc_down.push_back(result);
        nth_pad++;

        ana_helper::set_tdc_search_range(h_sum_tdc[1]);
        // -- SUM -----
        result = ana_helper::tdc_fit(h_htof_tdc[2][i], c_htof, nth_pad);
        tdc_sum.push_back(result);
        nth_pad++;
    }
    c_htof->Print(pdf_path);
    c_htof->Print(pdf_path + "]");
    delete c_htof;

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

    for (Int_t i = conf.htof_adc_exist_seg.first; i < conf.htof_adc_exist_seg.second+1; i++) {
        ch = i;
        tdc_p0_val.clear();
        tdc_p0_err.clear();

        const Int_t idx = i - conf.htof_adc_exist_seg.first;
        tdc_p0_val.push_back(tdc_up[idx].par[1]);
        tdc_p0_val.push_back(tdc_down[idx].par[1]);
        tdc_p0_val.push_back(tdc_sum[idx].par[1]);
        tdc_p0_err.push_back(tdc_up[idx].err[1]);
        tdc_p0_err.push_back(tdc_down[idx].err[1]);
        tdc_p0_err.push_back(tdc_sum[idx].err[1]);

        tree->Fill();
    }

    fout->cd();
    tree->Write();
    fout->Close();
}

Int_t main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <root file path> <particle>" << std::endl;
        std::cerr << "  particle: All (recommended), or Pi/K (forced to All)" << std::endl;
        return 1;
    }
    TString path = argv[1];
    TString particle = argv[2];
    if (particle != "All" && particle != "Pi" && particle != "K") {
        std::cerr << "Error: Unexpected particle name: " << particle
                  << " (expected All, Pi, or K)" << std::endl;
        return 1;
    }
    if (particle != "All") {
        std::cerr << "Warning: HTOF_HDPRM ignores particle='" << particle
                  << "'; always uses All-event TDC histograms" << std::endl;
    }

    conf.detector = "htof";
    conf.hdprm_mip_half_width_ratio = 0.1;
    conf.adc_ped_remove_nsigma = 15.0;
    analyze(path);
    return 0;
}
