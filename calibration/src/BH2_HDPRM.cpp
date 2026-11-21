// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <stdexcept>

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
#include "hodo_input.h"
#include "paths.h"
#include "params.h"

Config& conf = Config::getInstance();

void analyze(TString path, TString particle){
    const TString input_suffix = hodo_input::histogram_suffix(particle);

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
    conf.run_num = run_num;

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_path = Form("%s/run%05d_BH2_HDPRM_%s.root", root_base_dir.Data(), run_num, particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- tdc ----------
    TH1D *h_bh2_tdc[2][conf.num_of_ch.at("bh2")];
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++ ) h_bh2_tdc[0][i] = (TH1D*)f->Get(Form("BH2_TDC_seg%dU%s", i, input_suffix.Data()));
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++ ) h_bh2_tdc[1][i] = (TH1D*)f->Get(Form("BH2_TDC_seg%dD%s", i, input_suffix.Data()));
    // -- adc ----------
    TH1D *h_bh2_adc[2][conf.num_of_ch.at("bh2")];
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++ ) h_bh2_adc[0][i] = (TH1D*)f->Get(Form("BH2_ADC_seg%dU%s", i, input_suffix.Data()));
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++ ) h_bh2_adc[1][i] = (TH1D*)f->Get(Form("BH2_ADC_seg%dD%s", i, input_suffix.Data()));

    // Preserve the ADC input paired with this fit, before fitting rebins it.
    fout->cd();
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); ++i) {
        for (Int_t side = 0; side < 2; ++side) {
            h_bh2_adc[side][i]->Write();
        }
    }

    // -- set tdc range ----------
    TH1D *h_sum_tdc = (TH1D*)h_bh2_tdc[0][0]->Clone("h_sum_tdc");
    h_sum_tdc->Reset(); 
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++) {
        h_sum_tdc->Add(h_bh2_tdc[0][i]);
        h_sum_tdc->Add(h_bh2_tdc[1][i]);
    }
    ana_helper::set_tdc_search_range(h_sum_tdc);

    // -- set adc range ----------
    {
        Int_t seg = 6;
        auto c_dummy = new TCanvas("dummy", "", 1500, 1200);
        TString key = Form("bh2-%d-u", seg);
        conf.hdprm_mip_range_left = ana_helper::hdprm_mip_left_override(
            key.Data(),
            param::hdprm_params.count(key.Data()) ? param::hdprm_params.at(key.Data())[0] : -1.0);
        conf.hdprm_typical_value = ana_helper::adc_fit(h_bh2_adc[0][seg], c_dummy, 1);
        delete c_dummy;
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_BH2_HDPRM_%s.pdf", img_base_dir.Data(), run_num, particle.Data());

    // -- container -----
    std::vector<FitResult> adc_up;
    std::vector<FitResult> adc_down;
    ana_helper::MipFitState mip_states[2][15];
    std::vector<FitResult> tdc_up;
    std::vector<FitResult> tdc_down;

    Int_t n_rebin = 2;
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
        TString key;
        // -- UP -----
        result = ana_helper::tdc_fit(h_bh2_tdc[0][i], c_bh2, nth_pad);
        tdc_up.push_back(result);
        nth_pad++;

        key = Form("bh2-%d-u", i);
        conf.hdprm_mip_range_left = ana_helper::hdprm_mip_left_override(
            key.Data(),
            param::hdprm_params.count(key.Data()) ? param::hdprm_params.at(key.Data())[0] : -1.0);
        {
            ana_helper::MipFitSettings settings;
            bool has_settings = ana_helper::load_mip_settings(std::string(key.Data()) + ":" + particle.Data(), settings);
            result = ana_helper::adc_fit(h_bh2_adc[0][i], c_bh2, nth_pad, n_rebin,
                has_settings ? &settings : nullptr, &mip_states[0][i]);
            if (has_settings && !mip_states[0][i].success)
                throw std::runtime_error("manual MIP fit failed for " + std::string(key.Data()));
        }
        adc_up.push_back(result);
        nth_pad++;

        // -- DOWN -----
        result = ana_helper::tdc_fit(h_bh2_tdc[1][i], c_bh2, nth_pad);
        tdc_down.push_back(result);
        nth_pad++;

        key = Form("bh2-%d-d", i);
        conf.hdprm_mip_range_left = ana_helper::hdprm_mip_left_override(
            key.Data(),
            param::hdprm_params.count(key.Data()) ? param::hdprm_params.at(key.Data())[0] : -1.0);
        {
            ana_helper::MipFitSettings settings;
            bool has_settings = ana_helper::load_mip_settings(std::string(key.Data()) + ":" + particle.Data(), settings);
            result = ana_helper::adc_fit(h_bh2_adc[1][i], c_bh2, nth_pad, n_rebin,
                has_settings ? &settings : nullptr, &mip_states[1][i]);
            if (has_settings && !mip_states[1][i].success)
                throw std::runtime_error("manual MIP fit failed for " + std::string(key.Data()));
        }
        adc_down.push_back(result);
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
    std::vector<Double_t> adc_p0_val, adc_p1_val, tdc_p0_val; 
    std::vector<Double_t> adc_p0_err, adc_p1_err, tdc_p0_err; 
    std::vector<Int_t> tdc_fit_status, adc_ped_fit_status, adc_mip_fit_status;
    std::vector<Int_t> tdc_fit_valid, adc_ped_fit_valid, adc_mip_fit_valid;
    std::vector<Int_t> tdc_fit_covariance_status, adc_ped_fit_covariance_status, adc_mip_fit_covariance_status;
    std::vector<std::string> tdc_fit_method, adc_ped_fit_method, adc_mip_fit_method;
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("adc_p0_val", &adc_p0_val);
    tree->Branch("adc_p1_val", &adc_p1_val);
    tree->Branch("tdc_p0_val", &tdc_p0_val);
    tree->Branch("adc_p0_err", &adc_p0_err);
    tree->Branch("adc_p1_err", &adc_p1_err);
    tree->Branch("tdc_p0_err", &tdc_p0_err);
    tree->Branch("tdc_fit_status", &tdc_fit_status);
    tree->Branch("tdc_fit_valid", &tdc_fit_valid);
    tree->Branch("tdc_fit_method", &tdc_fit_method);
    tree->Branch("tdc_fit_covariance_status", &tdc_fit_covariance_status);
    tree->Branch("adc_ped_fit_status", &adc_ped_fit_status);
    tree->Branch("adc_ped_fit_valid", &adc_ped_fit_valid);
    tree->Branch("adc_ped_fit_method", &adc_ped_fit_method);
    tree->Branch("adc_ped_fit_covariance_status", &adc_ped_fit_covariance_status);
    tree->Branch("adc_mip_fit_status", &adc_mip_fit_status);
    tree->Branch("adc_mip_fit_valid", &adc_mip_fit_valid);
    tree->Branch("adc_mip_fit_method", &adc_mip_fit_method);
    tree->Branch("adc_mip_fit_covariance_status", &adc_mip_fit_covariance_status);
    
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++) {
        ch = i;
        adc_p0_val.clear(); adc_p1_val.clear(); tdc_p0_val.clear();
        adc_p0_err.clear(); adc_p1_err.clear(); tdc_p0_err.clear();
        tdc_fit_status.clear(); tdc_fit_valid.clear(); tdc_fit_method.clear();
        tdc_fit_covariance_status.clear();
        adc_ped_fit_status.clear(); adc_ped_fit_valid.clear(); adc_ped_fit_method.clear();
        adc_ped_fit_covariance_status.clear();
        adc_mip_fit_status.clear(); adc_mip_fit_valid.clear(); adc_mip_fit_method.clear();
        adc_mip_fit_covariance_status.clear();

        // -- pedestal -----
        adc_p0_val.push_back(adc_up[i].par[1]);
        adc_p0_val.push_back(adc_down[i].par[1]);
        adc_p0_err.push_back(adc_up[i].err[1]);
        adc_p0_err.push_back(adc_down[i].err[1]);
        for (const auto* fit : {&adc_up[i], &adc_down[i]}) {
            adc_ped_fit_status.push_back(fit->pedestal_fit_status);
            adc_ped_fit_valid.push_back(fit->pedestal_fit_valid);
            adc_ped_fit_method.push_back(fit->pedestal_fit_method);
            adc_ped_fit_covariance_status.push_back(fit->pedestal_covariance_status);
        }
    
        // -- mip -----
        adc_p1_val.push_back(adc_up[i].par[4]);
        adc_p1_val.push_back(adc_down[i].par[4]);
        adc_p1_err.push_back(adc_up[i].err[4]);
        adc_p1_err.push_back(adc_down[i].err[4]);
        for (const auto* fit : {&adc_up[i], &adc_down[i]}) {
            adc_mip_fit_status.push_back(fit->migrad_stats);
            adc_mip_fit_valid.push_back(fit->fit_valid);
            adc_mip_fit_method.push_back(fit->fit_method);
            adc_mip_fit_covariance_status.push_back(fit->covariance_status);
        }
    
        // -- tdc -----
        tdc_p0_val.push_back(tdc_up[i].par[1]);
        tdc_p0_val.push_back(tdc_down[i].par[1]);
        tdc_p0_err.push_back(tdc_up[i].err[1]);
        tdc_p0_err.push_back(tdc_down[i].err[1]);
        for (const auto* fit : {&tdc_up[i], &tdc_down[i]}) {
            tdc_fit_status.push_back(fit->migrad_stats);
            tdc_fit_valid.push_back(fit->fit_valid);
            tdc_fit_method.push_back(fit->fit_method);
            tdc_fit_covariance_status.push_back(fit->covariance_status);
        }
    
        tree->Fill();
    }
    
    fout->cd();
    tree->Write();
    TTree mip_tree("mip_fit", "MIP function parameters and final fit range");
    Int_t mip_seg, mip_side, mip_status, mip_covariance, mip_ndf;
    Bool_t mip_success;
    Double_t mip_parameters[4], mip_range[2], mip_chi_square;
    std::string mip_function;
    Double_t mip_statistics[10];
    Bool_t mip_statistics_valid;
    mip_tree.Branch("seg", &mip_seg);
    mip_tree.Branch("side", &mip_side);
    mip_tree.Branch("parameters", mip_parameters, "parameters[4]/D");
    mip_tree.Branch("range", mip_range, "range[2]/D");
    mip_tree.Branch("function", &mip_function);
    mip_tree.Branch("status", &mip_status);
    mip_tree.Branch("covariance", &mip_covariance);
    mip_tree.Branch("ndf", &mip_ndf);
    mip_tree.Branch("chi_square", &mip_chi_square);
    mip_tree.Branch("success", &mip_success);
    mip_tree.Branch("statistics", mip_statistics, "statistics[10]/D");
    mip_tree.Branch("statistics_valid", &mip_statistics_valid);
    for (mip_seg = 0; mip_seg < conf.num_of_ch.at("bh2"); ++mip_seg) {
        for (mip_side = 0; mip_side < 2; ++mip_side) {
            const auto& snapshot = mip_states[mip_side][mip_seg];
            std::copy(snapshot.parameters.begin(), snapshot.parameters.end(), mip_parameters);
            mip_range[0] = snapshot.lower;
            mip_range[1] = snapshot.upper;
            mip_function = snapshot.function;
            mip_status = snapshot.status;
            mip_covariance = snapshot.covariance;
            mip_ndf = snapshot.ndf;
            mip_chi_square = snapshot.chi_square;
            mip_success = snapshot.success;
            const auto& summary = snapshot.statistics;
            const std::array<double, 10> values{{summary.window_lower, summary.window_upper,
                summary.lower, summary.upper, summary.mean, summary.median, summary.window_median,
                summary.entries, summary.selected_entries, summary.error}};
            std::copy(values.begin(), values.end(), mip_statistics);
            mip_statistics_valid = summary.valid;
            mip_tree.Fill();
        }
    }
    mip_tree.Write();
    mip_tree.SetDirectory(nullptr);
    hodo_input::write_metadata(fout, particle);
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
    try {
        analyze(path, particle);
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
    return 0;
}
