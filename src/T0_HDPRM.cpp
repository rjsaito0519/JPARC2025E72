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

    // // +--------------------------+
    // // | prepare output root file |
    // // +--------------------------+
    // TString save_name;
    // Int_t dot_index = path.Last('.');
    // Int_t sla_index = path.Last('/');
    // for (Int_t i = sla_index+1; i < dot_index; i++) save_name += path[i];
    // TString output_path = Form("%s/root/acceptance_%s_%d.root", OUTPUT_DIR.Data(), save_name.Data(), focus_pdg_code);
    // if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    // TFile fout(output_path.Data(), "create");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- tdc ----------
    TH1D *h_t0_tdc[2][conf.num_of_ch.at("t0")];
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_tdc[0][i] = (TH1D*)f->Get(Form("T0_TDC_seg%dU_%s", i, particle.Data()));
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_tdc[1][i] = (TH1D*)f->Get(Form("T0_TDC_seg%dD_%s", i, particle.Data()));
    // -- adc ----------
    TH1D *h_t0_adc[2][conf.num_of_ch.at("t0")];
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_adc[0][i] = (TH1D*)f->Get(Form("T0_ADC_seg%dU_%s", i, particle.Data()));
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++ ) h_t0_adc[1][i] = (TH1D*)f->Get(Form("T0_ADC_seg%dD_%s", i, particle.Data()));

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString pdf_path = Form("%s/img/run%05d_t0_%s.pdf", OUTPUT_DIR.Data(), run_num, particle.Data());

    // -- adc up typical -----
    std::vector<std::vector<Double_t>> fit_result_ped_u;
    std::vector<std::vector<Double_t>> fit_result_mip_u;
    auto c_t0_adc_u = new TCanvas("t0_adc_u", "", 1500, 1200);
    c_t0_adc_u->Divide(cols, rows);
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("t0"); i++) {
        if (nth_pad > max_pads) {
            c_t0_adc_u->Print(pdf_path);
            c_t0_adc_u->Clear();
            c_t0_adc_u->Divide(cols, rows);
            nth_pad = 1;
        }
        // FitResult result1 = ana_helper::t0_tdc_fit(h_t0_tdc[0][i], c_t0_adc_u, nth_pad);
        // // fit_result_ped_u.push_back(par[0]);
        // // fit_result_mip_u.push_back(par[1]);
        // gPad->SetLogy(1);
        // nth_pad++;

        std::cout << i << std::endl;
        FitResult result = ana_helper::t0_adc_fit(h_t0_adc[0][i], c_t0_adc_u, nth_pad);
        // fit_result_ped_u.push_back(par[0]);
        // fit_result_mip_u.push_back(par[1]);
        gPad->SetLogy(1);
        nth_pad++;

    }
    c_t0_adc_u->Print(pdf_path);
    delete c_t0_adc_u;

    // // +-------+
    // // | Write |
    // // +-------+
    // // -- cal acceptance -----
    // h_acceptance->Divide( h_cos_theta_trig, h_cos_theta_raw, 1, 1 );
    // h_acceptance->GetYaxis()->SetRangeUser(0, 1.05);
    // h_acceptance->SetLineColor(kBlack);
    // h_acceptance->SetLineWidth(2);

    // h_acceptance_mp2->Divide( h_cos_theta_mp2, h_cos_theta_raw, 1, 1 );
    // h_acceptance_mp2->GetYaxis()->SetRangeUser(0, 1.05);
    // h_acceptance_mp2->SetLineColor(kBlue);
    // h_acceptance_mp2->SetLineWidth(2);

    // h_acceptance_htofp->Divide( h_cos_theta_htofp, h_cos_theta_raw, 1, 1 );
    // h_acceptance_htofp->GetYaxis()->SetRangeUser(0, 1.05);
    // h_acceptance_htofp->SetLineColor(kOrange);
    // h_acceptance_htofp->SetLineWidth(2);

    // // -- write -----
    // fout.cd();
    // h_cos_theta_raw->Write();
    // h_cos_theta_trig->Write();
    // h_cos_theta_mp2->Write();
    // h_cos_theta_htofp->Write();
    // h_acceptance->Write();
    // h_acceptance_mp2->Write();
    // h_acceptance_htofp->Write();
    // h_mom_dist->Write();

    // fout.Close(); 
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

    analyze(path, particle);
    return 0;
}