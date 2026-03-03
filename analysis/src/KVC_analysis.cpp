// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <cstring>

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
#include <TLatex.h>
#include <TF1.h>

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
    gStyle->SetTitleSize(0.06, "x"); 
    gStyle->SetTitleSize(0.06, "y"); 
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // +-----------+
    // | load file |
    // +-----------+
    TString root_file_path_hodo = Form("%s/run%05d_Hodo.root", DATA_DIR.Data(), run_num);
    auto *f_hodo = new TFile( root_file_path_hodo.Data() );
    if (!f_hodo || f_hodo->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path_hodo << std::endl;
        return;
    }
    TTreeReader reader_hodo("hodo", f_hodo);
    Int_t total_entry = reader_hodo.GetEntries();
    
    TTreeReaderValue<std::vector<Double_t>> kvc_hit_seg(reader_hodo, "kvc_hit_seg");
    TTreeReaderValue<std::vector<Double_t>> t1_hit_seg(reader_hodo, "t1_hit_seg");
    TTreeReaderValue<Double_t> btof(reader_hodo, "btof0");

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_path = Form("%s/run%05d_kvc.root", root_base_dir.Data(), run_num);
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    Int_t n_kvc = 8;
    if (conf.num_of_ch.count("kvc")) n_kvc = conf.num_of_ch.at("kvc");
    
    TH1D *h_kvc_npe[2][n_kvc]; 
    TH1D *h_kvc_npe_raw[2][n_kvc];
    TH1D *h_kvc_npe_all[2]; 
    TH1D *h_kvc_npe_all_raw[2];

    for (Int_t p = 0; p < 2; ++p) {
        TString p_str = (p == 0) ? "Pi" : "K";
        h_kvc_npe_all[p] = new TH1D(Form("KVC_Hit_Npe_all_offline_T1_%s", p_str.Data()), Form("KVC NPE all T1 %s", p_str.Data()), 500, -50, 450);
        h_kvc_npe_all_raw[p] = new TH1D(Form("KVC_Hit_Npe_all_offline_raw_%s", p_str.Data()), Form("KVC NPE all raw %s", p_str.Data()), 500, -50, 450);
        
        for (Int_t i = 0; i < n_kvc; ++i ) {
            h_kvc_npe[p][i] = (TH1D*)f_hodo->Get(Form("KVC_Hit_Npe_seg%dS_offline_T1_%s", i, p_str.Data()));
            h_kvc_npe_raw[p][i] = (TH1D*)f_hodo->Get(Form("KVC_Hit_Npe_seg%dS_offline_%s", i, p_str.Data()));
            
            if (h_kvc_npe[p][i]) {
                h_kvc_npe[p][i]->SetLineColor(kRed);
                h_kvc_npe_all[p]->Add(h_kvc_npe[p][i]);
            }
            if (h_kvc_npe_raw[p][i]) {
                h_kvc_npe_raw[p][i]->SetLineColor(kBlue);
                h_kvc_npe_all_raw[p]->Add(h_kvc_npe_raw[p][i]);
            }
        }
        h_kvc_npe_all[p]->SetLineColor(kRed);
        h_kvc_npe_all_raw[p]->SetLineColor(kBlue);
    }

    TH1D *h_btof[2]; 
    h_btof[0] = (TH1D*)f_hodo->Get("CBtof0_Pi");
    h_btof[1] = (TH1D*)f_hodo->Get("CBtof0_K");
    
    Double_t btof_threshold = -15.0;
    if (h_btof[0] && h_btof[1]) {
        btof_threshold = (h_btof[0]->GetBinCenter(h_btof[0]->GetMaximumBin()) + h_btof[1]->GetBinCenter(h_btof[1]->GetMaximumBin()))/2.0;
    }

    // +------------+
    // | efficiency |
    // +------------+
    std::vector<Int_t> n_hit(2, 0);
    std::vector<Int_t> n_hit_with_kvc(2, 0);
    
    Int_t evnum = 0;
    reader_hodo.Restart();
    while (reader_hodo.Next()){ displayProgressBar(++evnum, total_entry);
        Int_t particle = (*btof < btof_threshold) ? 1 : 0;
        if (t1_hit_seg->size() > 0) {
            n_hit[particle]++;
            if (kvc_hit_seg->size() > 0) n_hit_with_kvc[particle]++;
        }
    }

    std::cout << "\n--- Efficiency Results ---" << std::endl;
    for (Int_t p = 0; p < 2; ++p) {
        if (n_hit[p] > 0) {
            Double_t eff = static_cast<Double_t>(n_hit_with_kvc[p]) / n_hit[p];
            std::cout << Form("  %s Efficiency: %6.2f %% (%d/%d)", (p==0 ? "Pi" : "K"), eff*100.0, n_hit_with_kvc[p], n_hit[p]) << std::endl;
        }
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    TString img_data_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_KVC_NPE.pdf", img_data_dir.Data(), run_num);
    TCanvas *c = new TCanvas("c_npe", "", 1500, 1050);
    c->Divide(3, 3);
    c->Print(pdf_path + "[");

    fout->cd();
    TTree* tree = new TTree("tree", "KVC results");
    Int_t pid, is_t1, total_count, hit_count;
    Double_t constant, mean, sigma, eff;
    tree->Branch("pid", &pid, "pid/I"); 
    tree->Branch("is_t1", &is_t1, "is_t1/I"); // 0: raw, 1: t1
    tree->Branch("total_count", &total_count, "total_count/I");
    tree->Branch("hit_count", &hit_count, "hit_count/I");
    tree->Branch("constant", &constant, "constant/D");
    tree->Branch("mean", &mean, "mean/D");
    tree->Branch("sigma", &sigma, "sigma/D");
    tree->Branch("eff", &eff, "eff/D");

    for (Int_t p = 0; p < 2; ++p) {
        c->Clear();
        c->Divide(3, 3);
        for (Int_t ch = 0; ch < n_kvc; ++ch) {
            c->cd(ch+1);
            if (h_kvc_npe_raw[p][ch]) h_kvc_npe_raw[p][ch]->Draw();
            if (h_kvc_npe[p][ch]) h_kvc_npe[p][ch]->Draw("same");
        }
        c->cd(9);
        if (h_kvc_npe_all_raw[p]) {
            h_kvc_npe_all_raw[p]->Draw();
            if (h_kvc_npe_all[p]) h_kvc_npe_all[p]->Draw("same");

            // Raw fit
            h_kvc_npe_all_raw[p]->GetXaxis()->SetRangeUser(50, 400);
            Double_t m0 = h_kvc_npe_all_raw[p]->GetBinCenter(h_kvc_npe_all_raw[p]->GetMaximumBin());
            h_kvc_npe_all_raw[p]->GetXaxis()->SetRange();
            h_kvc_npe_all_raw[p]->Fit("gaus", "Q", "", m0 - 50, m0 + 50);
            TF1* f0 = h_kvc_npe_all_raw[p]->GetFunction("gaus");
            if (f0) {
                f0->SetLineColor(kBlue);
                pid = p; is_t1 = 0;
                total_count = n_hit[p]; hit_count = n_hit_with_kvc[p];
                eff = (total_count > 0) ? (Double_t)hit_count / total_count : 0.0;
                constant = f0->GetParameter(0); mean = f0->GetParameter(1); sigma = f0->GetParameter(2);
                tree->Fill();
            }
        }
        if (h_kvc_npe_all[p]) {
            // T1 fit
            h_kvc_npe_all[p]->GetXaxis()->SetRangeUser(50, 400);
            Double_t m1 = h_kvc_npe_all[p]->GetBinCenter(h_kvc_npe_all[p]->GetMaximumBin());
            h_kvc_npe_all[p]->GetXaxis()->SetRange();
            h_kvc_npe_all[p]->Fit("gaus", "Q", "", m1 - 50, m1 + 50);
            TF1* f1 = h_kvc_npe_all[p]->GetFunction("gaus");
            if (f1) {
                f1->SetLineColor(kRed);
                pid = p; is_t1 = 1;
                total_count = n_hit[p]; hit_count = n_hit_with_kvc[p];
                eff = (total_count > 0) ? (Double_t)hit_count / total_count : 0.0;
                constant = f1->GetParameter(0); mean = f1->GetParameter(1); sigma = f1->GetParameter(2);
                tree->Fill();
            }
        }
        
        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.06);
        if (h_kvc_npe_all_raw[p] && h_kvc_npe_all_raw[p]->GetFunction("gaus")) {
            tex.SetTextColor(kBlue); tex.DrawLatex(0.40, 0.85, Form("Raw Mean: %.1f", h_kvc_npe_all_raw[p]->GetFunction("gaus")->GetParameter(1)));
        }
        if (h_kvc_npe_all[p] && h_kvc_npe_all[p]->GetFunction("gaus")) {
            tex.SetTextColor(kRed); tex.DrawLatex(0.40, 0.78, Form("T1  Mean: %.1f", h_kvc_npe_all[p]->GetFunction("gaus")->GetParameter(1)));
        }
        c->Print(pdf_path);
    }
    c->Print(pdf_path + "]");

    // +-------+
    // | Write |
    // +-------+
    fout->cd();
    for (Int_t p = 0; p < 2; ++p) {
        if (h_kvc_npe_all[p]) h_kvc_npe_all[p]->Write();
        if (h_kvc_npe_all_raw[p]) h_kvc_npe_all_raw[p]->Write();
        for (Int_t i = 0; i < n_kvc; ++i) {
            if (h_kvc_npe[p][i]) h_kvc_npe[p][i]->Write();
            if (h_kvc_npe_raw[p][i]) h_kvc_npe_raw[p][i]->Write();
        }
    }
    tree->Write();
    fout->Close();
    delete c;
}

Int_t main(int argc, char** argv) {
    if (argc < 2) return 1;
    analyze(std::atoi(argv[1]));
    return 0;
}