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
#include <TMath.h>
#include <TLine.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "params.h"
#include "progress_bar.h"

Config& conf = Config::getInstance();
static const Double_t kvc_tdc_min = 680000.0;
static const Double_t kvc_tdc_max = 707470.0;
static const Double_t kvc_ratio_closeup_xmax = 45.;
static const Double_t kvc_ped_fit_min = -5.;
static const Double_t kvc_ped_fit_max = 30.;

// KVC gain (p1) per segment and channel [seg][a,b,c,d]. From HodoParam_run02140_Pi (CId=6, AorT=0, UorD=0,1,2,3).
static const std::vector<std::vector<Double_t>> kvc_gain = {
    {9.405720e+00, 8.864464e+00, 9.281763e+00, 9.215651e+00},  // seg 0
    {9.293611e+00, 1.002469e+01, 9.597398e+00, 9.984777e+00},  // seg 1
    {9.794817e+00, 9.032530e+00, 9.604148e+00, 9.305942e+00},  // seg 2
    {8.935510e+00, 8.566897e+00, 9.156487e+00, 9.693401e+00},  // seg 3
    {9.405720e+00, 8.864464e+00, 9.281763e+00, 9.215651e+00},  // seg 4
    {9.293611e+00, 1.002469e+01, 9.597398e+00, 9.984777e+00},  // seg 5
    {9.794817e+00, 9.032530e+00, 9.604148e+00, 9.305942e+00},  // seg 6
    {8.935510e+00, 8.566897e+00, 9.156487e+00, 9.693401e+00},  // seg 7
};

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
    TTreeReaderValue<std::vector<Double_t>> kvc_adc_a(reader_hodo, "kvc_adc_a");
    TTreeReaderValue<std::vector<Double_t>> kvc_adc_b(reader_hodo, "kvc_adc_b");
    TTreeReaderValue<std::vector<Double_t>> kvc_adc_c(reader_hodo, "kvc_adc_c");
    TTreeReaderValue<std::vector<Double_t>> kvc_adc_d(reader_hodo, "kvc_adc_d");    
    std::vector<TTreeReaderValue<std::vector<Double_t>>*> kvc_adc_vec = {
        &kvc_adc_a,
        &kvc_adc_b,
        &kvc_adc_c,
        &kvc_adc_d
    };
    TTreeReaderValue<std::vector<std::vector<Double_t>>> kvc_tdc(reader_hodo, "kvc_tdc_s");

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_path = Form("%s/run%05d_kvc_threshold.root", root_base_dir.Data(), run_num);
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    const Int_t num_seg = conf.num_of_ch.at("kvc");
    TH1D *h_nwt[num_seg];
    TH1D *h_npe_s[num_seg];
    for (Int_t seg = 0; seg < num_seg; ++seg) {
        h_npe_s[seg] = new TH1D(Form("h_npe_s_seg%d", seg), ";npe_{sum};count", conf.npe_bin_num, conf.npe_min, conf.npe_max);
        h_nwt[seg]   = new TH1D(Form("h_nwt_seg%d", seg), ";npe_{sum};count", conf.npe_bin_num, conf.npe_min, conf.npe_max);
    }

    // -- for pedestal fit ----------
    TH1D *h_indiv_adc[4][conf.num_of_ch.at("kvc")];
    for (Int_t i = 0; i < conf.num_of_ch.at("kvc"); ++i) {
        h_indiv_adc[0][i] = (TH1D*)f_hodo->Get(Form("KVC_ADC_seg%da", i));
        h_indiv_adc[1][i] = (TH1D*)f_hodo->Get(Form("KVC_ADC_seg%db", i));
        h_indiv_adc[2][i] = (TH1D*)f_hodo->Get(Form("KVC_ADC_seg%dc", i));
        h_indiv_adc[3][i] = (TH1D*)f_hodo->Get(Form("KVC_ADC_seg%dd", i));
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    TString img_data_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_KVC_threshold.pdf", img_data_dir.Data(), run_num);
    TCanvas *c = new TCanvas("c_npe", "", 1500, 1050);
    const Int_t rows = 2, cols = 2;
    const Int_t max_pads = rows * cols;
    c->Divide(cols, rows);
    c->Print(pdf_path + "["); // start

    auto maybe_new_page = [&](Int_t &pad) {
        if (pad > max_pads) {
            c->Print(pdf_path);
            c->Clear();
            c->Divide(cols, rows);
            pad = 1;
        }
    };

    Int_t nth_pad = 1;

    // suffix map for keys
    const std::vector<std::string> suffix = {"a", "b", "c", "d"};
    std::vector<std::vector<FitResult>> pedestal_results(num_seg);

    for (Int_t seg = 0; seg < conf.num_of_ch.at("kvc"); ++seg) {
        maybe_new_page(nth_pad);

        // -- ADC Fit (4 channels) --
        std::vector<FitResult> indiv_result;
        for (Int_t j = 0; j < 4; j++) {
            maybe_new_page(nth_pad);
            TString key = Form("kvc-%d-%s", seg, suffix[j].c_str());
            conf.hdprm_pedestal_range_right = param::hdprm_params.count(key.Data()) ? param::hdprm_params.at(key.Data())[0] : 2048.0;

            FitResult tmp_result = ana_helper::pedestal_fit(h_indiv_adc[j][seg], c, nth_pad);
            indiv_result.push_back(tmp_result);
            nth_pad++;
        }
        pedestal_results[seg] = indiv_result;
    }
    c->Print(pdf_path);
    c->Clear();
    c->Divide(cols, rows);
    nth_pad = 1;


    // +------------------------------------------+
    // | extract pedestals and fill ratio vs npe  |
    // +------------------------------------------+
    std::vector<std::vector<Double_t>> pedestal(num_seg, std::vector<Double_t>(4));
    for (Int_t seg = 0; seg < num_seg; ++seg)
        for (Int_t ch = 0; ch < 4; ++ch)
            pedestal[seg][ch] = pedestal_results[seg][ch].par[1];

    TH1D *h_ratio_npe[num_seg];
    std::vector<Double_t> threshold_at_05(num_seg, -1.0);

    reader_hodo.Restart();
    while (reader_hodo.Next()) {
        for (Int_t seg = 0; seg < num_seg; ++seg) {
            if (seg >= (Int_t)(*kvc_adc_a).size()) continue;
            Double_t adc_a = (*kvc_adc_a)[seg];
            Double_t adc_b = (*kvc_adc_b)[seg];
            Double_t adc_c = (*kvc_adc_c)[seg];
            Double_t adc_d = (*kvc_adc_d)[seg];
            Double_t g0 = kvc_gain[seg][0] > 1e-9 ? kvc_gain[seg][0] : 1.0;
            Double_t g1 = kvc_gain[seg][1] > 1e-9 ? kvc_gain[seg][1] : 1.0;
            Double_t g2 = kvc_gain[seg][2] > 1e-9 ? kvc_gain[seg][2] : 1.0;
            Double_t g3 = kvc_gain[seg][3] > 1e-9 ? kvc_gain[seg][3] : 1.0;
            Double_t npe_a = (adc_a - pedestal[seg][0]) / g0;
            Double_t npe_b = (adc_b - pedestal[seg][1]) / g1;
            Double_t npe_c = (adc_c - pedestal[seg][2]) / g2;
            Double_t npe_d = (adc_d - pedestal[seg][3]) / g3;
            Double_t npe_s = npe_a + npe_b + npe_c + npe_d;

            h_npe_s[seg]->Fill(npe_s);

            Bool_t nwt = kFALSE;
            if (seg < (Int_t)(*kvc_tdc).size()) {
                const std::vector<Double_t>& tdc_vec = (*kvc_tdc)[seg];
                for (size_t h = 0; h < tdc_vec.size(); ++h) {
                    if (tdc_vec[h] >= kvc_tdc_min && tdc_vec[h] <= kvc_tdc_max) {
                        nwt = kTRUE;
                        break;
                    }
                }
            }
            if (nwt) h_nwt[seg]->Fill(npe_s);
        }
    }

    for (Int_t seg = 0; seg < num_seg; ++seg) {
        h_ratio_npe[seg] = (TH1D*)h_nwt[seg]->Clone(Form("h_ratio_npe_seg%d", seg));
        h_ratio_npe[seg]->Divide(h_nwt[seg], h_npe_s[seg], 1., 1., "B");
        h_ratio_npe[seg]->SetTitle(";npe_{sum};ratio (nwt/all)");
    }

    for (Int_t seg = 0; seg < num_seg; ++seg) {
        TF1 *f_erf = new TF1("f_erf", "0.5*(1+TMath::Erf((x-[0])/[1]))", 0., conf.npe_max);
        f_erf->SetParameters(20., 5.);
        Int_t fit_status = h_ratio_npe[seg]->Fit(f_erf, "Q0", "", 0., 150.);
        if (fit_status == 0) threshold_at_05[seg] = f_erf->GetParameter(0);
    }

    // +----------------------------------------+
    // | pedestal Gauss fit -> n_sigma from 0.5 |
    // +----------------------------------------+
    std::vector<Double_t> ped_mean(num_seg, 0.), ped_sigma(num_seg, -1.), ped_amplitude(num_seg, 0.);
    std::vector<Double_t> n_sigma(num_seg, -1.0);
    for (Int_t seg = 0; seg < num_seg; ++seg) {
        TF1 *f_ped = new TF1("f_ped", "gaus", kvc_ped_fit_min, kvc_ped_fit_max);
        Double_t maxbin = h_npe_s[seg]->GetBinCenter(h_npe_s[seg]->GetMaximumBin());
        f_ped->SetParameters(h_npe_s[seg]->GetMaximum(), maxbin, 5.);
        Int_t ped_status = h_npe_s[seg]->Fit(f_ped, "Q0", "", kvc_ped_fit_min, kvc_ped_fit_max);
        if (ped_status == 0) {
            ped_amplitude[seg] = f_ped->GetParameter(0);
            ped_mean[seg] = f_ped->GetParameter(1);
            ped_sigma[seg] = f_ped->GetParameter(2);
            if (ped_sigma[seg] > 1e-9 && threshold_at_05[seg] > -0.5)
                n_sigma[seg] = (threshold_at_05[seg] - ped_mean[seg]) / ped_sigma[seg];
        }
    }

    // +---------------------------+
    // | draw npe_s+nwt, ratio     |
    // +---------------------------+
    nth_pad = 1;
    for (Int_t seg = 0; seg < num_seg; ++seg) {
        maybe_new_page(nth_pad);
        c->cd(nth_pad)->SetLogy(1);
        TH1D *h_npe_cu = (TH1D*)h_npe_s[seg]->Clone(Form("h_npe_cu_seg%d", seg));
        TH1D *h_nwt_cu = (TH1D*)h_nwt[seg]->Clone(Form("h_nwt_cu_seg%d", seg));
        h_npe_cu->GetXaxis()->SetRange(h_npe_cu->GetXaxis()->FindBin(0.), h_npe_cu->GetXaxis()->FindBin(kvc_ratio_closeup_xmax));
        h_nwt_cu->GetXaxis()->SetRange(h_nwt_cu->GetXaxis()->FindBin(0.), h_nwt_cu->GetXaxis()->FindBin(kvc_ratio_closeup_xmax));
        h_npe_cu->SetTitle(Form("KVC seg%d: npe_{sum} (all vs TDC, close-up);npe_{sum};count", seg));
        h_npe_cu->SetLineColor(kBlue);
        h_npe_cu->Draw();
        h_nwt_cu->SetLineColor(kRed);
        h_nwt_cu->Draw("same");
        TF1 *f_gaus = nullptr;
        if (ped_sigma[seg] > 1e-9) {
            f_gaus = new TF1(Form("f_gaus_ped_seg%d", seg), "gaus", 0., kvc_ratio_closeup_xmax);
            f_gaus->SetParameters(ped_amplitude[seg], ped_mean[seg], ped_sigma[seg]);
            f_gaus->SetLineColor(kOrange);
            f_gaus->SetNpx(500);
            f_gaus->Draw("same");
        }
        TLegend *leg = new TLegend(0.6, 0.75, 0.88, 0.88);
        leg->AddEntry(h_npe_cu, "npe_{sum} (all)", "l");
        leg->AddEntry(h_nwt_cu, "npe_{sum} (TDC)", "l");
        if (f_gaus) leg->AddEntry(f_gaus, "pedestal Gauss fit", "l");
        leg->Draw();
        nth_pad++;

        maybe_new_page(nth_pad);
        c->cd(nth_pad)->SetLogy(0);
        TF1 *f_sigmoid = new TF1("f_sigmoid", "0.5*(1+TMath::Erf((x-[0])/[1]))", 0., conf.npe_max);
        f_sigmoid->SetParameters(20., 5.);
        f_sigmoid->SetParNames("x0", "sigma");
        Int_t fit_status = h_ratio_npe[seg]->Fit(f_sigmoid, "Q0", "", 0., 150.);
        if (fit_status == 0) threshold_at_05[seg] = f_sigmoid->GetParameter(0);

        TH1D *h_closeup = (TH1D*)h_ratio_npe[seg]->Clone(Form("h_ratio_closeup_seg%d", seg));
        h_closeup->GetXaxis()->SetRange(h_closeup->GetXaxis()->FindBin(0.), h_closeup->GetXaxis()->FindBin(kvc_ratio_closeup_xmax));
        h_closeup->SetTitle(Form("KVC seg%d: ratio nwt/all (close-up);npe_{sum};ratio (nwt/all)", seg));
        h_closeup->Draw("E");
        if (fit_status == 0) {
            f_sigmoid->SetLineColor(kOrange);
            f_sigmoid->SetNpx(1000);
            f_sigmoid->Draw("same");
            TLine *vline = new TLine(threshold_at_05[seg], 0., threshold_at_05[seg], 1.);
            vline->SetLineColor(kRed);
            vline->SetLineStyle(2);
            vline->Draw();
        }
        TLatex *lt = new TLatex();
        lt->SetNDC(kTRUE);
        if (fit_status == 0)
            lt->DrawLatex(0.18, 0.82, Form("npe_{0.5} = %.2f", threshold_at_05[seg]));
        else
            lt->DrawLatex(0.18, 0.82, "npe_{0.5} = (fit failed)");
        if (n_sigma[seg] > -0.5)
            lt->DrawLatex(0.18, 0.75, Form("n_{#sigma} = %.2f", n_sigma[seg]));
        nth_pad++;
    }
    c->Print(pdf_path);
    c->Print(pdf_path + "]");
    delete c;

    // +-------+
    // | Write |
    // +-------+
    fout->cd();
    for (Int_t seg = 0; seg < num_seg; ++seg) {
        h_npe_s[seg]->Write();
        h_nwt[seg]->Write();
        h_ratio_npe[seg]->Write();
    }
    TTree *tree_thr = new TTree("threshold", "KVC ratio=0.5 threshold per segment");
    Int_t seg_out;
    Double_t thr_val;
    Double_t n_sigma_out;
    tree_thr->Branch("seg", &seg_out);
    tree_thr->Branch("npe_sum_at_05", &thr_val);
    tree_thr->Branch("n_sigma", &n_sigma_out);
    for (Int_t seg = 0; seg < num_seg; ++seg) {
        seg_out = seg;
        thr_val = threshold_at_05[seg];
        n_sigma_out = n_sigma[seg];
        tree_thr->Fill();
    }
    tree_thr->Write();
    fout->Close();

    Double_t sum_npe = 0., sum_nsigma = 0.;
    Int_t n_npe = 0, n_nsigma = 0;
    for (Int_t seg = 0; seg < num_seg; ++seg) {
        if (threshold_at_05[seg] > -0.5) { sum_npe += threshold_at_05[seg]; n_npe++; }
        if (n_sigma[seg] > -0.5) { sum_nsigma += n_sigma[seg]; n_nsigma++; }
    }
    if (n_npe > 0) std::cout << "KVC npe_{0.5} mean = " << (sum_npe / n_npe) << " (" << n_npe << " seg)" << std::endl;
    if (n_nsigma > 0) std::cout << "KVC n_sigma mean = " << (sum_nsigma / n_nsigma) << " (" << n_nsigma << " seg)" << std::endl;
}

Int_t main(int argc, char** argv) {
    if (argc < 2) return 1;
    analyze(std::atoi(argv[1]));
    return 0;
}