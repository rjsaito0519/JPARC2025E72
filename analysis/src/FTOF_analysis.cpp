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

    TTreeReaderValue<std::vector<Double_t>> cvc_hit_seg(reader_hodo, "cvc_hit_seg");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> cvc_time_u(reader_hodo, "cvc_time_u");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> cvc_time_d(reader_hodo, "cvc_time_d");

    TTreeReaderValue<std::vector<std::vector<Double_t>>> sfv_tdc(reader_hodo, "sfv_tdc_u");
    TTreeReaderValue<std::vector<std::vector<Double_t>>> sac3_tdc(reader_hodo, "sac3_tdc_u");

    TTreeReaderValue<Double_t> btof(reader_hodo, "btof0");
    TTreeReaderValue<Double_t> ftof(reader_hodo, "ftof0");

    // -- bcout-----
    TString root_file_path_bcout = Form("%s/dc_out_run%05d_Pi_resi.root", DATA_DIR.Data(), run_num);
    auto *f_bcout = new TFile( root_file_path_bcout.Data() );
    if (!f_bcout || f_bcout->IsZombie()) {
        std::cerr << "Error: Could not open file : " << root_file_path_bcout << std::endl;
        return;
    }
    TTreeReader reader_bcout("bcout", f_bcout);
    TTreeReaderValue<unsigned int> evnum_bcout(reader_bcout, "event_number");

    TTreeReaderValue<std::vector<Double_t>> x0(reader_bcout, "x0");
    TTreeReaderValue<std::vector<Double_t>> u0(reader_bcout, "u0");
    TTreeReaderValue<std::vector<Double_t>> y0(reader_bcout, "y0");
    TTreeReaderValue<std::vector<Double_t>> v0(reader_bcout, "v0");
    

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    // -- CVC -----    
    std::vector<HistPair> h_cvc_time_u;
    std::vector<HistPair> h_cvc_time_d;
    std::vector<HistPair> h_cvc_time_diff;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cvc"); ch++) {
        {
            TString name  = Form("CVCu_time_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d CVC UP(Time) ch%d;Time [ns];", run_num, ch + 1);
            h_cvc_time_u.emplace_back(name, title, conf.time_bin_num, conf.time_min, conf.time_max);
        }
        {
            TString name  = Form("CVCd_time_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d CVC Down(Time) ch%d;Time [ns];", run_num, ch + 1);
            h_cvc_time_d.emplace_back(name, title, conf.time_bin_num, conf.time_min, conf.time_max);
        }
        {
            TString name  = Form("CVCdiff_time_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d CVC Diff(Time) ch%d;Time [ns];", run_num, ch + 1);
            h_cvc_time_diff.emplace_back(name, title, 600, -15.0, 15.0);
        }
    }

    // -- SFV -----    
    std::vector<HistPair> h_sfv_tdc;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("sfv"); ch++) {
        TString name  = Form("SFV_tdc_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d SFV (TDC) ch%d;TDC;", run_num, ch + 1);
        h_sfv_tdc.emplace_back(name, title, conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    // -- SAC3 -----    
    std::vector<HistPair> h_sac3_tdc;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("sac3"); ch++) {
        TString name  = Form("SAC3_tdc_%d_%d", run_num, ch + 1);
        TString title = Form("run%05d SAC3 (TDC) ch%d;TDC;", run_num, ch + 1);
        h_sac3_tdc.emplace_back(name, title, conf.tdc_bin_num, conf.tdc_min, conf.tdc_max);
    }

    // -- ftof -----    
    std::vector<HistPair> h_ftof_pi;
    std::vector<HistPair> h_ftof_pbar;
    {
        TString name  = Form("FTOF_Pi_%d", run_num);
        TString title = Form("run%05d FTOF (#pi);Time [ns];", run_num);
        h_ftof_pi.emplace_back(name, title, conf.time_bin_num, conf.time_min, conf.time_max);
    }
    {
        TString name  = Form("FTOF_pbar_%d", run_num);
        TString title = Form("run%05d FTOF (pbar);Time [ns];", run_num);
        h_ftof_pbar.emplace_back(name, title, conf.time_bin_num, conf.time_min, conf.time_max);
    }

    // -- CVC profile -----
    std::vector<HistPair2D> h_cvc_seg_profile;
    std::vector<HistPair2D> h_cvc_time_profile;
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cvc")+1; ch++) {
        {
            TString name  = Form("CVC_profile_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d CVC seg%d;x;y;", run_num, ch + 1);
            h_cvc_seg_profile.emplace_back(name, title, 1000, -1000.0, 1000.0, 500, -500.0, 500.0);
        }
        {
            TString name  = Form("CVC_profile_time_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d CVC seg%d;time diff;y;", run_num, ch + 1);
            h_cvc_time_profile.emplace_back(name, title, 100, -10.0, 10.0, 500, -500.0, 500.0);
        }
    }


    // +------------------+
    // | Fill event (1st) |
    // +------------------+
    Int_t evnum = 0;
    reader_hodo.Restart();
    reader_bcout.Restart();
    while (reader_hodo.Next() && reader_bcout.Next()){ displayProgressBar(++evnum, total_entry);
        if (*evnum_hodo != *evnum_bcout) {
            std::cerr << "Error: event numbers are not mutched; hodo : " << *evnum_hodo << ", bcout : " << *evnum_bcout << std::endl;
            return;
        }

        // if (*btof0 > conf.btof_threshold) continue; 

        // // -- SFV -----
        // for (Int_t i = 0, n_i = (*sfv_tdc).size(); i < n_i; i++) {
        //     for (Int_t j = 0, n_j = (*sfv_tdc)[i].size(); j < n_j; j++) {
        //         h_sfv_tdc[i].raw->Fill((*sfv_tdc)[i][j]);
        //     }
        // }

        // -- SAC3 -----
        for (Int_t i = 0, n_i = (*sac3_tdc).size(); i < n_i; i++) {
            for (Int_t j = 0, n_j = (*sac3_tdc)[i].size(); j < n_j; j++) {
                h_sac3_tdc[i].raw->Fill((*sac3_tdc)[i][j]);
            }
        }
    }

    // +--------------+
    // | Fit and Plot |
    // +--------------+
    // メインフレームを作成
    TGMainFrame *main = new TGMainFrame(gClient->GetRoot(), 1000, 800);
    // Handle the window close event to terminate the application
    main->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    // タブウィジェットを作成
    TGTab *tab = new TGTab(main, 1000, 800);

    // -- BAC -----
    TCanvas *c_sac3 = ana_helper::add_tab(tab, "sac3");
    std::vector<FitResult> result_sac3;
    conf.detector = "sac3";
    for (Int_t ch = 0; ch < conf.num_of_ch.at("sac3"); ch++) {
        FitResult tmp_result = ana_helper::tdc_fit(h_sac3_tdc[ch].raw, c_sac3, ch+1);
        result_sac3.push_back(tmp_result);
    }
    
    // +------------------+
    // | Fill event (2nd) |
    // +------------------+
    evnum = 0;
    Double_t n_sigma = 5.0;
    Double_t sfv_tdc_min = 0.37E6;
    Double_t sfv_tdc_max = 0.42E6;
    reader_hodo.Restart();
    reader_bcout.Restart();
    while (reader_hodo.Next() && reader_bcout.Next()){ displayProgressBar(++evnum, total_entry);
        Bool_t is_btof_pi = true;
        if ((*btof) < -10.0) is_btof_pi = false;

        std::vector<Bool_t> has_hit_cvc(conf.num_of_ch.at("cvc"), false);
        for (Int_t i = 0, n_i = (*cvc_hit_seg).size(); i < n_i; i++) {
            Int_t index = static_cast<Int_t>((*cvc_hit_seg)[i]);
            if (0 <= index && index < conf.num_of_ch.at("cvc")) {
                Bool_t has_hit_cvc_ud[2] = {false, false};
                for (Int_t j = 0, n_j = (*cvc_time_u)[i].size(); j < n_j; j++) {
                    if ( 15.0 < (*cvc_time_u)[i][j] ) has_hit_cvc_ud[0] = true;
                }
                for (Int_t j = 0, n_j = (*cvc_time_d)[i].size(); j < n_j; j++) {
                    if ( 15.0 < (*cvc_time_d)[i][j] ) has_hit_cvc_ud[1] = true;
                }
                if (has_hit_cvc_ud[0] && has_hit_cvc_ud[1]) has_hit_cvc[index] = true;
            }
        }

        Bool_t has_hit_sfv = false;
        for (Int_t i = 0, n_i = (*sfv_tdc).size(); i < n_i; i++) {
            for (Int_t j = 0, n_j = (*sfv_tdc)[i].size(); j < n_j; j++) {
                Double_t lower = sfv_tdc_min;
                Double_t upper = sfv_tdc_max;
                if ( lower < (*sfv_tdc)[i][j] && (*sfv_tdc)[i][j] < upper ) has_hit_sfv = true;
            }
        }

        Bool_t has_hit_sac3 = false;
        for (Int_t i = 0, n_i = (*sac3_tdc).size(); i < n_i; i++) {
            for (Int_t j = 0, n_j = (*sac3_tdc)[i].size(); j < n_j; j++) {
                // Double_t lower = result_sac3[0].par[1] - 5.0*result_sac3[0].par[2];
                // Double_t upper = result_sac3[0].par[1] + 5.0*result_sac3[0].par[2];
                Double_t lower = 0.43E6;
                Double_t upper = 0.47E6;
                if ( lower < (*sac3_tdc)[i][j] && (*sac3_tdc)[i][j] < upper ) has_hit_sac3 = true;
            }
        }

        

        if (is_btof_pi && has_hit_sac3) {
            h_ftof_pi[0].raw->Fill(*ftof);
            if (has_hit_sfv) h_ftof_pi[0].trig->Fill(*ftof);
        } else if ( !is_btof_pi && !has_hit_sac3) {
            h_ftof_pbar[0].raw->Fill(*ftof);
            if (has_hit_sfv) h_ftof_pbar[0].trig->Fill(*ftof);
        }
 
        

        // -- tracking -----
        for (Int_t i = 0, n = (*x0).size(); i < n; i++) {
            Double_t x = (*x0)[i] + (*u0)[i]*conf.cvc_pos_z;
            Double_t y = (*y0)[i] + (*v0)[i]*conf.cvc_pos_z;
            for (Int_t ch = 0; ch < conf.num_of_ch.at("cvc"); ch++) {
                if (has_hit_cvc[ch]) {
                    h_cvc_seg_profile[ch].trig->Fill(x, y);
                    h_cvc_seg_profile[conf.num_of_ch.at("cvc")].trig->Fill(x, y);

                    for (Int_t i = 0, n_i = (*cvc_hit_seg).size(); i < n_i; i++) {
                        Int_t index = static_cast<Int_t>((*cvc_hit_seg)[i]);
                        if (0 <= index && index < conf.num_of_ch.at("cvc") && index == ch) {
                            if ( (*cvc_time_u)[i].size() > 0 && (*cvc_time_d)[i].size() > 0 ) {
                                Double_t diff = (*cvc_time_u)[i][0] - (*cvc_time_d)[i][0];
                                h_cvc_time_diff[ch].trig->Fill(diff);
                                h_cvc_time_profile[ch].trig->Fill(diff, y);
                                h_cvc_time_profile[conf.num_of_ch.at("cvc")].trig->Fill(diff, y);
                            }
                        }
                    }
                }
            }
        }

    }

    // -- FTOF -----
    TCanvas *c_ftof = ana_helper::add_tab(tab, "ftof");
    c_ftof->Divide(1, 2);
    c_ftof->cd(1)->SetLogy(1);
    h_ftof_pi[0].raw->SetLineColor(kBlue);
    h_ftof_pi[0].raw->Draw();
    h_ftof_pi[0].trig->SetLineColor(kRed);
    h_ftof_pi[0].trig->Draw();
    c_ftof->cd(2)->SetLogy(1);
    h_ftof_pbar[0].raw->SetLineColor(kBlue);
    h_ftof_pbar[0].raw->Draw();
    h_ftof_pbar[0].trig->SetLineColor(kRed);
    h_ftof_pbar[0].trig->Draw();

    TCanvas *c_cvc_diff = ana_helper::add_tab(tab, "diff");
    c_cvc_diff->Divide(3, 3);
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cvc"); ch++) {
        c_cvc_diff->cd(ch+1);
        h_cvc_time_diff[ch].trig->Draw();
    }

    TCanvas *c_cvc_seg_profile = ana_helper::add_tab(tab, "profile");
    c_cvc_seg_profile->Divide(3, 3);
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cvc")+1; ch++) {
        c_cvc_seg_profile->cd(ch+1);
        h_cvc_seg_profile[ch].trig->Draw("colz");
    }

    TCanvas *c_cvc_time_profile = ana_helper::add_tab(tab, "profile");
    c_cvc_time_profile->Divide(3, 3);
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cvc")+1; ch++) {
        c_cvc_time_profile->cd(ch+1);
        h_cvc_time_profile[ch].trig->Draw("colz");
    }

    // メインフレームにタブウィジェットを追加
    main->AddFrame(tab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    // ウィンドウを表示
    main->MapSubwindows();
    main->Resize(main->GetDefaultSize());
    main->MapWindow();
}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <run number>" << std::endl;
        return 1;
    }
    Int_t run_num = std::atoi(argv[1]);

    TApplication *theApp = new TApplication("App", &argc, argv);    
    analyze(run_num);
    theApp->Run();

    return 0;
}