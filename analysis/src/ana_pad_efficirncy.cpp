// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <set>
#include <cmath>
#include <memory>

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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPolyLine.h>
#include <TMath.h>
#include <Rtypes.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"

Config& conf = Config::getInstance();
const Int_t N_HIT_THRESHOLD = 15;

// DstTPCHitBcOutTracking 系: BcOut 直線との residual 閾値 [mm]（--mode=hit 時）
const Double_t RES_CUT_X_MM = 10.0;
const Double_t RES_CUT_Y_MM = 10.0;

namespace {

Bool_t good_res(Double_t v)
{
    return (Bool_t)(std::isfinite(v) && !TMath::IsNaN(v));
}

/** X と Y を別々に判定（どちらかが閾値外ならそのヒットは使わない） */
Bool_t pass_residual_xy(Double_t rx, Double_t ry)
{
    if (!good_res(rx) || TMath::Abs(rx) > RES_CUT_X_MM) return kFALSE;
    if (!good_res(ry) || TMath::Abs(ry) > RES_CUT_Y_MM) return kFALSE;
    return kTRUE;
}

enum InputMode { kModeRaw = 0, kModeHitBcOut = 1 };

static Bool_t parse_mode_args(int argc, char** argv, InputMode& out_mode, TString& err)
{
    out_mode = kModeRaw;
    err = "";
    Bool_t found = kFALSE;
    for (int i = 1; i < argc; i++) {
        TString a = argv[i];
        if (a == "--bcout") {
            found = kTRUE;
            out_mode = kModeHitBcOut;
        } else if (a.BeginsWith("--mode=")) {
            TString m = a(7, a.Length());
            m = m.Strip(TString::kBoth);
            if (m.IsNull()) {
                err = "empty value in --mode=";
                return kFALSE;
            }
            found = kTRUE;
            if (m == "bcout" || m == "hit" || m == "hitbcout") {
                out_mode = kModeHitBcOut;
            } else {
                err = Form("unknown --mode value: %s (use bcout)", m.Data());
                return kFALSE;
            }
        } else if (a == "--mode") {
            if (i + 1 >= argc) {
                err = "--mode needs a value (bcout)";
                return kFALSE;
            }
            TString m = argv[++i];
            m = m.Strip(TString::kBoth);
            found = kTRUE;
            if (m == "bcout" || m == "hit" || m == "hitbcout") {
                out_mode = kModeHitBcOut;
            } else {
                err = Form("unknown --mode value: %s (use bcout)", m.Data());
                return kFALSE;
            }
        }
    }
    return found;
}

}  // namespace

void analyze(TString path, InputMode mode, Bool_t mode_explicit)
{
    // +---------+
    // | setting |
    // +---------+
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0 / 256, 160.0 / 256, 44.0 / 256);

    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x");  // x軸のタイトルサイズ
    gStyle->SetTitleSize(0.06, "y");  // y軸のタイトルサイズ
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // +-----------+
    // | load file |
    // +-----------+
    auto* f = new TFile(path.Data());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << path << std::endl;
        return;
    }
    TTree* t = dynamic_cast<TTree*>(f->Get("tpc"));
    const InputMode mode_eff = mode_explicit ? mode : kModeRaw;
    const Bool_t has_layer_tpc = (Bool_t)(t && t->GetBranch("layerTpc"));

    if (mode_eff == kModeRaw) {
        if (!t || !has_layer_tpc) {
            std::cerr << "Error: tree \"tpc\" needs branch layerTpc when --mode is not set." << std::endl;
            delete f;
            return;
        }
    } else {
        const Bool_t has_cluster_layer = (Bool_t)(t && t->GetBranch("cluster_layer"));
        if (!t || (!has_layer_tpc && !has_cluster_layer) || !t->GetBranch("residual_x_hit") ||
            !t->GetBranch("residual_y_hit")) {
            std::cerr << "Error: --bcout requires residual_x_hit, residual_y_hit and "
                         "(layerTpc or cluster_layer)."
                      << std::endl;
            delete f;
            return;
        }
    }

    TTreeReader reader("tpc", f);
    Int_t total_entry = reader.GetEntries();

    TTreeReaderValue<UInt_t> run_number(reader, "run_number");

    std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> layer_tpc;
    std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> cluster_layer;
    std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> residual_x_hit;
    std::unique_ptr<TTreeReaderValue<std::vector<Double_t>>> residual_y_hit;

    if (mode_eff == kModeRaw || (mode_eff == kModeHitBcOut && has_layer_tpc)) {
        layer_tpc = std::make_unique<TTreeReaderValue<std::vector<Int_t>>>(reader, "layerTpc");
    }
    if (mode_eff == kModeHitBcOut && t->GetBranch("cluster_layer")) {
        cluster_layer = std::make_unique<TTreeReaderValue<std::vector<Int_t>>>(reader, "cluster_layer");
    }
    if (mode_eff == kModeHitBcOut) {
        residual_x_hit = std::make_unique<TTreeReaderValue<std::vector<Double_t>>>(reader, "residual_x_hit");
        residual_y_hit = std::make_unique<TTreeReaderValue<std::vector<Double_t>>>(reader, "residual_y_hit");
    }

    reader.SetEntry(0);
    Int_t run_num = static_cast<Int_t>(*run_number);

    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    TH1D* h_cobo_time[conf.num_of_ch.at("cobo")];
    TH1D* h_cobo_tdc[conf.num_of_ch.at("cobo")];
    for (Int_t ch = 0; ch < conf.num_of_ch.at("cobo"); ch++) {
        {
            TString name = Form("COBO_time_diff_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d COBO (diff) ch%d;Time [ns];", run_num, ch + 1);
            h_cobo_time[ch] = new TH1D(name, title, 100, 79.5, 80.5);
        }
        {
            TString name = Form("COBO_tdc_diff_%d_%d", run_num, ch + 1);
            TString title = Form("run%05d COBO (diff) ch%d;TDC [ch];", run_num, ch + 1);
            h_cobo_tdc[ch] = new TH1D(name, title, 1000, 81500.0, 82500.0);
        }
    }
    auto* h_cobo_1st_hit =
        new TH1D(Form("COBO_1st_hit_%d", run_num), Form("run%05d COBO 1st hit;TDC [ch];", run_num), 1000, 1.E6,
                 2.E6);

    // +------------+
    // | Fill event |
    // +------------+
    std::vector<Int_t> n_hit_wo_me(conf.NumOfLayersTPC, 0);
    std::vector<Int_t> n_hit_w_me(conf.NumOfLayersTPC, 0);

    reader.Restart();
    Int_t evnum = 0;
    while (reader.Next()) {
        displayProgressBar(++evnum, total_entry);
        std::set<Int_t> unique_layer;

        if (mode_eff == kModeRaw) {
            for (const auto it : **layer_tpc) unique_layer.insert(it);
        } else {
            const auto& layers = layer_tpc ? **layer_tpc : **cluster_layer;
            const auto& rx = **residual_x_hit;
            const auto& ry = **residual_y_hit;
            const size_t n = std::min({layers.size(), rx.size(), ry.size()});
            for (size_t i = 0; i < n; i++) {
                if (pass_residual_xy(rx[i], ry[i])) unique_layer.insert(layers[i]);
            }
        }

        for (Int_t layer = 0; layer < conf.NumOfLayersTPC; layer++) {
            Int_t n_hit_me = unique_layer.count(layer);
            Int_t tmp_n_hit_wo_me = static_cast<Int_t>(unique_layer.size()) - n_hit_me;
            if (tmp_n_hit_wo_me > N_HIT_THRESHOLD) {
                n_hit_wo_me[layer]++;
                if (n_hit_me) n_hit_w_me[layer]++;
            }
        }
    }

    delete f;

    // +----------------+
    // | cal efficiency |
    // +----------------+
    std::vector<Double_t> x, y, ex, ey;
    for (Int_t layer = 0; layer < conf.NumOfLayersTPC; layer++) {
        x.push_back(layer);
        ex.push_back(0.0);
        Double_t eff = 0.0;
        Double_t eff_err = 0.0;
        if (n_hit_wo_me[layer] > 0) {
            eff = static_cast<Double_t>(n_hit_w_me[layer]) / static_cast<Double_t>(n_hit_wo_me[layer]);
            eff_err = std::sqrt(eff * (1.0 - eff) / static_cast<Double_t>(n_hit_wo_me[layer]));
        }
        y.push_back(eff);
        ey.push_back(eff_err);
    }
    auto* g_eff = new TGraphErrors(conf.NumOfLayersTPC, x.data(), y.data(), ex.data(), ey.data());

    TString title;
    TString pdf_suffix;
    if (mode_eff == kModeRaw) {
        title = Form("Run %05d (Threshold: %d);Layer;Layer Efficiency", run_num, N_HIT_THRESHOLD);
        pdf_suffix = "";
    } else {
        title = Form("Run %05d (Threshold: %d, BcOut hit, |resX|<%.4g |resY|<%.4g mm);Layer;Layer Efficiency",
                     run_num, N_HIT_THRESHOLD, RES_CUT_X_MM, RES_CUT_Y_MM);
        pdf_suffix = Form("_bcout_hit_resX%g_resY%g", RES_CUT_X_MM, RES_CUT_Y_MM);
    }
    g_eff->SetTitle(title);
    g_eff->SetMarkerStyle(20);
    g_eff->SetMarkerSize(2.0);
    g_eff->SetMarkerColor(kBlue);
    g_eff->SetLineColor(kBlue);
    g_eff->SetMinimum(0.0);
    g_eff->SetMaximum(1.05);

    // +------+
    // | plot |
    // +------+
    Int_t nth_pad = 1;
    Int_t rows = 1, cols = 1;
    Int_t max_pads = rows * cols;

    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_pad_efficiency%s.pdf", img_base_dir.Data(), run_num, pdf_suffix.Data());

    auto* c_eff = new TCanvas("", "", 1500, 1200);
    c_eff->Divide(cols, rows);
    c_eff->Print(pdf_path + "[");  // start
    nth_pad = 1;
    for (Int_t i = 0; i < 1; i++) {
        if (nth_pad > max_pads) {
            c_eff->Print(pdf_path);
            c_eff->Clear();
            c_eff->Divide(cols, rows);
            nth_pad = 1;
        }

        c_eff->cd(nth_pad);
        g_eff->Draw("AP");
        nth_pad++;
    }

    c_eff->Print(pdf_path);
    c_eff->Print(pdf_path + "]");  // end
    delete c_eff;

    // // +--------------------------+
    // // | prepare output root file |
    // // +--------------------------+
    // TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    // TString output_path = Form("%s/run%05d_COBO_diff.root", root_base_dir.Data(), run_num);
    //
    // if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    // TFile* fout = new TFile(output_path.Data(), "RECREATE");
    //
    // // +-------+
    // // | Write |
    // // +-------+
    // TTree* tree = new TTree("tree", "");
    // Int_t ch;
    // std::vector<Double_t> diff;
    // tree->Branch("ch", &ch, "ch/I");
    // tree->Branch("diff", &diff);
    //
    // for (Int_t i = 0; i < conf.num_of_ch.at("cobo"); i++) {
    //     ch = i;
    //     diff = cobo_diff[i];
    //     tree->Fill();
    //
    //     h_cobo_tdc[i]->Write();
    //     h_cobo_time[i]->Write();
    // }
    // h_cobo_1st_hit->Write();
    //
    // fout->cd();
    // tree->Write();
    // fout->Close();
}

Int_t main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root file> [--bcout]\n"
                  << "  Without --bcout: use all hits from layerTpc (no residual-based selection).\n"
                  << "  --bcout: use DstTPCHitBcOutTracking style residual_x/y_hit against BcOut line.\n"
                  << "  (compat: --mode=bcout is also accepted)\n"
                  << "  Edit RES_CUT_X_MM / RES_CUT_Y_MM near the top of this file to change cuts."
                  << std::endl;
        return 1;
    }

    InputMode mode = kModeRaw;
    TString parse_err;
    const Bool_t mode_explicit = parse_mode_args(argc, argv, mode, parse_err);
    if (!parse_err.IsNull()) {
        std::cerr << "Error: " << parse_err << std::endl;
        return 1;
    }

    TString path;
    for (int i = 1; i < argc; i++) {
        TString a = argv[i];
        if (a == "--bcout") continue;
        if (a.BeginsWith("--mode=")) continue;
        if (a == "--mode") {
            ++i;
            continue;
        }
        if (a.BeginsWith("--")) {
            std::cerr << "Error: unknown option " << a << std::endl;
            return 1;
        }
        if (path.IsNull()) {
            path = a;
        } else {
            std::cerr << "Error: unexpected argument: " << a << std::endl;
            return 1;
        }
    }
    if (path.IsNull()) {
        std::cerr << "Error: no ROOT file path given." << std::endl;
        return 1;
    }

    analyze(path, mode, mode_explicit);
    return 0;
}
