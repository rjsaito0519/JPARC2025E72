// c++
#include <cstdlib>
#include <stdexcept>
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
#include "hodo_input.h"
#include "paths.h"

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
    TString output_path = Form("%s/run%05d_T0_Offset_%s.root", root_base_dir.Data(), run_num, particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // -- tdc ----------
    TH1D *h_t0_offset[conf.num_of_ch.at("bh2")];
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++ ) h_t0_offset[i] = (TH1D*)f->Get(Form("T0_seg%d_TimeOffset%s", i, input_suffix.Data()));

    fout->cd();
    for (Int_t segment = 0; segment < conf.num_of_ch.at("bh2"); ++segment)
        h_t0_offset[segment]->Write();

    // -- set tdc range ----------
    TH1D *h_sum_tdc = (TH1D*)h_t0_offset[0]->Clone("h_sum_tdc");
    h_sum_tdc->Reset(); 
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++) h_sum_tdc->Add(h_t0_offset[i]);
    ana_helper::set_tdc_search_range(h_sum_tdc);

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_T0_Offset_%s.pdf", img_base_dir.Data(), run_num, particle.Data());

    // -- container -----
    std::vector<FitResult> t0_offset;
    std::vector<ana_helper::MipFitState> snapshots;

    auto c_t0 = new TCanvas("t0", "", 1500, 1200);
    c_t0->Divide(cols, rows);
    c_t0->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++) {
        if (nth_pad > max_pads) {
            c_t0->Print(pdf_path);
            c_t0->Clear();
            c_t0->Divide(cols, rows);
            nth_pad = 1;
        }

        ana_helper::MipFitSettings settings;
        const char* override_path = std::getenv("E72_T0_FIT_OVERRIDE");
        const std::string key = "bh2-"+std::to_string(i)+"-t0:"+particle.Data();
        const bool manual = override_path && ana_helper::read_mip_settings(override_path, key, settings);
        // All can contain several populations. Seed the explicitly selected raw
        // T0 peak; manual Web settings still take precedence.
        bool selected_peak = false;
        if (hodo_input::all_histograms() && !manual) {
            const auto window = hodo_input::t0_window();
            auto* histogram = h_t0_offset[i];
            const int lower = histogram->FindBin(window.first);
            const int upper = histogram->FindBin(window.second);
            if (histogram->Integral(lower, upper) >= 100.) {
                histogram->GetXaxis()->SetRange(lower, upper);
                const double center = histogram->GetBinCenter(histogram->GetMaximumBin());
                const double half = std::min(.8, (window.second-window.first)*.5);
                settings.lower = std::max(window.first, center-half);
                settings.upper = std::min(window.second, center+half);
                settings.function = "gaussian";
                settings.range_mode = "auto";
                settings.parameters = {{histogram->Integral(lower, upper)*histogram->GetBinWidth(1)*3.,
                    center, .3, 0.}};
                histogram->GetXaxis()->UnZoom();
                selected_peak = true;
            }
        }
        ana_helper::MipFitState snapshot;
        FitResult result = ana_helper::t0_offset_fit(h_t0_offset[i], c_t0, nth_pad,
            manual || selected_peak ? &settings : nullptr, &snapshot);
        if ((manual || selected_peak) && !snapshot.success)
            throw std::runtime_error("selected T0 fit failed for "+key);
        snapshots.push_back(snapshot);
        t0_offset.push_back(result);
        nth_pad++;
    }
    c_t0->Print(pdf_path);
    c_t0->Print(pdf_path + "]"); // end
    delete c_t0;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> offset_p0_val; 
    std::vector<Double_t> offset_p0_err; 
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("offset_p0_val", &offset_p0_val);
    tree->Branch("offset_p0_err", &offset_p0_err);
    
    for (Int_t i = 0; i < conf.num_of_ch.at("bh2"); i++) {
        ch = i;
        offset_p0_val.clear();
        offset_p0_err.clear();
        
        offset_p0_val.push_back(t0_offset[i].par[1]);
        offset_p0_err.push_back(t0_offset[i].err[1]);
    
        tree->Fill();
    }
    
    fout->cd();
    tree->Write();
    TTree fit_tree("t0_fit", "Saved T0 Gaussian and fit range");
    Int_t segment, status, covariance, ndf;
    Bool_t success;
    Double_t parameters[3], range[2], chi_square;
    fit_tree.Branch("seg", &segment);
    fit_tree.Branch("parameters", parameters, "parameters[3]/D");
    fit_tree.Branch("range", range, "range[2]/D");
    fit_tree.Branch("status", &status);
    fit_tree.Branch("covariance", &covariance);
    fit_tree.Branch("ndf", &ndf);
    fit_tree.Branch("success", &success);
    fit_tree.Branch("chi_square", &chi_square);
    for (segment = 0; segment < conf.num_of_ch.at("bh2"); ++segment) {
        const auto& snapshot = snapshots[segment];
        std::copy_n(snapshot.parameters.begin(), 3, parameters);
        range[0] = snapshot.lower;
        range[1] = snapshot.upper;
        status = snapshot.status;
        covariance = snapshot.covariance;
        ndf = snapshot.ndf;
        success = snapshot.success;
        chi_square = snapshot.chi_square;
        fit_tree.Fill();
    }
    fit_tree.Write();
    fit_tree.SetDirectory(nullptr);
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