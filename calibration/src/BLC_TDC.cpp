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

void analyze(TString path, TString particle, Bool_t use_guess){    
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

    const char* treeName = nullptr;
    TString in_or_out = "";
    if (f->Get("bcout")) {
        treeName = "bcout";
        in_or_out = "2";
    }
    else if (f->Get("bcin")) {
        treeName = "bcin";
        in_or_out = "1";
    }
    else {
        std::cerr << "Error: neither 'bcout' nor 'bcin' tree exists.\n";
        return;
    }
    TTreeReader reader(treeName, f);
    TTreeReaderValue<unsigned int> run_number(reader, "run_number");
    reader.SetEntry(0);
    Int_t run_num = *run_number;

    // +--------------------------+
    // | prepare output root file |
    // +--------------------------+
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_path = Form("%s/run%05d_BLC%s_TDC_%s.root", root_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    // +-------------------+
    // | prepare histogram |
    // +-------------------+    
    // Halfwise TDC from 2D CTDC_vs_HitPat (segment vs channel)
    TH2D *h2_blca_tdc[conf.num_of_ch.at("blc")];
    TH2D *h2_blcb_tdc[conf.num_of_ch.at("blc")];
    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++ ) {
        h2_blca_tdc[i] = (TH2D*)f->Get(Form("BLC%sa_CTDC_vs_HitPat_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        h2_blcb_tdc[i] = (TH2D*)f->Get(Form("BLC%sb_CTDC_vs_HitPat_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        if (!h2_blca_tdc[i] || !h2_blcb_tdc[i]) {
            std::cerr << "Error: missing CTDC_vs_HitPat histogram for plane " << i
                      << " (BLC" << in_or_out << "a/b, particle=" << particle << ")" << std::endl;
            return;
        }
    }

    // +--------------+
    // | fit and plot |
    // +--------------+
    // -- prepare pdf -----
    Int_t nth_pad = 1;
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_BLC%s_TDC_%s.pdf", img_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());

    // -- container -----
    // index: 0=lo (wire 0-15), 1=hi (wire 16-31)
    std::vector<FitResult> blca_tdc_lo;
    std::vector<FitResult> blca_tdc_hi;
    std::vector<FitResult> blcb_tdc_lo;
    std::vector<FitResult> blcb_tdc_hi;

    auto c_blc = new TCanvas("blc", "", 1500, 1200);
    c_blc->Divide(cols, rows);
    c_blc->Print(pdf_path + "["); // start
    nth_pad = 1;
    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++) {
        // Determine bin ranges for wire 0-15 and 16-31 on the X axis
        TH2D* h2a = h2_blca_tdc[i];
        TH2D* h2b = h2_blcb_tdc[i];
        auto* xax = h2a->GetXaxis();

        // Use wire bin edges (wire integer centers are usually at 0..31 with edges -0.5..31.5)
        Int_t bLo1  = xax->FindBin(-0.5);
        Int_t bHi1  = xax->FindBin(15.5);
        Int_t bLo2  = xax->FindBin(16.5);
        Int_t bHi2  = xax->FindBin(31.5);
        Int_t bLoAll = xax->FindBin(-0.5);
        Int_t bHiAll = xax->FindBin(31.5);

        // Reserve pads in order: a_lo, a_hi, b_lo, b_hi
        auto get_pad = [&]() -> Int_t {
            if (nth_pad > max_pads) {
                c_blc->Print(pdf_path);
                c_blc->Clear();
                c_blc->Divide(cols, rows);
                nth_pad = 1;
            }
            return nth_pad++;
        };

        Int_t pad_a_lo = get_pad();
        Int_t pad_a_hi = get_pad();
        Int_t pad_b_lo = get_pad();
        Int_t pad_b_hi = get_pad();

        // Project all four halves first, then decide fit/copy to avoid Fit on empty hist
        TH1D* h_a_lo = h2a->ProjectionY(Form("BLC%sa_CTDC_plane%d_lo_%s", in_or_out.Data(), i, particle.Data()),
                                           bLo1, bHi1);
        TH1D* h_a_hi = h2a->ProjectionY(Form("BLC%sa_CTDC_plane%d_hi_%s", in_or_out.Data(), i, particle.Data()),
                                           bLo2, bHi2);
        TH1D* h_b_lo = h2b->ProjectionY(Form("BLC%sb_CTDC_plane%d_lo_%s", in_or_out.Data(), i, particle.Data()),
                                           bLo1, bHi1);
        TH1D* h_b_hi = h2b->ProjectionY(Form("BLC%sb_CTDC_plane%d_hi_%s", in_or_out.Data(), i, particle.Data()),
                                           bLo2, bHi2);
        TH1D* h_a_all = h2a->ProjectionY(Form("BLC%sa_CTDC_plane%d_all_%s", in_or_out.Data(), i, particle.Data()),
                                            bLoAll, bHiAll);
        TH1D* h_b_all = h2b->ProjectionY(Form("BLC%sb_CTDC_plane%d_all_%s", in_or_out.Data(), i, particle.Data()),
                                            bLoAll, bHiAll);

        auto empty_result = [&]() -> FitResult {
            FitResult r;
            r.par.assign(4, 0.0);
            r.err.assign(4, 0.0);
            r.chi_square = 0.0;
            r.ndf = 0;
            r.additional.assign(1, 0.0); // must exist: [0]
            return r;
        };

        // Default: each half uses its own init (nullptr seed).
        // --guess: borrow Erfc init / fit range from full-wire projection.
        ana_helper::DcTdcFitSeed seed_a;
        ana_helper::DcTdcFitSeed seed_b;
        const ana_helper::DcTdcFitSeed* pseed_a = nullptr;
        const ana_helper::DcTdcFitSeed* pseed_b = nullptr;
        if (use_guess) {
            seed_a = ana_helper::dc_tdc_seed(h_a_all);
            seed_b = ana_helper::dc_tdc_seed(h_b_all);
            pseed_a = seed_a.valid ? &seed_a : nullptr;
            pseed_b = seed_b.valid ? &seed_b : nullptr;
        }

        auto fit_halfwise = [&](TH1D* h_lo, TH1D* h_hi, TH1D* h_all,
                                Int_t pad_lo, Int_t pad_hi,
                                const ana_helper::DcTdcFitSeed* pseed,
                                FitResult& lo_fit, FitResult& hi_fit) {
            bool lo_ok = (h_lo && h_lo->GetEntries() > 0);
            bool hi_ok = (h_hi && h_hi->GetEntries() > 0);
            bool all_ok = (h_all && h_all->GetEntries() > 0);

            if (lo_ok && hi_ok) {
                lo_fit = ana_helper::dc_tdc_fit(h_lo, c_blc, pad_lo, pseed);
                hi_fit = ana_helper::dc_tdc_fit(h_hi, c_blc, pad_hi, pseed);
            } else if (lo_ok && !hi_ok) {
                lo_fit = ana_helper::dc_tdc_fit(h_lo, c_blc, pad_lo, pseed);
                hi_fit = lo_fit;
            } else if (!lo_ok && hi_ok) {
                hi_fit = ana_helper::dc_tdc_fit(h_hi, c_blc, pad_hi, pseed);
                lo_fit = hi_fit;
            } else if (all_ok) {
                FitResult all_fit = ana_helper::dc_tdc_fit(h_all, c_blc, pad_lo, nullptr);
                lo_fit = all_fit;
                hi_fit = all_fit;
            } else {
                lo_fit = empty_result();
                hi_fit = empty_result();
            }
        };

        FitResult a_lo_fit, a_hi_fit;
        fit_halfwise(h_a_lo, h_a_hi, h_a_all, pad_a_lo, pad_a_hi, pseed_a, a_lo_fit, a_hi_fit);

        FitResult b_lo_fit, b_hi_fit;
        fit_halfwise(h_b_lo, h_b_hi, h_b_all, pad_b_lo, pad_b_hi, pseed_b, b_lo_fit, b_hi_fit);

        blca_tdc_lo.push_back(a_lo_fit);
        blca_tdc_hi.push_back(a_hi_fit);
        blcb_tdc_lo.push_back(b_lo_fit);
        blcb_tdc_hi.push_back(b_hi_fit);
    }
    c_blc->Print(pdf_path);
    c_blc->Print(pdf_path + "]"); // end
    delete c_blc;

    // +-------+
    // | Write |
    // +-------+
    TTree* tree = new TTree("tree", "");
    Int_t ch;
    std::vector<Double_t> tdc_p0_val; 
    tree->Branch("ch", &ch, "ch/I");
    tree->Branch("tdc_p0_val", &tdc_p0_val);
    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++) {
        ch = i;
        tdc_p0_val.clear();
    
        // -- tdc (halfwise): [a_lo, a_hi, b_lo, b_hi] -----
        tdc_p0_val.push_back(blca_tdc_lo[i].additional[0]);
        tdc_p0_val.push_back(blca_tdc_hi[i].additional[0]);
        tdc_p0_val.push_back(blcb_tdc_lo[i].additional[0]);
        tdc_p0_val.push_back(blcb_tdc_hi[i].additional[0]);
    
        tree->Fill();
    }    
    
    fout->cd();
    tree->Write();
    fout->Close(); 
}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <root file path> <particle> [--guess]" << std::endl;
        return 1;
    }
    TString path = argv[1];
    TString particle = argv[2];
    if (particle != "Pi" && particle != "K") {
        std::cerr << "Error: Unexpected particle name: " << particle << std::endl;
        return 1;
    }

    Bool_t use_guess = kFALSE;
    for (Int_t i = 3; i < argc; ++i) {
        if (TString(argv[i]) == "--guess") {
            use_guess = kTRUE;
        } else {
            std::cerr << "Error: Unknown option: " << argv[i] << std::endl;
            std::cerr << "Usage: " << argv[0] << " <root file path> <particle> [--guess]" << std::endl;
            return 1;
        }
    }

    analyze(path, particle, use_guess);
    return 0;
}