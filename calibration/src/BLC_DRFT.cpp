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
#include <TDatime.h>
#include <TNamed.h>
#include <TGraph.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"

Config& conf = Config::getInstance();

namespace {

TGraph* clone_drift_graph(const TGraph* src, Int_t plane, const char* wire_range_suffix)
{
    if (!src || !wire_range_suffix || wire_range_suffix[0] == '\0') return nullptr;

    TString gname = Form("%s_Hit_DriftFunction_plane%d%s",
                         conf.detector.Data(), plane, wire_range_suffix);
    TString gtitle = Form("DriftFunction %s plane%d%s;Drift Time [ns];Drift Length [mm]",
                          conf.detector.Data(), plane, wire_range_suffix);
    TGraph* g = dynamic_cast<TGraph*>(src->Clone(gname));
    if (g) g->SetTitle(gtitle);
    return g;
}

void push_if(TGraph* g, std::vector<TGraph*>& out)
{
    if (g) out.push_back(g);
}

void draw_drift_hist(TH1D* h, TCanvas* c, Int_t pad)
{
    if (!h || !c || h->GetEntries() <= 0) return;
    c->cd(pad);
    h->Draw();
}

TGraph* draw_drift_graph(TH1D* h, TCanvas* c, Int_t pad, Int_t plane, const char* wire_range_suffix)
{
    if (!h || h->GetEntries() <= 0) return nullptr;
    return ana_helper::make_drift_function(h, c, pad, plane, wire_range_suffix);
}

// Legacy PDF layout per 2x2 page:
//   [lo DriftTime hist] [lo DriftFunction]
//   [hi DriftTime hist] [hi DriftFunction]
void make_halfwise_drift_page(TH1D* h_lo, TH1D* h_hi, TH1D* h_all,
                              TCanvas* c,
                              Int_t pad_lo_hist, Int_t pad_lo_graph,
                              Int_t pad_hi_hist, Int_t pad_hi_graph,
                              Int_t plane, std::vector<TGraph*>& out)
{
    const bool lo_ok  = (h_lo  && h_lo->GetEntries() > 0);
    const bool hi_ok  = (h_hi  && h_hi->GetEntries() > 0);
    const bool all_ok = (h_all && h_all->GetEntries() > 0);

    TGraph* g_lo = nullptr;
    TGraph* g_hi = nullptr;

    if (lo_ok && hi_ok) {
        draw_drift_hist(h_lo, c, pad_lo_hist);
        g_lo = draw_drift_graph(h_lo, c, pad_lo_graph, plane, "-w0-15");
        draw_drift_hist(h_hi, c, pad_hi_hist);
        g_hi = draw_drift_graph(h_hi, c, pad_hi_graph, plane, "-w16-31");
    } else if (lo_ok && !hi_ok) {
        draw_drift_hist(h_lo, c, pad_lo_hist);
        g_lo = draw_drift_graph(h_lo, c, pad_lo_graph, plane, "-w0-15");
        g_hi = clone_drift_graph(g_lo, plane, "-w16-31");
    } else if (!lo_ok && hi_ok) {
        draw_drift_hist(h_hi, c, pad_hi_hist);
        g_hi = draw_drift_graph(h_hi, c, pad_hi_graph, plane, "-w16-31");
        g_lo = clone_drift_graph(g_hi, plane, "-w0-15");
    } else if (all_ok) {
        draw_drift_hist(h_all, c, pad_lo_hist);
        g_lo = draw_drift_graph(h_all, c, pad_lo_graph, plane, "-w0-15");
        g_hi = clone_drift_graph(g_lo, plane, "-w16-31");
    }

    push_if(g_lo, out);
    push_if(g_hi, out);
}

void print_drift_page(TCanvas* c, const TString& pdf_path)
{
    c->Print(pdf_path);
    c->Clear();
    c->Divide(2, 2);
}

} // namespace

void analyze(TString path, TString particle, Bool_t debug){    
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
    TString output_path = Form("%s/param/DCDRFT/e72/DCDriftParam_run%05d_%s.root", ANALYZER_DIR.Data(), run_num, particle.Data());
    TFile* fout = nullptr;
    if (!debug) {
        fout = TFile::Open(output_path.Data(), "UPDATE");
        if (!fout || fout->IsZombie()) {
            delete fout; // just in case
            fout = TFile::Open(output_path.Data(), "RECREATE");
            std::cout << "Create new file: " << output_path << std::endl;
        } else {
            std::cout << "Update existing file: " << output_path << std::endl;
        }
    } else {
        std::cout << "[DEBUG] Skipping DCDRFT write: " << output_path << std::endl;
    }


    // +-------------------+
    // | prepare histogram |
    // +-------------------+
    TH1D *h_blca_drift[conf.num_of_ch.at("blc")];
    TH1D *h_blcb_drift[conf.num_of_ch.at("blc")];
    TH2D *h2_blca_drift[conf.num_of_ch.at("blc")];
    TH2D *h2_blcb_drift[conf.num_of_ch.at("blc")];
    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++ ) {
        h_blca_drift[i] = (TH1D*)f->Get(Form("BLC%sa_Hit_DriftTime_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        h_blcb_drift[i] = (TH1D*)f->Get(Form("BLC%sb_Hit_DriftTime_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        h2_blca_drift[i] = (TH2D*)f->Get(Form("BLC%sa_Hit_DriftTime_vs_HitPat_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        h2_blcb_drift[i] = (TH2D*)f->Get(Form("BLC%sb_Hit_DriftTime_vs_HitPat_plane%d_%s", in_or_out.Data(), i, particle.Data()));
        if (!h2_blca_drift[i] || !h2_blcb_drift[i]) {
            std::cerr << "Error: missing Hit_DriftTime_vs_HitPat histogram for plane " << i
                      << " (BLC" << in_or_out << "a/b, particle=" << particle << ")" << std::endl;
            return;
        }
    }
 
    // +--------------+
    // | fit and plot |
    // +--------------+
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_BLC%s_DRIFT_%s.pdf", img_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());

    std::vector<TGraph*> drift_graphs;

    auto c_blc = new TCanvas("blc", "", 1500, 1200);
    c_blc->Divide(2, 2);
    c_blc->Print(pdf_path + "["); // start

    for (Int_t i = 0; i < conf.num_of_ch.at("blc"); i++) {
        TH2D* h2a = h2_blca_drift[i];
        TH2D* h2b = h2_blcb_drift[i];
        auto* xax = h2a->GetXaxis();

        Int_t bLo1  = xax->FindBin(-0.5);
        Int_t bHi1  = xax->FindBin(15.5);
        Int_t bLo2  = xax->FindBin(16.5);
        Int_t bHi2  = xax->FindBin(31.5);
        Int_t bLoAll = xax->FindBin(-0.5);
        Int_t bHiAll = xax->FindBin(31.5);

        TH1D* h_a_lo = h2a->ProjectionY(Form("BLC%sa_Hit_DriftTime_plane%d_w0-15_%s", in_or_out.Data(), i, particle.Data()),
                                        bLo1, bHi1);
        TH1D* h_a_hi = h2a->ProjectionY(Form("BLC%sa_Hit_DriftTime_plane%d_w16-31_%s", in_or_out.Data(), i, particle.Data()),
                                        bLo2, bHi2);
        TH1D* h_b_lo = h2b->ProjectionY(Form("BLC%sb_Hit_DriftTime_plane%d_w0-15_%s", in_or_out.Data(), i, particle.Data()),
                                        bLo1, bHi1);
        TH1D* h_b_hi = h2b->ProjectionY(Form("BLC%sb_Hit_DriftTime_plane%d_w16-31_%s", in_or_out.Data(), i, particle.Data()),
                                        bLo2, bHi2);
        TH1D* h_a_all = h2a->ProjectionY(Form("BLC%sa_Hit_DriftTime_plane%d_all_%s", in_or_out.Data(), i, particle.Data()),
                                          bLoAll, bHiAll);
        TH1D* h_b_all = h2b->ProjectionY(Form("BLC%sb_Hit_DriftTime_plane%d_all_%s", in_or_out.Data(), i, particle.Data()),
                                          bLoAll, bHiAll);

        conf.detector = Form("BLC%sa", in_or_out.Data());
        if (h_blca_drift[i]) {
            push_if(ana_helper::make_drift_function(h_blca_drift[i], nullptr, 0, i, nullptr), drift_graphs);
        }

        // PDF page (2x2): a side — same pad order as legacy [lo hist][lo graph] / [hi hist][hi graph]
        c_blc->Clear();
        c_blc->Divide(2, 2);
        make_halfwise_drift_page(h_a_lo, h_a_hi, h_a_all, c_blc, 1, 2, 3, 4, i, drift_graphs);
        print_drift_page(c_blc, pdf_path);

        conf.detector = Form("BLC%sb", in_or_out.Data());
        if (h_blcb_drift[i]) {
            push_if(ana_helper::make_drift_function(h_blcb_drift[i], nullptr, 0, i, nullptr), drift_graphs);
        }

        // PDF page (2x2): b side
        make_halfwise_drift_page(h_b_lo, h_b_hi, h_b_all, c_blc, 1, 2, 3, 4, i, drift_graphs);
        print_drift_page(c_blc, pdf_path);
    }
    c_blc->Print(pdf_path + "]"); // end
    delete c_blc;

    // +-------+
    // | Write |
    // +-------+
    if (debug || !fout) {
        return;
    }

    fout->cd();

    TDatime now;
    TString datetime = now.AsString();
    TNamed("datetime", datetime.Data()).Write("", TObject::kOverwrite);

    TString ref_key = Form("reference_%s", in_or_out.Data());
    TNamed(ref_key.Data(), path.Data()).Write("", TObject::kOverwrite);

    for (auto* g : drift_graphs) {
        if (!g) continue;
        g->Write("", TObject::kOverwrite);
    }
    fout->Close();
}

Int_t main(int argc, char** argv) {

    // -- check argments -----
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <root file path> <particle> [--debug]" << std::endl;
        return 1;
    }
    TString path = argv[1];
    TString particle = argv[2];
    if (particle != "Pi" && particle != "K") {
        std::cerr << "Error: Unexpected particle name: " << particle << std::endl;
        return 1;
    }

    Bool_t debug = kFALSE;
    for (Int_t i = 3; i < argc; ++i) {
        if (TString(argv[i]) == "--debug") {
            debug = kTRUE;
        } else {
            std::cerr << "Error: Unknown option: " << argv[i] << std::endl;
            std::cerr << "Usage: " << argv[0] << " <root file path> <particle> [--debug]" << std::endl;
            return 1;
        }
    }

    analyze(path, particle, debug);
    return 0;
}
