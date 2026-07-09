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
#include <TPad.h>

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

// One canvas pad: DriftTime hist (left) + DriftFunction graph (right).
TGraph* draw_drift_half(TH1D* h, TCanvas* c, Int_t pad, Int_t plane, const char* wire_range_suffix)
{
    if (!h || h->GetEntries() <= 0) return nullptr;

    TGraph* g = ana_helper::make_drift_function(h, nullptr, 0, plane, wire_range_suffix);
    if (!c) return g;

    c->cd(pad);
    TString tag = Form("plane%d%s", plane, wire_range_suffix);
    auto* pad_hist = new TPad(Form("hist_%s", tag.Data()), "", 0.02, 0.02, 0.48, 0.98);
    auto* pad_graph = new TPad(Form("graph_%s", tag.Data()), "", 0.52, 0.02, 0.98, 0.98);
    pad_hist->Draw();
    pad_graph->Draw();

    pad_hist->cd();
    h->Draw();

    if (g) {
        pad_graph->cd();
        g->Draw("APC");
    }
    return g;
}

void make_halfwise_drift(TH1D* h_lo, TH1D* h_hi, TH1D* h_all,
                         TCanvas* c, Int_t pad_lo, Int_t pad_hi,
                         Int_t plane, std::vector<TGraph*>& out)
{
    const bool lo_ok  = (h_lo  && h_lo->GetEntries() > 0);
    const bool hi_ok  = (h_hi  && h_hi->GetEntries() > 0);
    const bool all_ok = (h_all && h_all->GetEntries() > 0);

    TGraph* g_lo = nullptr;
    TGraph* g_hi = nullptr;

    if (lo_ok && hi_ok) {
        g_lo = draw_drift_half(h_lo, c, pad_lo, plane, "-w0-15");
        g_hi = draw_drift_half(h_hi, c, pad_hi, plane, "-w16-31");
    } else if (lo_ok && !hi_ok) {
        g_lo = draw_drift_half(h_lo, c, pad_lo, plane, "-w0-15");
        g_hi = clone_drift_graph(g_lo, plane, "-w16-31");
    } else if (!lo_ok && hi_ok) {
        g_hi = draw_drift_half(h_hi, c, pad_hi, plane, "-w16-31");
        g_lo = clone_drift_graph(g_hi, plane, "-w0-15");
    } else if (all_ok) {
        g_lo = draw_drift_half(h_all, c, pad_lo, plane, "-w0-15");
        g_hi = clone_drift_graph(g_lo, plane, "-w16-31");
    }

    push_if(g_lo, out);
    push_if(g_hi, out);
}

} // namespace

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
    TFile* fout = TFile::Open(output_path.Data(), "UPDATE");
    if (!fout || fout->IsZombie()) {
        delete fout; // just in case
        fout = TFile::Open(output_path.Data(), "RECREATE");
        std::cout << "Create new file: " << output_path << std::endl;
    } else {
        std::cout << "Update existing file: " << output_path << std::endl;
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
    Int_t nth_pad = 1;
    // 2x2 per plane: a_lo, a_hi, b_lo, b_hi (each pad: hist | graph inside)
    Int_t rows = 2, cols = 2;
    Int_t max_pads = rows * cols;
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_BLC%s_DRIFT_%s.pdf", img_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());

    std::vector<TGraph*> drift_graphs;

    auto c_blc = new TCanvas("blc", "", 1500, 1200);
    c_blc->Divide(cols, rows);
    c_blc->Print(pdf_path + "["); // start
    nth_pad = 1;

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

        // --- a side: plane-level fallback TGraph (param only) + halfwise (default) ---
        conf.detector = Form("BLC%sa", in_or_out.Data());
        if (h_blca_drift[i]) {
            push_if(ana_helper::make_drift_function(h_blca_drift[i], nullptr, 0, i, nullptr), drift_graphs);
        }
        make_halfwise_drift(h_a_lo, h_a_hi, h_a_all, c_blc, pad_a_lo, pad_a_hi, i, drift_graphs);

        // --- b side ---
        conf.detector = Form("BLC%sb", in_or_out.Data());
        if (h_blcb_drift[i]) {
            push_if(ana_helper::make_drift_function(h_blcb_drift[i], nullptr, 0, i, nullptr), drift_graphs);
        }
        make_halfwise_drift(h_b_lo, h_b_hi, h_b_all, c_blc, pad_b_lo, pad_b_hi, i, drift_graphs);
    }
    c_blc->Print(pdf_path);
    c_blc->Print(pdf_path + "]"); // end
    delete c_blc;

    // +-------+
    // | Write |
    // +-------+
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
