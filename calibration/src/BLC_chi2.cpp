// c++
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TMath.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"

Config& conf = Config::getInstance();

static Double_t nu_from_nhit(TH1D* h_nhit)
{
    if (!h_nhit || h_nhit->GetEntries() <= 0) return -1.0;
    const Int_t imax = h_nhit->GetMaximumBin();
    return h_nhit->GetBinCenter(imax) - 4.0;
}

// Mean DCGEO Res [mm] for BLC1 (ids 1-16) or BLC2 (1001-1016).
static Double_t mean_res_from_dcgeo(Int_t run_num, const TString& particle, Bool_t is_blc1)
{
    const TString path = Form("%s/param/DCGEO/e72/DCGeomParam_run%05d_%s",
                              ANALYZER_DIR.Data(), run_num, particle.Data());
    std::ifstream ifs(path.Data());
    if (!ifs) {
        std::cerr << "Warning: cannot open DCGEO " << path << std::endl;
        return -1.0;
    }
    const Int_t id_lo = is_blc1 ? 1 : 1001;
    const Int_t id_hi = is_blc1 ? 16 : 1016;
    Double_t sum = 0.0;
    Int_t n = 0;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::vector<std::string> cols;
        std::string tok;
        while (iss >> tok) cols.push_back(tok);
        if (cols.size() < 10) continue;
        Int_t id = 0;
        try {
            id = std::stoi(cols[0]);
        } catch (...) {
            continue;
        }
        if (id < id_lo || id > id_hi) continue;
        try {
            sum += std::stod(cols[9]);
            ++n;
        } catch (...) {
            continue;
        }
    }
    if (n == 0) return -1.0;
    return sum / static_cast<Double_t>(n);
}

// Chamber-level residual = sum of plane residuals; pull axis = residual / Res.
static TH1D* make_chamber_pull(TFile* f, const TString& in_or_out, const TString& particle,
                               Double_t res_mm, const char* hname)
{
    if (res_mm <= 0.0) return nullptr;
    const Int_t nplane = conf.num_of_ch.at("blc");
    TH1D* h_sum = nullptr;
    for (Int_t i = 0; i < nplane; ++i) {
        for (const char* side : {"a", "b"}) {
            TH1D* h = (TH1D*)f->Get(Form("BLC%s%s_Track_Residual_plane%d_%s",
                                         in_or_out.Data(), side, i, particle.Data()));
            if (!h || h->GetEntries() <= 0) continue;
            if (!h_sum) {
                h_sum = (TH1D*)h->Clone("tmp_residual_sum");
                h_sum->SetDirectory(nullptr);
                h_sum->Reset();
            }
            h_sum->Add(h);
        }
    }
    if (!h_sum || h_sum->GetEntries() <= 0) {
        delete h_sum;
        return nullptr;
    }

    const Int_t nb = h_sum->GetNbinsX();
    const Double_t xmin = h_sum->GetXaxis()->GetXmin() / res_mm;
    const Double_t xmax = h_sum->GetXaxis()->GetXmax() / res_mm;
    TH1D* h_pull = new TH1D(hname,
                            Form("%s chamber pull;pull;count", hname),
                            nb, xmin, xmax);
    h_pull->SetDirectory(nullptr);
    for (Int_t b = 1; b <= nb; ++b) {
        h_pull->SetBinContent(b, h_sum->GetBinContent(b));
        h_pull->SetBinError(b, h_sum->GetBinError(b));
    }
    delete h_sum;
    return h_pull;
}

void analyze(TString path, TString particle)
{
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0 / 256, 160.0 / 256, 44.0 / 256);

    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x");
    gStyle->SetTitleSize(0.06, "y");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    auto* f = new TFile(path.Data());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file : " << path << std::endl;
        return;
    }

    const char* treeName = nullptr;
    TString in_or_out = "";
    TString bc_tag = "";
    Bool_t is_blc1 = kFALSE;
    if (f->Get("bcout")) {
        treeName = "bcout";
        in_or_out = "2";
        bc_tag = "BcOut";
        is_blc1 = kFALSE;
    } else if (f->Get("bcin")) {
        treeName = "bcin";
        in_or_out = "1";
        bc_tag = "BcIn";
        is_blc1 = kTRUE;
    } else {
        std::cerr << "Error: neither 'bcout' nor 'bcin' tree exists.\n";
        return;
    }

    TTreeReader reader(treeName, f);
    TTreeReaderValue<unsigned int> run_number(reader, "run_number");
    reader.SetEntry(0);
    Int_t run_num = *run_number;

    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_path = Form("%s/run%05d_BLC%s_chi2_%s.root",
                               root_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    TH1D* h_chi2 = (TH1D*)f->Get(Form("%sTrack_ChiSquare_%s", bc_tag.Data(), particle.Data()));
    TH1D* h_nhit = (TH1D*)f->Get(Form("%sTrack_NHit_%s", bc_tag.Data(), particle.Data()));

    if (!h_chi2) {
        std::cerr << "Error: ChiSquare hist not found for " << bc_tag << " " << particle << std::endl;
        fout->Close();
        return;
    }

    Double_t nu_eff = nu_from_nhit(h_nhit);
    if (nu_eff <= 2.0) {
        std::cerr << "Warning: invalid nu_eff=" << nu_eff << ", using nu=8 as fallback\n";
        nu_eff = 8.0;
    }

    const Double_t res_mm = mean_res_from_dcgeo(run_num, particle, is_blc1);
    TH1D* h_pull = make_chamber_pull(f, in_or_out, particle, res_mm,
                                     Form("h_pull_%s", bc_tag.Data()));

    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_BLC%s_chi2_%s.pdf",
                            img_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());

    // --- page 1: reduced χ² (closeup | wide) + fit & ideal s=1 ---
    auto c_chi2 = new TCanvas("c_chi2", "", 1600, 700);
    c_chi2->Print(pdf_path + "[");
    FitResult scale_fit = ana_helper::chi2_scale_fit(h_chi2, c_chi2, 0, nu_eff);
    c_chi2->Print(pdf_path);
    delete c_chi2;

    const Double_t scale_s = scale_fit.par.size() > 1 ? scale_fit.par[1] : 1.0;
    const Double_t sqrt_s = scale_fit.additional.size() > 0 ? scale_fit.additional[0]
                                                           : TMath::Sqrt(scale_s);
    const Double_t mean_chi2 = scale_fit.additional.size() > 2 ? scale_fit.additional[2]
                                                              : h_chi2->GetMean();
    const Int_t fit_ok = scale_fit.migrad_stats;

    // --- page 2: BcIn/BcOut chamber-level pull ---
    Double_t pull_sigma = 0.0;
    Int_t n_pull = 0;
    auto c_pull = new TCanvas("c_pull", "", 1000, 800);
    if (h_pull && h_pull->GetEntries() > 0) {
        FitResult pull_fit = ana_helper::residual_fit(h_pull, c_pull, 0);
        if (pull_fit.par.size() > 2 && pull_fit.par[2] > 0.0) {
            pull_sigma = pull_fit.par[2];
            n_pull = 1;
        }
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.55, 0.55, Form("Res = %.4f mm", res_mm));
        latex.DrawLatex(0.55, 0.49, Form("pull #sigma = %.3f", pull_sigma));
        latex.DrawLatex(0.55, 0.43, "(target #approx 1)");
    } else {
        c_pull->cd();
        TLatex latex;
        latex.SetNDC();
        latex.DrawLatex(0.2, 0.5, Form("%s pull unavailable (Res or residual missing)", bc_tag.Data()));
    }
    c_pull->Print(pdf_path);
    c_pull->Print(pdf_path + "]");
    delete c_pull;

    std::cout << Form("[BLC_chi2] %s %s: s=%.4f sqrt(s)=%.4f nu=%.1f mean=%.4f "
                      "Res=%.4f pull_sigma=%.4f fit_ok=%d\n",
                      bc_tag.Data(), particle.Data(), scale_s, sqrt_s, nu_eff,
                      mean_chi2, res_mm, pull_sigma, fit_ok);

    TTree* tree = new TTree("tree", "BLC chi2 scale / chamber pull summary");
    Double_t b_scale_s = scale_s;
    Double_t b_sqrt_s = sqrt_s;
    Double_t b_nu_eff = nu_eff;
    Double_t b_mean_chi2 = mean_chi2;
    Double_t b_pull_sigma = pull_sigma;
    Double_t b_res_mm = res_mm;
    Int_t b_fit_ok = fit_ok;
    Int_t b_n_pull = n_pull;
    tree->Branch("scale_s", &b_scale_s, "scale_s/D");
    tree->Branch("sqrt_s", &b_sqrt_s, "sqrt_s/D");
    tree->Branch("nu_eff", &b_nu_eff, "nu_eff/D");
    tree->Branch("mean_chi2", &b_mean_chi2, "mean_chi2/D");
    tree->Branch("pull_sigma", &b_pull_sigma, "pull_sigma/D");
    tree->Branch("res_mm", &b_res_mm, "res_mm/D");
    tree->Branch("fit_ok", &b_fit_ok, "fit_ok/I");
    tree->Branch("n_pull", &b_n_pull, "n_pull/I");
    tree->Fill();

    fout->cd();
    tree->Write();
    if (h_chi2) h_chi2->Write("h_chi2");
    if (h_nhit) h_nhit->Write("h_nhit");
    if (h_pull) h_pull->Write("h_pull");
    fout->Close();
    delete h_pull;

    std::cout << "Wrote " << output_path << std::endl;
    std::cout << "Wrote " << pdf_path << std::endl;
}

Int_t main(int argc, char** argv)
{
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
