// c++
#include <iostream>
#include <fstream>
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
    if (f->Get("bcout")) {
        treeName = "bcout";
        in_or_out = "2";
        bc_tag = "BcOut";
    } else if (f->Get("bcin")) {
        treeName = "bcin";
        in_or_out = "1";
        bc_tag = "BcIn";
    } else {
        std::cerr << "Error: neither 'bcout' nor 'bcin' tree exists.\n";
        return;
    }

    TTreeReader reader(treeName, f);
    TTreeReaderValue<unsigned int> run_number(reader, "run_number");
    reader.SetEntry(0);
    Int_t run_num = *run_number;

    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString output_path = Form("%s/run%05d_BLC%s_pull_%s.root",
                               root_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());
    if (std::ifstream(output_path.Data())) std::remove(output_path.Data());
    TFile* fout = new TFile(output_path.Data(), "RECREATE");

    const Int_t nplane = conf.num_of_ch.at("blc");
    TH1D* h_blca[nplane];
    TH1D* h_blcb[nplane];
    Int_t n_found = 0;
    for (Int_t i = 0; i < nplane; ++i) {
        h_blca[i] = (TH1D*)f->Get(Form("BLC%sa_Track_Pull_plane%d_%s",
                                       in_or_out.Data(), i, particle.Data()));
        h_blcb[i] = (TH1D*)f->Get(Form("BLC%sb_Track_Pull_plane%d_%s",
                                       in_or_out.Data(), i, particle.Data()));
        if (h_blca[i] && h_blca[i]->GetEntries() > 0) ++n_found;
        if (h_blcb[i] && h_blcb[i]->GetEntries() > 0) ++n_found;
    }
    if (n_found == 0) {
        std::cerr << "Error: no Track_Pull histograms found. "
                  << "Rebuild UserBc with #define FillPullExclusive 1.\n";
        fout->Close();
        return;
    }

    TH1D* h_sum = nullptr;
    for (Int_t i = 0; i < nplane; ++i) {
        for (TH1D* h : {h_blca[i], h_blcb[i]}) {
            if (!h || h->GetEntries() <= 0) continue;
            if (!h_sum) {
                h_sum = (TH1D*)h->Clone(Form("h_pull_sum_%s", bc_tag.Data()));
                h_sum->SetDirectory(nullptr);
                h_sum->Reset();
            }
            h_sum->Add(h);
        }
    }
    if (!h_sum || h_sum->GetEntries() <= 0) {
        std::cerr << "Error: empty Track_Pull sum histogram.\n";
        fout->Close();
        return;
    }
    h_sum->SetTitle(Form("%s exclusive pull;pull = (s_{hit} - s_{track}^{w/o hit}) / Res;count",
                         bc_tag.Data()));

    TH1D* h_chi2 = (TH1D*)f->Get(Form("%sTrack_ChiSquare_%s", bc_tag.Data(), particle.Data()));
    TH1D* h_nhit = (TH1D*)f->Get(Form("%sTrack_NHit_%s", bc_tag.Data(), particle.Data()));

    Double_t nu_eff = nu_from_nhit(h_nhit);
    if (nu_eff <= 2.0) {
        std::cerr << "Warning: invalid nu_eff=" << nu_eff << ", using nu=8 as fallback\n";
        nu_eff = 8.0;
    }

    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_BLC%s_pull_%s.pdf",
                            img_base_dir.Data(), run_num, in_or_out.Data(), particle.Data());

    // --- page 1: exclusive pull (update input) ---
    auto c_pull = new TCanvas("c_pull", "", 1000, 800);
    c_pull->Print(pdf_path + "[");
    FitResult fr = ana_helper::residual_fit(h_sum, c_pull, 0);
    Double_t pull_sigma = (fr.par.size() > 2) ? fr.par[2] : 0.0;
    Double_t pull_mean = (fr.par.size() > 1) ? fr.par[1] : 0.0;
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.55, 0.55, Form("pull #sigma = %.4f", pull_sigma));
    latex.DrawLatex(0.55, 0.49, Form("pull mean = %.4f", pull_mean));
    latex.DrawLatex(0.55, 0.43, "(target #sigma #approx 1)");
    latex.DrawLatex(0.55, 0.37, "(used for Res update)");
    c_pull->Print(pdf_path);
    delete c_pull;

    // --- page 2: same reduced-χ² scale analysis as BLC_chi2 (monitor only) ---
    Double_t scale_s = 1.0;
    Double_t sqrt_s = 1.0;
    Double_t mean_chi2 = -1.0;
    Int_t fit_ok = 0;
    auto c_chi2 = new TCanvas("c_chi2", "", 1600, 700);
    if (h_chi2 && h_chi2->GetEntries() > 0) {
        FitResult scale_fit = ana_helper::chi2_scale_fit(h_chi2, c_chi2, 0, nu_eff);
        scale_s = scale_fit.par.size() > 1 ? scale_fit.par[1] : 1.0;
        sqrt_s = scale_fit.additional.size() > 0 ? scale_fit.additional[0]
                                                 : TMath::Sqrt(scale_s);
        mean_chi2 = scale_fit.additional.size() > 2 ? scale_fit.additional[2]
                                                    : h_chi2->GetMean();
        fit_ok = scale_fit.migrad_stats;
        TLatex latex2;
        latex2.SetNDC();
        latex2.SetTextSize(0.035);
        latex2.DrawLatex(0.15, 0.92, "(#chi^{2} scale: monitor only, not used for Res update)");
    } else {
        c_chi2->cd();
        TLatex latex2;
        latex2.SetNDC();
        latex2.DrawLatex(0.2, 0.5, Form("%s ChiSquare hist unavailable", bc_tag.Data()));
        if (h_chi2) mean_chi2 = h_chi2->GetMean();
    }
    c_chi2->Print(pdf_path);
    c_chi2->Print(pdf_path + "]");
    delete c_chi2;

    std::cout << Form("[BLC_pull] %s %s: pull_sigma=%.4f pull_mean=%.4f "
                      "s=%.4f sqrt(s)=%.4f mean_chi2=%.4f nu=%.1f fit_ok=%d nplane_hist=%d\n",
                      bc_tag.Data(), particle.Data(), pull_sigma, pull_mean,
                      scale_s, sqrt_s, mean_chi2, nu_eff, fit_ok, n_found);

    TTree* tree = new TTree("tree", "BLC exclusive pull + chi2 monitor");
    Double_t b_pull_sigma = pull_sigma;
    Double_t b_pull_mean = pull_mean;
    Double_t b_scale_s = scale_s;
    Double_t b_sqrt_s = sqrt_s;
    Double_t b_nu_eff = nu_eff;
    Double_t b_mean_chi2 = mean_chi2;
    Int_t b_fit_ok = fit_ok;
    Int_t b_n_hist = n_found;
    tree->Branch("pull_sigma", &b_pull_sigma, "pull_sigma/D");
    tree->Branch("pull_mean", &b_pull_mean, "pull_mean/D");
    tree->Branch("scale_s", &b_scale_s, "scale_s/D");
    tree->Branch("sqrt_s", &b_sqrt_s, "sqrt_s/D");
    tree->Branch("nu_eff", &b_nu_eff, "nu_eff/D");
    tree->Branch("mean_chi2", &b_mean_chi2, "mean_chi2/D");
    tree->Branch("fit_ok", &b_fit_ok, "fit_ok/I");
    tree->Branch("n_hist", &b_n_hist, "n_hist/I");
    tree->Fill();

    fout->cd();
    tree->Write();
    h_sum->Write("h_pull");
    if (h_chi2) h_chi2->Write("h_chi2");
    if (h_nhit) h_nhit->Write("h_nhit");
    fout->Close();
    delete h_sum;

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
