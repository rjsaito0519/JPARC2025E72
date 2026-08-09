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

static const char* kPlaneLabels[] = {
    "U1", "UP1", "V1", "VP1", "U2", "UP2", "V2", "VP2"
};

struct PlanePull {
    bool ok = false;
    Double_t mean = 0.0;
    Double_t sigma = 0.0;
    Long64_t entries = 0;
};

static Double_t nu_from_nhit(TH1D* h_nhit)
{
    if (!h_nhit || h_nhit->GetEntries() <= 0) return -1.0;
    const Int_t imax = h_nhit->GetMaximumBin();
    return h_nhit->GetBinCenter(imax) - 4.0;
}

static PlanePull
fit_plane_pull(TH1D* h, TCanvas* c, Int_t pad)
{
    PlanePull out;
    if (!h || h->GetEntries() < 20) {
        c->cd(pad);
        TLatex tex;
        tex.SetNDC();
        tex.SetTextSize(0.06);
        tex.DrawLatex(0.2, 0.5, h ? "empty / low stats" : "NOT FOUND");
        return out;
    }
    out.entries = static_cast<Long64_t>(h->GetEntries());
    FitResult fr = ana_helper::residual_fit(h, c, pad);
    if (fr.par.size() > 2) {
        out.ok = true;
        out.mean = fr.par[1];
        out.sigma = fr.par[2];
    }
    c->cd(pad);
    TLatex tex;
    tex.SetNDC();
    tex.SetTextAlign(13);
    tex.SetTextSize(0.05);
    const Int_t ipl = pad - 1;
    if (ipl >= 0 && ipl < 8)
        tex.DrawLatex(0.15, 0.88, Form("%s", kPlaneLabels[ipl]));
    if (out.ok) {
        tex.SetTextSize(0.045);
        tex.DrawLatex(0.55, 0.88, Form("#mu=%.3f  #sigma=%.3f", out.mean, out.sigma));
    }
    return out;
}

static void
draw_plane_pull_page(TCanvas* c, const TString& pdf_path, const char* det_name,
                     TH1D** h_plane, Int_t nplane, std::vector<PlanePull>& results)
{
    c->Clear();
    c->Divide(4, 2);
    results.clear();
    results.resize(nplane);
    for (Int_t i = 0; i < nplane && i < 8; ++i) {
        results[i] = fit_plane_pull(h_plane[i], c, i + 1);
    }
    // title strip via latex on pad 1 overlay note
    c->cd();
    TLatex title;
    title.SetNDC();
    title.SetTextAlign(22);
    title.SetTextSize(0.03);
    title.DrawLatex(0.5, 0.985, Form("%s per-plane exclusive pull (grouping check)", det_name));
    c->Print(pdf_path);
}

static void
draw_plane_summary_page(TCanvas* c, const TString& pdf_path, const TString& bc_tag,
                        const char* name_a, const char* name_b,
                        const std::vector<PlanePull>& pa,
                        const std::vector<PlanePull>& pb,
                        Double_t sum_sigma, Double_t sum_mean)
{
    c->Clear();
    c->cd();
    TLatex tex;
    tex.SetNDC();
    tex.SetTextAlign(12);
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.08, 0.92, Form("%s per-plane pull summary (for Res grouping)", bc_tag.Data()));
    tex.SetTextSize(0.028);
    tex.DrawLatex(0.08, 0.86, Form("a+b sum (used for Res update now): #mu=%.4f  #sigma=%.4f",
                                   sum_mean, sum_sigma));
    tex.DrawLatex(0.08, 0.81, "Compare a vs b and plane-to-plane #sigma (~1) before choosing grouping.");

    tex.SetTextSize(0.032);
    tex.DrawLatex(0.08, 0.74, Form("%-6s", name_a));
    tex.DrawLatex(0.52, 0.74, Form("%-6s", name_b));
    tex.SetTextSize(0.026);
    tex.DrawLatex(0.08, 0.69, "plane    #mu      #sigma   n");
    tex.DrawLatex(0.52, 0.69, "plane    #mu      #sigma   n");

    Double_t y = 0.64;
    for (Int_t i = 0; i < 8; ++i) {
        const char* lab = (i < 8) ? kPlaneLabels[i] : Form("%d", i);
        if (i < static_cast<Int_t>(pa.size()) && pa[i].ok)
            tex.DrawLatex(0.08, y, Form("%-4s  %+6.3f  %6.3f  %lld",
                                        lab, pa[i].mean, pa[i].sigma,
                                        static_cast<long long>(pa[i].entries)));
        else
            tex.DrawLatex(0.08, y, Form("%-4s   n/a", lab));

        if (i < static_cast<Int_t>(pb.size()) && pb[i].ok)
            tex.DrawLatex(0.52, y, Form("%-4s  %+6.3f  %6.3f  %lld",
                                        lab, pb[i].mean, pb[i].sigma,
                                        static_cast<long long>(pb[i].entries)));
        else
            tex.DrawLatex(0.52, y, Form("%-4s   n/a", lab));
        y -= 0.055;
    }

    tex.SetTextSize(0.024);
    tex.DrawLatex(0.08, 0.12, "Grouping hint: similar #sigma within a block -> share Res; large a/b or U/V split -> finer Res.");
    tex.DrawLatex(0.08, 0.07, "Current update still uses a+b sum #sigma only (this page is diagnostic).");
    c->Print(pdf_path);
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

    const TString name_a = Form("BLC%sa", in_or_out.Data());
    const TString name_b = Form("BLC%sb", in_or_out.Data());

    // --- page 1: exclusive pull a+b sum (update input) ---
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
    latex.DrawLatex(0.55, 0.37, "(a+b sum; used for Res update)");
    c_pull->Print(pdf_path);
    delete c_pull;

    // --- pages: per-plane pull (a), (b) ---
    auto c_plane = new TCanvas("c_plane_pull", "", 1600, 1000);
    std::vector<PlanePull> pull_a;
    std::vector<PlanePull> pull_b;
    draw_plane_pull_page(c_plane, pdf_path, name_a.Data(), h_blca, nplane, pull_a);
    draw_plane_pull_page(c_plane, pdf_path, name_b.Data(), h_blcb, nplane, pull_b);

    // --- summary table ---
    draw_plane_summary_page(c_plane, pdf_path, bc_tag, name_a.Data(), name_b.Data(),
                            pull_a, pull_b, pull_sigma, pull_mean);
    delete c_plane;

    std::cout << Form("[BLC_pull] %s per-plane pull #sigma:\n", bc_tag.Data());
    for (Int_t i = 0; i < nplane && i < 8; ++i) {
        if (i < static_cast<Int_t>(pull_a.size()) && pull_a[i].ok)
            std::cout << Form("  %s-%s: mean=%+.4f sigma=%.4f\n",
                              name_a.Data(), kPlaneLabels[i],
                              pull_a[i].mean, pull_a[i].sigma);
        else
            std::cout << Form("  %s-%s: n/a\n", name_a.Data(), kPlaneLabels[i]);
        if (i < static_cast<Int_t>(pull_b.size()) && pull_b[i].ok)
            std::cout << Form("  %s-%s: mean=%+.4f sigma=%.4f\n",
                              name_b.Data(), kPlaneLabels[i],
                              pull_b[i].mean, pull_b[i].sigma);
        else
            std::cout << Form("  %s-%s: n/a\n", name_b.Data(), kPlaneLabels[i]);
    }

    // --- last page: same reduced-χ² scale analysis as BLC_chi2 (monitor only) ---
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
