#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    FitResult t0_adc_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        std::vector<Double_t> par, err;

        // -- pedestal -----
        Double_t ped_pos        = h->GetBinCenter(h->GetMaximumBin());
        Double_t ped_half_width = 5.0;
        std::pair<Double_t, Double_t> ped_n_sigma(2.0, 2.0);

        // -- first fit -----
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", ped_pos-ped_half_width, ped_pos+ped_half_width);
        f_prefit->SetParameter(1, ped_pos);
        f_prefit->SetParameter(2, ped_half_width);
        h->Fit(f_prefit, "0Q", "", ped_pos-ped_half_width, ped_pos+ped_half_width);
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));

        // -- second fit -----
        TF1 *f_fit_ped = new TF1( Form("ped_%s", h->GetName()), "gausn", par[1]-ped_n_sigma.first*par[2], par[1]+ped_n_sigma.second*par[2]);
        f_fit_ped->SetParameter(0, par[0]);
        f_fit_ped->SetParameter(1, par[1]);
        f_fit_ped->SetParameter(2, par[2]*0.9);
        f_fit_ped->SetLineColor(kOrange);
        f_fit_ped->SetLineWidth(2);
        f_fit_ped->SetNpx(1000);
        h->Fit(f_fit_ped, "0Q", "", par[1]-ped_n_sigma.first*par[2], par[1]+ped_n_sigma.second*par[2]);

        FitResult result;
        for (Int_t i = 0, n_par = f_fit_ped->GetNpar(); i < n_par; i++) {
            result.par.push_back(f_fit_ped->GetParameter(i));
            result.err.push_back(f_fit_ped->GetParError(i));
        }


        // -- mip -----
        par.clear(); err.clear();
        Double_t ped_mip_distance = 80.0;
        Double_t mip_pos          = ped_pos + ped_mip_distance;
        Double_t mip_half_width   = 20.0;
        std::pair<Double_t, Double_t> mip_n_sigma(1.5, 2.0);

        // -- first fit -----
        f_prefit->SetRange(mip_pos-mip_half_width, mip_pos+mip_half_width);
        f_prefit->SetParameter(1, mip_pos);
        f_prefit->SetParameter(2, mip_half_width);
        h->Fit(f_prefit, "0Q", "", mip_pos-mip_half_width, mip_pos+mip_half_width);
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
        delete f_prefit;

        // -- second fit -----
        TF1 *f_fit_mip_g = new TF1( Form("mip_gauss_%s", h->GetName()), "gausn", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        f_fit_mip_g->SetParameter(1, par[1]);
        f_fit_mip_g->SetParameter(2, par[2]*0.9);
        f_fit_mip_g->SetLineColor(kOrange);
        f_fit_mip_g->SetLineWidth(2.0);
        h->Fit(f_fit_mip_g, "0Q", "", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        Double_t p_value_g = TMath::Prob(f_fit_mip_g->GetChisquare(), f_fit_mip_g->GetNDF());

        TF1 *f_fit_mip_l = new TF1( Form("mip_landau_%s", h->GetName()), "landaun", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        f_fit_mip_l->SetParameter(1, par[1]);
        f_fit_mip_l->SetParameter(2, par[2]*0.9);
        f_fit_mip_l->SetLineColor(kOrange);
        f_fit_mip_l->SetLineWidth(2.0);
        h->Fit(f_fit_mip_l, "0Q", "", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        Double_t p_value_l = TMath::Prob(f_fit_mip_l->GetChisquare(), f_fit_mip_l->GetNDF());

        if (p_value_g > p_value_l) {
            for (Int_t i = 0, n_par = f_fit_mip_g->GetNpar(); i < n_par; i++) {
                result.par.push_back(f_fit_mip_g->GetParameter(i));
                result.err.push_back(f_fit_mip_g->GetParError(i));
            }

            // -- draw -----
            h->GetXaxis()->SetRangeUser(
                result.par[1]-10.0*result.par[2], 
                result.par[4]+ 5.0*result.par[5]
            );
            h->Draw();
            f_fit_mip_g->SetNpx(1000);
            f_fit_mip_g->Draw("same");
            delete f_fit_mip_l;
        } else {
            for (Int_t i = 0, n_par = f_fit_mip_l->GetNpar(); i < n_par; i++) {
                result.par.push_back(f_fit_mip_l->GetParameter(i));
                result.err.push_back(f_fit_mip_l->GetParError(i));
            }

            // -- draw -----
            h->GetXaxis()->SetRangeUser(
                result.par[1]-10.0*result.par[2], 
                result.par[4]+ 5.0*result.par[5]
            );
            h->Draw();
            f_fit_mip_l->SetNpx(1000);
            f_fit_mip_l->Draw("same");
            delete f_fit_mip_g;
        }
        f_fit_ped->Draw("same");
        TLine *line = new TLine(result.par[4], 0, result.par[4], h->GetMaximum());
        line->SetLineStyle(2); // 点線に設定
        line->SetLineColor(kRed); // 色を赤に設定
        c->Update();

        return result;
    }

    // // ____________________________________________________________________________________________
    // FitResult npe_gauss_fit(TH1D *h, TCanvas *c, Int_t n_c, Double_t n_sigma, Double_t cutoff_threshold) { // cutoff_thresholdは主にKVCの調整用
    //     c->cd(n_c);
    //     Config& conf = Config::getInstance();
    //     gPad->SetLogy(conf.log_flag);

    //     if (h->GetEntries() < 10) {
    //         h->Draw();
    //         FitResult result;
    //         return result;
    //     }

    //     std::vector<Double_t> par, err;
    //     h->GetXaxis()->SetRangeUser(cutoff_threshold, conf.npe_max);
    //     if ( h->GetMaximum() < 10) {
    //         h->Draw();
    //         FitResult result;
    //         return result;
    //     }
    //     Double_t peak_pos = h->GetMean();
    //     Double_t stdev = h->GetStdDev();

    //     // -- first fit -----
    //     Int_t n_iter = 3;
    //     par.insert(par.end(), {0.0, peak_pos, 2.0*stdev});
    //     for (Int_t dummy = 0; dummy < n_iter; dummy++) {
    //         Double_t fit_range_min = par[1]-n_sigma*par[2] > 5.0 ? par[1]-n_sigma*par[2] : 5.0;
    //         TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", fit_range_min, par[1]+n_sigma*par[2]);
    //         f_prefit->SetParameter(1, par[1]);
    //         f_prefit->SetParameter(2, par[2]*0.5);
    //         h->Fit(f_prefit, "0Q", "", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
    //         par.clear();
    //         for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
    //         delete f_prefit;            
    //     }

    //     Double_t fit_range_min = par[1]-n_sigma*par[2];
    //     Double_t fit_range_max = par[1]+n_sigma*par[2];
    //     TF1 *f_fit = new TF1(Form("gauss_%s", h->GetName()), "gausn", fit_range_min, fit_range_max);
    //     f_fit->SetParameters(par[0], par[1], par[2]);
    //     TString fit_option = "0Q";
    //     if ( h->GetBinContent( h->GetMaximumBin() ) < 50 ) fit_option += 'L';
    //     h->Fit(f_fit, fit_option.Data(), "", fit_range_min, fit_range_max);

    //     FitResult result;
    //     par.clear();
    //     Int_t n_par = f_fit->GetNpar();
    //     for (Int_t i = 0; i < n_par; i++) {
    //         par.push_back(f_fit->GetParameter(i));
    //         err.push_back(f_fit->GetParError(i));
    //     }
    //     Double_t chi2 = f_fit->GetChisquare();
    //     Double_t ndf  = f_fit->GetNDF();
    //     result.par = par;
    //     result.err = err;
    //     result.reduced_chi2 = (Double_t) chi2/ndf;

    //     // -- draw -----
    //     h->GetXaxis()->SetRangeUser(-3.0, result.par[1]+3.0*stdev);
    //     h->Draw();
    //     f_fit->SetLineColor(kOrange);
    //     f_fit->SetLineWidth(2);
    //     f_fit->SetNpx(1000);
    //     f_fit->Draw("same");

    //     auto *text = new TLatex();
    //     text->SetNDC();
    //     text->SetTextSize(0.1);
    //     text->DrawLatex(0.5, 0.8, Form("#lambda = %.2f", result.par[1]));
    //     text->Draw("same");

    //     return result;
    // }


    // // ____________________________________________________________________________________________
    // FitResult threshold_erf_fit(TH1D *h, TCanvas *c, Int_t n_c) {
    //     Config& conf = Config::getInstance();
    //     c->cd(n_c);
    //     std::vector<Double_t> par, err;

    //     // -- prepare smooth hist -----
    //     TH1D *h_clone = (TH1D*)h->Clone(Form("%s_clone", h->GetName()));
    //     h_clone->Smooth(10);

    //     // -- erf fit -----
    //     Double_t fit_range_min = conf.threshold_fit_range_min;
    //     Double_t fit_range_max = conf.threshold_fit_range_max;
    //     h_clone->GetXaxis()->SetRangeUser(fit_range_min, fit_range_max);
    //     TF1 *f_fit_erf = new TF1( Form("erf_fit_%s", h->GetName()), "[0]*TMath::Erf( (x-[1])/[2] ) + 0.5", fit_range_min, fit_range_max);
    //     f_fit_erf->FixParameter(0, 0.5);
    //     f_fit_erf->SetParameter(1, h_clone->GetBinCenter(h_clone->GetMaximumBin()) - 10.0);
    //     f_fit_erf->SetParameter(2, 10.0);
    //     f_fit_erf->SetLineColor(kOrange);
    //     f_fit_erf->SetLineWidth(2);
    //     f_fit_erf->SetNpx(1000);
    //     h_clone->Fit(f_fit_erf, "0", "", fit_range_min, fit_range_max);

    //     FitResult result;
    //     par.clear();
    //     Int_t n_par = f_fit_erf->GetNpar();
    //     for (Int_t i = 0; i < n_par; i++) {
    //         par.push_back(f_fit_erf->GetParameter(i));
    //         err.push_back(f_fit_erf->GetParError(i));
    //     }
    //     Double_t chi2 = f_fit_erf->GetChisquare();
    //     Double_t ndf  = f_fit_erf->GetNDF();
    //     result.par = par;
    //     result.err = err;
    //     result.reduced_chi2 = (Double_t) chi2/ndf;


    //     h->GetXaxis()->SetRangeUser(result.par[1] - 5.0*result.par[2], result.par[1] + 5.0*result.par[2]);
    //     h->Draw();
    //     h_clone->SetLineColor(kGreen);
    //     h_clone->Draw("same");
    //     f_fit_erf->Draw("same");

    //     TLine *line = new TLine(result.par[1], 0, result.par[1], 1.0);
    //     line->SetLineStyle(2);
    //     line->SetLineWidth(2);
    //     line->SetLineColor(kRed);
    //     line->Draw("same");

    //     // delete h_clone;
    //     return result;
    // }


}