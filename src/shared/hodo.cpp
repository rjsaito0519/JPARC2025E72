#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    void set_tdc_search_range(TH1D *h) {
        Config& conf = Config::getInstance();

        Double_t peak_pos = h->GetBinCenter(h->GetMaximumBin());
        Double_t stdev    = h->GetStdDev();
        std::pair<Double_t, Double_t> peak_n_sigma(5.0, 5.0);
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
        f_prefit->SetParameter(1, peak_pos);
        f_prefit->SetParameter(2, stdev*0.5);
        h->Fit(f_prefit, "0QEMR", "", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
        Double_t fit_center = f_prefit->GetParameter(1);
        Double_t fit_sigma  = f_prefit->GetParameter(2);
 
        Double_t range_n_sigma = 5.0;
        conf.tdc_search_range[conf.detector.Data()].first  = fit_center - range_n_sigma*fit_sigma;
        conf.tdc_search_range[conf.detector.Data()].second = fit_center + range_n_sigma*fit_sigma;
        
        delete f_prefit; 
    }

    // ____________________________________________________________________________________________
    FitResult tdc_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        gPad->SetLogy(1);
        std::vector<Double_t> par, err;
        TString fit_option = h->GetMaximum() > 500.0 ? "0QEMR" : "0QEMRL";
        FitResult result;

        h->GetXaxis()->SetRangeUser(
            conf.tdc_search_range[conf.detector.Data()].first, 
            conf.tdc_search_range[conf.detector.Data()].second
        );

        Double_t peak_pos = h->GetBinCenter(h->GetMaximumBin());
        Double_t stdev    = h->GetStdDev();
        std::pair<Double_t, Double_t> peak_n_sigma(2.0, 2.0);

        Double_t evnum_within_range = h->Integral(
            h->FindBin(conf.tdc_search_range[conf.detector.Data()].first),
            h->FindBin(conf.tdc_search_range[conf.detector.Data()].second)
        );
        
        if (evnum_within_range > 300.0) {
            // -- first fit -----
            TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
            f_prefit->SetParameter(1, peak_pos);
            f_prefit->SetParameter(2, stdev*0.9);
            h->Fit(f_prefit, fit_option.Data(), "", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;

            // -- second fit -----
            TF1 *f_fit = new TF1( Form("gaus_%s", h->GetName()), "gausn", par[1]-peak_n_sigma.first*par[2], par[1]+peak_n_sigma.second*par[2]);
            f_fit->SetParameter(0, par[0]);
            f_fit->SetParameter(1, par[1]);
            f_fit->SetParameter(2, par[2]*0.9);
            f_fit->SetLineColor(kOrange);
            f_fit->SetLineWidth(2);
            f_fit->SetNpx(1000);
            h->Fit(f_fit, fit_option.Data(), "", par[1]-peak_n_sigma.first*par[2], par[1]+peak_n_sigma.second*par[2]);

            // -- fill result -----
            for (Int_t i = 0, n_par = f_fit->GetNpar(); i < n_par; i++) {
                result.par.push_back(f_fit->GetParameter(i));
                result.err.push_back(f_fit->GetParError(i));
            }

            // -- draw figure -----
            h->GetXaxis()->SetRangeUser(
                result.par[1] - 15.0*result.par[2], 
                result.par[1] + 15.0*result.par[2]
            );
            h->Draw();
            f_fit->Draw("same");
        } else {
            Double_t err_value = 9999.0;
            // -- fill result -----
            result.par.push_back(evnum_within_range);
            result.err.push_back(err_value);

            result.par.push_back(h->GetMean());
            result.err.push_back(h->GetMeanError());
            
            result.par.push_back(h->GetStdDev());
            result.err.push_back(h->GetStdDevError());

            // -- draw figure -----
            h->GetXaxis()->SetRangeUser(
                result.par[1] - 15.0*result.par[2], 
                result.par[1] + 15.0*result.par[2]
            );
            h->Draw();
            TLatex* commnet = new TLatex();
            commnet->SetNDC();  // NDC座標（0〜1の正規化）を使う
            commnet->SetTextSize(0.04);
            commnet->DrawLatex(0.7, 0.85, "not fitting");
        }

        TLine *line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        line->SetLineStyle(2); // 点線に設定
        line->SetLineColor(kRed); // 色を赤に設定
        line->Draw("same");

        c->Update();

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult t0_adc_fit(TH1D *h, TCanvas *c, Int_t n_c, Double_t ped_mip_distance) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        gPad->SetLogy(1);
        std::vector<Double_t> par, err;

        // -- pedestal -----
        Double_t ped_pos        = h->GetBinCenter(h->GetMaximumBin());
        Double_t ped_half_width = 5.0;
        std::pair<Double_t, Double_t> ped_n_sigma(2.0, 2.0);

        // -- first fit -----
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", ped_pos-ped_half_width, ped_pos+ped_half_width);
        f_prefit->SetParameter(1, ped_pos);
        f_prefit->SetParameter(2, ped_half_width);
        h->Fit(f_prefit, "0QEMR", "", ped_pos-ped_half_width, ped_pos+ped_half_width);
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));

        // -- second fit -----
        TF1 *f_fit_ped = new TF1( Form("ped_%s", h->GetName()), "gausn", par[1]-ped_n_sigma.first*par[2], par[1]+ped_n_sigma.second*par[2]);
        f_fit_ped->SetParameter(0, par[0]);
        f_fit_ped->SetParameter(1, par[1]);
        f_fit_ped->SetParameter(2, par[2]*0.9);
        f_fit_ped->SetLineColor(kOrange);
        f_fit_ped->SetLineWidth(2);
        f_fit_ped->SetNpx(1000);
        h->Fit(f_fit_ped, "0QEMR", "", par[1]-ped_n_sigma.first*par[2], par[1]+ped_n_sigma.second*par[2]);

        FitResult result;
        for (Int_t i = 0, n_par = f_fit_ped->GetNpar(); i < n_par; i++) {
            result.par.push_back(f_fit_ped->GetParameter(i));
            result.err.push_back(f_fit_ped->GetParError(i));
        }


        // -- mip -----
        par.clear(); err.clear();
        Double_t mip_pos          = ped_pos + ped_mip_distance;
        Double_t mip_half_width   = 20.0;
        std::pair<Double_t, Double_t> mip_n_sigma(1.6, 2.0);

        // -- first fit -----
        f_prefit->SetRange(mip_pos-mip_half_width, mip_pos+mip_half_width);
        f_prefit->SetParameter(1, mip_pos);
        f_prefit->SetParameter(2, mip_half_width*0.9);
        h->Fit(f_prefit, "0QEMR", "", mip_pos-mip_half_width, mip_pos+mip_half_width);
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
        delete f_prefit;

        // -- second fit -----
        TF1 *f_fit_mip_g = new TF1( Form("mip_gauss_%s", h->GetName()), "gausn", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        f_fit_mip_g->SetParameter(1, par[1]);
        f_fit_mip_g->SetParameter(2, par[2]*0.9);
        f_fit_mip_g->SetLineColor(kOrange);
        f_fit_mip_g->SetLineWidth(2.0);
        h->Fit(f_fit_mip_g, "0QEMR", "", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        Double_t chi_square_g = f_fit_mip_g->GetChisquare();
        Double_t p_value_g = TMath::Prob(chi_square_g, f_fit_mip_g->GetNDF());

        TF1 *f_fit_mip_l = new TF1( Form("mip_landau_%s", h->GetName()), "landaun", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        f_fit_mip_l->SetParameter(1, par[1]);
        f_fit_mip_l->SetParameter(2, par[2]*0.9);
        f_fit_mip_l->SetLineColor(kOrange);
        f_fit_mip_l->SetLineWidth(2.0);
        h->Fit(f_fit_mip_l, "0QEMR", "", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
        Double_t chi_square_l = f_fit_mip_l->GetChisquare();
        Double_t p_value_l = TMath::Prob(chi_square_l, f_fit_mip_l->GetNDF());

        // // -- debug ------
        // std::cout << "gauss:  " << p_value_g << ", " << f_fit_mip_g->GetChisquare() << std::endl;
        // std::cout << "landau: " << p_value_l << ", " << f_fit_mip_l->GetChisquare() << std::endl;
        // Bool_t flag = p_value_g >= p_value_l;
        // std::cout << flag << std::endl;

        // if (p_value_g >= p_value_l) {
        if (chi_square_g <= chi_square_l) {
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

        TLine *ped_line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        ped_line->SetLineStyle(2); // 点線に設定
        ped_line->SetLineColor(kRed); // 色を赤に設定
        ped_line->Draw("same");
        TLine *mip_line = new TLine(result.par[4], 0, result.par[4], h->GetMaximum());
        mip_line->SetLineStyle(2); // 点線に設定
        mip_line->SetLineColor(kRed); // 色を赤に設定
        mip_line->Draw("same");

        c->Update();

        return result;
    }

    // ____________________________________________________________________________________________
    FitResult bht_tot_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        // gPad->SetLogy(1);
        std::vector<Double_t> par, err;
        TString fit_option = h->GetMaximum() > 500.0 ? "0QEMR" : "0QEMRL";

        Double_t peak_pos = h->GetBinCenter(h->GetMaximumBin());
        Double_t stdev    = h->GetStdDev();
        std::pair<Double_t, Double_t> peak_n_sigma(2.0, 2.0);

        // -- first fit -----
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
        f_prefit->SetParameter(1, peak_pos);
        f_prefit->SetParameter(2, stdev);
        h->Fit(f_prefit, fit_option.Data(), "", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));

        // -- second fit -----
        // TF1 *f_fit_g = new TF1( Form("tot_gauss_%s", h->GetName()), "gausn", par[1]-peak_n_sigma.first*par[2], par[1]+peak_n_sigma.second*par[2]);
        TF1 *f_fit_g = new TF1(
            Form("tot_gauss_%s", h->GetName()), 
            [](double *x, double *p) {
                return p[0] * TMath::Gaus(x[0], p[1], p[2], true) + p[3];
            },
            par[1]-peak_n_sigma.first*par[2],
            par[1]+peak_n_sigma.second*par[2],
            4
        );
        f_fit_g->SetParameter(0, par[0]);
        f_fit_g->SetParameter(1, par[1]);
        f_fit_g->SetParameter(2, par[2]*0.9);
        f_fit_g->SetParameter(3, 1.0);
        f_fit_g->SetParLimits(3, 0.0, 100000.0);        
        f_fit_g->SetLineColor(kOrange);
        f_fit_g->SetLineWidth(2.0);
        h->Fit(f_fit_g, fit_option.Data(), "", par[1]-peak_n_sigma.first*par[2], par[1]+peak_n_sigma.second*par[2]);
        Double_t chi_square_g = f_fit_g->GetChisquare();
        Double_t p_value_g = TMath::Prob(chi_square_g, f_fit_g->GetNDF());

        
        // -- first fit -----
        TF1 *f_prefit_l = new TF1("pre_fit_landau", "landaun", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
        f_prefit_l->SetParameter(1, peak_pos);
        f_prefit_l->SetParameter(2, stdev);
        h->Fit(f_prefit_l, fit_option.Data(), "", par[1]-peak_n_sigma.first*par[2], par[1]+peak_n_sigma.second*par[2]);
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit_l->GetParameter(i));
        delete f_prefit_l;

        // -- second fit -----
        // TF1 *f_fit_l = new TF1( Form("tot_landau_%s", h->GetName()), "landaun", par[1]-peak_n_sigma.first*par[2], par[1]+peak_n_sigma.second*par[2]);
        TF1 *f_fit_l = new TF1(
            Form("tot_landau_%s", h->GetName()), 
            [](double *x, double *p) {
                return p[0] * TMath::Landau(x[0], p[1], p[2], true) + p[3];
            },
            par[1]-peak_n_sigma.first*par[2],
            par[1]+peak_n_sigma.second*par[2],
            4
        );
        f_fit_l->SetParameter(0, par[3]);
        f_fit_l->SetParameter(1, par[4]);
        f_fit_l->SetParameter(2, par[5]*0.9);
        f_fit_l->SetParameter(3, 1.0);
        f_fit_l->SetParLimits(3, 0.0, 100000.0);
        f_fit_l->SetLineColor(kOrange);
        f_fit_l->SetLineWidth(2.0);
        h->Fit(f_fit_l, fit_option.Data(), "", par[1]-peak_n_sigma.first*par[2], par[1]+peak_n_sigma.second*par[2]);
        Double_t chi_square_l = f_fit_g->GetChisquare();
        Double_t p_value_l = TMath::Prob(chi_square_l, f_fit_l->GetNDF());

        // -- draw -----
        FitResult result;
        if (p_value_g >= p_value_l) {
        // if (chi_square_g <= chi_square_l) {
            for (Int_t i = 0, n_par = f_fit_g->GetNpar(); i < n_par; i++) {
                result.par.push_back(f_fit_g->GetParameter(i));
                result.err.push_back(f_fit_g->GetParError(i));
            }

            h->GetXaxis()->SetRangeUser(
                result.par[1]- 5.0*result.par[2], 
                result.par[1]+ 5.0*result.par[2]
            );
            h->Draw();
            f_fit_g->SetNpx(1000);
            f_fit_g->Draw("same");
            delete f_fit_l;
        } else {
            for (Int_t i = 0, n_par = f_fit_l->GetNpar(); i < n_par; i++) {
                result.par.push_back(f_fit_l->GetParameter(i));
                result.err.push_back(f_fit_l->GetParError(i));
            }

            h->GetXaxis()->SetRangeUser(
                result.par[1]- 5.0*result.par[2], 
                result.par[1]+ 10.0*result.par[2]
            );
            h->Draw();
            f_fit_l->SetNpx(1000);
            f_fit_l->Draw("same");
            delete f_fit_g;
        }
        TLine *line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        line->SetLineStyle(2); // 点線に設定
        line->SetLineColor(kRed); // 色を赤に設定
        line->Draw("same");

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
    //         h->Fit(f_prefit, "0QEMR", "", par[1]-n_sigma*par[2], par[1]+n_sigma*par[2]);
    //         par.clear();
    //         for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
    //         delete f_prefit;            
    //     }

    //     Double_t fit_range_min = par[1]-n_sigma*par[2];
    //     Double_t fit_range_max = par[1]+n_sigma*par[2];
    //     TF1 *f_fit = new TF1(Form("gauss_%s", h->GetName()), "gausn", fit_range_min, fit_range_max);
    //     f_fit->SetParameters(par[0], par[1], par[2]);
    //     TString fit_option = "0QEMR";
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