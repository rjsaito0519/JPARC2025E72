#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________
    void set_tdc_search_range(TH1D *h) {
        Config& conf = Config::getInstance();

        Double_t peak_pos = h->GetBinCenter(h->GetMaximumBin());
        Double_t stdev    = h->GetStdDev();
        stdev = 10.0*1000.0;
        std::pair<Double_t, Double_t> peak_n_sigma(5.0, 5.0);
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
        f_prefit->SetParameter(1, peak_pos);
        f_prefit->SetParameter(2, stdev*0.5);
        h->Fit(f_prefit, "0QEMR", "", peak_pos-peak_n_sigma.first*stdev, peak_pos+peak_n_sigma.second*stdev);
        Double_t fit_center = f_prefit->GetParameter(1);
        Double_t fit_sigma  = f_prefit->GetParameter(2);
 
        Double_t range_n_sigma = 15.0;
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
        Double_t width    = (conf.tdc_search_range[conf.detector.Data()].second - conf.tdc_search_range[conf.detector.Data()].first) / 3.0;
        std::pair<Double_t, Double_t> peak_n_sigma(2.0, 2.0);

        Double_t evnum_within_range = h->Integral(
            h->FindBin(conf.tdc_search_range[conf.detector.Data()].first),
            h->FindBin(conf.tdc_search_range[conf.detector.Data()].second)
        );
        
        if (evnum_within_range > 300.0) {
            // -- first fit -----
            Double_t lower = std::max(peak_pos-peak_n_sigma.first*width,  h->GetXaxis()->GetXmin());
            Double_t upper = std::min(peak_pos+peak_n_sigma.second*width, h->GetXaxis()->GetXmax());
            TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", lower, upper);
            f_prefit->SetParameter(1, peak_pos);
            f_prefit->SetParameter(2, width*0.9);
            h->Fit(f_prefit, fit_option.Data(), "", lower, upper);
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;

            // -- second fit -----
            lower = std::max(par[1]-peak_n_sigma.first*par[2],  h->GetXaxis()->GetXmin());
            upper = std::min(par[1]+peak_n_sigma.second*par[2], h->GetXaxis()->GetXmax());
            TF1 *f_fit = new TF1( Form("gaus_%s", h->GetName()), "gausn", lower, upper);
            f_fit->SetParameter(0, par[0]);
            f_fit->SetParameter(1, par[1]);
            f_fit->SetParameter(2, par[2]*0.9);
            f_fit->SetLineColor(kOrange);
            f_fit->SetLineWidth(2);
            f_fit->SetNpx(1000);
            h->Fit(f_fit, fit_option.Data(), "", lower, upper);

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
    FitResult adc_fit(TH1D *h, TCanvas *c, Int_t n_c, Int_t n_rebin) {
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
        Double_t mip_range_left = conf.hdprm_mip_range_left > 0.0 ? conf.hdprm_mip_range_left : result.par[1] + conf.adc_ped_remove_nsigma*result.par[2];
        h->GetXaxis()->SetRangeUser(
            mip_range_left, 
            h->GetXaxis()->GetXmax()
        );
        Double_t evnum_within_range = h->Integral(
            h->FindBin(mip_range_left),
            h->FindBin(h->GetXaxis()->GetXmax())
        );
        
        if (evnum_within_range > 300.) {
            Double_t mip_pos          = h->GetBinCenter(h->GetMaximumBin());
            // Double_t mip_half_width   = 20.0;
            Double_t mip_half_width   = h->GetStdDev();
            std::pair<Double_t, Double_t> mip_n_sigma(1.8, 2.2);
            h->GetXaxis()->UnZoom();
            
            // -- first fit -----
            f_prefit->SetRange(mip_pos-mip_half_width, mip_pos+mip_half_width);
            f_prefit->SetParameter(1, mip_pos);
            f_prefit->SetParLimits(1, mip_range_left, h->GetXaxis()->GetXmax());
            f_prefit->SetParameter(2, mip_half_width*0.9);
            h->Fit(f_prefit, "0QEMR", "", mip_range_left, mip_pos+mip_half_width);
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;

            // -- second fit -----
            TF1 *f_fit_mip_g = new TF1( Form("mip_gauss_%s", h->GetName()), "gausn", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
            f_fit_mip_g->SetParameter(1, par[1]);
            f_fit_mip_g->SetParLimits(1, mip_range_left, h->GetXaxis()->GetXmax());
            f_fit_mip_g->SetParameter(2, par[2]*0.9);
            f_fit_mip_g->SetLineColor(kOrange);
            f_fit_mip_g->SetLineWidth(2.0);
            h->Fit(f_fit_mip_g, "0QEMR", "", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
            Double_t chi_square_g = f_fit_mip_g->GetChisquare();
            Double_t p_value_g = TMath::Prob(chi_square_g, f_fit_mip_g->GetNDF());

            TF1 *f_fit_mip_l = new TF1( Form("mip_landau_%s", h->GetName()), "landaun", par[1]-mip_n_sigma.first*par[2], par[1]+mip_n_sigma.second*par[2]);
            f_fit_mip_l->SetParameter(1, par[1]);
            f_fit_mip_l->SetParLimits(1, mip_range_left, h->GetXaxis()->GetXmax());
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
                result.chi_square = chi_square_g;
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
                result.chi_square = chi_square_l;
                delete f_fit_mip_g;
            }
        } else {
            Double_t ped_to_mip = conf.hdprm_typical_value.par[4] - conf.hdprm_typical_value.par[1];
            h->GetXaxis()->SetRangeUser(
                mip_range_left, 
                result.par[1] + ped_to_mip + 5.0*conf.hdprm_typical_value.par[5]
            );

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
                result.par[1] - 10.0*result.par[2], 
                result.par[4] + 10.0*conf.hdprm_typical_value.par[5]
            );
            h->Draw();
            TLatex* commnet = new TLatex();
            commnet->SetNDC();  // NDC座標（0〜1の正規化）を使う
            commnet->SetTextSize(0.04);
            commnet->DrawLatex(0.7, 0.85, "not fitting");
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

    // ____________________________________________________________________________________________
    FitResult t0_offset_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        std::vector<Double_t> par, err;
        Double_t evnum_within_range = h->Integral(
            h->FindBin(h->GetXaxis()->GetXmin()),
            h->FindBin(h->GetXaxis()->GetXmax())
        );
        TString fit_option = evnum_within_range > 200.0 ? "0QEMR" : "0QEMRL";
        FitResult result;

        if (evnum_within_range > 200.0) {

            Double_t peak_pos = h->GetBinCenter(h->GetMaximumBin());
            Double_t width    = h->GetStdDev();
            std::pair<Double_t, Double_t> peak_n_sigma(2.0, 2.0);
            
            // -- first fit -----
            Double_t lower = std::max(peak_pos-peak_n_sigma.first*width,  h->GetXaxis()->GetXmin());
            Double_t upper = std::min(peak_pos+peak_n_sigma.second*width, h->GetXaxis()->GetXmax());
            TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", lower, upper);
            f_prefit->SetParameter(1, peak_pos);
            f_prefit->SetParameter(2, width*0.9);
            h->Fit(f_prefit, fit_option.Data(), "", lower, upper);
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;

            // -- second fit -----
            lower = std::max(par[1]-peak_n_sigma.first*par[2],  h->GetXaxis()->GetXmin());
            upper = std::min(par[1]+peak_n_sigma.second*par[2], h->GetXaxis()->GetXmax());
            TF1 *f_fit = new TF1( Form("gaus_%s", h->GetName()), "gausn", lower, upper);
            f_fit->SetParameter(0, par[0]);
            f_fit->SetParameter(1, par[1]);
            f_fit->SetParameter(2, par[2]*0.9);
            f_fit->SetLineColor(kOrange);
            f_fit->SetLineWidth(2);
            f_fit->SetNpx(1000);
            h->Fit(f_fit, fit_option.Data(), "", lower, upper);

            // -- fill result -----
            for (Int_t i = 0, n_par = f_fit->GetNpar(); i < n_par; i++) {
                result.par.push_back(f_fit->GetParameter(i));
                result.err.push_back(f_fit->GetParError(i));
            }

            // -- draw figure -----
            h->GetXaxis()->SetRangeUser(
                result.par[1] - 5.0*result.par[2], 
                result.par[1] + 5.0*result.par[2]
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
    std::pair<Double_t, Double_t> find_phc_range(TH1D* h, Double_t ratio) {
        Config& conf = Config::getInstance();
            
        // 1. x=1のときのbinとその値
        Int_t ref_bin = h->FindBin(1.0);
        Double_t ref_val = h->GetBinContent( ref_bin );
        Double_t threshold = ref_val * ratio;

        // 2. 右側方向（x増加方向）
        Int_t nbins = h->GetNbinsX();
        Int_t bin_right = -1;
        for (Int_t b = ref_bin + 1; b <= nbins; ++b) {
            if (h->GetBinContent(b) <= threshold) {
                bin_right = b;
                break;
            }
        }
        if (bin_right == -1) {
            bin_right = nbins;
        }

        // 3. 左側方向（x減少方向）
        Int_t bin_left = -1;
        for (Int_t b = ref_bin - 1; b >= 1; --b) {
            if (h->GetBinContent(b) <= threshold) {
                bin_left = b;
                break;
            }
        }
        if (bin_left == -1) {
            bin_left = 1;
        }

        return std::make_pair( 
            std::max(h->GetBinCenter(bin_left),  conf.phc_de_range_min),
            std::min(h->GetBinCenter(bin_right), conf.phc_de_range_max)
        );
    }
    
    // ____________________________________________________________________________________________
    FitResult phc_fit(TH2D *h, TCanvas *c, Int_t n_c) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        std::vector<Double_t> par, err;
        TString fit_option = "0QEMR";
        FitResult result;

        // -- make TProfile -----
        Int_t time_min = h->GetYaxis()->FindBin(conf.phc_time_window[conf.detector.Data()].first);
        Int_t time_max = h->GetYaxis()->FindBin(conf.phc_time_window[conf.detector.Data()].second);
        TProfile *pf   = h->ProfileX(Form("profile_%s", h->GetName()), time_min, time_max);
        pf->SetLineColor(kRed);

        // -- prepare parameter -----
        TH1D* de_proj  = h->ProjectionX(Form("projection_%s", h->GetName()), time_min, time_max);
        std::pair<Double_t, Double_t> fit_range = find_phc_range(de_proj, conf.phc_de_range_ratio[conf.detector.Data()]);        
        std::vector<std::vector<Double_t>> par_limits = {
            {0.001, 15.0},
            {-5.0, fit_range.first},
            {-1.0, 10.0}
        };

        // -- fit -----
        TF1 *f_fit = new TF1( Form("phc_%s", h->GetName()), "-[0]/TMath::Sqrt(TMath::Abs(x-[1]))+[2]", fit_range.first, fit_range.second);
        for (Int_t i = 0; i < 3; i++) {
            f_fit->SetParameter(i, (3.0*par_limits[i][0]+par_limits[i][1])/4.0);
            f_fit->SetParLimits(i, par_limits[i][0], par_limits[i][1]);
        }
        f_fit->SetParameter(1, fit_range.first - (fit_range.second - fit_range.first)*0.005);
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(1);
        pf->Fit(f_fit, fit_option.Data(), "", fit_range.first, fit_range.second);

        // -- fill result -----
        for (Int_t i = 0, n_par = f_fit->GetNpar(); i < n_par; i++) {
            result.par.push_back(f_fit->GetParameter(i));
            result.err.push_back(f_fit->GetParError(i));
        }

        // draw
        h->Draw("colz");
        pf->Draw("same");
        f_fit->Draw("same");

        Double_t chi_square = f_fit->GetChisquare();
        Double_t p_value = TMath::Prob(chi_square, f_fit->GetNDF());
        TLatex* commnet = new TLatex();
        commnet->SetNDC();  // NDC座標（0〜1の正規化）を使う
        commnet->SetTextSize(0.04);
        commnet->DrawLatex(0.7, 0.85, Form("%2f", chi_square/f_fit->GetNDF()));

        return result;
    }


}