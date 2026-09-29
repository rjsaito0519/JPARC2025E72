#include "ana_helper.h"
#include <TGraphErrors.h>

namespace ana_helper {

    // ____________________________________________________________________________________________
    void set_tdc_search_range(TH1D *h) {
        Config& conf = Config::getInstance();
        TString det_upper = conf.detector;
        det_upper.ToUpper();
        TString key = det_upper + "_TDC";

        if (conf.run_num > 0) {
            std::pair<Double_t, Double_t> user_range = get_user_param(key, conf.run_num);
            if (user_range.first != 0.0 || user_range.second != 0.0) {
                 if (user_range.first < user_range.second) {
                     conf.tdc_search_range[conf.detector.Data()].first  = user_range.first;
                     conf.tdc_search_range[conf.detector.Data()].second = user_range.second;
                     return;
                 }
            }
        }

        Double_t peak_pos = h->GetBinCenter(h->GetMaximumBin());
        Double_t stdev    = h->GetStdDev();
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
            std::max(conf.tdc_search_range[conf.detector.Data()].first,  h->GetXaxis()->GetXmin()), 
            std::min(conf.tdc_search_range[conf.detector.Data()].second, h->GetXaxis()->GetXmax())
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
            peak_n_sigma = std::make_pair(2.0, 2.0);
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
        h = (TH1D*)h->Rebin(n_rebin, h->GetName());

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
        
        if (evnum_within_range > 500. || conf.hdprm_typical_value.par.empty()) {
            Int_t mip_peak_bin     = h->GetMaximumBin();
            Double_t mip_pos      = h->GetBinCenter(mip_peak_bin);
            Double_t mip_peak_h   = h->GetBinContent(mip_peak_bin);
            std::pair<Double_t, Double_t> mip_n_sigma(1.8, 2.2);
            h->GetXaxis()->UnZoom();

            Int_t bin_left = mip_peak_bin;
            Int_t bin_left_lim = h->FindBin(mip_range_left);
            while (bin_left > bin_left_lim && h->GetBinContent(bin_left) > 0.5*mip_peak_h) --bin_left;
            Int_t bin_right = mip_peak_bin;
            while (bin_right < h->GetNbinsX() && h->GetBinContent(bin_right) > 0.3*mip_peak_h) ++bin_right;

            Double_t bin_w = h->GetBinWidth(mip_peak_bin);
            Double_t local_left  = mip_pos - h->GetBinCenter(bin_left);
            Double_t local_right = h->GetBinCenter(bin_right) - mip_pos;
            if (local_left  < bin_w) local_left  = bin_w;
            if (local_right < bin_w) local_right = bin_w;
            Double_t mip_half_width = 0.5*(local_left + local_right);

            Double_t prefit_min = mip_pos - local_left;
            if (prefit_min < mip_range_left) prefit_min = mip_range_left;
            Double_t prefit_max = mip_pos + local_right;

            // -- first fit -----
            f_prefit->SetRange(prefit_min, prefit_max);
            f_prefit->SetParameter(1, mip_pos);
            f_prefit->SetParLimits(1, mip_range_left, h->GetXaxis()->GetXmax());
            f_prefit->SetParameter(2, mip_half_width*conf.hdprm_mip_half_width_ratio);
            h->Fit(f_prefit, "0QEMR", "", prefit_min, prefit_max);
            for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));
            delete f_prefit;

            // -- second fit -----
            Double_t fit_mean  = par[1];
            Double_t fit_sigma = par[2];
            if (fit_sigma > 2.0*mip_half_width) fit_sigma = mip_half_width;
            if (fit_sigma < 0.3*mip_half_width) fit_sigma = mip_half_width;

            TF1 *f_fit_mip_g = new TF1( Form("mip_gauss_%s", h->GetName()), "gausn(0)+pol0(3)", prefit_min, prefit_max);
            f_fit_mip_g->SetParameter(0, par[0]);
            f_fit_mip_g->SetParameter(1, fit_mean);
            f_fit_mip_g->SetParameter(2, fit_sigma*0.9);
            f_fit_mip_g->SetParameter(3, 0.0);
            f_fit_mip_g->SetLineColor(kOrange);
            f_fit_mip_g->SetLineWidth(2.0);

            TF1 *f_fit_mip_l = new TF1( Form("mip_landau_%s", h->GetName()), "landaun(0)+pol0(3)", prefit_min, prefit_max);
            f_fit_mip_l->SetParameter(0, par[0]);
            f_fit_mip_l->SetParameter(1, fit_mean);
            f_fit_mip_l->SetParameter(2, fit_sigma*0.9);
            f_fit_mip_l->SetParameter(3, 0.0);
            f_fit_mip_l->SetLineColor(kOrange);
            f_fit_mip_l->SetLineWidth(2.0);

            Double_t fit_range_min = 0.0;
            Double_t fit_range_max = 0.0;
            Double_t chi_square_g = 0.0;
            Double_t chi_square_l = 0.0;
            Double_t p_value_g = 0.0;
            Double_t p_value_l = 0.0;

            for (Int_t iround = 0; iround < 2; ++iround) {
                fit_range_min = fit_mean - mip_n_sigma.first*fit_sigma;
                if (fit_range_min < mip_range_left) fit_range_min = mip_range_left;
                fit_range_max = fit_mean + mip_n_sigma.second*fit_sigma;

                Double_t mip_bg = h->GetBinContent(h->FindBin(fit_range_min));
                if (mip_bg < 0.0) mip_bg = 0.0;

                f_fit_mip_g->SetRange(fit_range_min, fit_range_max);
                f_fit_mip_g->SetParLimits(1, mip_range_left, h->GetXaxis()->GetXmax());
                f_fit_mip_g->SetParameter(3, mip_bg);
                f_fit_mip_g->SetParLimits(3, 0.0, h->GetMaximum());
                h->Fit(f_fit_mip_g, "0QEMR", "", fit_range_min, fit_range_max);
                chi_square_g = f_fit_mip_g->GetChisquare();
                p_value_g = TMath::Prob(chi_square_g, f_fit_mip_g->GetNDF());

                f_fit_mip_l->SetRange(fit_range_min, fit_range_max);
                f_fit_mip_l->SetParLimits(1, mip_range_left, h->GetXaxis()->GetXmax());
                f_fit_mip_l->SetParameter(3, mip_bg);
                f_fit_mip_l->SetParLimits(3, 0.0, h->GetMaximum());
                h->Fit(f_fit_mip_l, "0QEMR", "", fit_range_min, fit_range_max);
                chi_square_l = f_fit_mip_l->GetChisquare();
                p_value_l = TMath::Prob(chi_square_l, f_fit_mip_l->GetNDF());

                if (chi_square_g <= chi_square_l) {
                    fit_mean  = f_fit_mip_g->GetParameter(1);
                    fit_sigma = f_fit_mip_g->GetParameter(2);
                } else {
                    fit_mean  = f_fit_mip_l->GetParameter(1);
                    fit_sigma = f_fit_mip_l->GetParameter(2);
                }
            }

            // // -- debug ------
            // std::cout << "gauss:  " << p_value_g << ", " << f_fit_mip_g->GetChisquare() << std::endl;
            // std::cout << "landau: " << p_value_l << ", " << f_fit_mip_l->GetChisquare() << std::endl;
            // Bool_t flag = p_value_g >= p_value_l;
            // std::cout << flag << std::endl;

            // if (p_value_g >= p_value_l) {
            if (chi_square_g <= chi_square_l) {
                for (Int_t i = 0; i < 3; i++) {
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
                for (Int_t i = 0; i < 3; i++) {
                    result.par.push_back(f_fit_mip_l->GetParameter(i));
                    result.err.push_back(f_fit_mip_l->GetParError(i));
                }
                // Use actual maximum of the function as the representative peak (user suggestion)
                result.par[4] = f_fit_mip_l->GetMaximumX(fit_range_min, fit_range_max);

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
                result.par[1] + ped_to_mip + 10.0*conf.hdprm_typical_value.par[5]
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
    // HTOF ADC fit helpers (used only by htof_adc_fit below).

    // 3-bin running-mean smoothed copy of h rebinned by n (caller owns the result).
    static TH1D* htof_adc_coarse(TH1D *h, Int_t n, const char *name) {
        TH1D *hc  = (TH1D*)h->Rebin(n, name);
        TH1D *tmp = (TH1D*)hc->Clone(Form("%s_tmp", name));
        for (Int_t b = 2; b < hc->GetNbinsX(); ++b)
            hc->SetBinContent(b, (tmp->GetBinContent(b-1) + tmp->GetBinContent(b) + tmp->GetBinContent(b+1)) / 3.0);
        delete tmp;
        return hc;
    }

    // Bump search in [xmin, xmax] on the coarse histogram hc (optionally skipping a falling
    // edge at xmin up to its first local minimum), walk to 50% (left) / 30% (right) of the
    // maximum, then a bounded gaus fit on h and one refit within +-1 sigma.
    // Returns {mean, sigma}; ok=false if no usable peak.
    static Bool_t htof_adc_find_peak(TH1D *h, TH1D *hc, Double_t xmin, Double_t xmax, Bool_t skip_edge,
                                     Double_t &mean, Double_t &sigma) {
        mean = -1.0; sigma = -1.0;
        Int_t b0 = hc->FindBin(xmin) + 1, b1 = hc->FindBin(xmax);
        if (b1 - b0 < 3) return false;
        Int_t start = b0;
        if (skip_edge) {
            while (start < b1 && hc->GetBinContent(start+1) <= hc->GetBinContent(start)) ++start;
            if (start >= b1 - 1) start = b0; // monotonically falling: no bump found, fall back
        }
        Int_t pb = start;
        for (Int_t b = start; b <= b1; ++b) if (hc->GetBinContent(b) > hc->GetBinContent(pb)) pb = b;
        Double_t ph = hc->GetBinContent(pb);
        if (ph <= 0.0) return false;
        Int_t bl = pb; while (bl > start && hc->GetBinContent(bl-1) > 0.5*ph) --bl;
        Int_t br = pb; while (br < b1    && hc->GetBinContent(br+1) > 0.3*ph) ++br;
        Double_t lo = std::max(hc->GetBinLowEdge(bl), xmin);
        Double_t hi = hc->GetBinLowEdge(br+1);
        Double_t bw = hc->GetBinWidth(pb);
        if (hi - lo < 2.0*bw) { lo = std::max(hc->GetBinCenter(pb) - bw, xmin); hi = hc->GetBinCenter(pb) + bw; }

        TF1 g("htof_adc_seed_gaus", "gaus", lo, hi);
        g.SetParameters(ph * h->GetBinWidth(1) / bw, hc->GetBinCenter(pb), 0.5*(hi - lo));
        g.SetParLimits(1, lo, hi);
        g.SetParLimits(2, 0.5*bw, 2.0*(hi - lo));
        h->Fit(&g, "0QRB", "", lo, hi);
        Double_t m = g.GetParameter(1), s = g.GetParameter(2);
        Double_t lo2 = std::max(m - s, xmin), hi2 = m + s;
        if (hi2 - lo2 > 2.0*bw) {
            g.SetRange(lo2, hi2);
            g.SetParLimits(1, lo2, hi2);
            h->Fit(&g, "0QRB", "", lo2, hi2);
            m = g.GetParameter(1); s = g.GetParameter(2);
        }
        mean = m; sigma = s;
        return (s > 0.0 && m > xmin && m < xmax);
    }

    // ____________________________________________________________________________________________
    // HTOF ADC fit: pedestal from a species-unselected ("raw") per-segment/side histogram,
    // MIP peak from a species-tagged ("selected", e.g. pion) histogram of the same
    // segment/side.
    //
    // Rationale: HTOF's per-segment cluster match does not guarantee every individual
    // channel (U/D/S) saw real light, even for a correctly PID-tagged track - so a
    // residual pedestal-height bump can leak into the species-selected histogram, and
    // fitting the pedestal self-referentially from that same (selected) histogram is
    // unreliable. The pedestal must instead be read from the raw, unselected distribution.
    //
    // TENTATIVE MIP procedure:
    //  1) cut: first point right of the pedestal where raw (8-ch bins) drops below 1% of the
    //     pedestal peak, or its first local minimum (high-rate channels). A cut at
    //     ped + n*sigma would overshoot the MIP for broad (e.g. S) pedestals.
    //  2) raw MIP seed: bump search on coarse, smoothed raw above the cut. Most HTOF hits are
    //     non-proton, so the raw bump is a good seed (not used as the MIP value itself: it is
    //     systematically 10-30% lower than the pion-tagged peak on the edge segments).
    //  3) selected: rebin chosen from counts around the seed (n_rebin <= 0; ~100 counts/bin,
    //     2..32), peak search restricted to [seed-2s, seed+3s], then the final model
    //     (gausn+pol0 vs landaun+pol0, smaller chi2; MIP = gauss mean / landau maximum)
    //     in a single, fixed range (no iterative range update).
    //
    // Uses 3 canvas pads (n_c: raw/pedestal, n_c+1: selected/MIP, n_c+2: selected close-up
    // around the fit range, linear scale). par/err layout matches
    // adc_fit(): [0..2]=pedestal (amp,pos,sigma) from h_raw, [3..5]=MIP (amp,pos,sigma)
    // from h_selected; par[4] is the representative MIP position.
    // additional = {flag, rebin, cut, raw_seed_mean, raw_seed_sigma}; flag bits:
    //   1 = raw seed unusable, 2 = low statistics (<500) around the seed, 4 = ndf < 10,
    //   8 = MIP position at the fit-range edge / no fit, 16 = candidate two-peak structure,
    //   128 = MIP range / model taken from param::htof_adc_fit_hint[hint_key].
    //   (32, 64: beam-window weak-side refit, see htof_adc_fit_weak.)
    FitResult htof_adc_fit(TH1D *h_raw, TH1D *h_selected, TCanvas *c, Int_t n_c, Int_t n_rebin,
                           const std::string& hint_key) {
        const Double_t display_max = 2048.0; // HTOF ADC peaks never appear above ~2048 counts
        Config& conf = Config::getInstance();

        FitResult result;
        result.chi_square = 0.0;
        result.ndf = 0;
        result.migrad_stats = 0;
        Int_t flag = 0;

        // -- pedestal: fit from the raw (species-unselected) histogram -----
        Double_t ped_pos        = h_raw->GetBinCenter(h_raw->GetMaximumBin());
        Double_t ped_half_width = 5.0;
        TF1 f_prefit("pre_fit_gauss", "gausn", ped_pos-ped_half_width, ped_pos+ped_half_width);
        f_prefit.SetParameter(1, ped_pos);
        f_prefit.SetParameter(2, ped_half_width);
        h_raw->Fit(&f_prefit, "0QEMR", "", ped_pos-ped_half_width, ped_pos+ped_half_width);
        Double_t p1 = f_prefit.GetParameter(1), p2 = f_prefit.GetParameter(2);

        TF1 *f_fit_ped = new TF1(Form("ped_%s", h_raw->GetName()), "gausn", p1-2.0*p2, p1+2.0*p2);
        f_fit_ped->SetParameters(f_prefit.GetParameter(0), p1, p2*0.9);
        f_fit_ped->SetLineColor(kOrange);
        f_fit_ped->SetLineWidth(2);
        f_fit_ped->SetNpx(1000);
        h_raw->Fit(f_fit_ped, "0QEMR", "", p1-2.0*p2, p1+2.0*p2);
        for (Int_t i = 0; i < 3; i++) {
            result.par.push_back(f_fit_ped->GetParameter(i));
            result.err.push_back(f_fit_ped->GetParError(i));
        }
        Double_t ped = result.par[1];

        // -- cut: raw drops below 1% of the pedestal peak (or its first local minimum) -----
        Double_t cut;
        {
            TH1D *h8 = (TH1D*)h_raw->Rebin(8, Form("%s_cut8", h_raw->GetName()));
            Int_t b = h8->FindBin(ped);
            Double_t pmax = h8->GetBinContent(b);
            while (b < h8->GetNbinsX() && h8->GetBinContent(b) > 0.01*pmax
                   && !(h8->GetBinContent(b+1) > h8->GetBinContent(b) && h8->GetBinContent(b) < 0.5*pmax)) ++b;
            cut = h8->GetBinLowEdge(b+1);
            delete h8;
        }

        // -- raw MIP seed -----
        Double_t seed_m, seed_s;
        TH1D *h_raw_c = htof_adc_coarse(h_raw, 16, Form("%s_c16", h_raw->GetName()));
        Bool_t seed_ok = htof_adc_find_peak(h_raw, h_raw_c, cut, display_max, true, seed_m, seed_s);
        delete h_raw_c;
        // seed not usable if its width is comparable to its distance from the pedestal
        if (seed_ok && seed_s > 0.8*(seed_m - ped)) seed_ok = false;
        if (!seed_ok) flag |= 1;

        // -- selected: rebin and search window from the raw seed -----
        Double_t win_lo = cut, win_hi = display_max;
        if (seed_ok) {
            win_lo = std::max(cut, seed_m - 2.0*seed_s);
            win_hi = std::min(display_max, seed_m + 3.0*seed_s);
        }
        Double_t s_ref = seed_ok ? seed_s : 100.0;
        Double_t n_win = h_selected->Integral(h_selected->FindBin(std::max(cut, seed_m - s_ref)),
                                              h_selected->FindBin(seed_m + s_ref));
        if (n_win < 500.0) flag |= 2;
        if (n_rebin <= 0) {
            Double_t w_need = (n_win > 0.0) ? 2.0*s_ref*100.0/n_win : 32.0;
            n_rebin = 2;
            while (n_rebin < w_need && n_rebin < 32) n_rebin *= 2;
        }
        TH1D *h_sel = (TH1D*)h_selected->Rebin(n_rebin, Form("%s_rb", h_selected->GetName()));
        TH1D *h_sel_c = htof_adc_coarse(h_selected, std::max(16, n_rebin), Form("%s_c", h_selected->GetName()));

        Double_t sel_m, sel_s;
        Bool_t sel_ok = htof_adc_find_peak(h_sel, h_sel_c, win_lo, win_hi, false, sel_m, sel_s);

        // candidate two-peak structure: another local maximum >= 50% of the highest one in
        // [cut, display_max], separated by a dip below 70% of the lower of the two (not resolved here)
        {
            Int_t b0 = h_sel_c->FindBin(cut) + 1, b1 = h_sel_c->FindBin(display_max);
            std::vector<Int_t> mx;
            for (Int_t b = b0+1; b < b1; ++b)
                if (h_sel_c->GetBinContent(b) >= h_sel_c->GetBinContent(b-1) && h_sel_c->GetBinContent(b) > h_sel_c->GetBinContent(b+1)) mx.push_back(b);
            Double_t top = 0.0;
            for (Int_t b : mx) top = std::max(top, h_sel_c->GetBinContent(b));
            for (size_t i = 0; i < mx.size(); ++i) for (size_t j = i+1; j < mx.size(); ++j) {
                Double_t a = h_sel_c->GetBinContent(mx[i]), cc = h_sel_c->GetBinContent(mx[j]);
                if (a < 0.5*top || cc < 0.5*top) continue;
                Double_t dip = a;
                for (Int_t b = mx[i]; b <= mx[j]; ++b) dip = std::min(dip, h_sel_c->GetBinContent(b));
                if (dip < 0.7*std::min(a, cc)) flag |= 16;
            }
        }
        delete h_sel_c;

        // -- final MIP fit: fixed range, gauss vs landau (smaller chi2), or the range / model of
        //    param::htof_adc_fit_hint for this channel -----
        Int_t model = 0; // 0 both, 1 gauss, 2 landau
        const auto hint = param::htof_adc_fit_hint.find(hint_key);
        const Bool_t use_hint = (!hint_key.empty() && hint != param::htof_adc_fit_hint.end()
                                 && hint->second.size() >= 3);
        if (use_hint) {
            flag |= 128;
            model = static_cast<Int_t>(hint->second[2]);
            h_sel->GetXaxis()->SetRangeUser(hint->second[0], hint->second[1]);
            sel_m = h_sel->GetBinCenter(h_sel->GetMaximumBin());
            h_sel->GetXaxis()->SetRange(0, 0);
            sel_s = 0.25*(hint->second[1] - hint->second[0]);
            sel_ok = kTRUE;
        }
        TF1 *f_best = nullptr;
        Double_t fit_lo = 0.0, fit_hi = 0.0;
        if (sel_ok) {
            fit_lo = use_hint ? hint->second[0] : std::max(sel_m - 1.8*sel_s, cut);
            fit_hi = use_hint ? hint->second[1] : sel_m + 2.2*sel_s;
            Double_t bg  = std::max(0.0, h_sel->GetBinContent(h_sel->FindBin(fit_lo)));
            Double_t amp = h_sel->Integral(h_sel->FindBin(sel_m - sel_s), h_sel->FindBin(sel_m + sel_s), "width");

            TF1 *f_g = new TF1(Form("mip_gauss_%s", h_selected->GetName()), "gausn(0)+pol0(3)", fit_lo, fit_hi);
            f_g->SetParameters(amp, sel_m, sel_s*0.9, bg);
            f_g->SetParLimits(1, fit_lo, fit_hi);
            f_g->SetParLimits(3, 0.0, h_sel->GetMaximum());
            h_sel->Fit(f_g, "0QEMR", "", fit_lo, fit_hi);

            TF1 *f_l = new TF1(Form("mip_landau_%s", h_selected->GetName()), "landaun(0)+pol0(3)", fit_lo, fit_hi);
            f_l->SetParameters(amp, sel_m, sel_s*0.9, bg);
            f_l->SetParLimits(1, fit_lo, fit_hi);
            f_l->SetParLimits(3, 0.0, h_sel->GetMaximum());
            h_sel->Fit(f_l, "0QEMR", "", fit_lo, fit_hi);

            // Hint channels: the fit can stop in a local minimum (landau + const trading width
            // for background) depending on the initial width; start from several widths and keep
            // the smallest chi2.
            if (use_hint) {
                for (TF1 *f : {f_g, f_l}) {
                    std::vector<Double_t> best(f->GetNpar());
                    for (Int_t k = 0; k < f->GetNpar(); ++k) best[k] = f->GetParameter(k);
                    Double_t best_chi2 = f->GetChisquare();
                    for (Double_t frac : {0.1, 0.2, 0.3}) {
                        f->SetParameters(amp, sel_m, frac*(fit_hi - fit_lo), bg);
                        h_sel->Fit(f, "0QEMR", "", fit_lo, fit_hi);
                        if (f->GetChisquare() < best_chi2) {
                            best_chi2 = f->GetChisquare();
                            for (Int_t k = 0; k < f->GetNpar(); ++k) best[k] = f->GetParameter(k);
                        }
                    }
                    f->SetParameters(best.data());
                    h_sel->Fit(f, "0QEMR", "", fit_lo, fit_hi);
                }
            }

            const Bool_t take_gauss = (model == 1) || (model == 0 && f_g->GetChisquare() <= f_l->GetChisquare());
            if (take_gauss) { f_best = f_g; delete f_l; }
            else            { f_best = f_l; delete f_g; }
            for (Int_t i = 0; i < 3; i++) {
                result.par.push_back(f_best->GetParameter(i));
                result.err.push_back(f_best->GetParError(i));
            }
            // Landau: use actual maximum of the function as the representative peak
            if (f_best == f_l) result.par[4] = f_best->GetMaximumX(fit_lo, fit_hi);
            result.chi_square = f_best->GetChisquare();
            result.ndf        = f_best->GetNDF();
            if (result.ndf < 10) flag |= 4;
            Double_t pm = f_best->GetParameter(1);
            if (std::fabs(pm - fit_lo) < 1e-3*fit_hi || std::fabs(pm - fit_hi) < 1e-3*fit_hi) flag |= 8;
        } else {
            flag |= 8;
            Double_t err_value = 9999.0;
            result.par.push_back(n_win);             result.err.push_back(err_value);
            result.par.push_back(h_sel->GetMean());  result.err.push_back(h_sel->GetMeanError());
            result.par.push_back(h_sel->GetStdDev()); result.err.push_back(h_sel->GetStdDevError());
        }
        result.additional = { static_cast<Double_t>(flag), static_cast<Double_t>(n_rebin), cut, seed_m, seed_s };

        // -- draw raw/pedestal pad -----
        c->cd(n_c);
        gPad->SetLogy(1);
        h_raw->GetXaxis()->SetRangeUser(conf.adc_min, display_max);
        h_raw->Draw("HIST");
        f_fit_ped->Draw("same");
        Double_t y_raw = h_raw->GetMaximum();
        TLine *ped_line_raw = new TLine(ped, 0, ped, y_raw);
        ped_line_raw->SetLineStyle(2);
        ped_line_raw->SetLineColor(kRed);
        ped_line_raw->Draw("same");
        TLine *cut_line_raw = new TLine(cut, 0, cut, y_raw);
        cut_line_raw->SetLineStyle(3);
        cut_line_raw->SetLineColor(kBlue);
        cut_line_raw->Draw("same");
        if (seed_ok) {
            TLine *seed_line_raw = new TLine(seed_m, 0, seed_m, y_raw);
            seed_line_raw->SetLineColor(kGreen+2);
            seed_line_raw->SetLineWidth(2);
            seed_line_raw->Draw("same");
        }
        c->Update();

        // -- draw selected/MIP pad -----
        c->cd(n_c+1);
        gPad->SetLogy(1);
        h_sel->GetXaxis()->SetRangeUser(conf.adc_min, display_max);
        h_sel->Draw("HIST");
        if (f_best) {
            f_best->SetLineColor(kOrange);
            f_best->SetLineWidth(2);
            f_best->SetNpx(1000);
            f_best->Draw("same");
        } else {
            TLatex *comment = new TLatex();
            comment->SetNDC();
            comment->SetTextSize(0.04);
            comment->DrawLatex(0.6, 0.85, "not fitting");
        }
        Double_t y_sel = h_sel->GetMaximum();
        // Raw pedestal position (magenta) - lets you check by eye whether a peak visible
        // in the selected histogram near this line is leftover pedestal, not real MIP.
        TLine *ped_ref_line = new TLine(ped, 0, ped, y_sel);
        ped_ref_line->SetLineStyle(2);
        ped_ref_line->SetLineColor(kMagenta);
        ped_ref_line->Draw("same");
        TLine *cut_line_sel = new TLine(cut, 0, cut, y_sel);
        cut_line_sel->SetLineStyle(3);
        cut_line_sel->SetLineColor(kBlue);
        cut_line_sel->Draw("same");
        if (seed_ok) {
            TLine *seed_line_sel = new TLine(seed_m, 0, seed_m, y_sel);
            seed_line_sel->SetLineStyle(2);
            seed_line_sel->SetLineColor(kGreen+2);
            seed_line_sel->Draw("same");
        }
        TLine *mip_line_sel = new TLine(result.par[4], 0, result.par[4], y_sel);
        mip_line_sel->SetLineStyle(2);
        mip_line_sel->SetLineColor(kRed);
        mip_line_sel->Draw("same");
        TLatex *info = new TLatex();
        info->SetNDC();
        info->SetTextSize(0.04);
        info->DrawLatex(0.45, 0.80, Form("rebin %d, flag %d", n_rebin, flag));

        // -- draw selected/MIP close-up pad (linear, around the fit range) -----
        c->cd(n_c+2);
        gPad->SetLogy(0);
        Double_t z_lo, z_hi;
        if (f_best) {
            const Double_t w = fit_hi - fit_lo;
            z_lo = std::max(cut, fit_lo - 0.6*w);
            z_hi = std::min(display_max, fit_hi + 0.6*w);
        } else {
            z_lo = win_lo;
            z_hi = win_hi;
        }
        TH1D *h_zoom = (TH1D*)h_sel->Clone(Form("%s_zoom", h_sel->GetName()));
        h_zoom->SetTitle(Form("%s (close-up)", h_sel->GetTitle()));
        h_zoom->GetXaxis()->SetRangeUser(z_lo, z_hi);
        h_zoom->SetMinimum(0.0);
        h_zoom->Draw("HIST");
        if (f_best) f_best->Draw("same");
        Double_t y_zoom = h_zoom->GetMaximum();
        for (Double_t x : {fit_lo, fit_hi}) {
            if (!f_best) break;
            TLine *l = new TLine(x, 0, x, y_zoom);
            l->SetLineStyle(3);
            l->SetLineColor(kGray+2);
            l->Draw("same");
        }
        TLine *mip_line_zoom = new TLine(result.par[4], 0, result.par[4], y_zoom);
        mip_line_zoom->SetLineStyle(2);
        mip_line_zoom->SetLineColor(kRed);
        mip_line_zoom->Draw("same");
        TLatex *info_zoom = new TLatex();
        info_zoom->SetNDC();
        info_zoom->SetTextSize(0.05);
        info_zoom->DrawLatex(0.50, 0.82, Form("MIP %.1f", result.par[4]));
        if (result.ndf > 0) info_zoom->DrawLatex(0.50, 0.75, Form("#chi^{2}/ndf %.1f/%d", result.chi_square, result.ndf));

        c->Update();

        return result;
    }

    // ____________________________________________________________________________________________
    // HTOF ADC MIP fit for the weak side of the beam-window short segments (seg1 D, seg2 U, seg3 D,
    // seg4 U), TENTATIVE. Hits there concentrate next to the beam window and the far
    // readout end sees little light: for tracks whose strong side is MIP-like, the weak-side ADC is a
    // sharp spike just above the pedestal plus a broad MIP hump, and htof_adc_fit() picks the spike.
    // h_cond: weak-side ADC of the ADC sample restricted to tracks whose strong-side ADC is near its
    // MIP (selection done by the caller). Model: gaus (spike) + landau (hump) + const, likelihood fit
    // in [ped + 5, ped + 3 u], u = strong-side (MIP - ped) (initial values / limits scale with u).
    // base: htof_adc_fit() result of the same channel; its pedestal (par[0..2]) and cut are kept.
    // Returns par[3..5] = landau (amp, MPV, width) with par[4] replaced by the landau maximum (MIP).
    // additional = {flag, rebin, cut, spike mean, spike sigma}; flag bits:
    //   2 = low statistics (< 500 entries), 8 = MPV at its limit / no fit,
    //   32 = weak-side two-component fit used (always set), 64 = hump not separated from the spike
    //   (hump fraction < min_hump_frac: MIP uncertain), 128 = hump MPV range taken from
    //   param::htof_adc_fit_weak_hint[hint_key].
    // Redraws pads n_c (h_cond, log) and n_c+1 (close-up, linear) with both components.
    FitResult htof_adc_fit_weak(TH1D *h_cond, const FitResult& base, Double_t strong_unit,
                                TCanvas *c, Int_t n_c, Double_t min_hump_frac,
                                const std::string& hint_key) {
        FitResult result;
        result.par.assign(base.par.begin(), base.par.begin() + 3);
        result.err.assign(base.err.begin(), base.err.begin() + 3);
        result.chi_square = 0.0;
        result.ndf = 0;
        result.migrad_stats = 0;
        Int_t flag = 32;
        const Double_t ped = base.par[1];
        const Double_t cut = base.additional.size() > 2 ? base.additional[2] : ped;
        const Double_t u   = (strong_unit > 10.0) ? strong_unit : 100.0;

        const Double_t bw = h_cond->GetXaxis()->GetBinWidth(1);
        const Int_t n_rebin = std::max(1, static_cast<Int_t>(std::lround(4.0 / bw))); // ~4 ADC ch per bin
        TH1D *h = (TH1D*)h_cond->Rebin(n_rebin, Form("%s_rb", h_cond->GetName()));
        if (h_cond->GetEntries() < 500) flag |= 2;

        // spike seed: maximum in [ped+5, ped+80]
        h->GetXaxis()->SetRangeUser(ped + 5.0, ped + 80.0);
        const Double_t sp = h->GetBinCenter(h->GetMaximumBin());
        const Double_t sp_amp = std::max(1.0, h->GetMaximum());
        h->GetXaxis()->SetRange(0, 0);

        const Double_t lo = ped + 5.0, hi = ped + 3.0*u;
        // hump MPV range: automatic, or param::htof_adc_fit_weak_hint for this channel (flag 128)
        const auto hint = param::htof_adc_fit_weak_hint.find(hint_key);
        const Bool_t use_hint = (!hint_key.empty() && hint != param::htof_adc_fit_weak_hint.end()
                                 && hint->second.size() >= 2);
        const Double_t mpv_lo = use_hint ? hint->second[0] : sp + 30.0;
        const Double_t mpv_hi = use_hint ? hint->second[1] : ped + 2.0*u;
        if (use_hint) flag |= 128;
        const Double_t mpv0 = use_hint ? 0.5*(mpv_lo + mpv_hi) : ped + 0.6*u;
        TF1 *f = new TF1(Form("mip_weak_%s", h_cond->GetName()), "gaus(0)+landau(3)+[6]", lo, hi);
        f->SetParameters(sp_amp, sp, 12.0, 0.3*sp_amp, mpv0, 0.25*(mpv0 - ped), 0.0);
        f->SetParLimits(1, ped, ped + 80.0);
        f->SetParLimits(2, 2.0, 40.0);
        f->SetParLimits(3, 0.0, 10.0*sp_amp);
        f->SetParLimits(4, mpv_lo, mpv_hi);
        f->SetParLimits(5, 5.0, u);
        f->SetParLimits(6, 0.0, 0.05*sp_amp);
        f->SetNpx(1000);
        h->Fit(f, "0QRL");
        if (use_hint) {
            // start from several hump positions / widths and keep the best likelihood
            std::vector<Double_t> best(f->GetNpar());
            for (Int_t k = 0; k < f->GetNpar(); ++k) best[k] = f->GetParameter(k);
            Double_t best_chi2 = f->GetChisquare();
            for (Double_t mfrac : {0.3, 0.7}) for (Double_t wfrac : {0.1, 0.2, 0.35}) {
                f->SetParameters(sp_amp, sp, 12.0, 0.3*sp_amp, mpv_lo + mfrac*(mpv_hi - mpv_lo), wfrac*u, 0.0);
                h->Fit(f, "0QRL");
                if (f->GetChisquare() < best_chi2) {
                    best_chi2 = f->GetChisquare();
                    for (Int_t k = 0; k < f->GetNpar(); ++k) best[k] = f->GetParameter(k);
                }
            }
            f->SetParameters(best.data());
            h->Fit(f, "0QRL");
        }

        TF1 *f_s = new TF1(Form("%s_spike", f->GetName()), "gaus", lo, hi);
        f_s->SetParameters(f->GetParameter(0), f->GetParameter(1), f->GetParameter(2));
        TF1 *f_h = new TF1(Form("%s_hump", f->GetName()), "landau", lo, hi);
        f_h->SetParameters(f->GetParameter(3), f->GetParameter(4), f->GetParameter(5));
        const Double_t mip = f_h->GetMaximumX(lo, hi);
        const Double_t i_h = f_h->Integral(lo, hi), i_s = f_s->Integral(lo, hi);
        const Double_t hump_frac = (i_h + i_s > 0.0) ? i_h/(i_h + i_s) : 0.0;
        if (hump_frac < min_hump_frac) flag |= 64;
        const Double_t mpv = f->GetParameter(4);
        if (std::fabs(mpv - mpv_lo) < 1e-3*hi || std::fabs(mpv - mpv_hi) < 1e-3*hi) flag |= 8;

        result.par.push_back(f->GetParameter(3)); result.err.push_back(f->GetParError(3));
        result.par.push_back(mip);                result.err.push_back(f->GetParError(4));
        result.par.push_back(f->GetParameter(5)); result.err.push_back(f->GetParError(5));
        result.chi_square = f->GetChisquare();
        result.ndf        = f->GetNDF();
        result.additional = { static_cast<Double_t>(flag), static_cast<Double_t>(n_rebin), cut,
                              f->GetParameter(1), f->GetParameter(2) };

        // -- draw: selected pad (log) and close-up (linear) -----
        f->SetLineColor(kOrange);
        f->SetLineWidth(2);
        f_s->SetLineColor(kBlue);  f_s->SetLineStyle(2); f_s->SetNpx(1000);
        f_h->SetLineColor(kGreen+2); f_h->SetLineStyle(2); f_h->SetNpx(1000);
        for (Int_t k = 0; k < 2; k++) {
            c->cd(n_c + k);
            gPad->SetLogy(k == 0 ? 1 : 0);
            TH1D *hd = (TH1D*)h->Clone(Form("%s_%s", h->GetName(), k == 0 ? "log" : "zoom"));
            hd->SetTitle(Form("%s (weak side, strong side near MIP)%s", h_cond->GetTitle(), k == 0 ? "" : " (close-up)"));
            hd->GetXaxis()->SetRangeUser(k == 0 ? 0.0 : std::max(0.0, ped - 30.0), k == 0 ? 2048.0 : ped + 4.0*u);
            if (k == 1) hd->SetMinimum(0.0);
            hd->Draw("HIST");
            f->Draw("same");
            f_s->Draw("same");
            f_h->Draw("same");
            const Double_t y = hd->GetMaximum();
            for (auto [x, col] : {std::pair<Double_t, Int_t>{ped, kMagenta}, {mip, kRed}}) {
                TLine *l = new TLine(x, 0, x, y);
                l->SetLineStyle(2);
                l->SetLineColor(col);
                l->Draw("same");
            }
            TLatex *t = new TLatex();
            t->SetNDC();
            t->SetTextSize(k == 0 ? 0.04 : 0.05);
            t->DrawLatex(0.45, 0.82, Form("MIP %.1f (hump), flag %d", mip, flag));
            t->DrawLatex(0.45, 0.75, Form("spike %.0f, hump frac %.2f", f->GetParameter(1), hump_frac));
            if (k == 1 && result.ndf > 0) t->DrawLatex(0.45, 0.68, Form("#chi^{2}/ndf %.1f/%d", result.chi_square, result.ndf));
        }
        c->Update();
        return result;
    }

    // ____________________________________________________________________________________________
    // HTOF PHC (walk) fit, TENTATIVE. HTOF-only: the other counters use the shared phc_fit()
    // (fixed phc_time_window, profile over the whole window).
    //
    // The HTOF TOF-vs-dE band moves by several ns over the dE range and sits on a flat
    // background, so a fixed time window both cuts the low-dE part of the band and flattens
    // the profile with background. Instead:
    //  1) ridge: starting from the densest dE slice near 1 MIP, follow the TOF peak slice by
    //     slice (2 x-bins each) in both directions, searching within +-1 ns of the neighbour
    //     slice's peak (gaus fit, +-0.8 ns around the local maximum)
    //  2) fit range: contiguous slices around the seed whose band content (+-1 ns of the
    //     ridge) is >= dense_ratio (default 10%) of the maximum, clipped to
    //     [phc_de_range_min, phc_de_range_max]
    //  3) preliminary fit of the ridge points, then a TProfile built only from entries within
    //     +-w of that curve (w = max(2 sigma_median, 3 y-bins)), and the final fit on it.
    // Same function/par layout as phc_fit() and the analyzer's HodoPHCParam Type1:
    //   (time0 - time) = -p0/sqrt|dE - p1| + p2   (analyzer: ctime = time - p0/sqrt|dE-p1| + p2)
    // If ref is given (e.g. a fit of a higher-statistics sample), its parameters are the
    // initial values and those with bit i set in fix_mask are fixed to ref->par[i]; the
    // fit range and band are still determined from h itself.
    // additional = {fit_lo, fit_hi, band_half_width}.
    FitResult htof_phc_fit(TH2D *h, TCanvas *c, Int_t n_c, const FitResult *ref, Int_t fix_mask,
                           Double_t dense_ratio) {
        Config& conf = Config::getInstance();
        FitResult result;
        result.chi_square = 0.0;
        result.ndf = 0;
        result.migrad_stats = 0;

        const Int_t nx = h->GetNbinsX();
        const Double_t ybw = h->GetYaxis()->GetBinWidth(1);
        const Double_t y_axis_lo = h->GetYaxis()->GetXmin(), y_axis_hi = h->GetYaxis()->GetXmax();
        const Int_t slice = 2;
        const Int_t n_sl = nx / slice;

        struct RidgePt { Double_t x = 0, x_lo = 0, x_hi = 0, y = 0, ey = 0, s = 0, n = 0; Bool_t ok = false; };
        std::vector<RidgePt> pts(n_sl);

        // TOF peak of slice isl, searched in [y_c - y_half, y_c + y_half] (whole axis if y_half < 0)
        auto slice_fit = [&](Int_t isl, Double_t y_c, Double_t y_half) -> Bool_t {
            RidgePt &p = pts[isl];
            const Int_t bx0 = isl*slice + 1, bx1 = bx0 + slice - 1;
            p.x_lo = h->GetXaxis()->GetBinLowEdge(bx0);
            p.x_hi = h->GetXaxis()->GetBinUpEdge(bx1);
            p.x = 0.5*(p.x_lo + p.x_hi);
            TH1D *py = h->ProjectionY(Form("%s_sl%d", h->GetName(), isl), bx0, bx1);
            const Double_t lo = (y_half < 0) ? y_axis_lo : std::max(y_axis_lo, y_c - y_half);
            const Double_t hi = (y_half < 0) ? y_axis_hi : std::min(y_axis_hi, y_c + y_half);
            Int_t b_lo = py->FindBin(lo), b_hi = py->FindBin(hi);
            Int_t pb = b_lo;
            for (Int_t b = b_lo; b <= b_hi; ++b) {
                // 3-bin sum to be robust against single-bin fluctuations
                auto s3 = [&](Int_t k) { return py->GetBinContent(k-1) + py->GetBinContent(k) + py->GetBinContent(k+1); };
                if (s3(b) > s3(pb)) pb = b;
            }
            const Double_t pk = py->GetBinCenter(pb);
            const Double_t n_core = py->Integral(py->FindBin(pk - 1.0), py->FindBin(pk + 1.0));
            p.ok = false;
            if (n_core >= 30.0) {
                TF1 g("htof_phc_slice_gaus", "gaus", pk - 0.8, pk + 0.8);
                g.SetParameters(py->GetBinContent(pb), pk, 0.4);
                g.SetParLimits(1, pk - 0.8, pk + 0.8);
                g.SetParLimits(2, 0.5*ybw, 1.5);
                py->Fit(&g, "0QRB", "", pk - 0.8, pk + 0.8);
                p.y  = g.GetParameter(1);
                p.ey = std::max(g.GetParError(1), 1e-3);
                p.s  = g.GetParameter(2);
                p.n  = n_core;
                p.ok = true;
            }
            delete py;
            return p.ok;
        };

        // 1) seed slice: most entries with dE in [0.8, 1.3] MIP, then track outward
        Int_t seed = -1; Double_t seed_n = -1.0;
        for (Int_t isl = 0; isl < n_sl; ++isl) {
            const Double_t xc = h->GetXaxis()->GetBinUpEdge(isl*slice + 1);
            if (xc < 0.8 || xc > 1.3) continue;
            const Double_t n = h->Integral(isl*slice + 1, isl*slice + slice, 1, h->GetNbinsY());
            if (n > seed_n) { seed_n = n; seed = isl; }
        }
        Bool_t ridge_ok = (seed >= 0) && slice_fit(seed, 0.0, -1.0);
        if (ridge_ok) {
            for (Int_t dir : {+1, -1}) {
                Double_t y_prev = pts[seed].y;
                Int_t n_fail = 0;
                for (Int_t isl = seed + dir; isl >= 0 && isl < n_sl; isl += dir) {
                    if (slice_fit(isl, y_prev, 1.0)) { y_prev = pts[isl].y; n_fail = 0; }
                    else if (++n_fail >= 3) break;
                }
            }
        }

        // 2) fit range: contiguous dense slices around the seed
        Double_t fit_lo = conf.phc_de_range_min, fit_hi = conf.phc_de_range_max;
        std::vector<Double_t> sig;
        if (ridge_ok) {
            Double_t n_max = 0.0;
            for (const auto &p : pts) if (p.ok) n_max = std::max(n_max, p.n);
            auto dense = [&](Int_t isl) { return pts[isl].ok && pts[isl].n >= dense_ratio*n_max; };
            Int_t il = seed, ir = seed;
            while (il - 1 >= 0 && dense(il - 1)) --il;
            while (ir + 1 < n_sl && dense(ir + 1)) ++ir;
            fit_lo = std::max(pts[il].x_lo, conf.phc_de_range_min);
            fit_hi = std::min(pts[ir].x_hi, conf.phc_de_range_max);
            for (Int_t isl = il; isl <= ir; ++isl) if (pts[isl].ok) sig.push_back(pts[isl].s);
        }
        Double_t s_med = 0.4;
        if (!sig.empty()) { std::sort(sig.begin(), sig.end()); s_med = sig[sig.size()/2]; }
        const Double_t w = std::max(2.0*s_med, 3.0*ybw);

        // 3) preliminary fit of the ridge, then band-restricted TProfile and final fit
        TF1 *f_fit = new TF1(Form("phc_%s", h->GetName()), "-[0]/TMath::Sqrt(TMath::Abs(x-[1]))+[2]", fit_lo, fit_hi);
        const std::vector<std::vector<Double_t>> par_limits = {
            {0.001, 15.0},
            {-5.0, fit_lo},
            {-10.0, 10.0}
        };
        auto init_fit = [&]() {
            for (Int_t i = 0; i < 3; i++) {
                f_fit->SetParameter(i, (3.0*par_limits[i][0]+par_limits[i][1])/4.0);
                f_fit->SetParLimits(i, par_limits[i][0], par_limits[i][1]);
            }
            f_fit->SetParameter(1, fit_lo - (fit_hi - fit_lo)*0.005);
            if (ref && ref->par.size() >= 3) {
                for (Int_t i = 0; i < 3; i++) {
                    if (fix_mask & (1 << i)) {
                        f_fit->FixParameter(i, ref->par[i]);
                    } else {
                        // keep the initial value inside the limits of this fit
                        const Double_t v = std::min(std::max(ref->par[i], par_limits[i][0]), par_limits[i][1]);
                        f_fit->SetParameter(i, v);
                    }
                }
            }
        };
        init_fit();

        TGraphErrors *g_ridge = new TGraphErrors();
        for (const auto &p : pts) {
            if (!p.ok || p.x < fit_lo || p.x > fit_hi) continue;
            const Int_t k = g_ridge->GetN();
            g_ridge->SetPoint(k, p.x, p.y);
            g_ridge->SetPointError(k, 0.0, p.ey);
        }
        if (g_ridge->GetN() >= 4) g_ridge->Fit(f_fit, "0QR", "", fit_lo, fit_hi);

        TProfile *pf = new TProfile(Form("profile_%s", h->GetName()), "",
                                    nx, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        pf->SetDirectory(nullptr);
        for (Int_t bx = h->GetXaxis()->FindBin(fit_lo); bx <= h->GetXaxis()->FindBin(fit_hi); ++bx) {
            const Double_t xc = h->GetXaxis()->GetBinCenter(bx);
            const Double_t yr = f_fit->Eval(xc);
            for (Int_t by = 1; by <= h->GetNbinsY(); ++by) {
                const Double_t yc = h->GetYaxis()->GetBinCenter(by);
                const Double_t cnt = h->GetBinContent(bx, by);
                if (cnt > 0 && std::fabs(yc - yr) < w) pf->Fill(xc, yc, cnt);
            }
        }
        pf->SetLineColor(kRed);
        pf->SetMarkerColor(kRed);
        pf->SetMarkerStyle(20);
        pf->SetMarkerSize(0.5);
        pf->Fit(f_fit, "0QEMR", "", fit_lo, fit_hi);

        for (Int_t i = 0, n_par = f_fit->GetNpar(); i < n_par; i++) {
            result.par.push_back(f_fit->GetParameter(i));
            result.err.push_back(f_fit->GetParError(i));
        }
        result.chi_square = f_fit->GetChisquare();
        result.ndf        = f_fit->GetNDF();
        result.additional = { fit_lo, fit_hi, w };

        // -- draw -----
        c->cd(n_c);
        Double_t y_min = f_fit->Eval(fit_lo), y_max = f_fit->Eval(fit_hi);
        if (y_min > y_max) std::swap(y_min, y_max);
        h->GetYaxis()->SetRangeUser(std::max(y_axis_lo, y_min - 3.0), std::min(y_axis_hi, y_max + 3.0));
        h->Draw("colz");
        f_fit->SetLineColor(kOrange);
        f_fit->SetLineWidth(2);
        f_fit->SetNpx(500);
        for (Double_t sgn : {+1.0, -1.0}) {
            TF1 *f_band = new TF1(Form("%s_band%+d", f_fit->GetName(), (Int_t)sgn),
                                  Form("-[0]/TMath::Sqrt(TMath::Abs(x-[1]))+[2]%+f", sgn*w), fit_lo, fit_hi);
            for (Int_t i = 0; i < 3; i++) f_band->SetParameter(i, f_fit->GetParameter(i));
            f_band->SetLineColor(kGray+1);
            f_band->SetLineStyle(2);
            f_band->SetLineWidth(1);
            f_band->Draw("same");
        }
        pf->Draw("same");
        f_fit->Draw("same");
        for (Double_t x : {fit_lo, fit_hi}) {
            TLine *l = new TLine(x, std::max(y_axis_lo, y_min - 3.0), x, std::min(y_axis_hi, y_max + 3.0));
            l->SetLineStyle(3);
            l->SetLineColor(kBlack);
            l->Draw("same");
        }
        TLatex *comment = new TLatex();
        comment->SetNDC();
        comment->SetTextSize(0.045);
        comment->DrawLatex(0.20, 0.85, Form("#chi^{2}/ndf %.1f/%d, band #pm%.2f ns", result.chi_square, result.ndf, w));
        comment->DrawLatex(0.20, 0.79, Form("p0 %.3f%s p1 %.3f%s p2 %.3f%s",
                                            result.par[0], (ref && (fix_mask & 1)) ? "(fix)" : "",
                                            result.par[1], (ref && (fix_mask & 2)) ? "(fix)" : "",
                                            result.par[2], (ref && (fix_mask & 4)) ? "(fix)" : ""));
        c->Update();

        delete g_ridge;
        return result;
    }

    // ____________________________________________________________________________________________
    FitResult pedestal_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        gPad->SetLogy(1);
        std::vector<Double_t> par, err;

        h->GetXaxis()->SetRangeUser(
            h->GetXaxis()->GetXmin(),
            conf.hdprm_pedestal_range_right
        );

        // -- pedestal -----
        // Use a smoothed clone to find the robust maximum position (avoids narrow noise spikes)
        TH1D *h_tmp = (TH1D*)h->Clone("h_tmp_robust");
        h_tmp->Smooth(5); 
        Double_t ped_pos = h_tmp->GetBinCenter(h_tmp->GetMaximumBin());
        delete h_tmp;

        Double_t ped_half_width = 10.0;
        std::pair<Double_t, Double_t> ped_n_sigma(2.0, 2.0);

        // -- first fit (simple Gaussian for seeding) -----
        TF1 *f_prefit = new TF1("pre_fit_gauss", "gausn", ped_pos-ped_half_width, ped_pos+ped_half_width);
        f_prefit->SetParameter(1, ped_pos);
        f_prefit->SetParameter(2, ped_half_width*0.5);
        h->Fit(f_prefit, "0QEMR", "", ped_pos-ped_half_width, ped_pos+ped_half_width);
        for (Int_t i = 0; i < 3; i++) par.push_back(f_prefit->GetParameter(i));

        // -- candidate 1: Gaussian + background -----
        Double_t fit_min = par[1]-ped_n_sigma.first*par[2];
        Double_t fit_max = par[1]+ped_n_sigma.second*par[2];
        TF1 *f_gauss = new TF1(Form("ped_g_%s", h->GetName()), 
            [](double *x, double *p){ return p[0]*TMath::Gaus(x[0], p[1], p[2], true) + p[3]; }, 
            fit_min, fit_max, 4);
        f_gauss->SetParameters(par[0], par[1], par[2]*0.9, 1.0);
        f_gauss->SetParLimits(3, 0.0, h->GetMaximum());
        h->Fit(f_gauss, "0QEMR", "", fit_min, fit_max);

        // -- candidate 2: Landau + background -----
        TF1 *f_landau = new TF1(Form("ped_l_%s", h->GetName()), 
            [](double *x, double *p){ return p[0]*TMath::Landau(x[0], p[1], p[2], true) + p[3]; }, 
            fit_min, fit_max, 4);
        f_landau->SetParameters(par[0], par[1], par[2]*0.5, 1.0); // Landau sigma is usually narrower
        f_landau->SetParLimits(3, 0.0, h->GetMaximum());
        h->Fit(f_landau, "0QEMR", "", fit_min, fit_max);

        // -- selection -----
        TF1 *f_best = (f_gauss->GetChisquare() <= f_landau->GetChisquare()) ? f_gauss : f_landau;
        if (f_best == f_gauss) delete f_landau; else delete f_gauss;

        FitResult result;
        for (Int_t i = 0; i < 4; i++) {
            result.par.push_back(f_best->GetParameter(i));
            result.err.push_back(f_best->GetParError(i));
        }
        result.chi_square = f_best->GetChisquare();
        
        // Use actual maximum of the function as the representative peak (user suggestion)
        result.par[1] = f_best->GetMaximumX(fit_min, fit_max);

        // -- draw figure -----
        h->GetXaxis()->SetRangeUser(
            result.par[1] - 10.0*result.par[2], 
            result.par[1] + 10.0*result.par[2]
        );
        h->Draw();

        f_best->SetLineColor(kOrange);
        f_best->SetLineWidth(2.0);
        f_best->SetNpx(1000);
        f_best->Draw("same");
        
        TLine *ped_line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        ped_line->SetLineStyle(2);
        ped_line->SetLineColor(kRed);
        ped_line->Draw("same");

        delete f_prefit;
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

            // Use actual maximum of the function as the representative peak (user suggestion)
            result.par[1] = f_fit_l->GetMaximumX(result.par[1]- 5.0*result.par[2], result.par[1]+ 10.0*result.par[2]);

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

        h = (TH1D*)h->Rebin(3, h->GetName());
        h->GetXaxis()->SetRangeUser(
            conf.tdc_search_range[conf.detector.Data()].first, 
            conf.tdc_search_range[conf.detector.Data()].second
        );

        c->cd(n_c);
        std::vector<Double_t> par, err;
        Double_t max_bin_content = h->GetMaximum(); 
        TString fit_option = max_bin_content > 50.0 ? "0QEMR" : "0QEMRL";
        FitResult result;

        if (max_bin_content > 10.0) {

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
            Double_t evnum_within_range = h->Integral(
                h->FindBin(h->GetXaxis()->GetXmin()),
                h->FindBin(h->GetXaxis()->GetXmax())
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
            Double_t center = h->GetBinCenter(h->GetMaximumBin());
            Double_t sigma = h->GetStdDev();
            Double_t range_width = h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin();
            
            Double_t half_width = 5.0 * sigma;
            if (half_width < range_width * 0.01) half_width = range_width * 0.01;
            if (half_width > range_width * 0.2)  half_width = range_width * 0.2;

            h->GetXaxis()->SetRangeUser(center - half_width, center + half_width);
            h->Draw();
            TLatex* commnet = new TLatex();
            commnet->SetNDC();  // NDC座標（0〜1の正規化）を使う
            commnet->SetTextSize(0.04);
            commnet->DrawLatex(0.7, 0.85, "not fitting");
        }

        TLine *line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        line->SetLineStyle(2); // 点線に設定
        line->SetLineWidth(1.5); // 点線に設定
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