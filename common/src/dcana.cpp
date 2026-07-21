#include "ana_helper.h"

#include <cmath>
#include <iostream>

namespace ana_helper {

namespace {

constexpr Double_t kTdcFitWidth = 80.0;

std::pair<Double_t, Double_t> dc_tdc_default_range(TH1D *h)
{
    Double_t fit_range_min = h->GetXaxis()->GetBinCenter(h->GetMaximumBin()) - 5.0;
    return {fit_range_min, fit_range_min + kTdcFitWidth};
}

TF1 *make_dc_tdc_tf1(const char *name, Double_t fit_range_min, Double_t fit_range_max)
{
    return new TF1(name, "[0]*TMath::Erfc( (x-[1])/[2] ) + [3]", fit_range_min, fit_range_max);
}

void set_dc_tdc_init_params(TF1 *fit_f, TH1D *h, const DcTdcFitSeed *seed)
{
    fit_f->SetParameter(0, h->GetMaximum() / 2.0);
    if (seed && seed->valid) {
        fit_f->SetParameter(1, seed->par1);
        fit_f->SetParameter(2, seed->par2);
        fit_f->SetParameter(3, seed->par3);
    } else {
        fit_f->SetParameter(1, h->GetXaxis()->GetBinCenter(h->GetMaximumBin()) + 10.0);
        fit_f->SetParameter(2, 10.0);
        fit_f->SetParameter(3, 0.0);
    }
}

void clear_hist_fit_functions(TH1D *h)
{
    if (!h) return;
    TObject* obj = nullptr;
    while ((obj = h->GetListOfFunctions()->First())) {
        h->GetListOfFunctions()->Remove(obj);
        delete obj;
    }
}

Double_t dc_tdc_fixed_amplitude(TH1D *h, Double_t bg)
{
    Double_t signal_max = h->GetMaximum() - bg;
    if (signal_max < 1.0) signal_max = h->GetMaximum();
    return signal_max / 2.0;
}

TF1 *run_dc_tdc_fit(TH1D *h, TF1 *fit_f, const char *fit_option,
                    Double_t fit_range_min, Double_t fit_range_max)
{
    h->Fit(fit_f, fit_option, "", fit_range_min, fit_range_max);
    return fit_f;
}

Double_t dc_tdc_t0_from_fit(TF1 *fit_f, Double_t fit_range_max)
{
    std::vector<Double_t> par(4);
    for (Int_t i = 0; i < 4; i++) {
        par[i] = fit_f->GetParameter(i);
    }

    Double_t target_value_ratio = 0.005;
    Double_t target_value = 2.0 * par[0] * target_value_ratio;
    ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
    ROOT::Math::Functor1D erfc_func([=](Double_t x) {
        return par[0] * TMath::Erfc((x - par[1]) / par[2]) - target_value;
    });
    rootFinder.SetFunction(erfc_func, par[1], fit_range_max);
    return rootFinder.Solve() ? rootFinder.Root() : 0.0;
}

FitResult dc_tdc_fit_impl(TH1D *h, TCanvas *c, Int_t n_c, const DcTdcFitSeed *seed, Bool_t draw)
{
    FitResult result;
    if (!h || h->GetEntries() <= 0) {
        result.par.assign(4, 0.0);
        result.err.assign(4, 0.0);
        result.chi_square = 0.0;
        result.ndf = 0;
        result.additional.assign(1, 0.0);
        return result;
    }

    Double_t fit_range_min = 0.0;
    Double_t fit_range_max = 0.0;
    if (seed && seed->valid) {
        fit_range_min = seed->range_min;
        fit_range_max = seed->range_max;
    } else {
        std::tie(fit_range_min, fit_range_max) = dc_tdc_default_range(h);
    }

    if (draw) {
        c->cd(n_c);
    }

    TString fit_option = h->GetMaximum() > 500.0 ? "0QEMR" : "0QEMRL";

    // Stage 1: free 4-parameter fit to estimate shape / background.
    TF1 *prefit = make_dc_tdc_tf1(Form("dc_t0_%s_prefit", h->GetName()),
                                  fit_range_min, fit_range_max);
    set_dc_tdc_init_params(prefit, h, seed);
    run_dc_tdc_fit(h, prefit, fit_option.Data(), fit_range_min, fit_range_max);

    Double_t pre_par[4] = {0.0, 0.0, 10.0, 0.0};
    for (Int_t i = 0; i < 4; ++i) {
        pre_par[i] = prefit->GetParameter(i);
    }

    // Stage 2: fix amplitude at (peak - background) / 2, inherit other inits from stage 1.
    clear_hist_fit_functions(h);
    Double_t amp_fixed = dc_tdc_fixed_amplitude(h, pre_par[3]);
    TF1 *fit_f = make_dc_tdc_tf1(Form("dc_t0_%s", h->GetName()), fit_range_min, fit_range_max);
    fit_f->SetParameter(0, amp_fixed);
    fit_f->SetParameter(1, pre_par[1]);
    fit_f->SetParameter(2, pre_par[2]);
    fit_f->SetParameter(3, pre_par[3]);
    fit_f->FixParameter(0, amp_fixed);
    fit_f->SetLineColor(kOrange);
    fit_f->SetLineWidth(1);
    run_dc_tdc_fit(h, fit_f, fit_option.Data(), fit_range_min, fit_range_max);

    for (Int_t i = 0; i < 4; i++) {
        result.par.push_back(fit_f->GetParameter(i));
        result.err.push_back(fit_f->GetParError(i));
    }
    result.chi_square = fit_f->GetChisquare();
    result.ndf = fit_f->GetNDF();

    Double_t t0 = dc_tdc_t0_from_fit(fit_f, fit_range_max);
    result.additional.push_back(t0);

    if (draw) {
        Double_t xrange_min = fit_range_min - kTdcFitWidth / 2.0;
        Double_t xrange_max = fit_range_max;
        if (t0 > 0.0 && std::isfinite(t0)) {
            if (t0 < xrange_min) xrange_min = t0 - 10.0;
            if (t0 > xrange_max) xrange_max = t0 + 10.0;
        } else {
            std::cerr << "Warning: dc_tdc_fit could not determine t0 for "
                      << h->GetName() << std::endl;
        }

        h->GetXaxis()->SetRangeUser(xrange_min, xrange_max);
        h->Draw();
        fit_f->Draw("same");
        c->Update();

        if (t0 > 0.0 && std::isfinite(t0)) {
            TLine *line = new TLine(t0, gPad->GetUymin(), t0, gPad->GetUymax());
            line->SetLineStyle(2);
            line->SetLineColor(kRed);
            line->Draw("same");
        }
    }

    // fit_f is owned by the histogram after h->Fit(); do not delete here.
    return result;
}

} // namespace

    // ____________________________________________________________________________________________
    DcTdcFitSeed dc_tdc_seed(TH1D *h)
    {
        DcTdcFitSeed seed;
        if (!h || h->GetEntries() <= 0) {
            return seed;
        }

        FitResult seed_fit = dc_tdc_fit_impl(h, nullptr, 0, nullptr, kFALSE);
        if (seed_fit.par.size() < 4) {
            return seed;
        }

        std::tie(seed.range_min, seed.range_max) = dc_tdc_default_range(h);
        seed.par1 = seed_fit.par[1];
        seed.par2 = seed_fit.par[2];
        seed.par3 = seed_fit.par[3];
        seed.valid = kTRUE;
        return seed;
    }

    // ____________________________________________________________________________________________
    FitResult dc_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c, const DcTdcFitSeed *seed)
    {
        return dc_tdc_fit_impl(h, c, n_c, seed, kTRUE);
    }

    // ____________________________________________________________________________________________
    TGraph* make_drift_function(TH1D *h, TCanvas *c, Int_t n_c, Int_t plane,
                                const char* wire_range_suffix)
    {
        Config& conf = Config::getInstance();
        if (!h) return nullptr;
        if (c) c->cd(n_c);
    
        Double_t max_dt = conf.max_drift_time.at(conf.detector.Data()); // ns
        Double_t max_dl = conf.max_drift_length.at(conf.detector.Data()); // mm

        // 積分範囲を決める
        Int_t bin_min = h->FindBin(-1.0);    // -1 ns 付近
        Int_t bin_max = h->FindBin(max_dt);  // max_dt ns 付近

        Double_t entry = h->Integral(bin_min, bin_max);
        if (entry <= 0) {
            std::cerr << "Integral is zero for " << h->GetName() << std::endl;
            return nullptr;
        }

        std::vector<Double_t> dt;
        std::vector<Double_t> dl;
        dt.reserve(bin_max - bin_min + 1);
        dl.reserve(bin_max - bin_min + 1);

        for (Int_t bin = bin_min; bin <= bin_max; ++bin) {
            Double_t integral = h->Integral(bin_min, bin);
            Double_t t = h->GetBinCenter(bin);      // Drift Time
            Double_t L = max_dl * integral / entry; // Drift Length (max_dl に正規化)

            dt.push_back(t);
            dl.push_back(L);
        }

        Int_t npoints = dt.size();
        if (npoints == 0) {
            std::cerr << "No points for graph of " << h->GetName() << std::endl;
            return nullptr;
        }

        // TGraph のコンストラクタは double* を要求するので data() を渡す
        TGraph* g = new TGraph(npoints, dt.data(), dl.data());

        // 名前とタイトルを設定
        TString gname;
        if (wire_range_suffix && wire_range_suffix[0] != '\0') {
            gname = Form("%s_Hit_DriftFunction_plane%d%s", conf.detector.Data(), plane, wire_range_suffix);
        } else {
            gname = Form("%s_Hit_DriftFunction_plane%d", conf.detector.Data(), plane);
        }
        gname.ReplaceAll("DriftTime", "DriftFunction");  // 名前の一部を置き換え

        TString gtitle;
        if (wire_range_suffix && wire_range_suffix[0] != '\0') {
            gtitle.Form("DriftFunction %s plane%d%s;Drift Time [ns];Drift Length [mm]",
                        conf.detector.Data(), plane, wire_range_suffix);
        } else {
            gtitle.Form("DriftFunction %s plane%d;Drift Time [ns];Drift Length [mm]",
                        conf.detector.Data(), plane);
        }

        g->SetName(gname);
        g->SetTitle(gtitle);

        // 軸範囲などのスタイル
        g->GetXaxis()->SetTitleOffset(1.1);
        g->GetXaxis()->SetLimits(-0.1 * max_dt, 1.2 * max_dt);
        g->GetYaxis()->SetRangeUser(-0.1 * max_dl, 1.2 * max_dl);
        g->SetMarkerStyle(8);
        g->SetMarkerSize(0.4);

        if (c) g->Draw("APC");

        return g;
    }

    // ____________________________________________________________________________________________
    FitResult residual_fit(TH1D *h, TCanvas *c, Int_t n_c)
    {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        std::vector<Double_t> par, err;
        TString fit_option = h->GetMaximum() > 500.0 ? "0QEMR" : "0QEMRL";
        FitResult result;

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

        TLine *line = new TLine(result.par[1], 0, result.par[1], h->GetMaximum());
        line->SetLineStyle(2); // 点線に設定
        line->SetLineColor(kRed); // 色を赤に設定
        line->Draw("same");

        c->Update();

        return result;
    }


}
