#include "ana_helper.h"

namespace ana_helper {

    // ____________________________________________________________________________________________    
    FitResult dc_tdc_fit(TH1D *h, TCanvas *c, Int_t n_c) {
        Config& conf = Config::getInstance();

        c->cd(n_c);
        TString fit_option = h->GetMaximum() > 500.0 ? "0QEMR" : "0QEMRL";
        FitResult result;

        Double_t tdc_width = 50.0;
        Double_t fit_range_min = h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ) - 5.0;
        Double_t fit_range_max = fit_range_min + tdc_width;

        TF1 *fit_f = new TF1( Form("dc_t0_%s", h->GetName()), "[0]*TMath::Erfc( (x-[1])/[2] ) + [3]", fit_range_min, fit_range_max);
        fit_f->SetParameter(0, h->GetMaximum() / 2.0);
        fit_f->SetParameter(1, h->GetXaxis()->GetBinCenter( h->GetMaximumBin() ) + 10.0);
        fit_f->SetParameter(2, 10.0);
        fit_f->SetParameter(3, 0.0);
        fit_f->SetLineColor(kOrange);
        fit_f->SetLineWidth(1);
        h->Fit(fit_f, fit_option.Data(), "", fit_range_min, fit_range_max);

        for (Int_t i = 0; i < 4; i++) {
            result.par.push_back(fit_f->GetParameter(i));
            result.err.push_back(fit_f->GetParError(i));
        }
        result.chi_square = fit_f->GetChisquare();
        result.ndf        = fit_f->GetNDF();

        // 数値的に解を求めるためのRootFinderの設定, t0の値を調べる
        Double_t target_value_ratio = 0.005;
        Double_t target_value = 2.0*result.par[0] * target_value_ratio;
        ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
        ROOT::Math::Functor1D erfc_func([=](Double_t x) { return result.par[0]*TMath::Erfc( (x-result.par[1])/result.par[2] ) - target_value; });
        rootFinder.SetFunction(erfc_func, result.par[1], fit_range_max); // 探索範囲が第2, 3引数
        Bool_t is_success = rootFinder.Solve();
        Double_t t0 = is_success ? rootFinder.Root() : 0.0;
        result.additional.push_back(t0);

        // draw
        h->GetXaxis()->SetRangeUser(fit_range_min - tdc_width/2.0, fit_range_max);
        h->Draw();
        fit_f->Draw("same");

        // 縦の点線を引くためのラインオブジェクトの作成
        // TLine *line = new TLine(par[1]+par[2], 0, par[1]+par[2], h->GetMaximum());
        TLine *line = new TLine(t0, 0.0, t0, h->GetMaximum());
        line->SetLineStyle(2); // 点線に設定
        line->SetLineColor(kRed); // 色を赤に設定
        
        // ラインを描画
        line->Draw("same");

        return result;
    }

    // ____________________________________________________________________________________________    
    TGraph* make_drift_function(TH1D *h, TCanvas *c, Int_t n_c, Int_t plane)
    {
        Config& conf = Config::getInstance();
        if (!h) return nullptr;
        c->cd(n_c);
    
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
        TString gname = Form("%s_Hit_DriftFunction_plane%d", conf.detector.Data(), plane);
        gname.ReplaceAll("DriftTime", "DriftFunction");  // 名前の一部を置き換え

        TString gtitle;
        gtitle.Form("DriftFunction %s plane%d;Drift Time [ns];Drift Length [mm]", conf.detector.Data(), plane);

        g->SetName(gname);
        g->SetTitle(gtitle);

        // 軸範囲などのスタイル
        g->GetXaxis()->SetTitleOffset(1.1);
        g->GetXaxis()->SetLimits(-0.1 * max_dt, 1.2 * max_dt);
        g->GetYaxis()->SetRangeUser(-0.1 * max_dl, 1.2 * max_dl);
        g->SetMarkerStyle(8);
        g->SetMarkerSize(0.4);

        g->Draw("APC");

        return g;
    }

}
