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
        Double_t target_value_ratio = 0.0001;
        Double_t target_value = 2.0*result.par[0] * target_value_ratio;
        ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
        ROOT::Math::Functor1D erfc_func([=](Double_t x) { return result.par[0]*TMath::Erfc( (x-result.par[1])/result.par[2] ) - target_value; });
        rootFinder.SetFunction(erfc_func, result.par[1], fit_range_max); // 探索範囲が第2, 3引数
        Bool_t is_success = rootFinder.Solve();
        Double_t t0 = is_success ? rootFinder.Root() : 0.0;
        result.additional.push_back(t0);
        std::cout << t0 << ", " << result.par[1] - 2*result.par[2] << std::endl;

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
}
