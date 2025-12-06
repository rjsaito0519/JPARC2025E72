// hdprm_fit_tool.C
// ヒストの左右レンジを指定して、その範囲でフィットするインタラクティブツール
//
// 使い方例:
//
// root -l
// .L hdprm_fit_tool.C+
//
// hdprm::set_path(".../run00001_HTOF_HDPRM_Pi.root"); // 入力ROOT
// // 粒子はデフォルトで "Pi" なので省略可
// // hdprm::set_particle("K"); // K にしたい場合だけ設定
//
// hdprm::set_counter("htof-u-0"); // 形式: det-ud-seg or det-seg-ud どちらでもOK
// // hdprm::set_counter("htof-0-u");
// // hdprm::set_counter("bh2-u-1");
//
// hdprm::set_range(420, 650);
// hdprm::fit("auto"); // "gaus", "landau", "auto"

#include <iostream>
#include <cstdio>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TF1.h"
#include "TString.h"

namespace hdprm {

struct State {
    TString path;        // 入力 ROOT ファイル
    TString particle;    // "Pi" (デフォルト), "K" など
    TString counter;     // "htof-u-0" / "htof-0-u" / "bh2-u-1" ...
    TString histName;    // 実際のヒスト名 (HTOF_ADC_seg%dX_%s など)
    Double_t range_left;  // ユーザー指定: 左
    Double_t range_right; // ユーザー指定: 右 (0 or <= left のときは xmax を使う)
    Int_t    rebin;       // Rebin ファクタ
};

static State    g;
static bool     gInit   = false;
static TFile*   gFile   = nullptr;
static TCanvas* gCanvas = nullptr;

// 初期化（一度だけ呼ばれる）
void ensure_init()
{
    if (gInit) return;
    g.path        = "";
    g.particle    = "Pi"; // デフォルト粒子名
    g.counter     = "";
    g.histName    = "";
    g.range_left  = 0.0;
    g.range_right = 0.0;
    g.rebin       = 1;
    gInit         = true;
}

// ------------------------------------------------------
// counter 文字列パーサ
//   対応形式:
//     1) <det>-<ud>-<seg>  (例: "htof-u-0", "bh2-d-1")
//     2) <det>-<seg>-<ud>  (例: "htof-0-u")  ← 互換用
//   -> det, ud, seg に分解
// ------------------------------------------------------
struct DetectorKey {
    TString det;   // "htof", "bh2", ...
    TString ud;    // "u", "d", "s"
    int     seg;   // segment number
    bool    valid;
};

DetectorKey parse_detector_key(const TString& counter)
{
    ensure_init();

    DetectorKey res;
    res.valid = false;

    char det[32];
    char ud[32];
    int seg;

    // パターン1: <det>-<ud>-<seg>
    if (std::sscanf(counter.Data(), "%[^-]-%[^-]-%d", det, ud, &seg) == 3) {
        res.det   = det;
        res.ud    = ud;
        res.seg   = seg;
        res.valid = true;
        // ud は小文字・大文字に関係なく扱うので統一しておく
        res.ud.ToLower();
        return res;
    }

    // パターン2: <det>-<seg>-<ud> (旧形式互換)
    char ud2[32];
    int  seg2;
    if (std::sscanf(counter.Data(), "%[^-]-%d-%[^-]", det, &seg2, ud2) == 3) {
        res.det   = det;
        res.ud    = ud2;
        res.seg   = seg2;
        res.valid = true;
        res.ud.ToLower();
        return res;
    }

    std::cerr << "Error: counter format should be like 'htof-u-0' or 'htof-0-u', got: "
              << counter << std::endl;
    return res;
}

// ------------------------------------------------------
// detector/ud/seg -> ADC ヒスト名
//   ※自分のヒスト名ルールに合わせてここだけ必要に応じて書き換え
// ------------------------------------------------------
TString counter_to_hist_adc(const TString& counter, const TString& particle)
{
    DetectorKey key = parse_detector_key(counter);
    if (!key.valid) return "";

    TString ud = key.ud;
    ud.ToUpper(); // U / D / S

    if (key.det == "htof") {
        // 例: HTOF_ADC_seg%dU_%s
        return Form("HTOF_ADC_seg%d%s_%s", key.seg, ud.Data(), particle.Data());
    }
    else if (key.det == "bh2") {
        // 例: BH2_ADC_seg%dU_%s （ここは環境に合わせて調整）
        return Form("BH2_ADC_seg%d%s_%s", key.seg, ud.Data(), particle.Data());
    }
    else {
        // その他: 共通フォーマット "<DET>_ADC_seg%dX_%s" に従うならここを使う
        return Form("%s_ADC_seg%d%s_%s", key.det.Data(), key.seg, ud.Data(), particle.Data());
    }
}

// ------------------------------------------------------
// set_* 系
// ------------------------------------------------------
void set_path(const char* path)
{
    ensure_init();

    if (gFile) {
        gFile->Close();
        delete gFile;
        gFile = nullptr;
    }
    g.path = path;
    std::cout << "[hdprm] set_path: " << g.path << std::endl;

    gFile = TFile::Open(g.path, "READ");
    if (!gFile || gFile->IsZombie()) {
        std::cerr << "Error: cannot open file: " << g.path << std::endl;
        gFile = nullptr;
    }
}

// 粒子名を設定。デフォルトは "Pi" なので、Pi だけなら呼ばなくてOK
void set_particle(const char* particle)
{
    ensure_init();
    g.particle = particle;
    std::cout << "[hdprm] set_particle: " << g.particle << std::endl;
}

// カウンタ名を設定 ("htof-u-0" / "htof-0-u" / "bh2-u-1" など)
void set_counter(const char* counter)
{
    ensure_init();
    g.counter  = counter;
    g.histName = counter_to_hist_adc(g.counter, g.particle);
    std::cout << "[hdprm] set_counter: " << g.counter
              << " -> histName = " << g.histName << std::endl;
}

// 左端の範囲
void set_range_left(Double_t left)
{
    ensure_init();
    g.range_left = left;
    std::cout << "[hdprm] set_range_left: " << g.range_left << std::endl;
}

// 右端の範囲
void set_range_right(Double_t right)
{
    ensure_init();
    g.range_right = right;
    std::cout << "[hdprm] set_range_right: " << g.range_right << std::endl;
}

// 左右まとめて設定
void set_range(Double_t left, Double_t right)
{
    ensure_init();
    g.range_left  = left;
    g.range_right = right;
    std::cout << "[hdprm] set_range: [" << g.range_left
              << ", " << g.range_right << "]" << std::endl;
}

// Rebin
void set_rebin(Int_t r)
{
    ensure_init();
    g.rebin = (r > 0 ? r : 1);
    std::cout << "[hdprm] set_rebin: " << g.rebin << std::endl;
}

// ------------------------------------------------------
// 内部: fit 範囲内でのピーク位置探し
// ------------------------------------------------------
Double_t find_peak_x(TH1D* h, Double_t xleft, Double_t xright)
{
    Int_t bin_min = h->FindBin(xleft);
    Int_t bin_max = (xright > xleft) ? h->FindBin(xright) : h->GetNbinsX();

    Double_t max_val = -1.0;
    Int_t max_bin = bin_min;

    for (Int_t ib = bin_min; ib <= bin_max; ++ib) {
        Double_t c = h->GetBinContent(ib);
        if (c > max_val) {
            max_val = c;
            max_bin = ib;
        }
    }
    return h->GetBinCenter(max_bin);
}

// ------------------------------------------------------
// メイン: フィット実行
//   model = "gaus", "landau", "auto"
// ------------------------------------------------------
void fit(const char* model = "gaus", Bool_t logy = kTRUE)
{
    ensure_init();

    if (!gFile) {
        std::cerr << "Error: file not opened. Call hdprm::set_path() first." << std::endl;
        return;
    }
    if (g.histName == "") {
        std::cerr << "Error: histName not set. Call hdprm::set_counter() first." << std::endl;
        return;
    }

    TH1D* h_orig = dynamic_cast<TH1D*>(gFile->Get(g.histName));
    if (!h_orig) {
        std::cerr << "Error: histogram not found: " << g.histName << std::endl;
        return;
    }

    // Rebin 用に clone
    TH1D* h = (TH1D*)h_orig->Clone(Form("%s_reb", h_orig->GetName()));
    if (g.rebin > 1) {
        h = (TH1D*)h->Rebin(g.rebin, h->GetName());
    }

    Double_t xmin = h->GetXaxis()->GetXmin();
    Double_t xmax = h->GetXaxis()->GetXmax();

    Double_t fit_left  = (g.range_left  > xmin) ? g.range_left  : xmin;
    Double_t fit_right = (g.range_right > fit_left) ? g.range_right : xmax;

    // 範囲内でピーク位置を探して初期値に使う
    Double_t peakX = find_peak_x(h, fit_left, fit_right);
    Double_t width_guess = (fit_right - fit_left) / 6.0;
    if (width_guess <= 0) width_guess = h->GetStdDev();

    // キャンバス準備
    if (!gCanvas) {
        gCanvas = new TCanvas("c_hdprm_fit", "hdprm fit tool", 800, 600);
    }
    gCanvas->cd();
    if (logy) gPad->SetLogy(1);
    else      gPad->SetLogy(0);

    // 描画
    h->GetXaxis()->SetRangeUser(fit_left - (fit_right-fit_left)*0.2,
                                fit_right + (fit_right-fit_left)*0.2);
    h->Draw();

    // フィット関数（gaus/landau/auto）
    TF1* f_gaus   = nullptr;
    TF1* f_landau = nullptr;

    // Gauss
    f_gaus = new TF1("f_user_gaus", "gaus", fit_left, fit_right);
    f_gaus->SetParameter(0, h->GetMaximum()); // amplitude
    f_gaus->SetParameter(1, peakX);           // mean
    f_gaus->SetParameter(2, width_guess);     // sigma

    // Landau
    f_landau = new TF1("f_user_landau", "landaun", fit_left, fit_right);
    f_landau->SetParameter(0, h->GetMaximum()); // amplitude
    f_landau->SetParameter(1, peakX);           // MPV
    f_landau->SetParameter(2, width_guess);     // width

    TF1* f_chosen = nullptr;

    TString m(model);
    m.ToLower();

    if (m == "gaus") {
        h->Fit(f_gaus, "QEMR");
        f_chosen = f_gaus;
    } else if (m == "landau") {
        h->Fit(f_landau, "QEMR");
        f_chosen = f_landau;
    } else if (m == "auto") {
        // 両方フィットして χ²/ndf が小さい方を採用
        h->Fit(f_gaus, "QEMR");
        Double_t chi2_g = f_gaus->GetChisquare();
        Double_t ndf_g  = f_gaus->GetNDF();

        h->Fit(f_landau, "QEMR");
        Double_t chi2_l = f_landau->GetChisquare();
        Double_t ndf_l  = f_landau->GetNDF();

        Double_t val_g = (ndf_g > 0 ? chi2_g/ndf_g : 1e9);
        Double_t val_l = (ndf_l > 0 ? chi2_l/ndf_l : 1e9);

        if (val_g <= val_l) f_chosen = f_gaus;
        else                f_chosen = f_landau;
    } else {
        std::cerr << "Error: unknown model '" << model
                  << "', use 'gaus', 'landau', or 'auto'." << std::endl;
        delete f_gaus;
        delete f_landau;
        return;
    }

    // 塗りつぶし（指定レンジ）
    Double_t ymax = h->GetMaximum();

    TBox* box = new TBox(fit_left, 0, fit_right, ymax);
    box->SetFillColor(800); // 適当な色番号
    box->SetFillStyle(3004);
    box->Draw("same");

    // フィット関数を太線で
    f_chosen->SetLineWidth(2);
    f_chosen->Draw("same");

    // テキスト (counter, range, モデル名)
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.85,
                    Form("%s  [%.1f, %.1f]  (%s, %s)",
                         g.counter.Data(), fit_left, fit_right,
                         m.Data(), g.particle.Data()));

    gCanvas->Update();

    // 結果出力
    Double_t chi2  = f_chosen->GetChisquare();
    Double_t ndf   = f_chosen->GetNDF();
    Double_t p0    = f_chosen->GetParameter(0);
    Double_t p1    = f_chosen->GetParameter(1);
    Double_t p2    = f_chosen->GetParameter(2);
    Double_t e0    = f_chosen->GetParError(0);
    Double_t e1    = f_chosen->GetParError(1);
    Double_t e2    = f_chosen->GetParError(2);

    std::cout << "========================================" << std::endl;
    std::cout << " counter   : " << g.counter << std::endl;
    std::cout << " histName  : " << g.histName << std::endl;
    std::cout << " particle  : " << g.particle << std::endl;
    std::cout << " model     : " << m << std::endl;
    std::cout << " fit range : [" << fit_left << ", " << fit_right << "]" << std::endl;
    std::cout << " chi2/ndf  : " << chi2 << " / " << ndf
              << " = " << ( (ndf>0)?chi2/ndf:0 ) << std::endl;
    std::cout << " p0 (amp)  : " << p0 << " ± " << e0 << std::endl;
    if (m == "gaus") {
        std::cout << " p1 (mean) : " << p1 << " ± " << e1 << std::endl;
        std::cout << " p2 (sigma): " << p2 << " ± " << e2 << std::endl;
    } else {
        std::cout << " p1 (MPV)  : " << p1 << " ± " << e1 << std::endl;
        std::cout << " p2 (width): " << p2 << " ± " << e2 << std::endl;
    }
    std::cout << "========================================" << std::endl;
}

} // namespace hdprm
