// hdprm_fit_tool.C
// Hist の左右レンジを指定して、その範囲でフィットするインタラクティブツール
//
// 使い方例:
//
// root -l
// .L hdprm_fit_tool.C+
//
// set_path(".../run00001_HTOF_HDPRM_Pi.root"); // 入力ROOT
// // 粒子はデフォルト "Pi" なので省略可
//
// set_counter("htof-u-0");  // det-ud-seg or det-seg-ud
// // set_counter("htof-0-u");
// // set_counter("bh2-u-1");
//
// set_range(420, 650);
// fit_iter(3, "auto");   // 自動レンジ更新付きで最大3回まで回す
// // fit("auto");        // 1回だけ回す場合
//
// // 片側レンジ固定の例:
// // set_range(300, 900);
// // lock_left();         // 左端固定
// // fit_iter(5, "auto");

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
    TString particle;    // 粒子名（デフォルト "Pi"）
    TString counter;     // "htof-u-0" / "htof-0-u" / "bh2-u-1" など
    TString histName;    // 実際のヒスト名
    Double_t range_left;  // フィット左端
    Double_t range_right; // フィット右端
    Int_t    rebin;       // Rebin ファクタ

    // 直近の fit 結果（収束判定用）
    Double_t last_p1;
    Double_t last_p2;
    Bool_t   has_last_fit;

    // 範囲ロックフラグ
    Bool_t   lock_left;
    Bool_t   lock_right;
};

static State    g;
static bool     gInit   = false;
static TFile*   hdFile  = nullptr;   // ROOT の gFile マクロと衝突しないように
static TCanvas* gCanvas = nullptr;

//--------------------------------------------------
// 初期化
//--------------------------------------------------
void ensure_init()
{
    if (gInit) return;
    g.path        = "";
    g.particle    = "Pi";   // デフォルト粒子名
    g.counter     = "";
    g.histName    = "";
    g.range_left  = 0.0;
    g.range_right = 0.0;
    g.rebin       = 1;
    g.last_p1     = 0.0;
    g.last_p2     = 0.0;
    g.has_last_fit = kFALSE;
    g.lock_left   = kFALSE;
    g.lock_right  = kFALSE;
    gInit         = true;
}

//--------------------------------------------------
// カウンタキー (det, ud, seg)
//--------------------------------------------------
struct DetectorKey {
    TString det;   // "htof", "bh2", ...
    TString ud;    // "u", "d", "s"
    int     seg;   // segment number
    bool    valid;
};

// det-ud-seg または det-seg-ud をパース
DetectorKey parse_detector_key(const TString& counter)
{
    ensure_init();

    DetectorKey res;
    res.valid = false;

    char det[32];
    char ud[32];
    int seg;

    // pattern 1: det-ud-seg  (例: htof-u-0, bh2-d-1)
    if (std::sscanf(counter.Data(), "%[^-]-%[^-]-%d", det, ud, &seg) == 3) {
        res.det   = det;
        res.ud    = ud;
        res.seg   = seg;
        res.ud.ToLower();
        res.valid = true;
        return res;
    }

    // pattern 2: det-seg-ud  (例: htof-0-u)
    char ud2[32];
    int  seg2;
    if (std::sscanf(counter.Data(), "%[^-]-%d-%[^-]", det, &seg2, ud2) == 3) {
        res.det   = det;
        res.ud    = ud2;
        res.seg   = seg2;
        res.ud.ToLower();
        res.valid = true;
        return res;
    }

    std::cerr << "Error: counter format should be like 'htof-u-0' or 'htof-0-u', got: "
              << counter << std::endl;
    return res;
}

//--------------------------------------------------
// detector/ud/seg -> ADC ヒスト名
//   必要ならここを自分の環境に合わせて修正
//--------------------------------------------------
TString counter_to_hist_adc(const TString& counter, const TString& particle)
{
    DetectorKey k = parse_detector_key(counter);
    if (!k.valid) return "";

    TString ud = k.ud;
    ud.ToUpper(); // "u" -> "U"

    if (k.det == "htof") {
        // 例: HTOF_ADC_seg%dU_%s
        return Form("HTOF_ADC_seg%d%s_%s", k.seg, ud.Data(), particle.Data());
    }
    else if (k.det == "bh2") {
        // 例: BH2_ADC_seg%dU_%s
        return Form("BH2_ADC_seg%d%s_%s", k.seg, ud.Data(), particle.Data());
    }
    else {
        // その他は汎用形式 "<DET>_ADC_seg%dX_%s" と仮定
        return Form("%s_ADC_seg%d%s_%s", k.det.Data(), k.seg, ud.Data(), particle.Data());
    }
}

//--------------------------------------------------
// 設定系
//--------------------------------------------------
void set_path(const char* path)
{
    ensure_init();

    if (hdFile) {
        hdFile->Close();
        delete hdFile;
        hdFile = nullptr;
    }
    g.path = path;
    std::cout << "[hdprm] set_path: " << g.path << std::endl;

    hdFile = TFile::Open(g.path, "READ");
    if (!hdFile || hdFile->IsZombie()) {
        std::cerr << "Error: cannot open file: " << g.path << std::endl;
        hdFile = nullptr;
    }
}

// 粒子名を設定（デフォルト "Pi"）
void set_particle(const char* particle)
{
    ensure_init();
    g.particle = particle;
    std::cout << "[hdprm] set_particle: " << g.particle << std::endl;
}

// カウンタ名設定 ("htof-u-0", "htof-0-u", "bh2-u-1", ...)
void set_counter(const char* counter)
{
    ensure_init();
    g.counter  = counter;
    g.histName = counter_to_hist_adc(g.counter, g.particle);
    std::cout << "[hdprm] set_counter: " << g.counter
              << " -> histName = " << g.histName << std::endl;
}

void set_range_left(Double_t left)
{
    ensure_init();
    g.range_left = left;
    std::cout << "[hdprm] set_range_left: " << g.range_left << std::endl;
}

void set_range_right(Double_t right)
{
    ensure_init();
    g.range_right = right;
    std::cout << "[hdprm] set_range_right: " << g.range_right << std::endl;
}

void set_range(Double_t left, Double_t right)
{
    ensure_init();
    g.range_left  = left;
    g.range_right = right;
    std::cout << "[hdprm] set_range: [" << g.range_left
              << ", " << g.range_right << "]" << std::endl;
}

void set_rebin(Int_t r)
{
    ensure_init();
    g.rebin = (r > 0 ? r : 1);
    std::cout << "[hdprm] set_rebin: " << g.rebin << std::endl;
}

// 範囲ロック系 ------------------------------------
void lock_left(Bool_t on = kTRUE)
{
    ensure_init();
    g.lock_left = on;
    std::cout << "[hdprm] lock_left: " << (g.lock_left ? "ON" : "OFF") << std::endl;
}

void lock_right(Bool_t on = kTRUE)
{
    ensure_init();
    g.lock_right = on;
    std::cout << "[hdprm] lock_right: " << (g.lock_right ? "ON" : "OFF") << std::endl;
}

void unlock_all()
{
    ensure_init();
    g.lock_left  = kFALSE;
    g.lock_right = kFALSE;
    std::cout << "[hdprm] unlock_all" << std::endl;
}

//--------------------------------------------------
// 指定範囲内でのピーク位置を探す
//--------------------------------------------------
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

//--------------------------------------------------
// フィット本体
//   model = "gaus", "landau", "auto"
//   ・1回だけフィット
//   ・結果から次回用レンジを更新（ロック状態を考慮）
//   ・p1, p2 を State に保存（収束判定用）
//--------------------------------------------------
void fit(const char* model = "auto", Bool_t logy = kTRUE)
{
    ensure_init();

    if (!hdFile) {
        std::cerr << "Error: file not opened. Call set_path() first." << std::endl;
        return;
    }
    if (g.histName == "") {
        std::cerr << "Error: histName not set. Call set_counter() first." << std::endl;
        return;
    }

    TH1D* h_orig = dynamic_cast<TH1D*>(hdFile->Get(g.histName));
    if (!h_orig) {
        std::cerr << "Error: histogram not found: " << g.histName << std::endl;
        return;
    }

    // Rebin 用 clone
    TH1D* h = (TH1D*)h_orig->Clone(Form("%s_reb", h_orig->GetName()));
    if (g.rebin > 1) {
        h = (TH1D*)h->Rebin(g.rebin, h->GetName());
    }

    Double_t xmin = h->GetXaxis()->GetXmin();
    Double_t xmax = h->GetXaxis()->GetXmax();

    Double_t fit_left  = (g.range_left  > xmin) ? g.range_left  : xmin;
    Double_t fit_right = (g.range_right > fit_left) ? g.range_right : xmax;

    // 指定範囲内でピーク位置探索 (初期値用 & 可視化用)
    Double_t peakX = find_peak_x(h, fit_left, fit_right);
    Double_t width_guess = (fit_right - fit_left) / 6.0;
    if (width_guess <= 0) width_guess = h->GetStdDev();

    // Canvas
    if (!gCanvas) {
        gCanvas = new TCanvas("c_hdprm_fit", "hdprm fit tool", 800, 600);
    }
    gCanvas->cd();
    if (logy) gPad->SetLogy(1);
    else      gPad->SetLogy(0);

    // Draw hist
    h->GetXaxis()->SetRangeUser(fit_left - 0.2*(fit_right-fit_left),
                                fit_right + 0.2*(fit_right-fit_left));
    h->Draw();

    // Peak marker (赤い点線)
    TLine *peakLine = new TLine(peakX, 0, peakX, h->GetMaximum());
    peakLine->SetLineColor(kRed);
    peakLine->SetLineStyle(2);
    peakLine->SetLineWidth(2);
    peakLine->Draw("same");

    // Fit funcs
    TString m(model);
    m.ToLower();

    TF1* f_gaus   = new TF1("f_user_gaus",   "gaus",   fit_left, fit_right);
    TF1* f_landau = new TF1("f_user_landau", "landaun", fit_left, fit_right);

    f_gaus->SetParameters(h->GetMaximum(), peakX, width_guess);
    f_landau->SetParameters(h->GetMaximum(), peakX, width_guess);

    TF1* f_best = nullptr;

    if (m == "gaus") {
        h->Fit(f_gaus, "QEMR");
        f_best = f_gaus;
    }
    else if (m == "landau") {
        h->Fit(f_landau, "QEMR");
        f_best = f_landau;
    }
    else if (m == "auto") {
        h->Fit(f_gaus, "QEMR");
        Double_t chi2_g = f_gaus->GetChisquare();
        Double_t ndf_g  = f_gaus->GetNDF();

        h->Fit(f_landau, "QEMR");
        Double_t chi2_l = f_landau->GetChisquare();
        Double_t ndf_l  = f_landau->GetNDF();

        Double_t val_g = (ndf_g>0 ? chi2_g/ndf_g : 1e9);
        Double_t val_l = (ndf_l>0 ? chi2_l/ndf_l : 1e9);

        if (val_g <= val_l) {
            m = "gaus";
        } else {
            m = "landau";
        }
        f_best = (val_g <= val_l ? f_gaus : f_landau);
    }
    else {
        std::cerr << "Error: unknown model '" << model
                  << "', use 'gaus', 'landau', or 'auto'." << std::endl;
        return;
    }

    // フィット範囲塗りつぶし
    Double_t ymax = h->GetMaximum();
    TBox* box = new TBox(fit_left, 0, fit_right, ymax);
    box->SetFillStyle(3004);
    box->SetFillColor(kOrange);
    box->Draw("same");

    // フィット関数
    f_best->SetLineWidth(2);
    f_best->Draw("same");

    // ラベル
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.85,
                    Form("%s  [%.1f, %.1f]  (%s, %s)",
                         g.counter.Data(), fit_left, fit_right,
                         m.Data(), g.particle.Data()));

    gCanvas->Update();

    // 結果出力
    Double_t chi2  = f_best->GetChisquare();
    Double_t ndf   = f_best->GetNDF();
    Double_t p0    = f_best->GetParameter(0);
    Double_t p1    = f_best->GetParameter(1);
    Double_t p2    = f_best->GetParameter(2);
    Double_t e0    = f_best->GetParError(0);
    Double_t e1    = f_best->GetParError(1);
    Double_t e2    = f_best->GetParError(2);

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

    // ---- 次のイテレーション用にレンジを自動更新（ロック対応） ----
    Double_t new_left_calc  = p1 - 1.8 * p2;
    Double_t new_right_calc = p1 + 2.2 * p2;

    // ヒストの範囲をはみ出さないようにクリップ
    Double_t axis_min = h->GetXaxis()->GetXmin();
    Double_t axis_max = h->GetXaxis()->GetXmax();
    if (new_left_calc  < axis_min) new_left_calc  = axis_min;
    if (new_right_calc > axis_max) new_right_calc = axis_max;

    // ロック状態を考慮して反映
    Double_t new_left  = g.lock_left  ? g.range_left  : new_left_calc;
    Double_t new_right = g.lock_right ? g.range_right : new_right_calc;

    g.range_left  = new_left;
    g.range_right = new_right;

    // 直近の p1, p2 を保存（fit_iter の収束判定用）
    g.last_p1 = p1;
    g.last_p2 = p2;
    g.has_last_fit = kTRUE;

    // ログ＆コピペ用出力
    std::cout << Form("{ %.1f, %.1f, %.1f, %d }",
                      new_left, new_right, p1, m=="landau") << std::endl;
    std::cout << Form("set_range(%.1f, %.1f); fit(\"%s\")",
                      new_left, new_right, m.Data()) << std::endl;

    std::cout << "========================================" << std::endl;
}

//--------------------------------------------------
// イテレーション版 fit
//   n_iter: 最大イテレーション回数
//   model : "gaus", "landau", "auto"
//   logy  : yログ表示
//   内部で:
//     - 毎回 fit()
