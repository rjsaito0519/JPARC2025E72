//_____________________________________________________________________________
/**
 * DstTPCBcOutTracking の出力 ROOT から、
 * TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock (補正前) と
 * TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d (補正後) を描画し、PDF に出力する。
 *
 * オプションで TpcPhase.root も読み込み、RawClock ヒストグラムに
 * フィット曲線を重ね描きする。
 *
 * 2D ヒストが重いため、各ページを PNG 等のラスター形式で出力し、
 * ImageMagick の convert で PDF に埋め込む（軽量 PDF）。
 * convert がなければ PNG のみ出力し、手動で convert *.png out.pdf 可能。
 *
 * Usage:
 *   tpc_phase_plot <bcout.root> [output.pdf] [--phase <TpcPhase.root>] [--vdrift V]
 *
 * オプション:
 *   --phase <file>  : TpcPhase.root を指定してフィット曲線を重ね描き
 *   --vdrift V      : drift velocity [mm/ns] (デフォルト 0.055)
 *   --mode hit|trk  : ヒスト名の優先順（trk: TPCTrk_ResY_vs_ClockTime_... を先に、デフォルト hit）
 *
 * 出力 PDF 構成:
 *   Page 1: CoBo RawClock (補正前、2x4) + フィット曲線（--phase 指定時）
 *   Page 2: CoBo RawClock close-up（段差近傍、fit 品質確認）
 *   Page 3: CoBo 補正後 (2x4) - 補正効果の確認用
 *   Page 4: Asad RawClock (補正前、4x8)
 *   Page 5: Asad 補正後 (4x8) - 補正効果の確認用
 */
//_____________________________________________________________________________

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TPRegexp.h>

#include "paths.h"
#include "ana_helper.h"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

//_____________________________________________________________________________
static const Int_t NumOfSegCOBO = 8;
static const Int_t NumOfAsadTPC = 31;

// PNG 解像度（ピクセル）。ラスター PDF 用なので控えめに。
static const Int_t kCoboW = 1400, kCoboH = 760;
static const Int_t kAsadW = 1900, kAsadH = 980;

//_____________________________________________________________________________
// 検索するヒストグラム名のフォーマット候補 (優先順)
// フォーマット文字列には channel (%d または %02d) と suffix ("_RawClock" か empty) が入る
static const std::vector<std::string> HIST_FMTS_COBO_HIT = {
  "TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d%s",
  "TPCHit_ResY_vs_ClockTime_CoBo%d%s",
  "TPCCl_ResY_vs_ClockTime_CoBo%d%s",
  "TPCTrk_ResY_vs_ClockTime_CoBo%d%s",
};
static const std::vector<std::string> HIST_FMTS_COBO_TRK = {
  "TPCTrk_ResY_vs_ClockTime_CoBo%d%s",
  "TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d%s",
  "TPCHit_ResY_vs_ClockTime_CoBo%d%s",
  "TPCCl_ResY_vs_ClockTime_CoBo%d%s",
};

static const std::vector<std::string> HIST_FMTS_ASAD_HIT = {
  "TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d%s",
  "TPCHit_ResY_vs_ClockTime_Asad%02d%s",
  "TPCCl_ResY_vs_ClockTime_Asad%02d%s",
  "TPCTrk_ResY_vs_ClockTime_Asad%02d%s",
};
static const std::vector<std::string> HIST_FMTS_ASAD_TRK = {
  "TPCTrk_ResY_vs_ClockTime_Asad%02d%s",
  "TPC_ResidualY_BcOut_vs_ClockTime_Asad%02d%s",
  "TPCHit_ResY_vs_ClockTime_Asad%02d%s",
  "TPCCl_ResY_vs_ClockTime_Asad%02d%s",
};

static const std::vector<std::string>& HistFmtsCobo(const TString& mode)
{
  if (mode == "trk" || mode == "track" || mode == "TPCTrk")
    return HIST_FMTS_COBO_TRK;
  return HIST_FMTS_COBO_HIT;
}
static const std::vector<std::string>& HistFmtsAsad(const TString& mode)
{
  if (mode == "trk" || mode == "track" || mode == "TPCTrk")
    return HIST_FMTS_ASAD_TRK;
  return HIST_FMTS_ASAD_HIT;
}
//_____________________________________________________________________________
static TH2D* GetHist(TFile* f, const TString& hname)
{
  const TString paths[] = { hname, TString("hist/") + hname };
  for (const TString& path : paths) {
    TObject* obj = f->Get(path.Data());
    if (obj) {
      TH2D* h = dynamic_cast<TH2D*>(obj);
      if (h) return h;
    }
  }
  return nullptr;
}

//_____________________________________________________________________________
// TGraphErrors (clk, dclk[ns]) -> TGraph (clk, ResY[mm]); max_points<=0 は全点
static TGraph* MakeResyGraphFromDclock(TGraphErrors* gsrc, Double_t vdrift, Int_t max_points)
{
  if (!gsrc) return nullptr;
  const Int_t n = gsrc->GetN();
  if (n <= 0) return nullptr;
  Int_t step = 1;
  if (max_points > 0 && n > max_points)
    step = (n + max_points - 1) / max_points;
  Int_t nout = 0;
  for (Int_t i = 0; i < n; i += step) ++nout;
  TGraph* g = new TGraph(nout);
  Int_t j = 0;
  for (Int_t i = 0; i < n; i += step) {
    Double_t clk = 0.0, dclk = 0.0;
    gsrc->GetPoint(i, clk, dclk);
    g->SetPoint(j++, clk, -dclk * vdrift);
  }
  return g;
}

//_____________________________________________________________________________
// 式付き TF1 は内部が p0,p1,... 名のため GetParameters() 配列で読む
static Double_t Tf1Par(const TF1* f, Int_t i, Double_t def = 0.0)
{
  if (!f || i < 0 || i >= f->GetNpar()) return def;
  const Double_t* par = f->GetParameters();
  return par ? par[i] : def;
}

//_____________________________________________________________________________
// TF1 からステップ位置 t [ns] と幅 w [ns] を取得（par 名優先、なければレイアウト推定）
static void StepTwFromTf1(const TF1* f, Double_t& t_step, Double_t& w_step)
{
  t_step = 0.0;
  w_step = 0.01;
  if (!f) return;

  Bool_t have_t = kFALSE;
  Bool_t have_w = kFALSE;
  const Int_t np = f->GetNpar();
  for (Int_t i = 0; i < np; ++i) {
    TString nm = f->GetParName(i);
    nm.ToLower();
    if (nm == "t" || nm.Contains("t0")) {
      t_step = Tf1Par(f, i);
      have_t = kTRUE;
    }
    if (nm == "w" || nm.Contains("width")) {
      w_step = std::max(Tf1Par(f, i), 1e-6);
      have_w = kTRUE;
    }
  }

  if (np >= 4) {
    if (!have_t) t_step = Tf1Par(f, 2);
    if (!have_w) w_step = std::max(Tf1Par(f, 3), 1e-6);
  } else if (np == 3) {
    if (!have_t) t_step = Tf1Par(f, 1);
    if (!have_w) w_step = std::max(Tf1Par(f, 2), 1e-6);
  }
}

//_____________________________________________________________________________
// ROOT から読み戻した TF1 は C++ 関数ポインタが無効なのでパラメータから評価する。
static Bool_t IsBaselinePhaseTf1(const TF1* f)
{
  if (!f || f->GetNpar() != 4) return kFALSE;
  TString n = f->GetName();
  if (n.Contains("FitRaw")) return kTRUE;
  TString nm0 = f->GetParName(0);
  nm0.ToLower();
  if (nm0 == "c0" || nm0.Contains("base")) return kTRUE;
  // StepSum(N=1): p[0]≈1 は baseline ではない
  if (std::fabs(Tf1Par(f, 0) - 1.0) < 0.05) return kFALSE;
  return kFALSE;
}

static Bool_t PhaseTf1Usable(const TF1* f)
{
  if (!f) return kFALSE;
  const Int_t np = f->GetNpar();
  if (np == 3 || np == 4) return std::isfinite(Tf1Par(f, 0));
  const Int_t n = (Int_t)(Tf1Par(f, 0) + 0.5);
  return (n > 0 && np == 1 + 3 * n);
}

static Double_t EvalPhaseDclk(Double_t x, const TF1* f)
{
  if (!f) return 0.0;
  const Int_t np = f->GetNpar();

  if (np == 4 && IsBaselinePhaseTf1(f)) {
    const Double_t w = std::max(Tf1Par(f, 3), 1e-6);
    return Tf1Par(f, 0)
         + Tf1Par(f, 1) * TMath::Freq((x - Tf1Par(f, 2)) / w);
  }
  if (np == 3) {
    const Double_t w = std::max(Tf1Par(f, 2), 1e-6);
    return Tf1Par(f, 0) * TMath::Freq((x - Tf1Par(f, 1)) / w);
  }
  if (np >= 4) {
    const Int_t n = (Int_t)(Tf1Par(f, 0) + 0.5);
    if (n > 0 && np >= 1 + 3 * n) {
      Double_t s = 0.0;
      for (Int_t i = 0; i < n; ++i) {
        const Double_t A = Tf1Par(f, 1 + 3 * i);
        const Double_t t = Tf1Par(f, 2 + 3 * i);
        const Double_t w = std::max(Tf1Par(f, 3 + 3 * i), 1e-6);
        s += A * TMath::Freq((x - t) / w);
      }
      return s;
    }
  }
  return 0.0;
}

static Double_t EvalPhaseResy(Double_t x, const TF1* f, Double_t vdrift)
{
  return -EvalPhaseDclk(x, f) * vdrift;
}

//_____________________________________________________________________________
// ProfileFitUsed の X 範囲（実際の fit 窓）
static Bool_t GraphXRange(const TGraph* g, Double_t& xlo, Double_t& xhi)
{
  if (!g || g->GetN() < 1) return kFALSE;
  Double_t y0 = 0.0;
  g->GetPoint(0, xlo, y0);
  xhi = xlo;
  for (Int_t i = 1; i < g->GetN(); ++i) {
    Double_t x = 0.0, y = 0.0;
    g->GetPoint(i, x, y);
    xlo = std::min(xlo, x);
    xhi = std::max(xhi, x);
  }
  return (xhi > xlo);
}

//_____________________________________________________________________________
// DrawClone したオブジェクトをパッド最前面へ（pad が所有するので呼び出し側は g を delete 可）
static void DrawGraphCloneOnTop(TGraph* g, const char* opt)
{
  if (!g) return;
  TObject* cloned = g->DrawClone(opt);
  if (cloned && gPad) {
    gPad->GetListOfPrimitives()->Remove(cloned);
    gPad->GetListOfPrimitives()->Add(cloned);
  }
}

//_____________________________________________________________________________
// TGraph::Eval は x が定義域外だと 0 を返す → 端点でクランプしてパルス状の偽段差を防ぐ
static Double_t EvalPhaseGraphDclk(const TGraph* g, Double_t x)
{
  if (!g || g->GetN() < 1) return 0.0;
  Double_t x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;
  g->GetPoint(0, x0, y0);
  g->GetPoint(g->GetN() - 1, x1, y1);
  if (x <= x0) return y0;
  if (x >= x1) return y1;
  return g->Eval(x);
}

//_____________________________________________________________________________
// TpcPhase_Cobo TGraph (Δclock) を ResY 赤線として描画
static TGraph* MakeResyGraphFromPhase(TGraph* g_phase, Double_t vdrift,
                                      Double_t x1, Double_t x2, Int_t nseg = 300)
{
  if (!g_phase || g_phase->GetN() < 2 || vdrift <= 0.0 || !(x2 > x1))
    return nullptr;
  TGraph* g = new TGraph(nseg);
  for (Int_t i = 0; i < nseg; ++i) {
    const Double_t x = x1 + (x2 - x1) * (Double_t)i / (Double_t)(nseg - 1);
    g->SetPoint(i, x, -EvalPhaseGraphDclk(g_phase, x) * vdrift);
  }
  return g;
}

//_____________________________________________________________________________
// TF1 (Δclock[ns]) を ResY[mm] 曲線として描画（パラメータ評価、Eval 非使用）。
static void DrawPhaseTf1Resy(TF1* tf1, Double_t vdrift,
                             Double_t x1, Double_t x2,
                             Color_t color, Width_t width, Style_t style,
                             Int_t nseg = 300)
{
  if (!tf1 || vdrift <= 0.0 || !(x2 > x1))
    return;
  TGraph* g = new TGraph(nseg);
  for (Int_t i = 0; i < nseg; ++i) {
    const Double_t x = x1 + (x2 - x1) * (Double_t)i / (Double_t)(nseg - 1);
    g->SetPoint(i, x, EvalPhaseResy(x, tf1, vdrift));
  }
  g->SetLineColor(color);
  g->SetLineWidth(width);
  g->SetLineStyle(style);
  DrawGraphCloneOnTop(g, "L same");
  delete g;
}

// TpcPhase_CoboFallback を一度だけ読み込み（SetBranchAddress 競合回避）
struct CoboFallbackRow {
  Int_t cobo = -1;
  Int_t nstep_fit = 0;
  Double_t p0 = 0.0;
  Double_t p1 = 0.0;
  Double_t p2 = 1.0;
  Double_t t0_init = 0.0;
  Double_t amp_init = 0.0;
};

static std::vector<CoboFallbackRow> g_cobo_fb_rows;
static TFile* g_cobo_fb_src = nullptr;

static void LoadCoboFallbackRows(TFile* pf)
{
  if (!pf) return;
  if (pf == g_cobo_fb_src && !g_cobo_fb_rows.empty()) return;
  g_cobo_fb_src = pf;
  g_cobo_fb_rows.clear();

  TTree* tree = (TTree*)pf->Get("TpcPhase_CoboFallback");
  if (!tree) return;

  Int_t cb = 0, nstep = 0;
  Double_t p0 = 0.0, p1 = 0.0, p2 = 1.0;
  Double_t t0i = 0.0, ampi = 0.0;
  tree->ResetBranchAddresses();
  tree->SetBranchAddress("cobo", &cb);
  tree->SetBranchAddress("nstep_fit", &nstep);
  tree->SetBranchAddress("p0", &p0);
  tree->SetBranchAddress("p1", &p1);
  tree->SetBranchAddress("p2", &p2);
  const Bool_t has_t0 = (tree->GetBranch("t0_init") != nullptr);
  const Bool_t has_amp = (tree->GetBranch("amp_init") != nullptr);
  if (has_t0) tree->SetBranchAddress("t0_init", &t0i);
  if (has_amp) tree->SetBranchAddress("amp_init", &ampi);

  const Long64_t nent = tree->GetEntries();
  g_cobo_fb_rows.reserve((size_t)nent);
  for (Long64_t i = 0; i < nent; ++i) {
    tree->GetEntry(i);
    CoboFallbackRow row;
    row.cobo = cb;
    row.nstep_fit = nstep;
    row.p0 = p0;
    row.p1 = p1;
    row.p2 = p2;
    row.t0_init = has_t0 ? t0i : 0.0;
    row.amp_init = has_amp ? ampi : 0.0;
    g_cobo_fb_rows.push_back(row);
  }
}

static const CoboFallbackRow* FindCoboFallbackRow(Int_t cobo)
{
  for (const auto& row : g_cobo_fb_rows) {
    if (row.cobo == cobo) return &row;
  }
  return nullptr;
}

//_____________________________________________________________________________
// 段差位置 t・幅 w・振幅 A を phase ファイルから取得（TF1 無しでも可）
static Bool_t StepTwFromPhaseFile(TFile* pf, Int_t cobo,
                                  Double_t& t_step, Double_t& w_step,
                                  Double_t& amp_dclk)
{
  t_step = 0.0;
  w_step = 0.005;
  amp_dclk = 0.0;
  if (!pf) return kFALSE;

  LoadCoboFallbackRows(pf);
  const CoboFallbackRow* row = FindCoboFallbackRow(cobo);
  if (row && row->nstep_fit == 1) {
    amp_dclk = row->p0;
    t_step = row->p1;
    w_step = std::max(row->p2, 1e-6);
    return kTRUE;
  }

  TF1* f_fold = (TF1*)pf->Get(Form("TpcPhase_Fit_Cobo%d", cobo));
  if (f_fold && PhaseTf1Usable(f_fold)) {
    StepTwFromTf1(f_fold, t_step, w_step);
    amp_dclk = Tf1Par(f_fold, 0);
    return kTRUE;
  }
  TF1* f_raw = (TF1*)pf->Get(Form("TpcPhase_FitRaw_Cobo%d", cobo));
  if (f_raw && PhaseTf1Usable(f_raw)) {
    StepTwFromTf1(f_raw, t_step, w_step);
    amp_dclk = IsBaselinePhaseTf1(f_raw) ? Tf1Par(f_raw, 1) : Tf1Par(f_raw, 0);
    return kTRUE;
  }

  TGraph* g = (TGraph*)pf->Get(Form("TpcPhase_Cobo%d", cobo));
  if (!g || g->GetN() < 4) return kFALSE;
  Double_t best_slope = 0.0;
  for (Int_t i = 1; i < g->GetN(); ++i) {
    Double_t x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;
    g->GetPoint(i - 1, x0, y0);
    g->GetPoint(i, x1, y1);
    const Double_t dx = x1 - x0;
    if (dx <= 0.0) continue;
    const Double_t slope = std::fabs((y1 - y0) / dx);
    if (slope > best_slope) {
      best_slope = slope;
      t_step = 0.5 * (x0 + x1);
      amp_dclk = y1 - y0;
    }
  }
  return (best_slope > 0.0);
}

// params.h 由来の t0 初期値（TpcPhase_CoboFallback.t0_init）
static Bool_t T0InitFromPhaseFile(TFile* pf, Int_t cobo, Double_t& t0_init)
{
  t0_init = 0.0;
  if (!pf) return kFALSE;
  LoadCoboFallbackRows(pf);
  const CoboFallbackRow* row = FindCoboFallbackRow(cobo);
  if (!row) return kFALSE;
  t0_init = row->t0_init;
  return (std::fabs(t0_init) > 1e-6);
}

static void DrawT0InitVLine(Double_t t0_init, Double_t y1, Double_t y2)
{
  if (std::fabs(t0_init) < 1e-6) return;
  TLine* l = new TLine(t0_init, y1, t0_init, y2);
  l->SetLineColor(kMagenta + 1);
  l->SetLineStyle(7);
  l->SetLineWidth(2);
  l->Draw("same");
}

// fit 窓 X 範囲: ProfileFitUsed 優先、無ければ t0 ± 20 ns
static const Double_t kPhaseStepFitHalfWindowNs = 20.0;

static Bool_t FitWindowXRange(TFile* pf, Int_t cobo, Double_t& xlo, Double_t& xhi)
{
  if (!pf) return kFALSE;
  TGraphErrors* g_fitused =
      (TGraphErrors*)pf->Get(Form("TpcPhase_ProfileFitUsed_Cobo%d", cobo));
  if (GraphXRange(g_fitused, xlo, xhi)) return kTRUE;
  Double_t t0 = 0.0, w = 0.005, amp = 0.0;
  if (!StepTwFromPhaseFile(pf, cobo, t0, w, amp)) return kFALSE;
  xlo = t0 - kPhaseStepFitHalfWindowNs;
  xhi = t0 + kPhaseStepFitHalfWindowNs;
  return (xhi > xlo);
}

// fit に使った X 範囲の縦線
static void DrawFitRangeVLines(Double_t fxmin, Double_t fxmax,
                               Double_t y1, Double_t y2)
{
  if (!(fxmax > fxmin)) return;
  TLine* l_xmin = new TLine(fxmin, y1, fxmin, y2);
  TLine* l_xmax = new TLine(fxmax, y1, fxmax, y2);
  l_xmin->SetLineColor(kBlue);
  l_xmax->SetLineColor(kBlue);
  l_xmin->SetLineStyle(3);
  l_xmax->SetLineStyle(3);
  l_xmin->SetLineWidth(2);
  l_xmax->SetLineWidth(2);
  l_xmin->Draw("same");
  l_xmax->Draw("same");
}

//_____________________________________________________________________________
static void SetupStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleSize(0.06, "t");
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(50);
}

//_____________________________________________________________________________
static void DrawCoboPage(TFile* f, TCanvas* c, Bool_t corrected,
                         const std::vector<std::string>& hist_fmts_cobo,
                         TFile* phase_file = nullptr, Double_t vdrift = 0.055)
{
  c->Clear();
  c->Divide(4, 2, 0.001, 0.001);

  const char* suffix = corrected ? "_Corrected" : "";
  const char* label = corrected ? "Corrected" : "Raw";

  for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
    // 候補フォーマットを順に試す
    TH2D* h = nullptr;
    TString found_name = "";
    
    // suffix: RawClock or empty (corrected)
    const char* hist_suffix = !corrected ? "_RawClock" : "";

    for (const auto& fmt : hist_fmts_cobo) {
      TString hname = Form(fmt.c_str(), cobo, hist_suffix);
      h = GetHist(f, hname);
      if (h) {
        found_name = hname;
        break;
      }
    }

    c->cd(cobo + 1);
    gPad->SetRightMargin(0.18);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);
    if (h) {
      h->SetTitle(Form("CoBo %d (%s);Clock Time [ns];Residual Y [mm]", cobo, label));
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->GetXaxis()->SetLabelSize(0.045);
      h->GetYaxis()->SetLabelSize(0.045);
      h->Draw("colz");
      
      // フィット曲線と補助線を重ね描き（RawClockかつphase_file指定時のみ）
      if (!corrected && phase_file) {
        TGraph* g_phase = (TGraph*)phase_file->Get(Form("TpcPhase_Cobo%d", cobo));
        TF1* f_raw  = (TF1*)phase_file->Get(Form("TpcPhase_FitRaw_Cobo%d", cobo));
        TF1* f_fold = (TF1*)phase_file->Get(Form("TpcPhase_Fit_Cobo%d", cobo));
        TF1* f_fitplot = nullptr;
        if (f_raw && PhaseTf1Usable(f_raw))
          f_fitplot = f_raw;
        else if (f_fold && PhaseTf1Usable(f_fold))
          f_fitplot = f_fold;

        const Double_t x_plot_lo = h->GetXaxis()->GetXmin();
        const Double_t x_plot_hi = h->GetXaxis()->GetXmax();

        // 全 profile 点をオレンジ丸で描画（fit 窓外も含む）
        TGraphErrors* g_profile = (TGraphErrors*)phase_file->Get(Form("TpcPhase_Profile_Cobo%d", cobo));
        if (g_profile) {
          TGraph* g_pts = MakeResyGraphFromDclock(g_profile, vdrift, 0);
          if (g_pts) {
            g_pts->SetMarkerStyle(24);
            g_pts->SetMarkerSize(0.75);
            g_pts->SetMarkerColor(kOrange + 7);
            g_pts->SetLineColor(kOrange + 7);
            g_pts->DrawClone("P same");
            delete g_pts;
          }
        }

        Double_t t0_fit = 0.0, w_fit = 0.005, amp_dummy = 0.0;
        const Bool_t have_t0 = StepTwFromPhaseFile(phase_file, cobo, t0_fit, w_fit, amp_dummy);
        const Double_t y_min = h->GetYaxis()->GetXmin();
        const Double_t y_max = h->GetYaxis()->GetXmax();
        Double_t t0_init = 0.0;
        if (T0InitFromPhaseFile(phase_file, cobo, t0_init))
          DrawT0InitVLine(t0_init, y_min, y_max);
        if (have_t0) {
          TLine* l_t0 = new TLine(t0_fit, y_min, t0_fit, y_max);
          l_t0->SetLineColor(kGreen + 2);
          l_t0->SetLineStyle(2);
          l_t0->Draw("same");
        }

        // 赤 fit 線: ヒスト表示範囲だけ（TF1 定義域 ±60 全体は描かない）
        if (g_phase && g_phase->GetN() > 1) {
          TGraph* g_red = MakeResyGraphFromPhase(g_phase, vdrift, x_plot_lo, x_plot_hi, 400);
          if (g_red) {
            g_red->SetLineColor(kRed);
            g_red->SetLineWidth(3);
            g_red->SetLineStyle(1);
            DrawGraphCloneOnTop(g_red, "L same");
            delete g_red;
          }
        } else if (f_fitplot) {
          DrawPhaseTf1Resy(f_fitplot, vdrift, x_plot_lo, x_plot_hi, kRed, 3, 1);
        }

        // fit 窓の縦線（赤線の上に描いて見えるようにする）
        Double_t fxmin = 0.0, fxmax = 0.0;
        if (FitWindowXRange(phase_file, cobo, fxmin, fxmax))
          DrawFitRangeVLines(fxmin, fxmax, y_min, y_max);
      }
    } else {
      TH2D* dummy = new TH2D(Form("dummy_cobo%d_%s", cobo, label),
                             Form("CoBo %d (%s) - not found", cobo, label),
                             1, -60, 50, 1, -20, 20);
      dummy->Draw();
    }
  }
  c->Update();
}

//_____________________________________________________________________________
static void DrawCoboCloseupPage(TFile* f, TCanvas* c,
                                const std::vector<std::string>& hist_fmts_cobo,
                                TFile* phase_file, Double_t vdrift = 0.055,
                                Double_t half_window_ns = 6.0)
{
  c->Clear();
  c->Divide(4, 2, 0.001, 0.001);

  for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
    TH2D* h = nullptr;
    const char* hist_suffix = "_RawClock";
    for (const auto& fmt : hist_fmts_cobo) {
      TString hname = Form(fmt.c_str(), cobo, hist_suffix);
      h = GetHist(f, hname);
      if (h) break;
    }
    c->cd(cobo + 1);
    gPad->SetRightMargin(0.18);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);

    if (!h || !phase_file) {
      TH2D* dummy = new TH2D(Form("dummy_closeup_cobo%d", cobo),
                             Form("CoBo %d close-up - not found", cobo),
                             1, -10, 10, 1, -20, 20);
      dummy->Draw();
      continue;
    }

    TGraph* g_phase = (TGraph*)phase_file->Get(Form("TpcPhase_Cobo%d", cobo));
    TGraphErrors* g_fitused =
        (TGraphErrors*)phase_file->Get(Form("TpcPhase_ProfileFitUsed_Cobo%d", cobo));

    Double_t t0 = 0.0, w_step = 0.005, amp_dclk = 0.0;
    StepTwFromPhaseFile(phase_file, cobo, t0, w_step, amp_dclk);

    const Double_t half_win =
        std::max(half_window_ns, std::max(6.0, 80.0 * w_step));
    const Double_t hxmin = h->GetXaxis()->GetXmin();
    const Double_t hxmax = h->GetXaxis()->GetXmax();
    Double_t x1 = t0 - half_win;
    Double_t x2 = t0 + half_win;
    x1 = std::max(x1, hxmin);
    x2 = std::min(x2, hxmax);
    if (x2 <= x1) {
      x1 = hxmin;
      x2 = hxmax;
    }

    Double_t y_lo = h->GetYaxis()->GetXmin();
    Double_t y_hi = h->GetYaxis()->GetXmax();
    std::vector<Double_t> y_probe;
    y_probe.reserve(512);
    if (g_fitused) {
      for (Int_t i = 0; i < g_fitused->GetN(); ++i) {
        Double_t x = 0.0, dclk = 0.0;
        g_fitused->GetPoint(i, x, dclk);
        if (x < x1 || x > x2) continue;
        y_probe.push_back(-dclk * vdrift);
      }
    }
    if (g_phase && g_phase->GetN() > 1) {
      for (Int_t i = 0; i < 200; ++i) {
        const Double_t x = x1 + (x2 - x1) * (Double_t)i / 199.0;
        y_probe.push_back(-EvalPhaseGraphDclk(g_phase, x) * vdrift);
      }
    }
    const Double_t amp_resy = std::fabs(amp_dclk) * vdrift;
    if (amp_resy > 0.2 && g_phase && g_phase->GetN() > 1) {
      const Double_t y_mid = -EvalPhaseGraphDclk(g_phase, t0) * vdrift;
      y_lo = y_mid - 0.65 * amp_resy;
      y_hi = y_mid + 0.65 * amp_resy;
    }
    if (!y_probe.empty()) {
      const auto mm = std::minmax_element(y_probe.begin(), y_probe.end());
      const Double_t span = std::max(0.8, *mm.second - *mm.first);
      y_lo = std::min(y_lo, *mm.first - 0.20 * span);
      y_hi = std::max(y_hi, *mm.second + 0.20 * span);
    }

    TH2D* h_zoom = (TH2D*)h->Clone(Form("h_closeup_cobo%d", cobo));
    h_zoom->SetDirectory(nullptr);
    h_zoom->SetTitle(Form("CoBo %d close-up;Clock Time [ns];Residual Y [mm]", cobo));
    h_zoom->GetXaxis()->SetRangeUser(x1, x2);
    h_zoom->GetYaxis()->SetRangeUser(y_lo, y_hi);
    h_zoom->DrawClone("colz");

    TGraphErrors* g_profile = (TGraphErrors*)phase_file->Get(Form("TpcPhase_Profile_Cobo%d", cobo));
    if (g_profile) {
      TGraph* g_pts = new TGraph();
      for (Int_t i = 0; i < g_profile->GetN(); ++i) {
        Double_t x = 0.0, dclk = 0.0;
        g_profile->GetPoint(i, x, dclk);
        if (x < x1 || x > x2) continue;
        const Int_t j = g_pts->GetN();
        g_pts->SetPoint(j, x, -dclk * vdrift);
      }
      if (g_pts->GetN() > 0) {
        g_pts->SetMarkerStyle(24);
        g_pts->SetMarkerSize(1.0);
        g_pts->SetMarkerColor(kOrange + 7);
        g_pts->SetLineColor(kOrange + 7);
        g_pts->DrawClone("P same");
      }
      delete g_pts;
    }

    Double_t t0_init = 0.0;
    if (T0InitFromPhaseFile(phase_file, cobo, t0_init))
      DrawT0InitVLine(t0_init, y_lo, y_hi);

    TLine* lt0 = new TLine(t0, y_lo, t0, y_hi);
    lt0->SetLineColor(kGreen + 2);
    lt0->SetLineStyle(2);
    lt0->Draw("same");

    // 赤 fit 線
    if (g_phase && g_phase->GetN() > 1) {
      TGraph* g_red = MakeResyGraphFromPhase(g_phase, vdrift, x1, x2, 300);
      if (g_red) {
        g_red->SetLineColor(kRed);
        g_red->SetLineWidth(3);
        g_red->SetLineStyle(1);
        DrawGraphCloneOnTop(g_red, "L same");
        delete g_red;
      }
    } else {
      TF1* f_fold = (TF1*)phase_file->Get(Form("TpcPhase_Fit_Cobo%d", cobo));
      if (f_fold && PhaseTf1Usable(f_fold))
        DrawPhaseTf1Resy(f_fold, vdrift, x1, x2, kRed, 3, 1);
    }

    // fit 窓の縦線
    Double_t fxmin = 0.0, fxmax = 0.0;
    if (FitWindowXRange(phase_file, cobo, fxmin, fxmax))
      DrawFitRangeVLines(fxmin, fxmax, y_lo, y_hi);

    delete h_zoom;
  }
  c->Update();
}

//_____________________________________________________________________________
static void DrawAsadPage(TFile* f, TCanvas* c, Bool_t corrected,
                         const std::vector<std::string>& hist_fmts_asad)
{
  c->Clear();
  c->Divide(8, 4, 0.001, 0.001);

  const char* suffix = corrected ? "_Corrected" : "";
  const char* label = corrected ? "Corrected" : "Raw";

  for (Int_t asad = 0; asad < NumOfAsadTPC; ++asad) {
    // 候補フォーマットを順に試す
    TH2D* h = nullptr;
    TString found_name = "";

    // suffix: RawClock or empty (corrected)
    const char* hist_suffix = !corrected ? "_RawClock" : "";

    for (const auto& fmt : hist_fmts_asad) {
      TString hname = Form(fmt.c_str(), asad, hist_suffix);
      h = GetHist(f, hname);
      if (h) {
        found_name = hname;
        break;
      }
    }
    
    c->cd(asad + 1);
    gPad->SetRightMargin(0.18);
    gPad->SetLeftMargin(0.11);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.10);
    if (h) {
      h->SetTitle(Form("Asad %d (%s)", asad, label));
      h->GetXaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelSize(0.05);
      h->Draw("colz");
    } else {
      TH2D* dummy = new TH2D(Form("dummy_asad%d_%s", asad, label),
                             Form("Asad %d - not found", asad),
                             1, -60, 50, 1, -20, 20);
      dummy->Draw();
    }
  }
  c->Update();
}

//_____________________________________________________________________________
static int run_convert(const std::vector<std::string>& pngs, const std::string& pdf)
{
  std::ostringstream cmd;
  cmd << "convert";
  for (const auto& p : pngs) cmd << " \"" << p << "\"";
  cmd << " \"" << pdf << "\"";
  return std::system(cmd.str().c_str());
}

//_____________________________________________________________________________
int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <bcout.root> [output.pdf] [--phase <TpcPhase.root>] [--vdrift V] [--mode hit|trk]\n"
              << "  Draws 2D hists as PNG, embeds into PDF via 'convert' (ImageMagick).\n"
              << "  Options:\n"
              << "    --phase <file> : overlay fit curves from TpcPhase.root\n"
              << "    --vdrift V     : drift velocity [mm/ns] (default: 0.055)\n"
              << "    --mode hit|trk : histogram name priority (default: hit)\n"
              << "  If convert is missing, PNGs are kept; merge with:\n"
              << "    convert <base>_page*.png <output.pdf>\n";
    return 1;
  }

  std::string inpath = argv[1];
  std::string outpath;
  std::string phase_path;
  Double_t vdrift = 0.055;
  TString mode = "hit";

  // パース引数
  Int_t argidx = 2;
  while (argidx < argc) {
    std::string arg = argv[argidx];
    if (arg == "--phase" && argidx + 1 < argc) {
      phase_path = argv[++argidx];
      ++argidx;
    } else if (arg == "--mode" && argidx + 1 < argc) {
      mode = argv[++argidx];
      ++argidx;
    } else if (arg == "--vdrift" && argidx + 1 < argc) {
      vdrift = std::atof(argv[++argidx]);
      if (vdrift <= 0.0) vdrift = 0.055;
      ++argidx;
    } else if (arg.find("--") == 0) {
      ++argidx;
    } else {
      outpath = arg;
      ++argidx;
      break;
    }
  }

  if (outpath.empty()) {
    Int_t run_number = -1;
    TFile* f_tmp = TFile::Open(inpath.c_str());
    if (f_tmp && !f_tmp->IsZombie()) {
      TTree *tree = (TTree*)f_tmp->Get("tpc");
      if (tree) {
        try {
          TTreeReader reader(tree);
          TTreeReaderValue<UInt_t> rv(reader, "run_number");
          if (reader.Next()) {
            run_number = *rv;
            std::cout << "Auto-detected run number from TTree (TTreeReader): " << run_number << std::endl;
          }
        } catch (...) {
          // Ignore errors and try fallback
        }
      }
      f_tmp->Close();
      delete f_tmp;
    }

    if (run_number < 0) {
      TString filename = inpath.c_str();
      TPRegexp run_regex("run([0-9]+)");
      TString run_str = filename(run_regex);
      if (!run_str.IsNull()) {
        run_str.ReplaceAll("run", "");
        run_number = run_str.Atoi();
        std::cout << "Auto-detected run number from filename: " << run_number << std::endl;
      }
    }

    fs::path p(inpath);
    std::string base = p.stem().string();
    if (run_number >= 0) {
      TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
      outpath = std::string(img_dir.Data()) + "/" + base + "_TPC_Phase.pdf";
    } else {
      outpath = (p.parent_path() / (base + "_TPC_Phase.pdf")).string();
    }
  }

  fs::path out(outpath);
  fs::path outdir = out.parent_path();
  std::string outbase = outdir.empty()
    ? out.stem().string()
    : (outdir / out.stem()).string();

  TFile* f = TFile::Open(inpath.c_str());
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << inpath << std::endl;
    return 1;
  }
  
  // Phase ファイルを開く（オプション）
  TFile* f_phase = nullptr;
  if (!phase_path.empty()) {
    f_phase = TFile::Open(phase_path.c_str());
    if (!f_phase || f_phase->IsZombie()) {
      std::cerr << "Warning: cannot open phase file " << phase_path << ", skipping fit overlay\n";
      f_phase = nullptr;
    } else {
      std::cout << "Phase file loaded: " << phase_path << " (vdrift=" << vdrift << " mm/ns)\n";
    }
  }

  const std::vector<std::string>& hist_fmts_cobo = HistFmtsCobo(mode);
  const std::vector<std::string>& hist_fmts_asad = HistFmtsAsad(mode);
  std::cout << "Histogram mode: " << mode << std::endl;

  SetupStyle();
  gROOT->SetBatch(kTRUE);

  fs::path tmpdir = fs::temp_directory_path() / ("tpc_phase_plot_" + std::to_string(getpid()));
  std::error_code ec;
  if (!fs::create_directories(tmpdir, ec)) {
    std::cerr << "Error: cannot create temp dir " << tmpdir << std::endl;
    f->Close();
    delete f;
    return 1;
  }

  TCanvas* c_cobo = new TCanvas("c_cobo", "CoBo", kCoboW, kCoboH);
  TCanvas* c_asad = new TCanvas("c_asad", "Asad", kAsadW, kAsadH);

  std::vector<std::string> pngs;
  auto save = [&](const char* tag) {
    std::string png = (tmpdir / (std::string("page_") + tag + ".png")).string();
    pngs.push_back(png);
  };

  // Page 1: CoBo RawClock (補正前) + フィット曲線
  DrawCoboPage(f, c_cobo, false, hist_fmts_cobo, f_phase, vdrift);
  save("cobo_raw");
  c_cobo->Print(pngs.back().c_str());

  // Page 2: CoBo RawClock close-up（段差近傍の fit 品質確認）
  if (f_phase) {
    DrawCoboCloseupPage(f, c_cobo, hist_fmts_cobo, f_phase, vdrift, 6.0);
    save("cobo_raw_closeup");
    c_cobo->Print(pngs.back().c_str());
  }

  // Page 3: CoBo 補正後（補正効果の確認用）
  DrawCoboPage(f, c_cobo, true, hist_fmts_cobo, nullptr, vdrift);
  save("cobo_corrected");
  c_cobo->Print(pngs.back().c_str());

  // Page 4: Asad RawClock (補正前)
  DrawAsadPage(f, c_asad, false, hist_fmts_asad);
  save("asad_raw");
  c_asad->Print(pngs.back().c_str());

  // Page 5: Asad 補正後（補正効果の確認用）
  DrawAsadPage(f, c_asad, true, hist_fmts_asad);
  save("asad_corrected");
  c_asad->Print(pngs.back().c_str());

  if (f_phase) {
    f_phase->Close();
    delete f_phase;
    f_phase = nullptr;
  }
  f->Close();
  delete f;
  f = nullptr;
  delete c_cobo;
  c_cobo = nullptr;
  delete c_asad;
  c_asad = nullptr;

  int ret = run_convert(pngs, outpath);
  for (const auto& p : pngs)
    fs::remove(p, ec);

  if (ret != 0) {
    std::cerr << "Warning: 'convert' failed (ImageMagick?). Saving PNGs next to output.\n";
    if (!outdir.empty()) fs::create_directories(outdir, ec);
    TCanvas* c_cobo2 = new TCanvas("c_cobo", "CoBo", kCoboW, kCoboH);
    TCanvas* c_asad2 = new TCanvas("c_asad", "Asad", kAsadW, kAsadH);
    TFile* f2 = TFile::Open(inpath.c_str());
    TFile* f_phase2 = nullptr;
    if (!phase_path.empty()) {
      f_phase2 = TFile::Open(phase_path.c_str());
    }
    if (f2 && !f2->IsZombie()) {
      DrawCoboPage(f2, c_cobo2, false, hist_fmts_cobo, f_phase2, vdrift);
      c_cobo2->Print((outbase + "_page1_cobo_raw.png").c_str());
      if (f_phase2) {
        DrawCoboCloseupPage(f2, c_cobo2, hist_fmts_cobo, f_phase2, vdrift, 6.0);
        c_cobo2->Print((outbase + "_page2_cobo_raw_closeup.png").c_str());
      }
      DrawAsadPage(f2, c_asad2, false, hist_fmts_asad);
      c_asad2->Print((outbase + "_page3_asad_raw.png").c_str());
      f2->Close();
      delete f2;
    }
    if (f_phase2) {
      f_phase2->Close();
      delete f_phase2;
    }
    delete c_cobo2;
    delete c_asad2;
    std::cerr << "PNGs: " << outbase << "_page{1,2}_*.png\n"
              << "Merge: convert " << outbase << "_page*.png " << outpath << "\n";
    fs::remove_all(tmpdir, ec);
    std::exit(0);
  }

  fs::remove(tmpdir, ec);
  std::cout << "Saved: " << outpath << " (raster PDF via convert)\n";
  // ROOT 終了時のデストラクタで落ちることがあるため正常終了は即抜け
  std::exit(0);
}
