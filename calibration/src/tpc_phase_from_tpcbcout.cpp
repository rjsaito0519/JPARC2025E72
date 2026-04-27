//_____________________________________________________________________________
/**
 * DstTPCBcOutTracking で作った ROOT ファイルから、
 * TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d の Profile を取り、
 * Δclock(clk) = -mean_ResY(clk) / Vdrift を計算し、
 * TPCParamMan が読む補正用 ROOT ファイルを出力する。
 *
 * TPCParamMan は TpcPhase_Cobo%d という名前の TGraph を期待する。
 * (TF1 でフィット後、細かく Eval して TGraph に変換して保存)
 *
 * さらに TTree "TpcPhase_CoboFallback" を保存する:
 *   --fit-step 1 のとき各 CoBo の 1 段ステップ振幅・位置・幅を (p0,p1,p2) として記録し、
 *   TPCPRM の kCobo 行（aty=3）フォールバック式 PhaseShift と一致させるのに使う。
 *   nstep_fit: 0=flat 初期, 1=1 段フィット成功（p0,p1,p2 有効）, >1=多段（p は無効）, -1=フィット失敗
 *
 * Usage:
 *   tpc_phase_from_tpcbcout <tpcbcout.root> <TpcPhase.root> [--smooth N] [--vdrift V] [--fit-step [N]] [--graph-points N]
 *
 * 動作:
 *   - 生ヒスト (TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock) から Profile を取り、フィットして TpcPhase.root を出力
 *   - TF1 でフィット後、細かく評価して TGraph に変換
 *
 * オプション:
 *   --smooth N        : 移動平均の半窓幅（デフォルト 0 = スムージングなし）
 *   --vdrift V        : drift velocity [mm/ns]（デフォルト 0.055）
 *   --fit-step [N]     : ステップ関数 Freq の N 重ね合わせでフィット（N=0 で全ゼロの初期ファイル作成、省略時 N=1）
 *   --min-entries N   : Profile で使う最小エントリ数（デフォルト 5、統計が少ないビンをスキップ）
 *   --graph-points N  : TGraph の点数（デフォルト 10000）
 *   --mode hit|trk    : 入力 2D ヒストの優先順（trk: TPCTrk_ResY_vs_ClockTime_... を先に、デフォルト hit）
 *
 * 例:
 *   # 初期ファイル（全ゼロ、tpcbcout 不要）
 *   tpc_phase_from_tpcbcout param/TPCPHASE/TpcPhase_02601.root --fit-step 0
 *   # フィット
 *   tpc_phase_from_tpcbcout tpcbcout_run02601.root param/TPCPHASE/TpcPhase_02601.root --fit-step 2
 *
 * Ref: HistTools.cc TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock,
 *      DstTPCBcOutTracking.cc, TPCParamMan.cc TpcPhase_Cobo%d
 */
//_____________________________________________________________________________

#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TKey.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>

#include "params.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

//_____________________________________________________________________________
// HistTools.cc / DetectorID.hh と合わせる
static const Int_t NumOfSegCOBO = 8;

// StepSum フィットに使う clock 窓（t0 推定値の周り）[ns]
static const Double_t kPhaseStepFitHalfWindowNs = 40.0;
// 各 clock ビンの Y 射影に対する Gauss 窓の半幅の下限（±40 ns を ResY に換算）[mm]
static const Double_t kPhaseSliceGaussHalfWidthNsEquiv = 40.0;

// TPCParamMan が期待する TF1 名（補正用）
static const char* PHASE_NAME_FMT = "TpcPhase_Cobo%d";
// フィット関数の結果を保存する名前（参照用）
static const char* FIT_NAME_FMT = "TpcPhase_Fit_Cobo%d";

// 検索するヒストグラム名のフォーマット候補 (優先順)。--mode trk では TPCTrk を先に試す。
static const std::vector<std::string> HIST_FMTS_HIT = {
  "TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock",
  "TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock",
  "TPCCl_ResY_vs_ClockTime_CoBo%d_RawClock",
  "TPCTrk_ResY_vs_ClockTime_CoBo%d_RawClock",
};
static const std::vector<std::string> HIST_FMTS_TRK = {
  "TPCTrk_ResY_vs_ClockTime_CoBo%d_RawClock",
  "TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock",
  "TPCHit_ResY_vs_ClockTime_CoBo%d_RawClock",
  "TPCCl_ResY_vs_ClockTime_CoBo%d_RawClock",
};

// ResY [mm] -> Δclock [ns] の変換で使う drift velocity [mm/ns]
static const Double_t DEFAULT_VDRIFT = 0.055;

// デフォルトのスムージング半窓幅（0 = スムージングなし）
static const Int_t DEFAULT_SMOOTH = 0;

// TGraph の点数（TF1 を細かく Eval する）
static const Int_t DEFAULT_GRAPH_POINTS = 10000;

// ステップ関数の重ね合わせ: dclk(x) = sum_i A_i * Freq((x - t_i) / w_i)
// パラメータ: p[0]=N, p[1+3*i]=A_i, p[2+3*i]=t_i, p[3+3*i]=w_i
static Double_t StepSumFunc(Double_t* x, Double_t* p)
{
  const Double_t clk = x[0];
  const Int_t n = (Int_t)(p[0] + 0.5);
  if (n <= 0) return 0.0;
  
  Double_t s = 0.0;
  for (Int_t i = 0; i < n; ++i) {
    Double_t A = p[1 + 3 * i];
    Double_t t = p[2 + 3 * i];
    Double_t w = p[3 + 3 * i];
    if (w <= 0.0) w = 1e-6;
    s += A * TMath::Freq((clk - t) / w);
  }
  return s;
}

// ステップ初期値推定用
struct StepGuess {
  Double_t c_left;   // 左プラトー Δclock
  Double_t c_right;  // 右プラトー Δclock
  Double_t t0;       // ステップ位置 [ns]
};

// 前方宣言: ResY[mm] -> Δclock[ns] 変換
static Double_t ResYToDeltaClock(Double_t mean_resy_mm, Double_t vdrift_mm_ns);

// Y 射影 (ProjectionY) からステップの位置と高さを推定
static StepGuess
EstimateStepFromProjectionY(const TH2D* h, Double_t vdrift)
{
  StepGuess g{0.0, 0.0, 0.0};
  if (!h || vdrift <= 0.0) return g;

  // X 全体をいくつかのセグメントに分け、それぞれで Y 射影をとって
  // 「低いプラトー」「高いプラトー」の代表値と、その切り替わり位置を探す。
  const Double_t x_min = h->GetXaxis()->GetXmin();
  const Double_t x_max = h->GetXaxis()->GetXmax();
  const Double_t xrange = x_max - x_min;
  if (xrange <= 0.0) return g;

  // X を等間隔セグメントに分割して、それぞれで「ピーク付近のみ」を使って Y の代表値を取る
  const Int_t n_seg = 10;                 // X 分割数（あまり細かくしない）
  const Int_t min_entries_seg = 50;       // 1 セグメントあたり必要な最小エントリ数

  std::vector<Double_t> seg_x_center;
  std::vector<Double_t> seg_dclk;         // 各セグメントの Δclock 代表値

  for (Int_t i = 0; i < n_seg; ++i) {
    const Double_t x1 = x_min + xrange * (double)i / n_seg;
    const Double_t x2 = x_min + xrange * (double)(i + 1) / n_seg;

    const Int_t bin1 = h->GetXaxis()->FindFixBin(x1);
    const Int_t bin2 = h->GetXaxis()->FindFixBin(x2);
    if (bin2 < bin1) continue;

    TString proj_name = Form("tpc_phase_segY_%d", i);
    TH1D* py = h->ProjectionY(proj_name, bin1, bin2);
    if (!py) continue;

    const Double_t entries = py->GetEntries();
    if (entries < min_entries_seg) {
      delete py;
      continue;
    }

    // Y 軸ヒストから「最大ビンの 1/5 以上」のビンだけで重み付き平均をとり、
    // プラトーの中心に近い値を代表値とする（上下クラスタの中間を避ける狙い）
    const Int_t ny = py->GetNbinsX();
    Double_t max_bin_content = 0.0;
    for (Int_t iy = 1; iy <= ny; ++iy) {
      const Double_t w = py->GetBinContent(iy);
      if (w > max_bin_content) max_bin_content = w;
    }
    const Double_t threshold = max_bin_content / 5.0;

    Double_t sum = 0.0;
    Double_t sumw = 0.0;
    for (Int_t iy = 1; iy <= ny; ++iy) {
      const Double_t yc = py->GetXaxis()->GetBinCenter(iy);
      const Double_t w  = py->GetBinContent(iy);
      if (w >= threshold) {
        sum  += w * yc;
        sumw += w;
      }
    }

    if (sumw > 0.0) {
      const Double_t mean_y = sum / sumw;
      const Double_t dclk   = ResYToDeltaClock(mean_y, vdrift);
      const Double_t xc     = 0.5 * (x1 + x2);
      seg_x_center.push_back(xc);
      seg_dclk.push_back(dclk);
    }

    delete py;
  }

  if (seg_x_center.size() < 2) {
    // セグメントが十分に取れなかった場合は何もせず 0 のまま返す
    return g;
  }

  // Δclock のレンジが小さすぎる場合もステップがないとみなす
  Double_t dmin = *std::min_element(seg_dclk.begin(), seg_dclk.end());
  Double_t dmax = *std::max_element(seg_dclk.begin(), seg_dclk.end());
  if (dmax - dmin < 1.0) { // 1ns 未満の変化ならステップ無しとみなす
    return g;
  }

  // 2 プラトーへのざっくりな分類: 中間値で左右クラスタを分ける
  const Double_t dmid = 0.5 * (dmin + dmax);
  std::vector<int> label(seg_dclk.size(), 0);
  for (size_t i = 0; i < seg_dclk.size(); ++i) {
    label[i] = (seg_dclk[i] > dmid) ? 1 : 0;
  }

  // ラベルが切り替わる位置（候補が複数ある場合は Δ の差が最も大きいところ）を探す
  Double_t best_x_mid = 0.0;
  Double_t best_jump  = 0.0;
  for (size_t i = 0; i + 1 < seg_x_center.size(); ++i) {
    if (label[i] == label[i + 1]) continue;
    const Double_t jump = std::fabs(seg_dclk[i + 1] - seg_dclk[i]);
    if (jump > best_jump) {
      best_jump  = jump;
      best_x_mid = 0.5 * (seg_x_center[i] + seg_x_center[i + 1]);
    }
  }

  if (best_jump <= 0.0) {
    // ラベル変化が見つからなければ素直に中点を使う
    best_x_mid = 0.5 * (x_min + x_max);
  }

  // 左右プラトー値を、境界より左／右にあるセグメントから平均して求める
  Double_t sum_left  = 0.0;
  Double_t sum_right = 0.0;
  Int_t    n_left    = 0;
  Int_t    n_right   = 0;
  for (size_t i = 0; i < seg_x_center.size(); ++i) {
    if (seg_x_center[i] <= best_x_mid) {
      sum_left += seg_dclk[i];
      ++n_left;
    } else {
      sum_right += seg_dclk[i];
      ++n_right;
    }
  }

  if (n_left > 0 && n_right > 0) {
    g.c_left  = sum_left  / n_left;
    g.c_right = sum_right / n_right;
  } else {
    // どちらか一方しか取れなかった場合は全体平均で埋めておく
    const Double_t dmean = 0.5 * (dmin + dmax);
    g.c_left  = dmean;
    g.c_right = dmean;
  }

  // ステップ位置の初期値: 物理的レンジ ±35ns 内にクリップ
  g.t0 = std::clamp(best_x_mid, -35.0, 35.0);

  return g;
}

//_____________________________________________________________________________
// 1D Profile (delta_clock vs x_center) からステップ位置をロバストに推定
// 前後ウィンドウの平均差が最大となる位置をエッジとして検出する。
static Double_t EstimateStepPositionFromProfile(const std::vector<Double_t>& x_center,
                                                const std::vector<Double_t>& delta_clock)
{
  const Int_t n = (Int_t)x_center.size();
  if (n < 7 || delta_clock.size() != x_center.size()) return 0.0;

  const Int_t W = 3; // 前後3点ずつ
  Double_t best_score = 0.0;
  Int_t best_idx = -1;

  for (Int_t i = W; i < n - W; ++i) {
    Double_t sumL = 0.0, sumR = 0.0;
    for (Int_t j = i - W; j < i; ++j) {
      sumL += delta_clock[j];
    }
    for (Int_t j = i + 1; j <= i + W; ++j) {
      sumR += delta_clock[j];
    }
    const Double_t meanL = sumL / (Double_t)W;
    const Double_t meanR = sumR / (Double_t)W;
    const Double_t score = std::fabs(meanR - meanL);
    if (score > best_score) {
      best_score = score;
      best_idx = i;
    }
  }

  if (best_idx < 0) return 0.0;

  // さらに鋭い段差を捉えるために、best_idx 近傍で隣接ペアの差分を調べて
  // 一番大きく飛んでいる 2 点の「中点」をステップ位置とみなす。
  Int_t k_start = std::max(best_idx - 2, 0);
  Int_t k_end   = std::min(best_idx + 2, n - 2); // (k, k+1) なので n-2 まで

  Double_t best_jump = 0.0;
  Int_t best_pair_k = best_idx;
  for (Int_t k = k_start; k <= k_end; ++k) {
    const Double_t jump = std::fabs(delta_clock[k + 1] - delta_clock[k]);
    if (jump > best_jump) {
      best_jump = jump;
      best_pair_k = k;
    }
  }

  // 隣接 2 点の X 座標の中点を採用
  const Double_t x1 = x_center[best_pair_k];
  const Double_t x2 = x_center[best_pair_k + 1];
  return 0.5 * (x1 + x2);
}

//_____________________________________________________________________________
// N段ステップ関数でフィッティング
static TF1* FitStepSum(Int_t nstep, std::vector<Double_t>& x_center,
                       std::vector<Double_t>& delta_clock,
                       Double_t xmin, Double_t xmax,
                       const std::vector<Double_t>* err_delta = nullptr,
                       const StepGuess* guess = nullptr,
                       Double_t run_amp_ref_ns = 0.0,
                       Double_t amp_override_ns = 0.0,
                       Double_t step_width_ns = 0.005,
                       Bool_t free_step_width = kFALSE,
                       TF1** initial_out = nullptr)
{
  if (nstep <= 0 || x_center.empty() || delta_clock.size() != x_center.size())
    return nullptr;

  const Int_t n = (Int_t)x_center.size();
  
  // データのY範囲を取得
  Double_t ymin = *std::min_element(delta_clock.begin(), delta_clock.end());
  Double_t ymax = *std::max_element(delta_clock.begin(), delta_clock.end());
  if (ymax <= ymin) ymax = ymin + 1.0;

  // ============================================================
  // ステップ位置の検出: 勾配 × √(Y変化量) でスコアリング
  // ============================================================
  std::vector<Double_t> step_positions;
  std::vector<Double_t> step_amplitudes;
  
  if (n < nstep + 1) {
    // データ点数が不足 → 等間隔配置
    for (Int_t i = 0; i < nstep; ++i) {
      step_positions.push_back(xmin + (xmax - xmin) * (i + 1.0) / (nstep + 1.0));
      step_amplitudes.push_back((ymax - ymin) / nstep);
    }
    std::cout << "    Few data points, using evenly spaced steps\n";
  } else {
    // 隣接点間の勾配とY変化量を計算
    std::vector<std::pair<Double_t, Double_t>> score_x_pairs;  // (score, x)
    const Double_t x_range = xmax - xmin;
    const Double_t x_safe_min = xmin + x_range * 0.1;  // 端10%を除外
    const Double_t x_safe_max = xmax - x_range * 0.1;
    const Double_t y_threshold = (ymax - ymin) * 0.05;  // 5%未満の変化を除外
    
    for (Int_t i = 0; i < n - 1; ++i) {
      Double_t dx = x_center[i + 1] - x_center[i];
      Double_t dy = delta_clock[i + 1] - delta_clock[i];
      Double_t x_mid = (x_center[i] + x_center[i + 1]) / 2.0;
      
      // 端領域除外、小変化除外
      if (x_mid < x_safe_min || x_mid > x_safe_max) continue;
      if (std::abs(dy) < y_threshold) continue;
      
      Double_t gradient = (dx > 0) ? std::abs(dy / dx) : 0.0;
      Double_t score = gradient * std::sqrt(std::abs(dy));
      score_x_pairs.push_back({score, x_mid});
    }
    
    // スコア降順でソート
    std::sort(score_x_pairs.begin(), score_x_pairs.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    std::cout << "    Detecting steps (edges excluded):\n";
    
    // 上位スコアの位置をステップとする（近すぎる位置は除外）
    const Double_t min_separation = x_range / (nstep * 3.0);
    for (const auto& [score, x_pos] : score_x_pairs) {
      if ((Int_t)step_positions.size() >= nstep) break;
      
      Bool_t too_close = false;
      for (const auto& pos : step_positions) {
        if (std::abs(x_pos - pos) < min_separation) {
          too_close = true;
          break;
        }
      }
      
      if (!too_close) {
        step_positions.push_back(x_pos);
        
        // 振幅推定: 前後2nsの範囲で平均を取る
        Double_t y_before = 0.0, y_after = 0.0;
        Int_t n_before = 0, n_after = 0;
        for (Int_t j = 0; j < n; ++j) {
          if (x_center[j] < x_pos - 2.0) {
            y_before += delta_clock[j];
            n_before++;
          } else if (x_center[j] > x_pos + 2.0) {
            y_after += delta_clock[j];
            n_after++;
          }
        }
        
        Double_t amp = (n_before > 0 && n_after > 0) 
                       ? (y_after / n_after - y_before / n_before)
                       : (ymax - ymin) / nstep;
        step_amplitudes.push_back(amp);
        std::cout << "      Step at x=" << x_pos << " (score=" << score << ")\n";
      }
    }
    
    // 不足分は等間隔で補完
    while ((Int_t)step_positions.size() < nstep) {
      Double_t pos = xmin + x_range * (step_positions.size() + 1.0) / (nstep + 1.0);
      step_positions.push_back(pos);
      step_amplitudes.push_back((ymax - ymin) / nstep);
      std::cout << "      Step at x=" << pos << " (fallback)\n";
    }
    
    std::sort(step_positions.begin(), step_positions.end());
  }
  
  // ============================================================
  // フィット関数の設定
  // ============================================================
  const Int_t npar = 1 + 3 * nstep;  // N, A_i, t_i, w_i
  TF1* f = new TF1("step_sum", StepSumFunc, xmin, xmax, npar);
  
  // パラメータ0: ステップ数（固定）
  f->SetParName(0, "N");
  f->FixParameter(0, (Double_t)nstep);
  
  // パラメータ2以降: 各ステップの振幅、位置、幅
  const Double_t fixed_width = (step_width_ns > 0.0) ? step_width_ns : 0.005;  // 固定幅 [ns]
  // 位置調整範囲 [ns]。ステップ位置は物理的に |ClockTime|≲35ns の範囲内にある前提。
  const Double_t position_range = 35.0;
  const Double_t min_step_amp = 2.0;      // 振幅の絶対値の下限 [ns]
  const Double_t amp_rel_limit = 0.25;    // 振幅微調整の相対許容幅 (±25%)
  
  for (Int_t i = 0; i < nstep; ++i) {
    f->SetParName(1 + 3 * i, Form("A%d", i));
    f->SetParName(2 + 3 * i, Form("t%d", i));
    f->SetParName(3 + 3 * i, Form("w%d", i));
    
    // 振幅の初期値:
    //  - 1段目: run 全体の代表振幅 run_amp_ref_ns を中心に、Y 射影からの符号・微調整のみ許可
    //  - それ以外: 検出した step_amplitudes を使用
    Double_t A_init = step_amplitudes[i];
    if (guess && i == 0) {
      A_init = guess->c_right - guess->c_left;
    }
    Double_t sign = (A_init >= 0.0) ? 1.0 : -1.0;
    Double_t absA = std::fabs(A_init);
    if (absA < min_step_amp) absA = min_step_amp;

    if (i == 0 && std::fabs(amp_override_ns) > 0.0) {
      // run-cobo ごとの明示的な振幅指定（Δclock[ns]）を最優先
      Double_t ref_abs = std::max(std::fabs(amp_override_ns), min_step_amp);
      Double_t A_center = (amp_override_ns >= 0.0) ? ref_abs : -ref_abs;
      A_init = A_center;
      Double_t A_min = A_center * (1.0 - amp_rel_limit);
      Double_t A_max = A_center * (1.0 + amp_rel_limit);
      if (A_min > A_max) std::swap(A_min, A_max);
      f->SetParLimits(1 + 3 * i, A_min, A_max);
    } else if (i == 0 && run_amp_ref_ns > 0.0) {
      // run ごとの基本振幅を優先し、±amp_rel_limit の範囲だけ動けるように制限
      Double_t ref_abs = std::max(run_amp_ref_ns, min_step_amp);
      Double_t A_center = sign * ref_abs;
      A_init = A_center;
      Double_t A_min = A_center * (1.0 - amp_rel_limit);
      Double_t A_max = A_center * (1.0 + amp_rel_limit);
      if (A_min > A_max) std::swap(A_min, A_max);
      f->SetParLimits(1 + 3 * i, A_min, A_max);
    } else {
      // 通常の下限のみ
      A_init = sign * absA;
    }

    f->SetParameter(1 + 3 * i, A_init);   // 振幅（フィッティング）

    // 位置の初期値と制約: 候補位置をベースにしつつ、物理的には ±35ns に制限
    Double_t t0 = step_positions[i];
    if (guess && i == 0) {
      t0 = guess->t0;
    }
    if (t0 < -position_range) t0 = -position_range;
    if (t0 >  position_range) t0 =  position_range;
    f->SetParameter(2 + 3 * i, t0);    // 位置（初期値）
    f->SetParLimits(2 + 3 * i, -position_range, position_range);
    if (free_step_width) {
      f->SetParameter(3 + 3 * i, fixed_width);
      // 幅を可変にする場合は、過度に縛らず緩めのレンジで探索させる
      f->SetParLimits(3 + 3 * i, 0.0005, 0.50);
    } else {
      f->FixParameter(3 + 3 * i, fixed_width);        // 幅（固定）
    }
  }
  
  std::cout << "    Fit settings: width=" << fixed_width
            << " ns (" << (free_step_width ? "free" : "fixed")
            << "), position range=±" << position_range << " ns\n";

  // 初期パラメータのスナップショットを保存（デバッグ・可視化用）
  if (initial_out) {
    *initial_out = (TF1*)f->Clone("step_sum_initial");
  }

  // フィッティング実行
  Int_t fit_status = -1;
  if (err_delta && err_delta->size() == (size_t)n) {
    TGraphErrors gr(n, x_center.data(), delta_clock.data(), nullptr, err_delta->data());
    fit_status = gr.Fit(f, "QRN0MS");
  } else {
    TGraph gr(n, x_center.data(), delta_clock.data());
    fit_status = gr.Fit(f, "QRN0MS");
  }
  
  if (fit_status != 0) {
    std::cerr << "Warning: Fit failed (status=" << fit_status << ")\n";
    std::cerr << "  Range: x=[" << xmin << ", " << xmax << "], y=[" << ymin << ", " << ymax << "]\n";
    std::cerr << "  Points=" << n << ", Steps=" << nstep << "\n";
  } else {
    const Double_t chi2 = f->GetChisquare();
    const Int_t ndf = f->GetNDF();
    const Double_t chi2ndf = (ndf > 0) ? (chi2 / (Double_t)ndf) : -1.0;
    const Double_t prob = (ndf > 0) ? TMath::Prob(chi2, ndf) : 0.0;
    std::cout << "    Fit successful!\n";
    std::cout << "    Fit quality: chi2=" << chi2
              << ", ndf=" << ndf
              << ", chi2/ndf=" << chi2ndf
              << ", prob=" << prob << "\n";
    // フィット結果の概要を出力（ステップ位置・振幅など）
    std::cout << "    Fitted parameters:\n";
    for (Int_t i = 0; i < nstep; ++i) {
      const Double_t A = f->GetParameter(1 + 3 * i);
      const Double_t t = f->GetParameter(2 + 3 * i);
      const Double_t w = f->GetParameter(3 + 3 * i);
      const Double_t eA = f->GetParError(1 + 3 * i);
      const Double_t et = f->GetParError(2 + 3 * i);
      const Double_t ew = f->GetParError(3 + 3 * i);
      std::cout << "      Step " << i
                << ": A=" << A << " +- " << eA << " [ns]"
                << ", t=" << t << " +- " << et << " [ns]"
                << ", w=" << w << " +- " << ew << " [ns]\n";
    }
  }
 
  for (size_t i = 0; i < x_center.size(); ++i)
    delta_clock[i] = f->Eval(x_center[i]);
  return f;
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
// ProfileX: 各 x ビンで Y の加重平均を計算。mean_y に加え err_y に誤差（加重 RMS/sqrt(N)）を返す
static void ProfileX(const TH2D* h, std::vector<Double_t>& x_center,
                    std::vector<Double_t>& mean_y, std::vector<Double_t>& err_y,
                    Int_t min_entries = 5)
{
  x_center.clear();
  mean_y.clear();
  err_y.clear();
  Int_t nxbins = h->GetNbinsX();
  Int_t nybins = h->GetNbinsY();
  
  for (Int_t ix = 1; ix <= nxbins; ++ix) {
    Double_t xc = h->GetXaxis()->GetBinCenter(ix);
    
    // まず各Yビンのエントリを調べて、最大値を見つける
    Double_t max_bin_content = 0.0;
    for (Int_t iy = 1; iy <= nybins; ++iy) {
      Double_t w = h->GetBinContent(ix, iy);
      if (w > max_bin_content) max_bin_content = w;
    }
    
    // 閾値: 最大エントリの1/5以下は除外
    Double_t threshold = max_bin_content / 5.0;
    
    // 閾値以上のビンのみで加重平均を計算
    Double_t sum = 0.0;
    Double_t sumw = 0.0;
    for (Int_t iy = 1; iy <= nybins; ++iy) {
      Double_t yc = h->GetYaxis()->GetBinCenter(iy);
      Double_t w = h->GetBinContent(ix, iy);
      if (w >= threshold) {  // 閾値以上のみ
        sum += w * yc;
        sumw += w;
      }
    }
    
    if (sumw < min_entries) continue;  // 統計が少ないビンはスキップ
    
    Double_t mean = sum / sumw;
    
    // RMS計算も閾値以上のビンのみ
    Double_t sum_sq = 0.0;
    for (Int_t iy = 1; iy <= nybins; ++iy) {
      Double_t yc = h->GetYaxis()->GetBinCenter(iy);
      Double_t w = h->GetBinContent(ix, iy);
      if (w >= threshold) {
        sum_sq += w * (yc - mean) * (yc - mean);
      }
    }
    
    Double_t rms = (sum_sq > 0 && sumw > 1) ? TMath::Sqrt(sum_sq / sumw) : 0.0;
    Double_t err = (sumw > 0) ? rms / TMath::Sqrt(sumw) : 0.0;
    if (err <= 0.0) err = 1e-6;  // ゼロ誤差はフィットで問題になるので最小値を与える
    
    x_center.push_back(xc);
    mean_y.push_back(mean);
    err_y.push_back(err);
  }
}

//_____________________________________________________________________________
static Double_t ResYToDeltaClock(Double_t mean_resy_mm, Double_t vdrift_mm_ns)
{
  return -mean_resy_mm / vdrift_mm_ns;
}

//_____________________________________________________________________________
// スライス統計が少ないとき Y 軸をまとめてからガウスフィット（ビン統計のガタつき低減）
// 返り値: フィットに使うヒスト（rebin 無しなら py 本体）。*out_rebinned が true なら
//         呼び出し側は (rebinned を delete) と py の両方を解放すること。
static TH1D* RebinnedProjectionForGaussFit(TH1D* py, Int_t cobo, Int_t k, Bool_t* out_rebinned)
{
  *out_rebinned = kFALSE;
  if (!py || py->GetNbinsX() < 8)
    return py;
  const Double_t sum = py->Integral();
  const Int_t nbx = py->GetNbinsX();
  Int_t rb = 1;
  if (sum < 130.0 && nbx >= 16)
    rb = 2;
  if (sum < 60.0 && nbx >= 32)
    rb = 4;
  if (rb <= 1)
    return py;
  TH1D* h = (TH1D*)py->Rebin(rb, Form("tpc_phase_reb_c%d_%d", cobo, k));
  if (h && h->GetNbinsX() >= 4) {
    *out_rebinned = kTRUE;
    return h;
  }
  if (h)
    delete h;
  return py;
}

//_____________________________________________________________________________
// 移動平均で y をスムージング（in-place）。半窓幅 half_window（0 なら何もしない）
// err_y が非 null なら、移動平均後の誤差も更新（err_smoothed = sqrt(sum(err^2))/n）
static void SmoothMovingAverage(std::vector<Double_t>& y, Int_t half_window,
                                 std::vector<Double_t>* err_y = nullptr)
{
  if (half_window <= 0 || y.empty()) return;
  const Int_t n = (Int_t)y.size();
  std::vector<Double_t> out(n);
  std::vector<Double_t> out_err;
  if (err_y && err_y->size() == (size_t)n) out_err.resize(n);
  for (Int_t i = 0; i < n; ++i) {
    Int_t lo = std::max(0, i - half_window);
    Int_t hi = std::min(n - 1, i + half_window);
    Double_t sum = 0.0;
    Double_t sum_sq_err = 0.0;
    Int_t count = 0;
    for (Int_t j = lo; j <= hi; ++j) {
      sum += y[j];
      if (!out_err.empty()) sum_sq_err += (*err_y)[j] * (*err_y)[j];
      ++count;
    }
    out[i] = (count > 0) ? sum / count : y[i];
    if (!out_err.empty())
      out_err[i] = (count > 0) ? TMath::Sqrt(sum_sq_err) / count : (*err_y)[i];
  }
  y = std::move(out);
  if (!out_err.empty()) *err_y = std::move(out_err);
}

//_____________________________________________________________________________
int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <TpcPhase.root> [tpcbcout.root] [--fit-step N] [--smooth N] [--vdrift V] ...\n"
              << "\n"
              << "  --fit-step 0         : create initial flat-zero file (tpcbcout.root not needed)\n"
              << "  --fit-step N (N>0)   : fit with N step functions, need tpcbcout.root\n"
              << "  1 arg:  TpcPhase.root only (for --fit-step 0)\n"
              << "  2 args: tpcbcout.root TpcPhase.root (for --fit-step N>0)\n"
              << "\n"
              << "  --smooth N       : moving average half-window (default: 0 = no smoothing)\n"
              << "  --vdrift V      : drift velocity [mm/ns] (default: 0.055)\n"
              << "  --step-width W  : fixed step width [ns] (smaller => steeper, default: 0.005)\n"
              << "  --free          : free step width in fit (default: fixed)\n"
              << "  --min-entries N : min entries per bin for profile (default: 5)\n"
              << "  --graph-points N: number of points in output TGraph (default: 10000)\n"
              << "  --mode hit|trk  : input 2D hist priority (trk: TPCTrk_ResY_vs_ClockTime_... first)\n";
    return 1;
  }

  std::string tpcbcout_path;
  std::string phase_path;
  Double_t vdrift = DEFAULT_VDRIFT;
  Int_t smooth_half_window = DEFAULT_SMOOTH;
  Int_t fit_nstep = 0;
  Int_t min_entries = 5;
  Int_t graph_points = DEFAULT_GRAPH_POINTS;
  Double_t step_width_ns = 0.005;
  Bool_t free_step_width = kFALSE;
  std::string mode = "hit";

  // 非オプション引数を収集。1個=phaseのみ(--fit-step 0用)、2個=tpcbcout,phase(従来互換)
  std::vector<std::string> positionals;
  for (Int_t i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--mode" && i + 1 < argc) {
      mode = argv[++i];
    } else if (a == "--vdrift" && i + 1 < argc) {
      vdrift = std::atof(argv[++i]);
      if (vdrift <= 0.0) vdrift = DEFAULT_VDRIFT;
    } else if (a == "--smooth" && i + 1 < argc) {
      smooth_half_window = std::atoi(argv[++i]);
      if (smooth_half_window < 0) smooth_half_window = 0;
    } else if (a == "--step-width" && i + 1 < argc) {
      step_width_ns = std::atof(argv[++i]);
      if (!(step_width_ns > 0.0)) step_width_ns = 0.005;
    } else if (a == "--free") {
      free_step_width = kTRUE;
    } else if (a == "--fit-step") {
      fit_nstep = 1;  // デフォルト1個
      if (i + 1 < argc) {
        char* end = nullptr;
        long n = std::strtol(argv[i + 1], &end, 10);
        if (end != argv[i + 1] && n >= 0) {
          fit_nstep = (Int_t)n;
          ++i;
        }
      }
    } else if (a == "--min-entries" && i + 1 < argc) {
      min_entries = std::atoi(argv[++i]);
      if (min_entries < 1) min_entries = 1;
    } else if (a == "--graph-points" && i + 1 < argc) {
      graph_points = std::atoi(argv[++i]);
      if (graph_points < 100) graph_points = DEFAULT_GRAPH_POINTS;
    } else if (a.find("-") != 0) {
      positionals.push_back(a);
    }
  }
  if (positionals.size() >= 1) phase_path = positionals.back();  // 最後が常に phase
  if (positionals.size() >= 2) tpcbcout_path = positionals[0];   // 2個なら 1番目が tpcbcout

  if (phase_path.empty()) {
    std::cerr << "Error: need <TpcPhase.root>\n";
    return 1;
  }

  const std::vector<std::string>& hist_fmts =
      (mode == "trk" || mode == "track" || mode == "TPCTrk") ? HIST_FMTS_TRK : HIST_FMTS_HIT;

  if (tpcbcout_path.empty() && fit_nstep != 0) {
    std::cerr << "Error: need <tpcbcout.root> for fitting (--fit-step N)\n";
    return 1;
  }

  // 入力ファイルを絶対パスに
  namespace fs = std::filesystem;
  if (!tpcbcout_path.empty()) tpcbcout_path = fs::absolute(tpcbcout_path).string();
  phase_path = fs::absolute(phase_path).string();

  // run 番号（連続した数字列）を tpcbcout ファイル名から抽出
  std::string run_id;
  if (!tpcbcout_path.empty()) {
    const std::string fname = fs::path(tpcbcout_path).filename().string();
    for (size_t i = 0; i < fname.size(); ++i) {
      if (std::isdigit(static_cast<unsigned char>(fname[i]))) {
        size_t j = i;
        while (j < fname.size() && std::isdigit(static_cast<unsigned char>(fname[j]))) {
          ++j;
        }
        run_id = fname.substr(i, j - i);
        break;
      }
    }
  }

  // 安全チェック: 出力ファイルが既存のDstファイルっぽい場合は警告
  if (fs::exists(phase_path)) {
    std::string phase_lower = phase_path;
    std::transform(phase_lower.begin(), phase_lower.end(), phase_lower.begin(), ::tolower);
    
    // "dst" を含み、かつ "phase" を含まない場合は危険
    if (phase_lower.find("dst") != std::string::npos && 
        phase_lower.find("phase") == std::string::npos) {
      std::cerr << "\n!!! WARNING !!!" << std::endl;
      std::cerr << "Output file looks like a Dst file (contains 'dst' but not 'phase'):" << std::endl;
      std::cerr << "  " << phase_path << std::endl;
      std::cerr << "This will OVERWRITE the existing file!" << std::endl;
      std::cerr << "Did you swap input/output arguments?" << std::endl;
      std::cerr << "\nCorrect usage:" << std::endl;
      std::cerr << "  " << argv[0] << " <input_tpcbcout.root> <output_TpcPhase.root> --fit-step N" << std::endl;
      std::cerr << "\nType 'yes' to continue anyway, or anything else to abort: " << std::flush;
      
      std::string response;
      std::getline(std::cin, response);
      if (response != "yes") {
        std::cerr << "Aborted." << std::endl;
        return 1;
      }
    }
  }

  // --fit-step 0: 全ゼロの初期ファイル作成（tpcbcout 不要）
  Bool_t flat_zero_mode = (fit_nstep == 0);

  // パラメータ表示
  std::cout << "========================================" << std::endl;
  std::cout << "TPC Phase Correction File Generator" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Output file      : " << phase_path << std::endl;
  if (!tpcbcout_path.empty()) {
    std::cout << "Input file       : " << tpcbcout_path << std::endl;
  }
  std::cout << "Mode             : " << (flat_zero_mode ? "Flat Zero (no fitting)" : "Fitting") << std::endl;
  std::cout << "Histogram mode   : " << mode << std::endl;
  if (!flat_zero_mode) {
    std::cout << "Fit steps        : " << fit_nstep << std::endl;
    std::cout << "Step width [ns]  : " << step_width_ns
              << (free_step_width ? " (initial, free)" : " (fixed)") << std::endl;
    std::cout << "Free width fit   : " << (free_step_width ? "ON" : "OFF") << std::endl;
    std::cout << "Smooth window    : " << smooth_half_window << " (half-width)" << std::endl;
    std::cout << "Drift velocity   : " << vdrift << " mm/ns" << std::endl;
    std::cout << "Min entries/bin  : " << min_entries << std::endl;
  }
  std::cout << "Graph points     : " << graph_points << std::endl;
  std::cout << "========================================" << std::endl;

  TFile* fin = nullptr;
  if (!flat_zero_mode && !tpcbcout_path.empty()) {
    fin = TFile::Open(tpcbcout_path.c_str());
    if (!fin || fin->IsZombie()) {
      std::cerr << "Error: cannot open " << tpcbcout_path << std::endl;
      return 1;
    }
    // params.h の tpc_phase_step_init は "runNNNNN-cobo" キー。複合ファイルではファイル名が不正確なことがあるので、
    // tpc ツリーの先頭イベントの run_number を優先する（無ければ従来どおりファイル名由来の run_id）。
    TTree* tr = (TTree*)fin->Get("tpc");
    if (tr && tr->GetBranch("run_number")) {
      try {
        TTreeReader rd(tr);
        TTreeReaderValue<UInt_t> rv(rd, "run_number");
        if (rd.Next()) {
          run_id = TString::Format("%05u", *rv).Data();
          std::cout << "run_id for params.h (from tpc/run_number, first entry): " << run_id << std::endl;
        }
      } catch (...) {
      }
    }
  }

  // 一時ファイルに出力（あとでリネーム）
  std::string tmp_path = phase_path + ".tmp";
  TFile* fout = TFile::Open(tmp_path.c_str(), "recreate");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: cannot create " << tmp_path << std::endl;
    if (fin) { fin->Close(); delete fin; }
    return 1;
  }

  fout->cd();
  Int_t fb_cobo = 0;
  Int_t fb_nstep_fit = 0;
  Double_t fb_p0 = 0.0, fb_p1 = 0.0, fb_p2 = 1.0;
  TTree* t_cobo_fb = new TTree("TpcPhase_CoboFallback",
                                "TPCPRM kCobo fallback; p0,p1,p2 for PhaseShift when nstep_fit==1");
  t_cobo_fb->Branch("cobo", &fb_cobo, "cobo/I");
  t_cobo_fb->Branch("nstep_fit", &fb_nstep_fit, "nstep_fit/I");
  t_cobo_fb->Branch("p0", &fb_p0, "p0/D");
  t_cobo_fb->Branch("p1", &fb_p1, "p1/D");
  t_cobo_fb->Branch("p2", &fb_p2, "p2/D");

  // PDF出力は将来の実装のため、現在は無効化

  Int_t n_created = 0;

  //==========================================================================
  // run ごとの代表振幅 (Δclock) を Y 射影ベースで推定
  //==========================================================================
  Double_t run_amp_ref_ns = 0.0;
  std::vector<Double_t> amp_samples;
  std::vector<StepGuess> step_guesses(NumOfSegCOBO);
  std::vector<bool> has_guess(NumOfSegCOBO, false);
  std::vector<bool> has_manual_param(NumOfSegCOBO, false);  // params.h で手入力されたかどうか

  if (!flat_zero_mode && fin) {
    for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
      TH2D* raw = nullptr;
      for (const auto& fmt : hist_fmts) {
        TString s = Form(fmt.c_str(), cobo);
        raw = GetHist(fin, s);
        if (raw) break;
      }
      if (!raw) continue;
      StepGuess g = EstimateStepFromProjectionY(raw, vdrift);

      // params.h に run-cobo 指定の t0 / 高さ初期値があればそれを優先
      if (!run_id.empty()) {
        const std::string key_run_cobo = run_id + "-" + std::to_string(cobo);
        auto it = param::tpc_phase_step_init.find(key_run_cobo);
        if (it == param::tpc_phase_step_init.end()) {
          const std::string key_default = "default-" + std::to_string(cobo);
          it = param::tpc_phase_step_init.find(key_default);
        }
        if (it != param::tpc_phase_step_init.end() && !it->second.empty()) {
          has_manual_param[cobo] = true;
          const auto& v = it->second;
          // v[0]: t0[ns]（有効なら）
          if (v.size() >= 1 && std::isfinite(v[0])) {
            g.t0 = std::clamp(v[0], -35.0, 35.0);
          }
          // v[1]: 高さ[mm]（上 − 下）。有効なら c_left/c_right を上書き。
          if (v.size() >= 2 && std::isfinite(v[1]) && vdrift > 0.0) {
            // 現在の Δclock 振幅（ns）からおおよその下側 Δclock を維持しつつ高さだけ合わせる
            const Double_t current_amp_ns = g.c_right - g.c_left;
            const Double_t current_mid_ns = 0.5 * (g.c_right + g.c_left);
            const Double_t target_amp_ns  = -v[1] / vdrift; // ResY[mm] -> Δclock[ns]
            g.c_left  = current_mid_ns - 0.5 * target_amp_ns;
            g.c_right = current_mid_ns + 0.5 * target_amp_ns;
          }
        }
      }

      step_guesses[cobo] = g;
      has_guess[cobo] = true;
      Double_t amp = std::fabs(g.c_right - g.c_left);
      if (amp > 0.0) amp_samples.push_back(amp);
    }
    if (!amp_samples.empty()) {
      run_amp_ref_ns = std::accumulate(amp_samples.begin(), amp_samples.end(), 0.0) /
                       (Double_t)amp_samples.size();
      std::cout << "Estimated run-level step amplitude (|Δclock|) ≈ "
                << run_amp_ref_ns << " ns" << std::endl;
    }
  }

  std::cout << "\nProcessing CoBo:" << std::endl;
  for (Int_t cobo = 0; cobo < NumOfSegCOBO; ++cobo) {
    // --fit-step 0: 全ゼロのフラットな TGraph を保存
    if (flat_zero_mode) {
      fout->cd();
      TGraph* g_flat = new TGraph(graph_points);
      g_flat->SetName(Form(PHASE_NAME_FMT, cobo));
      g_flat->SetTitle(Form("TPC Phase CoBo %d (flat zero);Clock Time [ns];#Delta clock [ns]", cobo));
      Double_t dx = 100.0 / (graph_points - 1);  // -50 to 50 ns
      for (Int_t i = 0; i < graph_points; ++i) {
        Double_t x = -50.0 + i * dx;
        g_flat->SetPoint(i, x, 0.0);
      }
      g_flat->Write();
      delete g_flat;
      fb_cobo = cobo;
      fb_nstep_fit = 0;
      fb_p0 = 0.0;
      fb_p1 = 0.0;
      fb_p2 = 1.0;
      t_cobo_fb->Fill();
      std::cout << "  CoBo " << cobo << ": Flat zero created" << std::endl;
      ++n_created;
      continue;
    }

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "  CoBo " << cobo << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "  CoBo " << cobo << std::endl;
    TString hname = "";
    TH2D* raw = nullptr;
    
    // 候補フォーマットを順に試す
    for (const auto& fmt : hist_fmts) {
      TString s = Form(fmt.c_str(), cobo);
      raw = GetHist(fin, s);
      if (raw) {
        hname = s;
        break;
      }
    }

    if (!raw) {
      std::cout << "    Histogram not found (tried all known formats), skipped" << std::endl;
      continue;
    }
    std::cout << "    Found " << hname << ". Processing..." << std::flush;
    TH2D* h = (TH2D*)raw->Clone(Form("%s_tmp", hname.Data()));

    // 事前に推定した Y 射影ベースのステップ初期値（なければゼロ初期値）
    StepGuess guess{0.0, 0.0, 0.0};
    if (cobo >= 0 && cobo < NumOfSegCOBO && has_guess[cobo]) {
      guess = step_guesses[cobo];
    }

    std::vector<Double_t> x_center, mean_y, err_y;
    ProfileX(h, x_center, mean_y, err_y, min_entries);

    Int_t npts = (Int_t)x_center.size();
    if (npts < 1) {
      std::cout << " No valid points, skipped" << std::endl;
      delete h;
      continue;
    }

    std::vector<Double_t> delta_clock(npts);
    std::vector<Double_t> err_delta(npts);
    for (Int_t i = 0; i < npts; ++i) {
      delta_clock[i] = ResYToDeltaClock(mean_y[i], vdrift);
      err_delta[i] = (vdrift > 0) ? err_y[i] / vdrift : 1e-6;
    }
    SmoothMovingAverage(delta_clock, smooth_half_window, &err_delta);

    // --- Save profile points for later visualization ---
    // Profile points after smoothing: (clk=x_center) vs (dclk=delta_clock)
    fout->cd();
    if (!x_center.empty() && delta_clock.size() == x_center.size()) {
      TGraphErrors* g_profile = new TGraphErrors((Int_t)x_center.size());
      g_profile->SetName(Form("TpcPhase_Profile_Cobo%d", cobo));
      g_profile->SetTitle(Form("TPC Phase Profile CoBo %d;Clock Time [ns];#Delta clock [ns]", cobo));
      for (Int_t i = 0; i < (Int_t)x_center.size(); ++i) {
        g_profile->SetPoint(i, x_center[i], delta_clock[i]);
        const Double_t ey = (i < (Int_t)err_delta.size()) ? err_delta[i] : 0.0;
        g_profile->SetPointError(i, 0.0, ey);
      }
      g_profile->Write();
      delete g_profile;
    }

    // ステップ位置の 1D プロファイルベース推定（delta_clock を利用）
    // params.h で明示的に t0 が与えられている CoBo は、その値を尊重し補正しない。
    if (!has_manual_param[cobo]) {
      Double_t t0_profile = EstimateStepPositionFromProfile(x_center, delta_clock);
      if (t0_profile != 0.0) {
        // Y 射影ベースの t0 と組み合わせる（単純平均）
        if (std::abs(guess.t0) > 0.0) {
          guess.t0 = 0.5 * (guess.t0 + t0_profile);
        } else {
          guess.t0 = t0_profile;
        }
      }
    }

    // ステップでフィット（X 範囲を絶対座標 [-40, +40] ns に固定）
    Double_t xmin_all = *std::min_element(x_center.begin(), x_center.end());
    Double_t xmax_all = *std::max_element(x_center.begin(), x_center.end());
    Double_t x_range  = xmax_all - xmin_all;
    if (x_range <= 0.0) {
      std::cout << "    Invalid x-range, skipped" << std::endl;
      delete h;
      continue;
    }

    Double_t fit_xmin = xmin_all;
    Double_t fit_xmax = xmax_all;

    const Double_t half_window = kPhaseStepFitHalfWindowNs;  // [ns]
    fit_xmin = std::max(xmin_all, -half_window);
    fit_xmax = std::min(xmax_all, +half_window);
    if (fit_xmax - fit_xmin < 10.0) {
      // 入力レンジが狭すぎるときだけ最小限のフォールバック
      fit_xmin = xmin_all + 0.10 * x_range;
      fit_xmax = xmax_all - 0.10 * x_range;
    }

    // 実際にフィットに使う点だけ抽出
    std::vector<Double_t> x_fit, dclk_fit, err_fit;
    std::vector<Int_t> idx_fit;  // x_center の元インデックス
    for (Int_t i = 0; i < npts; ++i) {
      const Double_t x = x_center[i];
      if (x < fit_xmin || x > fit_xmax) continue;
      x_fit.push_back(x);
      dclk_fit.push_back(delta_clock[i]);
      err_fit.push_back(err_delta[i]);
      idx_fit.push_back(i);
    }

    // --- gauss + constant をコア点（x_fit）だけで y スライスに対して当て直す ---
    // 目的:
    //   delta_clock(x) 作成時の「左右非対称な分布」由来の mean_ResY 偏りを、各スライスの
    //   gauss(peak) + pol0(定数ベースライン) で吸収し、StepSumFunc のフィットを安定化する。
    //
    // 注意:
    //   ProfileX は別途保存・可視化されるので、ここで更新するのは StepSumFunc に投入する
    //   dclk_fit/err_fit のみ。
    {
      const Double_t gaussian_fit_nsigma = 2.0;
      const Double_t min_step_gauss_sigma = 0.01;
      const Double_t kMinPeakRatioGauss = 1.05; // ほぼ平坦スライスは避ける
      // rebin 後はビンが粗く peak 比の判定が緩むので、しきい値をわずかに上げる
      const Double_t kMinPeakRatioGaussRebinned = 1.08;

      for (Int_t k = 0; k < (Int_t)x_fit.size(); ++k) {
        const Int_t i = idx_fit.empty() ? k : idx_fit[k];
        const Double_t x = x_center[i];
        const Int_t binx = h->GetXaxis()->FindFixBin(x);

        // x ビン 1個だけの Projection を作ってガウス + 定数背景を当て直す
        TH1D* py = h->ProjectionY(Form("tpc_phase_py_c%d_%d", cobo, k), binx, binx);
        if (!py) continue;
        if (py->GetEntries() <= 0) {
          delete py;
          continue;
        }

        Bool_t used_rebin = kFALSE;
        TH1D* py_fit = RebinnedProjectionForGaussFit(py, cobo, k, &used_rebin);
        const Double_t peak_ratio_cut =
            used_rebin ? kMinPeakRatioGaussRebinned : kMinPeakRatioGauss;

        // ほぼ平坦スライスは避ける（=無駄な fit を減らす）
        const Double_t max_content = py_fit->GetMaximum();
        const Double_t mean_content_bins =
            (py_fit->GetNbinsX() > 0) ? (py_fit->Integral() / py_fit->GetNbinsX()) : 0.0;
        if (mean_content_bins <= 0.0 || max_content < mean_content_bins * peak_ratio_cut) {
          if (used_rebin)
            delete py_fit;
          delete py;
          continue;
        }

        const Int_t max_bin = py_fit->GetMaximumBin();
        const Double_t peak_resy = py_fit->GetBinCenter(max_bin);
        Double_t rms_resy = py_fit->GetRMS();
        if (rms_resy < 0.1) rms_resy = 0.5;
        if (rms_resy > 5.0) rms_resy = 2.5;
        // rebin 後は RMS がやや小さめに出ることがあるので、フィット窓の下限を確保
        if (used_rebin)
          rms_resy = TMath::Max(rms_resy, 0.25);

        const Double_t half_resy_mm =
            TMath::Max(gaussian_fit_nsigma * rms_resy, kPhaseSliceGaussHalfWidthNsEquiv * vdrift);
        const Double_t fit_min_resy = peak_resy - half_resy_mm;
        const Double_t fit_max_resy = peak_resy + half_resy_mm;

        if (!(fit_min_resy < fit_max_resy)) {
          if (used_rebin)
            delete py_fit;
          delete py;
          continue;
        }

        TF1 sliceFit(Form("tpc_phase_sliceFit_c%d_%d", cobo, k),
                     "gaus(0)+pol0(3)", fit_min_resy, fit_max_resy);
        sliceFit.SetParameters(max_content, peak_resy,
                                std::max(rms_resy * 0.5, min_step_gauss_sigma),
                                0.0);
        sliceFit.SetParNames("Constant", "Mean", "Sigma", "Background");
        sliceFit.SetParLimits(2, min_step_gauss_sigma, 20.0);

        const Int_t fit_status = py_fit->Fit(&sliceFit, "Q");
        if (fit_status == 0) {
          const Double_t mean_resy = sliceFit.GetParameter(1);
          Double_t mean_err = sliceFit.GetParError(1);
          if (used_rebin && mean_err > 0.0)
            mean_err *= 1.15; // 相関をざっくり反映した誤差の下限寄せ
          if (vdrift > 0.0 && std::isfinite(mean_resy) && std::isfinite(mean_err) && mean_err > 0.0) {
            dclk_fit[k] = ResYToDeltaClock(mean_resy, vdrift);
            err_fit[k] = mean_err / vdrift;
          }
        }
        if (used_rebin)
          delete py_fit;
        delete py;
      }

      // コア点列に沿った移動平均（プロット・StepSum への入力のガタつき低減）
      static const Int_t kSmoothAfterGaussHalfWin = 1;
      SmoothMovingAverage(dclk_fit, kSmoothAfterGaussHalfWin, &err_fit);
    }

    // Save points that will be used in fitting (core points).
    fout->cd();
    if (!x_fit.empty() && x_fit.size() == dclk_fit.size()) {
      TGraphErrors* g_fitused = new TGraphErrors((Int_t)x_fit.size());
      g_fitused->SetName(Form("TpcPhase_ProfileFitUsed_Cobo%d", cobo));
      g_fitused->SetTitle(Form("TPC Phase FitUsed CoBo %d;Clock Time [ns];#Delta clock [ns]", cobo));
      for (Int_t i = 0; i < (Int_t)x_fit.size(); ++i) {
        g_fitused->SetPoint(i, x_fit[i], dclk_fit[i]);
        const Double_t ey = (i < (Int_t)err_fit.size()) ? err_fit[i] : 0.0;
        g_fitused->SetPointError(i, 0.0, ey);
      }
      g_fitused->Write();
      delete g_fitused;
    }

    if (x_fit.size() < 5) {
      std::cout << "    Too few points in core range, skipped" << std::endl;
      delete h;
      continue;
    }

    Double_t xmin = *std::min_element(x_fit.begin(), x_fit.end());
    Double_t xmax = *std::max_element(x_fit.begin(), x_fit.end());

    // run-cobo ごとの振幅初期値（ResidualY[mm] → Δclock[ns] に変換）
    Double_t amp_override_ns = 0.0;
    if (!run_id.empty()) {
      const std::string key_run_cobo = run_id + "-" + std::to_string(cobo);
      auto it = param::tpc_phase_step_init.find(key_run_cobo);
      if (it == param::tpc_phase_step_init.end()) {
        const std::string key_default = "default-" + std::to_string(cobo);
        it = param::tpc_phase_step_init.find(key_default);
      }
      if (it != param::tpc_phase_step_init.end() && it->second.size() >= 2 &&
          std::isfinite(it->second[1]) && vdrift > 0.0) {
        const Double_t height_mm = it->second[1];
        amp_override_ns = -height_mm / vdrift;
      }
    }

    TF1* init_func = nullptr;
    TF1* fit_func = FitStepSum(fit_nstep, x_fit, dclk_fit,
                               xmin, xmax, &err_fit, &guess,
                               run_amp_ref_ns, amp_override_ns, step_width_ns, free_step_width, &init_func);

    fout->cd();
    fb_cobo = cobo;
    if (!fit_func) {
      fb_nstep_fit = -1;
      fb_p0 = fb_p1 = 0.0;
      fb_p2 = 1.0;
      t_cobo_fb->Fill();
    } else {
      // TF1 を細かく Eval して TGraph に変換
      TGraph* g_save = new TGraph(graph_points);
      g_save->SetName(Form(PHASE_NAME_FMT, cobo));
      g_save->SetTitle(Form("TPC Phase CoBo %d;Clock Time [ns];#Delta clock [ns]", cobo));
      Double_t dx = (xmax_all - xmin_all) / (graph_points - 1);
      for (Int_t i = 0; i < graph_points; ++i) {
        Double_t x = xmin_all + i * dx;
        Double_t y = fit_func->Eval(x);
        g_save->SetPoint(i, x, y);
      }
      g_save->Write();

      if (fit_nstep == 1) {
        fb_nstep_fit = 1;
        fb_p0 = fit_func->GetParameter(1);
        fb_p1 = fit_func->GetParameter(2);
        fb_p2 = fit_func->GetParameter(3);
      } else {
        fb_nstep_fit = fit_nstep;
        fb_p0 = fb_p1 = 0.0;
        fb_p2 = 1.0;
      }
      t_cobo_fb->Fill();

      // フィット関数自体も参照用に保存
      if (fit_nstep > 0) {
        TF1* f_fit = (TF1*)fit_func->Clone(Form(FIT_NAME_FMT, cobo));
        f_fit->SetTitle(Form("TPC Phase Fit CoBo %d;Clock Time [ns];#Delta clock [ns]", cobo));
        f_fit->Write();
        if (init_func) {
          TF1* f_init = (TF1*)init_func->Clone(Form("TpcPhase_Init_Cobo%d", cobo));
          f_init->SetTitle(Form("TPC Phase Init CoBo %d;Clock Time [ns];#Delta clock [ns]", cobo));
          f_init->Write();
          delete f_init;
        }
      }

      delete g_save;
      if (init_func) delete init_func;
      delete fit_func;
      std::cout << " Done (" << npts << " profile points -> " << graph_points << " graph points)" << std::endl;
    }
    delete h;
    ++n_created;
  }

  if (fin) { fin->Close(); delete fin; }
  fout->cd();
  t_cobo_fb->Write();
  delete t_cobo_fb;
  fout->Close();
  delete fout;

  // tmp → 本番にリネーム
  fs::rename(tmp_path, phase_path);

  // PDF生成は将来の実装として保留
  // 現在はROOTファイルを直接開いて結果を確認してください

  std::cout << "\n========================================" << std::endl;
  std::cout << "Result: Created " << n_created << " CoBo phase corrections" << std::endl;
  std::cout << "Output: " << phase_path << std::endl;
  std::cout << "========================================" << std::endl;
  return 0;
}
