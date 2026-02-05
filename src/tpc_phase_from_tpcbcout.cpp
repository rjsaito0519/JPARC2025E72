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
#include <TKey.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

//_____________________________________________________________________________
// HistTools.cc / DetectorID.hh と合わせる
static const Int_t NumOfSegCOBO = 8;

// TPCParamMan が期待する TF1 名（補正用）
static const char* PHASE_NAME_FMT = "TpcPhase_Cobo%d";
// フィット関数の結果を保存する名前（参照用）
static const char* FIT_NAME_FMT = "TpcPhase_Fit_Cobo%d";

// ResY [mm] -> Δclock [ns] の変換で使う drift velocity [mm/ns]
static const Double_t DEFAULT_VDRIFT = 0.055;

// デフォルトのスムージング半窓幅（0 = スムージングなし）
static const Int_t DEFAULT_SMOOTH = 0;

// TGraph の点数（TF1 を細かく Eval する）
static const Int_t DEFAULT_GRAPH_POINTS = 10000;

// ステップ関数の重ね合わせ: dclk(x) = C + sum_i A_i * Freq((x - t_i) / w_i)
// パラメータ: p[0]=N, p[1]=C (定数項), p[2+3*i]=A_i, p[3+3*i]=t_i, p[4+3*i]=w_i
static Double_t StepSumFunc(Double_t* x, Double_t* p)
{
  const Double_t clk = x[0];
  const Int_t n = (Int_t)(p[0] + 0.5);
  if (n <= 0) return p[1];  // 定数項のみ
  
  Double_t s = p[1];  // 定数項（オフセット）
  for (Int_t i = 0; i < n; ++i) {
    Double_t A = p[2 + 3 * i];
    Double_t t = p[3 + 3 * i];
    Double_t w = p[4 + 3 * i];
    if (w <= 0.0) w = 1e-6;
    s += A * TMath::Freq((clk - t) / w);
  }
  return s;
}

//_____________________________________________________________________________
// N段ステップ関数でフィッティング
static TF1* FitStepSum(Int_t nstep, std::vector<Double_t>& x_center,
                       std::vector<Double_t>& delta_clock,
                       Double_t xmin, Double_t xmax,
                       const std::vector<Double_t>* err_delta = nullptr)
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
  const Int_t npar = 2 + 3 * nstep;  // N, C, A_i, t_i, w_i
  TF1* f = new TF1("step_sum", StepSumFunc, xmin, xmax, npar);
  
  // パラメータ0: ステップ数（固定）
  f->SetParName(0, "N");
  f->FixParameter(0, (Double_t)nstep);
  
  // パラメータ1: 定数項オフセット（最初のステップ前の中央値）
  std::vector<Double_t> y_before_first;
  for (size_t j = 0; j < x_center.size(); ++j) {
    if (x_center[j] < step_positions[0]) {
      y_before_first.push_back(delta_clock[j]);
    }
  }
  Double_t offset = y_before_first.empty() ? ymin 
                    : y_before_first[y_before_first.size() / 2];
  
  f->SetParName(1, "C");
  f->SetParameter(1, offset);
  f->SetParLimits(1, ymin - (ymax - ymin), ymax + (ymax - ymin));
  
  // パラメータ2以降: 各ステップの振幅、位置、幅
  const Double_t fixed_width = 0.01;      // 固定幅 [ns]: 超急峻
  const Double_t position_range = 5.0;    // 位置調整範囲 [ns]
  
  for (Int_t i = 0; i < nstep; ++i) {
    f->SetParName(2 + 3 * i, Form("A%d", i));
    f->SetParName(3 + 3 * i, Form("t%d", i));
    f->SetParName(4 + 3 * i, Form("w%d", i));
    
    f->SetParameter(2 + 3 * i, step_amplitudes[i]);  // 振幅（フィッティング）
    f->SetParameter(3 + 3 * i, step_positions[i]);    // 位置（微調整可）
    f->SetParLimits(3 + 3 * i, step_positions[i] - position_range, 
                                step_positions[i] + position_range);
    f->FixParameter(4 + 3 * i, fixed_width);          // 幅（固定）
  }
  
  std::cout << "    Fit settings: width=" << fixed_width 
            << " ns (fixed), position range=±" << position_range << " ns\n";

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
    std::cout << "    Fit successful!\n";
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
              << "  --min-entries N : min entries per bin for profile (default: 5)\n"
              << "  --graph-points N: number of points in output TGraph (default: 10000)\n";
    return 1;
  }

  std::string tpcbcout_path;
  std::string phase_path;
  Double_t vdrift = DEFAULT_VDRIFT;
  Int_t smooth_half_window = DEFAULT_SMOOTH;
  Int_t fit_nstep = 0;
  Int_t min_entries = 5;
  Int_t graph_points = DEFAULT_GRAPH_POINTS;

  // 非オプション引数を収集。1個=phaseのみ(--fit-step 0用)、2個=tpcbcout,phase(従来互換)
  std::vector<std::string> positionals;
  for (Int_t i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--vdrift" && i + 1 < argc) {
      vdrift = std::atof(argv[++i]);
      if (vdrift <= 0.0) vdrift = DEFAULT_VDRIFT;
    } else if (a == "--smooth" && i + 1 < argc) {
      smooth_half_window = std::atoi(argv[++i]);
      if (smooth_half_window < 0) smooth_half_window = 0;
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
  if (tpcbcout_path.empty() && fit_nstep != 0) {
    std::cerr << "Error: need <tpcbcout.root> for fitting (--fit-step N)\n";
    return 1;
  }

  // 入力ファイルを絶対パスに
  namespace fs = std::filesystem;
  if (!tpcbcout_path.empty()) tpcbcout_path = fs::absolute(tpcbcout_path).string();
  phase_path = fs::absolute(phase_path).string();

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
  if (!flat_zero_mode) {
    std::cout << "Fit steps        : " << fit_nstep << std::endl;
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
  }

  // 一時ファイルに出力（あとでリネーム）
  std::string tmp_path = phase_path + ".tmp";
  TFile* fout = TFile::Open(tmp_path.c_str(), "recreate");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: cannot create " << tmp_path << std::endl;
    if (fin) { fin->Close(); delete fin; }
    return 1;
  }

  static const char* hist_fmt = "TPC_ResidualY_BcOut_vs_ClockTime_CoBo%d_RawClock";

  // PDF出力は将来の実装のため、現在は無効化

  Int_t n_created = 0;
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
      std::cout << "  CoBo " << cobo << ": Flat zero created" << std::endl;
      ++n_created;
      continue;
    }

    TString hname = Form(hist_fmt, cobo);
    TH2D* raw = GetHist(fin, hname);
    if (!raw) {
      std::cout << "  CoBo " << cobo << ": Histogram not found, skipped" << std::endl;
      continue;
    }
    std::cout << "  CoBo " << cobo << ": Processing..." << std::flush;
    TH2D* h = (TH2D*)raw->Clone(Form("%s_tmp", hname.Data()));

    std::vector<Double_t> x_center, mean_y, err_y;
    ProfileX(h, x_center, mean_y, err_y, min_entries);
    delete h;

    Int_t npts = (Int_t)x_center.size();
    if (npts < 1) {
      std::cout << " No valid points, skipped" << std::endl;
      continue;
    }

    std::vector<Double_t> delta_clock(npts);
    std::vector<Double_t> err_delta(npts);
    for (Int_t i = 0; i < npts; ++i) {
      delta_clock[i] = ResYToDeltaClock(mean_y[i], vdrift);
      err_delta[i] = (vdrift > 0) ? err_y[i] / vdrift : 1e-6;
    }
    SmoothMovingAverage(delta_clock, smooth_half_window, &err_delta);

    // ステップでフィット
    Double_t xmin = *std::min_element(x_center.begin(), x_center.end());
    Double_t xmax = *std::max_element(x_center.begin(), x_center.end());
    TF1* fit_func = FitStepSum(fit_nstep, x_center, delta_clock, xmin, xmax, &err_delta);

    fout->cd();
    if (fit_func) {
      // TF1 を細かく Eval して TGraph に変換
      TGraph* g_save = new TGraph(graph_points);
      g_save->SetName(Form(PHASE_NAME_FMT, cobo));
      g_save->SetTitle(Form("TPC Phase CoBo %d;Clock Time [ns];#Delta clock [ns]", cobo));
      Double_t dx = (xmax - xmin) / (graph_points - 1);
      for (Int_t i = 0; i < graph_points; ++i) {
        Double_t x = xmin + i * dx;
        Double_t y = fit_func->Eval(x);
        g_save->SetPoint(i, x, y);
      }
      g_save->Write();
      
      // フィット関数自体も参照用に保存
      if (fit_nstep > 0) {
        TF1* f_fit = (TF1*)fit_func->Clone(Form(FIT_NAME_FMT, cobo));
        f_fit->SetTitle(Form("TPC Phase Fit CoBo %d;Clock Time [ns];#Delta clock [ns]", cobo));
        f_fit->Write();
      }
      
      delete g_save;
      delete fit_func;
      std::cout << " Done (" << npts << " profile points -> " << graph_points << " graph points)" << std::endl;
    }
    ++n_created;
  }

  if (fin) { fin->Close(); delete fin; }
  fout->Close(); delete fout;

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
