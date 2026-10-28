#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TLine.h>
#include <TMath.h>
#include <TObject.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TPRegexp.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "ana_helper.h"
#include "paths.h"
#include "progress_bar.h"
#include "TPCPadHelper.hh"
#include "tpc_gain_tree.h"
#include "tpc_pad_fill.h"

namespace fs = std::filesystem;

struct TpcParamEntry {
  Int_t layer{};
  Int_t row{};
  Int_t aty{};
  std::vector<Double_t> p;
  TString raw_line;
  Bool_t is_comment{kTRUE};
};

// Only target_mpv, threshold and replace_mode are set from the command line; the rest are fixed.
struct GainFitConfig {
  Double_t target_mpv = 200.0;
  Int_t threshold = 20;
  Double_t stat_range = 500.0;
  Double_t min_peak_ratio = 1.05;
  Double_t fit_nsigma = 2.0;
  Double_t local_peak_half = 50.0;
  Int_t local_peak_sep = 3;
  Double_t min_width = 5.0;
  Double_t max_chi2_ndf = 20.0;
  Int_t min_fit_ndf = 2;
  Double_t min_fit_prob = 0.0;
  Double_t mpv_min = 100.0;
  Double_t mpv_max = 0.0;
  Bool_t replace_mode = kFALSE;
};

static Double_t LocalWeightedRmsAroundPeak(const TH1D* h, Double_t peak_x, Double_t half_width);

struct GainFitResult {
  Bool_t ok = kFALSE;
  Double_t mpv = 0.0;
  Double_t width = 0.0;
  Double_t scale = 1.0;
  Double_t chi2ndf = 0.0;
  Double_t fit_min = 0.0;
  Double_t fit_max = 0.0;
};

// 各ビンを中心に ±sep の窓で局所最大。x_min 未満は候補から除外。
static Int_t DominantLocalMaximumBin(const TH1D* h, Int_t sep, Double_t x_min = -1e30) {
  if (!h) return 0;
  if (sep < 1) {
    Int_t best = h->GetMaximumBin();
    for (Int_t i = 1; i <= h->GetNbinsX(); ++i) {
      if (h->GetBinCenter(i) < x_min) continue;
      if (h->GetBinContent(i) > h->GetBinContent(best)) best = i;
    }
    return best;
  }
  const Int_t nb = h->GetNbinsX();
  Double_t best_y = -1.0;
  Int_t best_bin = h->GetMaximumBin();
  for (Int_t i = 1; i <= nb; ++i) {
    if (h->GetBinCenter(i) < x_min) continue;
    const Double_t y = h->GetBinContent(i);
    if (y <= 0.0) continue;
    Bool_t is_peak = kTRUE;
    const Int_t j0 = std::max(1, i - sep);
    const Int_t j1 = std::min(nb, i + sep);
    for (Int_t j = j0; j <= j1; ++j) {
      if (j != i && h->GetBinContent(j) > y) {
        is_peak = kFALSE;
        break;
      }
    }
    if (!is_peak) continue;
    const Double_t xi = TMath::Abs(h->GetBinCenter(i));
    const Double_t xb = TMath::Abs(h->GetBinCenter(best_bin));
    if (y > best_y + 1e-12) {
      best_y = y;
      best_bin = i;
    } else if (TMath::Abs(y - best_y) <= 1e-12 && xi < xb) {
      best_bin = i;
    }
  }
  if (best_y < 0.0 || h->GetBinCenter(best_bin) < x_min) {
    for (Int_t i = 1; i <= nb; ++i) {
      if (h->GetBinCenter(i) < x_min) continue;
      if (h->GetBinContent(i) > best_y) {
        best_y = h->GetBinContent(i);
        best_bin = i;
      }
    }
  }
  return best_bin;
}

static Int_t SelectGainPeakBin(const TH1D* h, const GainFitConfig& cfg) {
  if (!h) return 0;
  const Double_t x_floor = (cfg.mpv_min > 0.0) ? cfg.mpv_min : h->GetXaxis()->GetXmin();
  if (cfg.local_peak_sep >= 1) {
    return DominantLocalMaximumBin(h, cfg.local_peak_sep, x_floor);
  }
  Int_t best_bin = 0;
  Double_t best_y = -1.0;
  for (Int_t i = 1; i <= h->GetNbinsX(); ++i) {
    if (h->GetBinCenter(i) < x_floor) continue;
    const Double_t y = h->GetBinContent(i);
    if (y > best_y) {
      best_y = y;
      best_bin = i;
    }
  }
  return best_bin > 0 ? best_bin : h->GetMaximumBin();
}

static Bool_t ComputeGainFitWindow(const TH1D* h_fit, const GainFitConfig& cfg, Int_t peak_bin,
                                   Double_t mpv_center_hint, Double_t& mpv_init, Double_t& width_init,
                                   Double_t& fit_min, Double_t& fit_max) {
  if (!h_fit || peak_bin < 1) return kFALSE;
  const Double_t peak_x = h_fit->GetBinCenter(peak_bin);
  const Double_t peak_y = h_fit->GetBinContent(peak_bin);
  const Double_t n_entries = h_fit->GetEntries();
  const Double_t xmin = h_fit->GetXaxis()->GetXmin();
  const Double_t xmax = h_fit->GetXaxis()->GetXmax();
  const Double_t x_floor = (cfg.mpv_min > 0.0) ? cfg.mpv_min : xmin;

  const Double_t rms_ref_x =
      (mpv_center_hint > 0.0 && std::isfinite(mpv_center_hint)) ? mpv_center_hint : peak_x;
  Double_t rms_local = LocalWeightedRmsAroundPeak(h_fit, rms_ref_x, cfg.local_peak_half);
  const Double_t rms_global = h_fit->GetRMS();
  if (!(rms_local > 0.0)) rms_local = rms_global;
  if (!(rms_local > 0.0)) rms_local = 15.0;
  // 全幅 RMS で窓を広げると layer ヒストで実質 0–800 に近づくため、局所 RMS を主に使う。
  Double_t rms_win = std::max(rms_local, cfg.min_width);
  const Bool_t low_stat = (n_entries < 1500.0) || (peak_y < 25.0);
  if (rms_global > 0.0) {
    rms_win = std::max(rms_win, (low_stat ? 0.55 : 0.35) * rms_global);
  } else if (low_stat) {
    rms_win = std::max(rms_win, 45.0);
  }

  mpv_init = (mpv_center_hint > 0.0 && std::isfinite(mpv_center_hint)) ? mpv_center_hint : peak_x;
  if (low_stat && mpv_center_hint > 0.0 && std::isfinite(mpv_center_hint)) {
    // 低統計 pad は layer seed と peak の加重平均で窓中心を安定化
    mpv_init = 0.55 * mpv_center_hint + 0.45 * peak_x;
  }
  width_init = std::max(cfg.min_width, std::max(0.5, 0.45 * rms_local));
  if (cfg.mpv_min > 0.0) mpv_init = std::max(mpv_init, cfg.mpv_min);
  if (cfg.mpv_max > 0.0) mpv_init = std::min(mpv_init, cfg.mpv_max);
  mpv_init = std::min(std::max(mpv_init, xmin), xmax);

  const Double_t sigma_core = std::max(cfg.fit_nsigma * rms_win, 2.6 * width_init);
  const Double_t left_sigma = std::max(12.0, 0.75 * sigma_core);
  const Double_t right_sigma = std::max(35.0, (low_stat ? 3.0 : 2.5) * sigma_core);
  fit_min = mpv_init - left_sigma;
  fit_max = mpv_init + right_sigma;
  fit_min = std::max(fit_min, x_floor);
  fit_min = std::max(fit_min, xmin);
  fit_max = std::min(fit_max, xmax);
  if (cfg.mpv_max > 0.0) fit_max = std::min(fit_max, cfg.mpv_max + 1.0 * right_sigma);

  Double_t min_span = low_stat ? 95.0 : 55.0;
  if (peak_y < 15.0) min_span = std::max(min_span, 110.0);
  if (low_stat && rms_global > 0.0) {
    min_span = std::max(min_span, std::min(200.0, 1.10 * rms_global));
  }
  const Double_t max_span = std::max(min_span, std::min(220.0, 5.5 * rms_win));
  if (fit_max - fit_min > max_span) {
    fit_min = std::max(mpv_init - 0.32 * max_span, x_floor);
    fit_max = std::min(mpv_init + 0.68 * max_span, xmax);
  }
  if (fit_max - fit_min < min_span) {
    fit_min = std::max(mpv_init - 0.35 * min_span, x_floor);
    fit_max = std::min(mpv_init + 0.65 * min_span, xmax);
  }
  return fit_max > fit_min;
}

static Double_t LocalWeightedRmsAroundPeak(const TH1D* h, Double_t peak_x, Double_t half_width) {
  if (!h || half_width <= 0.0) return -1.0;
  Double_t sumw = 0.0, sumxw = 0.0, sumx2w = 0.0;
  const Int_t nb = h->GetNbinsX();
  for (Int_t i = 1; i <= nb; ++i) {
    const Double_t x = h->GetBinCenter(i);
    if (TMath::Abs(x - peak_x) > half_width) continue;
    const Double_t w = h->GetBinContent(i);
    if (w <= 0.0) continue;
    sumw += w;
    sumxw += w * x;
    sumx2w += w * x * x;
  }
  if (sumw < 1.0) return -1.0;
  const Double_t mean = sumxw / sumw;
  const Double_t var = sumx2w / sumw - mean * mean;
  if (!(var > 0.0) || !std::isfinite(var)) return -1.0;
  const Double_t rms = TMath::Sqrt(var);
  return std::isfinite(rms) ? rms : -1.0;
}

static Bool_t GainFitNeedsRefit(Double_t peak_x, Double_t mpv, Double_t fit_min, Double_t fit_max, Double_t peak_y,
                                Double_t n_entries) {
  if (!(mpv > 0.0) || !std::isfinite(mpv)) return kFALSE;
  // 低統計では再 fit で窓がさらに崩れやすいのでスキップ
  if (n_entries < 800.0 || peak_y < 20.0) return kFALSE;
  const Double_t span = fit_max - fit_min;
  if (!(span > 0.0)) return kFALSE;
  const Double_t pull_thr = std::max(12.0, 0.10 * std::max(peak_x, mpv));
  if (TMath::Abs(mpv - peak_x) > pull_thr) return kTRUE;
  if (peak_x < fit_min || peak_x > fit_max) return kTRUE;
  const Double_t edge = 0.12 * span;
  if (mpv <= fit_min + edge || mpv >= fit_max - edge) return kTRUE;
  return kFALSE;
}

// 1 回目 fit 後: 窓を peak 中心に組み直して再 fit（MPV が窓端・peak 外れのとき）。
static Bool_t ApplyGainRefitPass(TH1D* h_fit, TF1*& fFit, const char* fit_name, const GainFitConfig& cfg,
                                 Int_t peak_bin, Double_t peak_x, Double_t peak_y, Double_t mean_content,
                                 Double_t& fit_min, Double_t& fit_max, Double_t& width_init, Double_t& mpv,
                                 Double_t& width, Int_t& ndf, Double_t& chi2ndf) {
  if (!GainFitNeedsRefit(peak_x, mpv, fit_min, fit_max, peak_y, h_fit->GetEntries())) return kFALSE;

  const Double_t refit_center = peak_x;
  Double_t mpv_reinit = refit_center;
  width_init = std::max(cfg.min_width, width);
  if (!ComputeGainFitWindow(h_fit, cfg, peak_bin, refit_center, mpv_reinit, width_init, fit_min, fit_max)) {
    return kFALSE;
  }

  delete fFit;
  fFit = new TF1(fit_name, "landau(0)+pol0(3)", fit_min, fit_max);
  fFit->SetParNames("Norm", "MPV", "Width", "Background");
  fFit->SetNpx(1000);
  const Double_t fit_span = fit_max - fit_min;
  width_init = std::min(std::max(width_init, 0.08 * fit_span), 0.30 * fit_span);

  fFit->SetRange(fit_min, fit_max);
  fFit->SetParameters(std::max(peak_y, 1.0), mpv_reinit, width_init, std::max(0.0, mean_content * 0.2));
  fFit->SetParLimits(0, 0.0, peak_y * 100.0 + 1.0);
  fFit->SetParLimits(1, fit_min, fit_max);
  fFit->SetParLimits(2, 0.35, std::max(1.0, 0.35 * fit_span));
  fFit->SetParLimits(3, 0.0, peak_y * 10.0 + 1.0);
  if (h_fit->Fit(fFit, "Q0") != 0) return kFALSE;

  mpv = fFit->GetParameter(1);
  width = fFit->GetParameter(2);
  ndf = fFit->GetNDF();
  const Double_t chi2 = fFit->GetChisquare();
  chi2ndf = (ndf > 0) ? (chi2 / static_cast<Double_t>(ndf)) : 0.0;
  return kTRUE;
}

static Bool_t FitLandauGain(TH1D* h_fit, Int_t layer, Int_t row, const GainFitConfig& cfg, GainFitResult& res,
                            TF1** fFit_out = nullptr, Bool_t strict_quality = kTRUE) {
  res = GainFitResult{};
  if (!h_fit || h_fit->GetEntries() <= 0) return kFALSE;

  const Int_t bmin_stat = h_fit->FindBin(-cfg.stat_range);
  const Int_t bmax_stat = h_fit->FindBin(cfg.stat_range);
  const Double_t n_stat = h_fit->Integral(bmin_stat, bmax_stat);
  if (strict_quality && cfg.threshold > 0 && n_stat < cfg.threshold) return kFALSE;

  const Int_t peak_bin = SelectGainPeakBin(h_fit, cfg);
  const Double_t peak_x = h_fit->GetBinCenter(peak_bin);
  const Double_t peak_y = h_fit->GetBinContent(peak_bin);
  const Int_t n_bins_window = bmax_stat - bmin_stat + 1;
  const Double_t mean_content = (n_bins_window > 0) ? (n_stat / n_bins_window) : 0.0;
  if (strict_quality && (mean_content <= 0.0 || peak_y < mean_content * cfg.min_peak_ratio)) return kFALSE;

  Double_t mpv_init = peak_x;
  Double_t width_init = cfg.min_width;
  Double_t fit_min = 0.0;
  Double_t fit_max = 0.0;
  Double_t mpv_hint = -1.0;
  if (!ComputeGainFitWindow(h_fit, cfg, peak_bin, mpv_hint, mpv_init, width_init, fit_min, fit_max)) return kFALSE;

  const TString fit_name =
      (row >= 0) ? Form("fLandau_L%02d_R%03d", layer, row) : Form("fLandau_L%02d_layer", layer);
  TF1* fFit = new TF1(fit_name.Data(), "landau(0)+pol0(3)", fit_min, fit_max);
  fFit->SetParNames("Norm", "MPV", "Width", "Background");
  fFit->SetNpx(1000);
  const Double_t fit_span = fit_max - fit_min;
  width_init = std::min(std::max(width_init, 0.08 * fit_span), 0.30 * fit_span);

  auto fit_once = [&](Double_t mpv_seed) {
    fFit->SetRange(fit_min, fit_max);
    fFit->SetParameters(std::max(peak_y, 1.0), mpv_seed, width_init, std::max(0.0, mean_content * 0.2));
    fFit->SetParLimits(0, 0.0, peak_y * 100.0 + 1.0);
    fFit->SetParLimits(1, fit_min, fit_max);
    fFit->SetParLimits(2, 0.35, std::max(1.0, 0.35 * fit_span));
    fFit->SetParLimits(3, 0.0, peak_y * 10.0 + 1.0);
    return h_fit->Fit(fFit, "Q0");
  };

  if (fit_once(mpv_init) != 0) {
    delete fFit;
    return kFALSE;
  }

  Double_t mpv = fFit->GetParameter(1);
  Double_t width = fFit->GetParameter(2);
  Int_t ndf = fFit->GetNDF();
  Double_t chi2 = fFit->GetChisquare();
  Double_t chi2ndf = (ndf > 0) ? (chi2 / static_cast<Double_t>(ndf)) : 0.0;
  ApplyGainRefitPass(h_fit, fFit, fit_name.Data(), cfg, peak_bin, peak_x, peak_y, mean_content, fit_min,
                     fit_max, width_init, mpv, width, ndf, chi2ndf);
  chi2 = fFit->GetChisquare();
  chi2ndf = (ndf > 0) ? (chi2 / static_cast<Double_t>(ndf)) : 0.0;
  Double_t fit_prob = (ndf > 0) ? TMath::Prob(chi2, ndf) : 0.0;

  Bool_t quality_ok = kTRUE;
  if (!(mpv > 0.0) || !std::isfinite(mpv)) quality_ok = kFALSE;
  if (strict_quality) {
    if (cfg.mpv_min > 0.0 && mpv < cfg.mpv_min) quality_ok = kFALSE;
    if (cfg.mpv_max > 0.0 && mpv > cfg.mpv_max) quality_ok = kFALSE;
    if (width < cfg.min_width) quality_ok = kFALSE;
    if (cfg.min_fit_ndf > 0 && ndf < cfg.min_fit_ndf) quality_ok = kFALSE;
    if (cfg.max_chi2_ndf > 0.0 && ndf > 0 && chi2ndf > cfg.max_chi2_ndf) quality_ok = kFALSE;
    if (cfg.min_fit_prob > 0.0 && ndf > 0 && fit_prob < cfg.min_fit_prob) quality_ok = kFALSE;
  }
  if (!quality_ok) {
    delete fFit;
    return kFALSE;
  }

  Double_t scale = cfg.target_mpv / mpv;
  if (!std::isfinite(scale)) scale = 1.0;
  if (strict_quality && (scale < (1.0 / 3.0) || scale > 10.0)) scale = 1.0;

  res.ok = kTRUE;
  res.mpv = mpv;
  res.width = width;
  res.scale = scale;
  res.chi2ndf = chi2ndf;
  res.fit_min = fit_min;
  res.fit_max = fit_max;
  if (fFit_out) {
    *fFit_out = fFit;
  } else {
    delete fFit;
  }
  return kTRUE;
}

// QA: ヒストは全 x レンジ。赤線は [fit_min, fit_max] だけ（TF1::SetRange は描画に効かないので定義域を切る）。
static void DrawGainCalibPanel(TH1D* h, TF1* fFit, Bool_t fit_ok, Double_t fit_min, Double_t fit_max,
                               Double_t mpv, std::vector<TObject*>& page_owned) {
  if (!h) return;
  if (h->GetListOfFunctions()) h->GetListOfFunctions()->Clear();

  h->GetXaxis()->SetRangeUser(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  h->Draw();
  if (gPad) gPad->Update();

  if (!fit_ok || !fFit || !(fit_max > fit_min)) return;

  TF1* fDraw = new TF1(Form("%s_drawpanel", fFit->GetName()), "landau(0)+pol0(3)", fit_min, fit_max);
  for (Int_t p = 0; p < fFit->GetNpar(); ++p) {
    fDraw->SetParameter(p, fFit->GetParameter(p));
  }
  fDraw->SetLineColor(kRed);
  fDraw->SetLineWidth(1);
  fDraw->SetNpx(400);
  fDraw->Draw("same");
  page_owned.push_back(fDraw);

  if (mpv > 0.0 && std::isfinite(mpv) && gPad) {
    TLine* v = new TLine(mpv, gPad->GetUymin(), mpv, gPad->GetUymax());
    v->SetLineStyle(2);
    v->SetLineColor(kGreen + 2);
    v->SetLineWidth(1);
    v->Draw("same");
    page_owned.push_back(v);
  }
}

static void PrintUsage(const char* prog) {
  std::cerr << "Usage: " << prog
            << " <input.root> <run> [--threshold N] [--target-mpv M] [--rebin N] [--de-source all|pion] "
               "[--replace] [--ave] [--debug]\n"
               "       " << prog << " <helix.root> <run> --tree [--threshold N] [--ave] [--debug]\n"
               "  <run>        run number (auto-detected from the tpc tree or file name if omitted)\n"
               "  --threshold  minimum entries in a pad histogram to fit (default 20; --tree: 80)\n"
               "  --target-mpv dE MPV the gains are scaled to (default 200)\n"
               "  --rebin      rebin factor of the dE histograms before fitting (default 1)\n"
               "  --de-source  pion: pion-tagged dE only (default), all: every cluster\n"
               "  --tree       pion-beam Helix DST: pion sample from the tpc tree, dE per path normalized to\n"
               "               1 GeV/c pions, pads equalized within the inner / outer layers, inner layers matched\n"
               "               to the outer ones on the track dE/dx (see tpc_gain_tree.h)\n"
               "  --replace    write target/MPV as the gain instead of multiplying the current gain\n"
               "  --ave        fill unfitted pads from row neighbours / (layer, AsAd) median of fitted pads\n"
               "  --debug      do not write the parameter file"
            << std::endl;
}

int main(Int_t argc, char** argv) {
  if (argc < 2) {
    PrintUsage(argv[0]);
    return EXIT_FAILURE;
  }

  TString input_root = argv[1];
  Int_t run_number = -1;

  Double_t target_mpv = 200.0;
  Int_t threshold = 20;
  Int_t rebin = 1;
  TString de_source = "pion";
  Bool_t debug_mode = kFALSE;
  Bool_t replace_mode = kFALSE;
  Bool_t use_neighbor_mean_fill = kFALSE;  // --ave: fill unfitted pads (tpc_fill::FillUnfittedPads)
  Bool_t use_tree = kFALSE;
  Bool_t threshold_given = kFALSE;
  Bool_t hist_option_given = kFALSE;  // an option of the histogram mode only

  Int_t arg_idx = 2;
  if (arg_idx < argc && TString(argv[arg_idx]).IsDigit()) {
    run_number = std::atoi(argv[arg_idx++]);
  }

  for (Int_t i = arg_idx; i < argc; ++i) {
    const TString arg = argv[i];
    const Bool_t takes_value =
        (arg == "--target-mpv" || arg == "--threshold" || arg == "--rebin" || arg == "--de-source");
    if (takes_value && i + 1 >= argc) {
      std::cerr << "Error: option " << arg << " requires a value." << std::endl;
      PrintUsage(argv[0]);
      return EXIT_FAILURE;
    }
    if (arg == "--target-mpv") {
      target_mpv = std::atof(argv[++i]);
      hist_option_given = kTRUE;
    } else if (arg == "--threshold") {
      threshold = std::atoi(argv[++i]);
      threshold_given = kTRUE;
    } else if (arg == "--rebin") {
      rebin = std::max(1, std::atoi(argv[++i]));
      hist_option_given = kTRUE;
    } else if (arg == "--de-source") {
      de_source = argv[++i];
      hist_option_given = kTRUE;
    } else if (arg == "--debug") {
      debug_mode = kTRUE;
    } else if (arg == "--replace") {
      replace_mode = kTRUE;
      hist_option_given = kTRUE;
    } else if (arg == "--ave") {
      use_neighbor_mean_fill = kTRUE;
    } else if (arg == "--tree") {
      use_tree = kTRUE;
    } else {
      std::cerr << "Error: unknown option " << arg << std::endl;
      PrintUsage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  if (use_tree && hist_option_given) {
    std::cerr << "Error: --target-mpv, --rebin, --de-source and --replace do not apply to --tree." << std::endl;
    PrintUsage(argv[0]);
    return EXIT_FAILURE;
  }
  if (use_tree && !threshold_given) threshold = 80;

  TFile* f = TFile::Open(input_root.Data());
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << input_root << std::endl;
    return EXIT_FAILURE;
  }

  if (run_number < 0) {
    TTree* tree = (TTree*)f->Get("tpc");
    if (tree) {
      try {
        TTreeReader reader(tree);
        TTreeReaderValue<UInt_t> rv(reader, "run_number");
        if (reader.Next()) run_number = *rv;
      } catch (...) {
      }
    }
    if (run_number < 0) {
      TString filename = input_root;
      TPRegexp run_regex("run([0-9]+)");
      TString run_str = filename(run_regex);
      if (!run_str.IsNull()) {
        run_str.ReplaceAll("run", "");
        run_number = run_str.Atoi();
      }
    }
  }
  if (run_number < 0) {
    std::cerr << "Error: run number not provided and could not be auto-detected." << std::endl;
    f->Close();
    return EXIT_FAILURE;
  }

  const TString base_dir = ANALYZER_DIR + "/param/TPCPRM";
  const TString e72_dir = base_dir + "/e72";
  const TString run_param_path = Form("%s/TPCParam_e72_run%05d", e72_dir.Data(), run_number);
  const TString template_path = ANALYZER_DIR + "/param/TPCPRM/TPCParam_example";

  TString input_param_path = run_param_path;
  if (!fs::exists(run_param_path.Data())) {
    input_param_path = template_path;
    std::cout << "Input param file not found in e72, using template: " << template_path << std::endl;
  } else {
    std::cout << "Using existing param file for update: " << input_param_path << std::endl;
  }
  if (!fs::exists(input_param_path.Data())) {
    std::cerr << "Error: template/input param file does not exist: " << input_param_path << std::endl;
    f->Close();
    return EXIT_FAILURE;
  }

  std::vector<TpcParamEntry> entries;
  std::ifstream ifs(input_param_path.Data());
  if (!ifs.is_open()) {
    std::cerr << "Error: cannot open param file " << input_param_path << std::endl;
    f->Close();
    return EXIT_FAILURE;
  }

  std::string line;
  bool in_tpcphase_stamp_block = false;
  while (std::getline(ifs, line)) {
    if (line.find("# --- begin TPCPHASE pointer") != std::string::npos) {
      in_tpcphase_stamp_block = true;
      continue;
    }
    if (in_tpcphase_stamp_block) {
      if (line.find("# --- end TPCPHASE pointer ---") != std::string::npos) {
        in_tpcphase_stamp_block = false;
      }
      continue;
    }

    TpcParamEntry entry;
    entry.raw_line = line;
    std::string trimmed = line;
    size_t first = trimmed.find_first_not_of(" \t");
    if (first == std::string::npos) {
      entry.is_comment = kTRUE;
    } else {
      trimmed.erase(0, first);
      if (trimmed[0] == '#') {
        entry.is_comment = kTRUE;
      } else {
        std::stringstream ss(trimmed);
        if (!(ss >> entry.layer >> entry.row >> entry.aty)) {
          entry.is_comment = kTRUE;
        } else {
          entry.is_comment = kFALSE;
          Double_t val = 0.0;
          while (ss >> val) entry.p.push_back(val);
        }
      }
    }
    entries.push_back(entry);
  }
  ifs.close();

  std::map<std::pair<Int_t, Int_t>, Double_t> old_gain;
  for (const auto& e : entries) {
    if (!e.is_comment && e.aty == 0 && e.p.size() >= 2) {
      old_gain[{e.layer, e.row}] = e.p[1];
    }
  }

  // センターフレーム上のパッドは gain 校正対象外（デフォルト p0=0, gain=0）
  for (auto& e : entries) {
    if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
    if (!tpc::IsPadOnCenterFrame(e.layer, e.row)) continue;
    if (e.p.size() >= 1) e.p[0] = 0.0;
    e.p[1] = 0.0;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
  const TString pdf_path =
      use_tree ? img_dir + "/tpc_gain_calib_tree.pdf" : img_dir + Form("/tpc_gain_calib_%s.pdf", de_source.Data());

  GainFitConfig fit_cfg;
  fit_cfg.target_mpv = target_mpv;
  fit_cfg.threshold = threshold;
  fit_cfg.replace_mode = replace_mode;
  TCanvas* canvas = new TCanvas("c_gain", "c_gain", 1200, 1000);
  canvas->Print(pdf_path + "[");
  const Int_t nx = 4, ny = 4;
  Int_t ipad = 0;
  std::vector<TObject*> page_owned_draw;

  TH1D* hGainScale = new TH1D("hTPC_GainScale", "TPC gain scale (target/MPV);target/MPV;pads / bin", 120, 0.4, 1.8);
  TH1D* hMPV = new TH1D("hTPC_Gain_MPV", Form("TPC dE MPV (Run %05d);MPV;pads / bin", run_number), 120, 0, 1200);
  TH2Poly* hGainMap = new TH2Poly("hTPC_GainMap",
                                  Form("TPC gain map after calib (Run %05d);Z [mm];X [mm]", run_number), -300, 300,
                                  -300, 300);
  TH2Poly* hDeltaGainMap =
      new TH2Poly("hTPC_DeltaGainMap",
                  Form("TPC #Delta gain (new-old, fitted+ave-filled pads) (Run %05d);Z [mm];X [mm]", run_number), -300,
                  300, -300, 300);
  hGainScale->SetDirectory(nullptr);
  hMPV->SetDirectory(nullptr);
  hGainMap->SetDirectory(nullptr);
  hDeltaGainMap->SetDirectory(nullptr);
  tpc::InitializeHistograms(hGainMap);
  tpc::InitializeHistograms(hDeltaGainMap);

  Int_t n_updated = 0;
  Int_t n_filled_by_neighbor = 0;
  Int_t n_filled_by_asad = 0;
  Int_t n_total_pads = 0;
  std::set<std::pair<Int_t, Int_t>> pads_updated;
  std::set<std::pair<Int_t, Int_t>> pads_outlier_scale;  // scale forced to 1; not used as --ave sources
  std::set<std::pair<Int_t, Int_t>> pads_filled_by_neighbor;
  std::set<std::pair<Int_t, Int_t>> pads_filled_by_asad;
  std::map<Int_t, std::pair<Double_t, Int_t>> layer_mpv_stats;
  std::map<Int_t, std::pair<Double_t, Double_t>> layer_prefit_seed;

  if (use_tree) {
    tpc_gain_tree::Result tres;
    if (!tpc_gain_tree::Run(f, threshold, tres)) {
      f->Close();
      delete hGainScale;
      delete hMPV;
      delete hGainMap;
      delete hDeltaGainMap;
      delete canvas;
      return EXIT_FAILURE;
    }
    Int_t n_fit_inner = 0, n_fit_outer = 0;
    for (const auto& kv : tres.pad_mpv) (tpc_gain_tree::IsInner(kv.first.first) ? n_fit_inner : n_fit_outer)++;
    std::cout << Form("Tree sample: dE/dx scale %.1f, tracks %lld read / %lld in the pion window, %lld hits used",
                      tres.dedx_scale, tres.n_tracks_read, tres.n_tracks_pass, tres.n_hits_used)
              << std::endl;
    std::cout << Form("Pads fitted: inner %d (median MPV %.1f), outer %d (median MPV %.1f)", n_fit_inner,
                      tres.median_mpv_inner, n_fit_outer, tres.median_mpv_outer)
              << std::endl;
    std::cout << Form("Inner/outer truncated mean on %lld tracks: current gain %.3f, pad factors %.3f -> inner scale "
                      "%.4f -> %.3f",
                      tres.n_ratio_tracks, tres.ratio_before, tres.ratio_pads, tres.inner_scale, tres.ratio_after)
              << std::endl;

    for (auto& entry : entries) {
      if (entry.is_comment || entry.aty != 0 || entry.p.size() < 2) continue;
      if (tpc::IsPadOnCenterFrame(entry.layer, entry.row)) continue;
      n_total_pads++;
      const auto key = std::make_pair(entry.layer, entry.row);
      const Bool_t fitted = tres.pad_mpv.count(key) > 0;
      const Double_t factor = fitted ? tres.factor.at(key)
                                     : (tpc_gain_tree::IsInner(entry.layer) ? tres.inner_scale : 1.0);
      const Double_t old = entry.p[1];
      entry.p[1] = old * factor;
      if (fitted) {
        pads_updated.insert(key);
        hGainScale->Fill(factor);
        hMPV->Fill(tres.pad_mpv.at(key));
        std::cout << Form("Updated L=%2d R=%3d: MPV=%8.3f factor=%8.4f gain %10.6f -> %10.6f", entry.layer, entry.row,
                          tres.pad_mpv.at(key), factor, old, entry.p[1])
                  << std::endl;
      }
    }
    n_updated = static_cast<Int_t>(pads_updated.size());

    // QA: summary, pad MPV vs layer, inner/outer ratio, example pads
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    {
      TPaveText* pt = new TPaveText(0.02, 0.02, 0.98, 0.98, "NDC");
      pt->SetTextAlign(12);
      pt->SetFillColor(0);
      pt->AddText(Form("Run %05d, --tree (pion sample)", run_number));
      pt->AddText(Form("dE/dx scale / Bethe(pion): %.1f", tres.dedx_scale));
      pt->AddText(Form("tracks: %lld read, %lld in pion window", tres.n_tracks_read, tres.n_tracks_pass));
      pt->AddText(Form("hits used: %lld, threshold %d", tres.n_hits_used, threshold));
      pt->AddText(Form("fitted pads: inner %d, outer %d", n_fit_inner, n_fit_outer));
      pt->AddText(Form("median pad MPV: inner %.1f, outer %.1f", tres.median_mpv_inner, tres.median_mpv_outer));
      pt->AddText(Form("inner/outer trunc. mean: %.3f -> %.3f -> %.3f", tres.ratio_before, tres.ratio_pads,
                       tres.ratio_after));
      pt->AddText(Form("inner scale: %.4f", tres.inner_scale));
      pt->Draw();
      page_owned_draw.push_back(pt);
    }
    canvas->cd(2);
    TH2D* hMpvLayer = new TH2D("hGainTree_MpvLayer", "pad MPV (current gain);layer;MPV", 32, -0.5, 31.5, 100, 0, 400);
    hMpvLayer->SetDirectory(nullptr);
    for (const auto& kv : tres.pad_mpv) hMpvLayer->Fill(kv.first.first, kv.second);
    hMpvLayer->Draw("COLZ");
    page_owned_draw.push_back(hMpvLayer);
    canvas->cd(3);
    TH1D* hRb = new TH1D("hGainTree_RatioBefore", "inner/outer truncated mean (same track);ratio;tracks", 100, 0.4, 2.0);
    TH1D* hRa = new TH1D("hGainTree_RatioAfter", "", 100, 0.4, 2.0);
    hRb->SetDirectory(nullptr);
    hRa->SetDirectory(nullptr);
    for (Double_t r : tres.ratio_before_v) hRb->Fill(r);
    for (Double_t r : tres.ratio_after_v) hRa->Fill(r);
    hRa->SetLineColor(kRed);
    hRb->Draw("HIST");
    hRa->Draw("HIST SAME");
    page_owned_draw.push_back(hRb);
    page_owned_draw.push_back(hRa);
    canvas->cd(4);
    TH1D* hF = new TH1D("hGainTree_Factor", "factor of fitted pads;new/old gain;pads", 100, 0.4, 1.6);
    hF->SetDirectory(nullptr);
    for (const auto& kv : tres.pad_mpv) hF->Fill(tres.factor.at(kv.first));
    hF->Draw("HIST");
    page_owned_draw.push_back(hF);
    canvas->Print(pdf_path);
    for (TObject* obj : page_owned_draw) delete obj;
    page_owned_draw.clear();

    canvas->Clear();
    canvas->Divide(nx, ny);
    Int_t ipanel = 0;
    for (Int_t layer = 0; layer < tpc::NumOfLayersTPC && ipanel < nx * ny; layer += 2) {
      // the fitted pad with the most hits in this layer
      const TH1D* best = nullptr;
      Double_t best_mpv = 0.0;
      for (const auto& kv : tres.pad_mpv) {
        if (kv.first.first != layer) continue;
        const TH1D* h = tres.pad_hist.at(kv.first);
        if (!best || h->GetEntries() > best->GetEntries()) {
          best = h;
          best_mpv = kv.second;
        }
      }
      if (!best) continue;
      canvas->cd(++ipanel);
      TH1D* hc = (TH1D*)best->Clone(Form("%s_draw", best->GetName()));
      hc->SetDirectory(nullptr);
      hc->Draw("HIST");
      gPad->Update();
      TLine* v = new TLine(best_mpv, gPad->GetUymin(), best_mpv, gPad->GetUymax());
      v->SetLineColor(kRed);
      v->Draw();
      page_owned_draw.push_back(hc);
      page_owned_draw.push_back(v);
    }
    canvas->Print(pdf_path);
    for (TObject* obj : page_owned_draw) delete obj;
    page_owned_draw.clear();
  } else {
  // Layer 合算 dE ヒストを先に fit して、pad fit の初期値 seed を作る。
  // Width は合算ヒストの方が太く出やすいので少し縮めて使う。
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    TH1D* h_layer = nullptr;
    {
      const TString hname_layer =
          (de_source == "pion")
              ? Form("TPCCl_dE_Pion_Layer%02d", layer)
              : Form("TPCCl_dE_Layer%02d", layer);
      f->GetObject(hname_layer.Data(), h_layer);
      if (!h_layer) {
        const TString hname_hist = "hist/" + hname_layer;
        f->GetObject(hname_hist.Data(), h_layer);
      }
    }
    if (!h_layer || h_layer->GetEntries() <= 0) continue;

    TH1D* h_layer_fit = h_layer;
    TH1D* h_layer_rebinned = nullptr;
    if (rebin > 1) {
      h_layer_rebinned = (TH1D*)h_layer->Clone(Form("hTPCCl_dE_layer_rebin_L%02d", layer));
      h_layer_rebinned->SetDirectory(nullptr);
      h_layer_rebinned->Rebin(rebin);
      h_layer_fit = h_layer_rebinned;
    }

    const Int_t peak_bin = SelectGainPeakBin(h_layer_fit, fit_cfg);
    const Double_t peak_x = h_layer_fit->GetBinCenter(peak_bin);
    const Double_t peak_y = h_layer_fit->GetBinContent(peak_bin);
    Double_t mpv_init_layer = peak_x;
    Double_t width_init_layer = fit_cfg.min_width;
    Double_t fit_min = 0.0;
    Double_t fit_max = 0.0;
    if (!ComputeGainFitWindow(h_layer_fit, fit_cfg, peak_bin, -1.0, mpv_init_layer, width_init_layer, fit_min,
                              fit_max)) {
      if (h_layer_rebinned) delete h_layer_rebinned;
      continue;
    }

    TF1* fLayer = new TF1(Form("fLayerLandau_L%02d", layer), "landau(0)+pol0(3)", fit_min, fit_max);
    const Double_t span = fit_max - fit_min;
    width_init_layer = std::max(width_init_layer, std::max(1.0, 0.18 * span));
    fLayer->SetParameters(std::max(peak_y, 1.0), mpv_init_layer, width_init_layer, 0.0);
    fLayer->SetParLimits(1, fit_min, fit_max);
    fLayer->SetParLimits(2, 0.35, std::max(1.5, 0.45 * span));
    const Int_t st = h_layer_fit->Fit(fLayer, "Q");
    if (st == 0) {
      const Double_t mpv_seed = fLayer->GetParameter(1);
      const Double_t width_seed = std::max(0.35, 0.70 * fLayer->GetParameter(2));
      if (std::isfinite(mpv_seed) && mpv_seed > 0.0 && std::isfinite(width_seed) && width_seed > 0.0) {
        layer_prefit_seed[layer] = {mpv_seed, width_seed};
      }
    }
    delete fLayer;
    if (h_layer_rebinned) delete h_layer_rebinned;
  }

  for (auto& entry : entries) {
    if (entry.is_comment || entry.aty != 0 || entry.p.size() < 2) continue;
    if (tpc::IsPadOnCenterFrame(entry.layer, entry.row)) continue;
    n_total_pads++;

    TH1D* h = nullptr;
    {
      const TString hname =
          (de_source == "pion")
              ? Form("TPCCl_dE_Pion_Layer%02d_Row%03d", entry.layer, entry.row)
              : Form("TPCCl_dE_Layer%02d_Row%03d", entry.layer, entry.row);
      f->GetObject(hname.Data(), h);
      if (!h) {
        const TString hname_hist = "hist/" + hname;
        f->GetObject(hname_hist.Data(), h);
      }
    }
    if (!h || h->GetEntries() <= 0) continue;

    TH1D* h_fit = h;
    TH1D* h_rebinned = nullptr;
    bool keep_h_rebinned_for_page = false;
    if (rebin > 1) {
      h_rebinned = (TH1D*)h->Clone(Form("hTPCCl_dE_rebin_L%02d_R%03d", entry.layer, entry.row));
      h_rebinned->SetDirectory(nullptr);
      h_rebinned->Rebin(rebin);
      h_fit = h_rebinned;
    }

    const Int_t bmin_stat = h_fit->FindBin(-fit_cfg.stat_range);
    const Int_t bmax_stat = h_fit->FindBin(fit_cfg.stat_range);
    const Double_t n_stat = h_fit->Integral(bmin_stat, bmax_stat);
    if (fit_cfg.threshold > 0 && n_stat < fit_cfg.threshold) {
      if (h_rebinned) delete h_rebinned;
      continue;
    }

    const Int_t peak_bin = SelectGainPeakBin(h_fit, fit_cfg);
    const Double_t peak_x = h_fit->GetBinCenter(peak_bin);
    const Double_t peak_y = h_fit->GetBinContent(peak_bin);
    const Int_t n_bins_window = bmax_stat - bmin_stat + 1;
    const Double_t mean_content = (n_bins_window > 0) ? (n_stat / n_bins_window) : 0.0;
    if (mean_content <= 0.0 || peak_y < mean_content * fit_cfg.min_peak_ratio) continue;

    Double_t mpv_hint = -1.0;
    Double_t width_init = fit_cfg.min_width;
    auto it_prefit = layer_prefit_seed.find(entry.layer);
    if (it_prefit != layer_prefit_seed.end()) {
      mpv_hint = it_prefit->second.first;
      width_init = std::max(fit_cfg.min_width, it_prefit->second.second);
      const Double_t xmin = h_fit->GetXaxis()->GetXmin();
      const Double_t xmax = h_fit->GetXaxis()->GetXmax();
      if (fit_cfg.mpv_min > 0.0) mpv_hint = std::max(mpv_hint, fit_cfg.mpv_min);
      if (fit_cfg.mpv_max > 0.0) mpv_hint = std::min(mpv_hint, fit_cfg.mpv_max);
      mpv_hint = std::min(std::max(mpv_hint, xmin), xmax);
    }
    Double_t mpv_init_main = peak_x;
    Double_t fit_min = 0.0;
    Double_t fit_max = 0.0;
    if (!ComputeGainFitWindow(h_fit, fit_cfg, peak_bin, mpv_hint, mpv_init_main, width_init, fit_min, fit_max)) {
      if (h_rebinned) delete h_rebinned;
      continue;
    }
    if (fit_max <= fit_min) {
      if (h_rebinned) delete h_rebinned;
      continue;
    }

    TF1* fFit = new TF1(Form("fLandau_L%02d_R%03d", entry.layer, entry.row), "landau(0)+pol0(3)", fit_min, fit_max);
    fFit->SetParNames("Norm", "MPV", "Width", "Background");
    fFit->SetNpx(1000);
    const Double_t fit_span = fit_max - fit_min;
    width_init = std::min(std::max(width_init, 0.08 * fit_span), 0.30 * fit_span);

    auto fit_once = [&](Double_t mpv_init) {
      fFit->SetRange(fit_min, fit_max);
      fFit->SetParameters(std::max(peak_y, 1.0), mpv_init, width_init, std::max(0.0, mean_content * 0.2));
      fFit->SetParLimits(0, 0.0, peak_y * 100.0 + 1.0);
      fFit->SetParLimits(1, fit_min, fit_max);
      // 上限を抑えて「広すぎる解」へ逃げるのを防ぐ（必要なら min-width で最終採否）。
      fFit->SetParLimits(2, 0.35, std::max(1.0, 0.35 * fit_span));
      fFit->SetParLimits(3, 0.0, peak_y * 10.0 + 1.0);
      return h_fit->Fit(fFit, "Q0");
    };

    Int_t fit_status = fit_once(mpv_init_main);
    if (fit_status == 0) {
      Double_t mpv = fFit->GetParameter(1);
      Double_t width = fFit->GetParameter(2);
      Int_t ndf = fFit->GetNDF();
      Double_t chi2 = fFit->GetChisquare();
      Double_t chi2ndf = (ndf > 0) ? (chi2 / static_cast<Double_t>(ndf)) : 0.0;
      const TString fit_pad_name = Form("fLandau_L%02d_R%03d", entry.layer, entry.row);
      ApplyGainRefitPass(h_fit, fFit, fit_pad_name.Data(), fit_cfg, peak_bin, peak_x, peak_y, mean_content,
                         fit_min, fit_max, width_init, mpv, width, ndf, chi2ndf);
      chi2 = fFit->GetChisquare();
      chi2ndf = (ndf > 0) ? (chi2 / static_cast<Double_t>(ndf)) : 0.0;
      Double_t fit_prob = (ndf > 0) ? TMath::Prob(chi2, ndf) : 0.0;

      Bool_t quality_ok = kTRUE;
      if (!(mpv > 0.0) || !std::isfinite(mpv)) quality_ok = kFALSE;
      if (fit_cfg.mpv_min > 0.0 && mpv < fit_cfg.mpv_min) quality_ok = kFALSE;
      if (fit_cfg.mpv_max > 0.0 && mpv > fit_cfg.mpv_max) quality_ok = kFALSE;
      if (width < fit_cfg.min_width) quality_ok = kFALSE;
      if (fit_cfg.min_fit_ndf > 0 && ndf < fit_cfg.min_fit_ndf) quality_ok = kFALSE;
      if (fit_cfg.max_chi2_ndf > 0.0 && ndf > 0 && chi2ndf > fit_cfg.max_chi2_ndf) quality_ok = kFALSE;
      if (fit_cfg.min_fit_prob > 0.0 && ndf > 0 && fit_prob < fit_cfg.min_fit_prob) quality_ok = kFALSE;

      if (quality_ok) {
        Double_t scale = fit_cfg.target_mpv / mpv;
        Bool_t outlier_scale = (!std::isfinite(scale) || scale < (1.0 / 3.0) || scale > 10.0);

        // 外れ値になりそうなら、同 layer の平均 MPV を初期値にして 1 回だけ再 fit。
        if (outlier_scale) {
          auto it = layer_mpv_stats.find(entry.layer);
          if (it != layer_mpv_stats.end() && it->second.second > 0) {
            const Double_t layer_mpv_mean = it->second.first / static_cast<Double_t>(it->second.second);
            if (fit_once(layer_mpv_mean) == 0) {
              mpv = fFit->GetParameter(1);
              width = fFit->GetParameter(2);
              ndf = fFit->GetNDF();
              chi2 = fFit->GetChisquare();
              chi2ndf = (ndf > 0) ? (chi2 / static_cast<Double_t>(ndf)) : 0.0;
              fit_prob = (ndf > 0) ? TMath::Prob(chi2, ndf) : 0.0;
              scale = fit_cfg.target_mpv / mpv;
              outlier_scale = (!std::isfinite(scale) || scale < (1.0 / 3.0) || scale > 10.0);
            }
          }
        }
        if (outlier_scale) scale = 1.0;

        const Double_t old = entry.p[1];
        entry.p[1] = replace_mode ? scale : (old * scale);
        n_updated++;
        const auto pad_key = std::make_pair(entry.layer, entry.row);
        pads_updated.insert(pad_key);
        if (outlier_scale) pads_outlier_scale.insert(pad_key);
        hGainScale->Fill(scale);
        hMPV->Fill(mpv);
        if (std::isfinite(mpv) && mpv > 0.0) {
          auto& st = layer_mpv_stats[entry.layer];
          st.first += mpv;
          st.second += 1;
        }
        std::cout << Form("Updated L=%2d R=%3d: MPV=%8.3f width=%6.3f chi2/ndf=%.3f scale=%8.4f (%s) gain %10.6f -> %10.6f",
                          entry.layer, entry.row, mpv, width, chi2ndf, scale,
                          replace_mode ? "replace" : "multiply", old, entry.p[1])
                  << std::endl;

        if (ipad == 0) {
          canvas->Clear();
          canvas->Divide(nx, ny);
        }
        canvas->cd(++ipad);
        if (h_rebinned) {
          // Drawn objects must stay alive until canvas->Print() flushes this page.
          page_owned_draw.push_back(h_rebinned);
          keep_h_rebinned_for_page = true;
        }
        DrawGainCalibPanel(h_fit, fFit, kTRUE, fit_min, fit_max, mpv, page_owned_draw);
        if (ipad >= nx * ny) {
          canvas->Print(pdf_path);
          for (TObject* obj : page_owned_draw) delete obj;
          page_owned_draw.clear();
          ipad = 0;
        }
      }
    }
    delete fFit;
    if (h_rebinned && !keep_h_rebinned_for_page) delete h_rebinned;
  }

  if (ipad > 0) {
    canvas->Print(pdf_path);
    for (TObject* obj : page_owned_draw) delete obj;
    page_owned_draw.clear();
  }

  }  // end histogram mode

  // --ave: fill pads not updated by the fit from fitted pads (tpc_fill::FillUnfittedPads): mean of row+-1 in
  // the same AsAd, otherwise median of the (layer, AsAd) group. Sources are fitted pads whose scale was not
  // forced and whose gain is within 1/2 to 2 times the group median.
  Int_t n_left_unfilled = 0;
  if (use_neighbor_mean_fill) {
    std::map<tpc_fill::PadKey, Double_t> fill_values;
    for (const auto& e : entries) {
      if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
      fill_values[{e.layer, e.row}] = e.p[1];
    }
    auto gain_consistent = [](Double_t v, const tpc_fill::GroupStat& g) -> Bool_t {
      return v > 0.0 && v >= 0.5 * g.median && v <= 2.0 * g.median;
    };
    const tpc_fill::FillResult fill =
        tpc_fill::FillUnfittedPads(fill_values, pads_updated, pads_outlier_scale, gain_consistent);

    for (auto& e : entries) {
      if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
      const auto key = std::make_pair(e.layer, e.row);
      const auto it = fill.value.find(key);
      if (it == fill.value.end()) continue;
      const Double_t old_g = e.p[1];
      e.p[1] = it->second;
      if (fill.kind.at(key) == tpc_fill::FillKind::kNeighbor) {
        pads_filled_by_neighbor.insert(key);
        std::cout << Form("  neighbor-mean fill (gain): L=%2d R=%3d gain: %10.6f -> %10.6f", e.layer, e.row,
                          old_g, e.p[1])
                  << std::endl;
      } else {
        pads_filled_by_asad.insert(key);
        std::cout << Form("  asad-median fill (gain): L=%2d R=%3d (ASAD=%2d) gain: %10.6f -> %10.6f", e.layer,
                          e.row, tpc::GetASADId(e.layer, e.row), old_g, e.p[1])
                  << std::endl;
      }
    }
    n_filled_by_neighbor = fill.n_neighbor;
    n_filled_by_asad = fill.n_asad_median;
    n_left_unfilled = fill.n_left;
  }

  Double_t max_abs_delta_gain = 0.0;
  for (const auto& e : entries) {
    if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
    const auto key = std::make_pair(e.layer, e.row);
    const auto it_old = old_gain.find(key);
    if (it_old == old_gain.end()) continue;
    const Double_t old = it_old->second;
    const Double_t now = e.p[1];
    if (!std::isfinite(now)) continue;
    const Int_t bin_idx = tpc::GetPadId(e.layer, e.row) + 1;
    hGainMap->SetBinContent(bin_idx, now);
    if (pads_updated.count(key) || pads_filled_by_neighbor.count(key) || pads_filled_by_asad.count(key)) {
      const Double_t delta = now - old;
      hDeltaGainMap->SetBinContent(bin_idx, delta);
      max_abs_delta_gain = std::max(max_abs_delta_gain, std::abs(delta));
    }
  }
  if (max_abs_delta_gain > 0.0) {
    hDeltaGainMap->SetMinimum(-max_abs_delta_gain);
    hDeltaGainMap->SetMaximum(+max_abs_delta_gain);
  }

  // 外れ gain でスケールが潰れないよう、pad map の z 上限を固定（データ自体はそのまま）
  constexpr Double_t kGainMapZMax = 3.0;
  constexpr Double_t kGainMapZMin = 0.0;
  hGainMap->SetMinimum(kGainMapZMin);
  hGainMap->SetMaximum(kGainMapZMax);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hMPV->Draw("HIST");
  canvas->cd(2);
  hGainScale->Draw("HIST");
  canvas->cd(3);
  hGainMap->Draw("COLZ");
  canvas->cd(4);
  hDeltaGainMap->Draw("COLZ");
  canvas->Print(pdf_path);
  canvas->Print(pdf_path + "]");

  if (!debug_mode) {
    if (!fs::exists(e72_dir.Data())) fs::create_directories(e72_dir.Data());
    if (fs::exists(run_param_path.Data())) {
      fs::copy(run_param_path.Data(), (run_param_path + ".bak").Data(), fs::copy_options::overwrite_existing);
    }
    if (!fs::exists(run_param_path.Data())) {
      fs::copy(input_param_path.Data(), run_param_path.Data(), fs::copy_options::overwrite_existing);
    }

    std::ofstream ofs(run_param_path.Data());
    if (!ofs.is_open()) {
      std::cerr << "Error: cannot open output file " << run_param_path << std::endl;
      if (f) f->Close();
      delete hGainScale;
      delete hMPV;
      delete canvas;
      return EXIT_FAILURE;
    }

    for (const auto& entry : entries) {
      if (entry.is_comment) {
        ofs << entry.raw_line.Data() << "\n";
      } else {
        ofs << entry.layer << "\t" << entry.row << "\t" << entry.aty;
        for (size_t i = 0; i < entry.p.size(); ++i) {
          ofs << "\t" << std::fixed << std::setprecision(8) << entry.p[i];
        }
        ofs << "\n";
      }
    }
    ofs.close();
  }

  std::cout << "Done! Updated " << n_updated << " out of " << n_total_pads << " aty==0 pads." << std::endl;
  if (use_neighbor_mean_fill) {
    std::cout << "Unfitted pads filled (gain): neighbor-mean " << n_filled_by_neighbor << ", asad-median "
              << n_filled_by_asad << ", left unfilled " << n_left_unfilled << std::endl;
  }
  if (debug_mode) {
    std::cout << "Debug mode: parameter file was NOT updated." << std::endl;
  } else {
    std::cout << "Result saved to: " << run_param_path << std::endl;
  }
  std::cout << "QA PDF saved to: " << pdf_path << std::endl;

  if (f) f->Close();
  delete hGainScale;
  delete hMPV;
  delete hGainMap;
  delete hDeltaGainMap;
  delete canvas;
  return EXIT_SUCCESS;
}
