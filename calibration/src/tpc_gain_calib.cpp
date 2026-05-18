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
#include <TH2Poly.h>
#include <TLine.h>
#include <TMath.h>
#include <TObject.h>
#include <TPad.h>
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

namespace fs = std::filesystem;

struct TpcParamEntry {
  Int_t layer{};
  Int_t row{};
  Int_t aty{};
  std::vector<Double_t> p;
  TString raw_line;
  Bool_t is_comment{kTRUE};
};

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

static Bool_t FillTreeHists(TFile* f, const TString& tree_name, Int_t clsize_min, Int_t clsize_max,
                            Double_t min_abs_cos_theta, Bool_t pion_tracks_only, Bool_t fill_per_pad,
                            std::vector<TH1D*>& h_layers,
                            std::map<std::pair<Int_t, Int_t>, TH1D*>& h_pads, Long64_t& n_hits_filled) {
  h_layers.clear();
  h_layers.resize(tpc::NumOfLayersTPC, nullptr);
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    h_layers[layer] = new TH1D(Form("hTPCCl_dE_tree_L%02d", layer),
                               Form("Cluster dE (tree) Layer %02d;Cluster dE;Counts", layer), 400, 0., 800.);
    h_layers[layer]->SetDirectory(nullptr);
  }

  TTree* tree = (TTree*)f->Get(tree_name.Data());
  if (!tree) {
    std::cerr << "Error: tree '" << tree_name << "' not found in " << f->GetName() << std::endl;
    return kFALSE;
  }

  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> rv_ntTpc(reader, "ntTpc");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_clsize(reader, "track_cluster_size");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_theta_diff(reader, "theta_diff");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_clde(reader, "track_cluster_de");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> rv_hitlayer(reader, "hitlayer");
  std::unique_ptr<TTreeReaderValue<std::vector<std::vector<Double_t>>>> rv_hitrow;
  if (fill_per_pad) {
    rv_hitrow =
        std::make_unique<TTreeReaderValue<std::vector<std::vector<Double_t>>>>(reader, "track_cluster_row_center");
  }

  std::unique_ptr<TTreeReaderValue<std::vector<Int_t>>> rv_pid;
  if (pion_tracks_only) {
    rv_pid = std::make_unique<TTreeReaderValue<std::vector<Int_t>>>(reader, "pid");
  }

  const Long64_t n_entries = tree->GetEntries();
  if (n_entries <= 0) {
    std::cerr << "Error: tree '" << tree_name << "' has no entries." << std::endl;
    return kFALSE;
  }

  n_hits_filled = 0;
  Long64_t iev = 0;
  while (reader.Next()) {
    displayProgressBar(static_cast<Int_t>(++iev), static_cast<Int_t>(n_entries));
    const Int_t nt = *rv_ntTpc;
    const auto& clde_ev = *rv_clde;
    const auto& clsize_ev = *rv_clsize;
    const auto& theta_ev = *rv_theta_diff;
    const auto& hitlayer_ev = *rv_hitlayer;
    const auto* hitrow_ev = fill_per_pad ? &(**rv_hitrow) : nullptr;
    for (Int_t it = 0; it < nt; ++it) {
      if (pion_tracks_only && rv_pid) {
        const auto& pid_ev = **rv_pid;
        if (it >= static_cast<Int_t>(pid_ev.size()) || !(pid_ev[it] & 0x1)) continue;
      }
      if (it >= static_cast<Int_t>(clde_ev.size())) continue;
      const auto& v_de = clde_ev[it];
      const auto& v_cs = clsize_ev[it];
      const auto& v_td = theta_ev[it];
      const auto& v_ly = hitlayer_ev[it];
      const auto* v_rw = (fill_per_pad && it < static_cast<Int_t>(hitrow_ev->size())) ? &(*hitrow_ev)[it] : nullptr;
      const Int_t nh = static_cast<Int_t>(v_de.size());
      for (Int_t ih = 0; ih < nh; ++ih) {
        if (ih >= static_cast<Int_t>(v_cs.size()) || ih >= static_cast<Int_t>(v_td.size()) ||
            ih >= static_cast<Int_t>(v_ly.size()))
          continue;
        const Int_t cs = static_cast<Int_t>(v_cs[ih]);
        if (cs < clsize_min || cs > clsize_max) continue;
        if (TMath::Abs(TMath::Cos(v_td[ih])) < min_abs_cos_theta) continue;
        const Int_t layer = static_cast<Int_t>(v_ly[ih]);
        if (layer < 0 || layer >= tpc::NumOfLayersTPC) continue;
        h_layers[layer]->Fill(v_de[ih]);
        if (fill_per_pad && v_rw) {
          if (ih >= static_cast<Int_t>(v_rw->size())) continue;
          const Int_t row = static_cast<Int_t>((*v_rw)[ih]);
          const Int_t n_row = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
          if (row < 0 || row >= n_row) continue;
          const auto pad_key = std::make_pair(layer, row);
          auto it_pad = h_pads.find(pad_key);
          if (it_pad == h_pads.end()) {
            TH1D* h_pad = new TH1D(Form("hTPCCl_dE_tree_L%02d_R%03d", layer, row),
                                   Form("Cluster dE (tree) L%02d R%03d;Cluster dE;Counts", layer, row), 400, 0.,
                                   800.);
            h_pad->SetDirectory(nullptr);
            h_pads[pad_key] = h_pad;
          }
          h_pads[pad_key]->Fill(v_de[ih]);
        }
        ++n_hits_filled;
      }
    }
  }

  std::cout << "Tree fill: " << n_hits_filled << " hits passed cuts (" << iev << " events) into "
            << (fill_per_pad ? "pad+layer" : "layer") << " histograms." << std::endl;
  if (fill_per_pad) {
    std::cout << "  " << h_pads.size() << " pads with at least one hit." << std::endl;
  }
  return kTRUE;
}

int main(Int_t argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <input.root> [run_number] [--target-mpv M] [--threshold N] [--stat-range R] "
                 "[--min-peak-ratio X] [--fit-nsigma N] [--local-peak-half W] [--local-peak-sep B] "
                 "[--min-width W] "
                 "[--max-chi2-ndf X] [--min-fit-ndf N] [--min-fit-prob P] [--rebin N] "
                 "[--mpv-min X] [--mpv-max X] [--de-source all|pion] [--tree] [--layer] "
                 "[--tree-name NAME] [--clsize-min N] [--clsize-max N] [--min-abs-cos-theta X] "
                 "[--ave] [--debug] [--replace] [--mode hit|trk]"
              << std::endl;
    return EXIT_FAILURE;
  }

  TString input_root = argv[1];
  Int_t run_number = -1;

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
  Int_t rebin = 1;
  Double_t mpv_min = 100.0;
  Double_t mpv_max = 0.0;  // <=0 means disabled
  TString de_source = "all";
  Bool_t debug_mode = kFALSE;
  Bool_t replace_mode = kFALSE;
  Bool_t use_asad_mean_fill = kFALSE;
  Bool_t use_tree = kFALSE;
  Bool_t layer_fit = kFALSE;
  TString tree_name = "tpc";
  Int_t clsize_min = 1;
  Int_t clsize_max = 1;
  Double_t min_abs_cos_theta = 0.95;
  TString mode = "hit";

  Int_t arg_idx = 2;
  if (arg_idx < argc && TString(argv[arg_idx]).IsDigit()) {
    run_number = std::atoi(argv[arg_idx++]);
  }

  for (Int_t i = arg_idx; i < argc; ++i) {
    TString arg = argv[i];
    if (arg == "--target-mpv" && i + 1 < argc) {
      target_mpv = std::atof(argv[++i]);
    } else if (arg == "--threshold" && i + 1 < argc) {
      threshold = std::atoi(argv[++i]);
    } else if (arg == "--stat-range" && i + 1 < argc) {
      stat_range = std::atof(argv[++i]);
    } else if (arg == "--min-peak-ratio" && i + 1 < argc) {
      min_peak_ratio = std::atof(argv[++i]);
    } else if (arg == "--fit-nsigma" && i + 1 < argc) {
      fit_nsigma = std::atof(argv[++i]);
    } else if (arg == "--local-peak-half" && i + 1 < argc) {
      local_peak_half = std::atof(argv[++i]);
    } else if (arg == "--local-peak-sep" && i + 1 < argc) {
      local_peak_sep = std::max(0, std::atoi(argv[++i]));
    } else if (arg == "--min-width" && i + 1 < argc) {
      min_width = std::atof(argv[++i]);
    } else if (arg == "--max-chi2-ndf" && i + 1 < argc) {
      max_chi2_ndf = std::atof(argv[++i]);
    } else if (arg == "--min-fit-ndf" && i + 1 < argc) {
      min_fit_ndf = std::atoi(argv[++i]);
    } else if (arg == "--min-fit-prob" && i + 1 < argc) {
      min_fit_prob = std::atof(argv[++i]);
    } else if (arg == "--rebin" && i + 1 < argc) {
      rebin = std::max(1, std::atoi(argv[++i]));
    } else if (arg == "--mpv-min" && i + 1 < argc) {
      mpv_min = std::atof(argv[++i]);
    } else if (arg == "--mpv-max" && i + 1 < argc) {
      mpv_max = std::atof(argv[++i]);
    } else if (arg == "--de-source" && i + 1 < argc) {
      de_source = argv[++i];
    } else if (arg == "--debug") {
      debug_mode = kTRUE;
    } else if (arg == "--replace") {
      replace_mode = kTRUE;
    } else if (arg == "--mode" && i + 1 < argc) {
      mode = argv[++i];
    } else if (arg == "--ave" || arg == "--asad-avg") {
      use_asad_mean_fill = kTRUE;
    } else if (arg == "--tree") {
      use_tree = kTRUE;
    } else if (arg == "--layer") {
      layer_fit = kTRUE;
    } else if (arg == "--tree-name" && i + 1 < argc) {
      tree_name = argv[++i];
    } else if (arg == "--clsize-min" && i + 1 < argc) {
      clsize_min = std::atoi(argv[++i]);
    } else if (arg == "--clsize-max" && i + 1 < argc) {
      clsize_max = std::atoi(argv[++i]);
    } else if (arg == "--min-abs-cos-theta" && i + 1 < argc) {
      min_abs_cos_theta = std::atof(argv[++i]);
    }
  }

  const Bool_t tree_layer_mode = use_tree && layer_fit;
  const Bool_t tree_pad_mode = use_tree && !layer_fit;

  (void)mode;

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
  const TString template_path = ANALYZER_DIR + "/param/TPCPRM/TPCParam_e72_run02594";

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

  // センターフレーム上のパッドは gain 校正対象外（デフォルト p0=0, gain=1）
  for (auto& e : entries) {
    if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
    if (!tpc::IsPadOnCenterFrame(e.layer, e.row)) continue;
    if (e.p.size() >= 1) e.p[0] = 0.0;
    e.p[1] = 1.0;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  const TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
  TString pdf_path = img_dir + Form("/tpc_gain_calib_%s.pdf", de_source.Data());
  if (tree_layer_mode) {
    pdf_path = img_dir + Form("/tpc_gain_calib_%s_tree_layer.pdf", de_source.Data());
  } else if (tree_pad_mode) {
    pdf_path = img_dir + Form("/tpc_gain_calib_%s_tree.pdf", de_source.Data());
  }

  GainFitConfig fit_cfg;
  fit_cfg.target_mpv = target_mpv;
  fit_cfg.threshold = threshold;
  fit_cfg.stat_range = stat_range;
  fit_cfg.min_peak_ratio = min_peak_ratio;
  fit_cfg.fit_nsigma = fit_nsigma;
  fit_cfg.local_peak_half = local_peak_half;
  fit_cfg.local_peak_sep = local_peak_sep;
  fit_cfg.min_width = min_width;
  fit_cfg.max_chi2_ndf = max_chi2_ndf;
  fit_cfg.min_fit_ndf = min_fit_ndf;
  fit_cfg.min_fit_prob = min_fit_prob;
  fit_cfg.mpv_min = mpv_min;
  fit_cfg.mpv_max = mpv_max;
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
  Int_t n_filled_by_asad = 0;
  Int_t n_total_pads = 0;
  std::set<std::pair<Int_t, Int_t>> pads_updated;
  std::set<std::pair<Int_t, Int_t>> pads_filled_by_asad;
  std::map<Int_t, std::pair<Double_t, Int_t>> layer_mpv_stats;
  std::map<Int_t, std::pair<Double_t, Double_t>> layer_prefit_seed;
  std::vector<TH1D*> h_tree_layers;
  std::map<std::pair<Int_t, Int_t>, TH1D*> h_tree_pads;

  if (tree_layer_mode) {
    const Bool_t pion_tracks_only = (de_source == "pion");
    std::vector<TH1D*> h_layers;
    std::map<std::pair<Int_t, Int_t>, TH1D*> h_pads_unused;
    Long64_t n_hits_filled = 0;
    if (!FillTreeHists(f, tree_name, clsize_min, clsize_max, min_abs_cos_theta, pion_tracks_only, kFALSE, h_layers,
                       h_pads_unused, n_hits_filled)) {
      f->Close();
      delete hGainScale;
      delete hMPV;
      delete hGainMap;
      delete hDeltaGainMap;
      delete canvas;
      return EXIT_FAILURE;
    }

    for (const auto& e : entries) {
      if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
      if (tpc::IsPadOnCenterFrame(e.layer, e.row)) continue;
      n_total_pads++;
    }

    for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
      TH1D* h_layer = h_layers[layer];
      if (!h_layer || h_layer->GetEntries() <= 0) {
        std::cout << Form("Layer %02d (tree): no entries after cuts (hist skipped)", layer) << std::endl;
        continue;
      }

      TH1D* h_layer_fit = h_layer;
      TH1D* h_layer_rebinned = nullptr;
      if (rebin > 1) {
        h_layer_rebinned = (TH1D*)h_layer->Clone(Form("hTPCCl_dE_tree_rebin_L%02d", layer));
        h_layer_rebinned->SetDirectory(nullptr);
        h_layer_rebinned->Rebin(rebin);
        h_layer_fit = h_layer_rebinned;
      }

      GainFitResult fres;
      TF1* fFit = nullptr;
      const Bool_t fit_ok =
          FitLandauGain(h_layer_fit, layer, -1, fit_cfg, fres, &fFit, kFALSE /* tree: loose quality */);

      if (fit_ok) {
        std::cout << Form("Layer %02d (tree): MPV=%8.3f width=%6.3f chi2/ndf=%.3f scale=%8.4f n=%.0f", layer,
                          fres.mpv, fres.width, fres.chi2ndf, fres.scale, h_layer_fit->GetEntries())
                  << std::endl;
        hGainScale->Fill(fres.scale);
        hMPV->Fill(fres.mpv);
        for (auto& entry : entries) {
          if (entry.is_comment || entry.aty != 0 || entry.p.size() < 2) continue;
          if (entry.layer != layer) continue;
          if (tpc::IsPadOnCenterFrame(entry.layer, entry.row)) continue;
          const Double_t old = entry.p[1];
          entry.p[1] = replace_mode ? fres.scale : (old * fres.scale);
          pads_updated.insert({entry.layer, entry.row});
        }
      } else {
        std::cout << Form("Layer %02d (tree): fit failed (hist only) n=%.0f", layer, h_layer_fit->GetEntries())
                  << std::endl;
      }

      if (ipad == 0) {
        canvas->Clear();
        canvas->Divide(nx, ny);
      }
      canvas->cd(++ipad);
      h_layer_fit->SetTitle(Form("Layer %02d  n=%.0f%s", layer, h_layer_fit->GetEntries(),
                                 fit_ok ? "" : "  (fit failed)"));
      if (h_layer_rebinned) {
        page_owned_draw.push_back(h_layer_rebinned);
      }
      DrawGainCalibPanel(h_layer_fit, fFit, fit_ok, fres.fit_min, fres.fit_max, fres.mpv, page_owned_draw);
      if (fFit) delete fFit;
      if (ipad >= nx * ny) {
        canvas->Print(pdf_path);
        for (TObject* obj : page_owned_draw) delete obj;
        page_owned_draw.clear();
        ipad = 0;
      }
    }

    n_updated = static_cast<Int_t>(pads_updated.size());

    for (TH1D* h : h_layers) {
      if (h) delete h;
    }
    h_layers.clear();

    if (ipad > 0) {
      canvas->Print(pdf_path);
      for (TObject* obj : page_owned_draw) delete obj;
      page_owned_draw.clear();
    }
  } else {
  if (tree_pad_mode) {
    const Bool_t pion_tracks_only = (de_source == "pion");
    Long64_t n_hits_filled = 0;
    if (!FillTreeHists(f, tree_name, clsize_min, clsize_max, min_abs_cos_theta, pion_tracks_only, kTRUE,
                       h_tree_layers, h_tree_pads, n_hits_filled)) {
      f->Close();
      delete hGainScale;
      delete hMPV;
      delete hGainMap;
      delete hDeltaGainMap;
      delete canvas;
      return EXIT_FAILURE;
    }
  }

  // Layer 合算 dE ヒストを先に fit して、pad fit の初期値 seed を作る。
  // Width は合算ヒストの方が太く出やすいので少し縮めて使う。
  for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
    TH1D* h_layer = nullptr;
    if (tree_pad_mode) {
      if (layer >= static_cast<Int_t>(h_tree_layers.size())) continue;
      h_layer = h_tree_layers[layer];
    } else {
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
    Double_t width_init_layer = min_width;
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
    if (tree_pad_mode) {
      const auto it_pad = h_tree_pads.find({entry.layer, entry.row});
      if (it_pad == h_tree_pads.end()) continue;
      h = it_pad->second;
    } else {
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

    const Int_t bmin_stat = h_fit->FindBin(-stat_range);
    const Int_t bmax_stat = h_fit->FindBin(stat_range);
    const Double_t n_stat = h_fit->Integral(bmin_stat, bmax_stat);
    if (threshold > 0 && n_stat < threshold) {
      if (h_rebinned) delete h_rebinned;
      continue;
    }

    const Int_t peak_bin = SelectGainPeakBin(h_fit, fit_cfg);
    const Double_t peak_x = h_fit->GetBinCenter(peak_bin);
    const Double_t peak_y = h_fit->GetBinContent(peak_bin);
    const Int_t n_bins_window = bmax_stat - bmin_stat + 1;
    const Double_t mean_content = (n_bins_window > 0) ? (n_stat / n_bins_window) : 0.0;
    if (mean_content <= 0.0 || peak_y < mean_content * min_peak_ratio) continue;

    Double_t mpv_hint = -1.0;
    Double_t width_init = min_width;
    auto it_prefit = layer_prefit_seed.find(entry.layer);
    if (it_prefit != layer_prefit_seed.end()) {
      mpv_hint = it_prefit->second.first;
      width_init = std::max(min_width, it_prefit->second.second);
      const Double_t xmin = h_fit->GetXaxis()->GetXmin();
      const Double_t xmax = h_fit->GetXaxis()->GetXmax();
      if (mpv_min > 0.0) mpv_hint = std::max(mpv_hint, mpv_min);
      if (mpv_max > 0.0) mpv_hint = std::min(mpv_hint, mpv_max);
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
      if (mpv_min > 0.0 && mpv < mpv_min) quality_ok = kFALSE;
      if (mpv_max > 0.0 && mpv > mpv_max) quality_ok = kFALSE;
      if (width < min_width) quality_ok = kFALSE;
      if (min_fit_ndf > 0 && ndf < min_fit_ndf) quality_ok = kFALSE;
      if (max_chi2_ndf > 0.0 && ndf > 0 && chi2ndf > max_chi2_ndf) quality_ok = kFALSE;
      if (min_fit_prob > 0.0 && ndf > 0 && fit_prob < min_fit_prob) quality_ok = kFALSE;

      if (quality_ok) {
        Double_t scale = target_mpv / mpv;
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
              scale = target_mpv / mpv;
              outlier_scale = (!std::isfinite(scale) || scale < (1.0 / 3.0) || scale > 10.0);
            }
          }
        }
        if (outlier_scale) scale = 1.0;

        const Double_t old = entry.p[1];
        entry.p[1] = replace_mode ? scale : (old * scale);
        n_updated++;
        pads_updated.insert({entry.layer, entry.row});
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

  if (tree_pad_mode) {
    for (auto& kv : h_tree_pads) delete kv.second;
    h_tree_pads.clear();
    for (TH1D* hl : h_tree_layers) delete hl;
    h_tree_layers.clear();
  }
  }  // end hist / tree-pad mode

  // 未更新パッド補完（同一レイヤーかつ同一 ASAD の、Landau 更新済みパッドの gain 平均）
  if (use_asad_mean_fill) {
    std::map<std::pair<Int_t, Int_t>, std::pair<Double_t, Int_t>> lasad_stats;  // (layer, asad) -> (sum gain, count)
    for (const auto& e : entries) {
      if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
      if (tpc::IsPadOnCenterFrame(e.layer, e.row)) continue;
      const auto key = std::make_pair(e.layer, e.row);
      if (!pads_updated.count(key)) continue;
      const Int_t asad = tpc::GetASADId(e.layer, e.row);
      auto& st = lasad_stats[{e.layer, asad}];
      st.first += e.p[1];
      st.second += 1;
    }
    for (auto& e : entries) {
      if (e.is_comment || e.aty != 0 || e.p.size() < 2) continue;
      if (tpc::IsPadOnCenterFrame(e.layer, e.row)) continue;
      const auto key = std::make_pair(e.layer, e.row);
      if (pads_updated.count(key)) continue;
      const Int_t asad = tpc::GetASADId(e.layer, e.row);
      const auto it = lasad_stats.find({e.layer, asad});
      if (it == lasad_stats.end() || it->second.second <= 0) continue;
      const Double_t old_g = e.p[1];
      const Double_t mean_gain = it->second.first / static_cast<Double_t>(it->second.second);
      e.p[1] = mean_gain;
      pads_filled_by_asad.insert(key);
      n_filled_by_asad++;
      std::cout << Form("  layer+ASAD-mean fill (gain): L=%2d R=%3d (ASAD=%2d, n=%3d) gain: %10.6f -> %10.6f",
                        e.layer, e.row, asad, it->second.second, old_g, e.p[1])
                << std::endl;
    }
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
    if (pads_updated.count(key) || pads_filled_by_asad.count(key)) {
      const Double_t delta = now - old;
      hDeltaGainMap->SetBinContent(bin_idx, delta);
      max_abs_delta_gain = std::max(max_abs_delta_gain, std::abs(delta));
    }
  }
  if (max_abs_delta_gain > 0.0) {
    hDeltaGainMap->SetMinimum(-max_abs_delta_gain);
    hDeltaGainMap->SetMaximum(+max_abs_delta_gain);
  }

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
  if (use_asad_mean_fill) {
    std::cout << "layer+ASAD-mean filled pads (gain): " << n_filled_by_asad << std::endl;
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
