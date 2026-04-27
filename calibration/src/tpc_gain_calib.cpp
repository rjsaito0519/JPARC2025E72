#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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
#include <TLine.h>
#include <TMath.h>
#include <TPRegexp.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "ana_helper.h"
#include "paths.h"
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

int main(Int_t argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <input.root> [run_number] [--target-mpv M] [--threshold N] [--stat-range R] "
                 "[--min-peak-ratio X] [--fit-nsigma N] [--local-peak-half W] [--min-width W] "
                 "[--max-chi2-ndf X] [--min-fit-ndf N] [--min-fit-prob P] [--rebin N] "
                 "[--mpv-min X] [--mpv-max X] [--de-source all|pion] [--debug] [--replace] [--mode hit|trk]"
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
    }
  }

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
  const TString pdf_path = img_dir + Form("/tpc_gain_calib_%s.pdf", de_source.Data());
  TCanvas* canvas = new TCanvas("c_gain", "c_gain", 1200, 1000);
  canvas->Print(pdf_path + "[");
  const Int_t nx = 4, ny = 4;
  Int_t ipad = 0;
  std::vector<TH1D*> page_owned_hists;

  TH1D* hGainScale = new TH1D("hTPC_GainScale", "TPC gain scale (target/MPV);target/MPV;pads / bin", 120, 0.4, 1.8);
  TH1D* hMPV = new TH1D("hTPC_Gain_MPV", Form("TPC dE MPV (Run %05d);MPV;pads / bin", run_number), 120, 0, 1200);
  TH2D* hGainMap = new TH2D("hTPC_GainMap",
                            Form("TPC gain map after calib (Run %05d);layer;row", run_number),
                            tpc::NumOfLayersTPC, -0.5, tpc::NumOfLayersTPC - 0.5, 260, -0.5, 259.5);
  TH2D* hDeltaGainMap = new TH2D("hTPC_DeltaGainMap",
                                 Form("TPC #Delta gain map (new-old) (Run %05d);layer;row", run_number),
                                 tpc::NumOfLayersTPC, -0.5, tpc::NumOfLayersTPC - 0.5, 260, -0.5, 259.5);
  hGainScale->SetDirectory(nullptr);
  hMPV->SetDirectory(nullptr);
  hGainMap->SetDirectory(nullptr);
  hDeltaGainMap->SetDirectory(nullptr);

  Int_t n_updated = 0;
  Int_t n_total_pads = 0;
  std::set<std::pair<Int_t, Int_t>> pads_updated;
  std::map<Int_t, std::pair<Double_t, Int_t>> layer_mpv_stats;
  std::map<Int_t, std::pair<Double_t, Double_t>> layer_prefit_seed;

  // Layer 合算 dE ヒストを先に fit して、pad fit の初期値 seed を作る。
  // Width は合算ヒストの方が太く出やすいので少し縮めて使う。
  for (Int_t layer = 0; layer < 32; ++layer) {
    const TString hname_layer =
        (de_source == "pion")
            ? Form("TPCCl_dE_Pion_Layer%02d", layer)
            : Form("TPCCl_dE_Layer%02d", layer);
    TH1D* h_layer = nullptr;
    f->GetObject(hname_layer.Data(), h_layer);
    if (!h_layer) {
      const TString hname_hist = "hist/" + hname_layer;
      f->GetObject(hname_hist.Data(), h_layer);
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

    const Int_t max_bin = h_layer_fit->GetMaximumBin();
    const Double_t peak_x = h_layer_fit->GetBinCenter(max_bin);
    const Double_t peak_y = h_layer_fit->GetBinContent(max_bin);
    Double_t rms = h_layer_fit->GetRMS();
    if (!(rms > 0.0)) rms = 20.0;

    Double_t fit_min = peak_x - 2.5 * rms;
    Double_t fit_max = peak_x + 2.5 * rms;
    const Double_t xmin = h_layer_fit->GetXaxis()->GetXmin();
    const Double_t xmax = h_layer_fit->GetXaxis()->GetXmax();
    fit_min = std::max(fit_min, xmin);
    fit_max = std::min(fit_max, xmax);
    if (fit_max - fit_min < 40.0) {
      fit_min = std::max(peak_x - 20.0, xmin);
      fit_max = std::min(peak_x + 20.0, xmax);
    }
    if (fit_max <= fit_min) {
      if (h_layer_rebinned) delete h_layer_rebinned;
      continue;
    }

    TF1* fLayer = new TF1(Form("fLayerLandau_L%02d", layer), "landau(0)+pol0(3)", fit_min, fit_max);
    const Double_t span = fit_max - fit_min;
    const Double_t width_init = std::max(1.0, 0.18 * span);
    fLayer->SetParameters(std::max(peak_y, 1.0), peak_x, width_init, 0.0);
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

    const TString hname =
        (de_source == "pion")
            ? Form("TPCCl_dE_Pion_Layer%02d_Row%03d", entry.layer, entry.row)
            : Form("TPCCl_dE_Layer%02d_Row%03d", entry.layer, entry.row);
    TH1D* h = nullptr;
    f->GetObject(hname.Data(), h);
    if (!h) {
      const TString hname_hist = "hist/" + hname;
      f->GetObject(hname_hist.Data(), h);
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

    const Int_t max_bin = h_fit->GetMaximumBin();
    const Double_t peak_x = h_fit->GetBinCenter(max_bin);
    const Double_t peak_y = h_fit->GetBinContent(max_bin);
    const Int_t n_bins_window = bmax_stat - bmin_stat + 1;
    const Double_t mean_content = (n_bins_window > 0) ? (n_stat / n_bins_window) : 0.0;
    if (mean_content <= 0.0 || peak_y < mean_content * min_peak_ratio) continue;

    Double_t rms_local = LocalWeightedRmsAroundPeak(h_fit, peak_x, local_peak_half);
    const Double_t rms_global = h_fit->GetRMS();
    if (!(rms_local > 0.0)) rms_local = rms_global;
    if (!(rms_local > 0.0)) rms_local = 15.0;
    // 狭すぎる窓を避けるため、local を主にしつつ global 成分を少し混ぜる。
    Double_t rms_for_window = rms_local;
    if (rms_global > 0.0) rms_for_window = std::max(rms_local, 0.65 * rms_global);

    Double_t mpv_init_main = peak_x;
    const Double_t xmin = h_fit->GetXaxis()->GetXmin();
    const Double_t xmax = h_fit->GetXaxis()->GetXmax();
    // 細いピークに追従しやすいよう、初期 width は local RMS ベースで与える。
    const Double_t width_seed = std::max(0.5, 0.45 * rms_local);
    Double_t width_init = std::max(min_width, width_seed);
    auto it_prefit = layer_prefit_seed.find(entry.layer);
    if (it_prefit != layer_prefit_seed.end()) {
      mpv_init_main = it_prefit->second.first;
      width_init = std::max(min_width, it_prefit->second.second);
    }

    // MPV 初期値のレンジ制約（物理的にあり得る範囲に寄せる）
    if (mpv_min > 0.0) mpv_init_main = std::max(mpv_init_main, mpv_min);
    if (mpv_max > 0.0) mpv_init_main = std::min(mpv_init_main, mpv_max);
    mpv_init_main = std::min(std::max(mpv_init_main, xmin), xmax);

    // fit 範囲は「初期MPVの近傍」に限定
    const Double_t sigma_for_window = std::max(fit_nsigma * rms_for_window, 2.6 * width_init);
    Double_t fit_min = mpv_init_main - std::max(10.0, 1.0 * sigma_for_window);
    Double_t fit_max = mpv_init_main + std::max(18.0, 2.0 * sigma_for_window);
    fit_min = std::max(fit_min, xmin);
    fit_max = std::min(fit_max, xmax);
    if (mpv_min > 0.0) fit_min = std::max(fit_min, mpv_min - 0.5 * sigma_for_window);
    if (mpv_max > 0.0) fit_max = std::min(fit_max, mpv_max + 1.0 * sigma_for_window);
    if (fit_max - fit_min < 28.0) {
      fit_min = std::max(mpv_init_main - 14.0, xmin);
      fit_max = std::min(mpv_init_main + 14.0, xmax);
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
      return h_fit->Fit(fFit, "Q");
    };

    Int_t fit_status = fit_once(mpv_init_main);
    if (fit_status == 0) {
      Double_t mpv = fFit->GetParameter(1);
      Double_t width = fFit->GetParameter(2);
      Int_t ndf = fFit->GetNDF();
      Double_t chi2 = fFit->GetChisquare();
      Double_t chi2ndf = (ndf > 0) ? (chi2 / static_cast<Double_t>(ndf)) : 0.0;
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
        const Double_t draw_margin = std::max(2.0, 0.25 * fit_span);
        const Double_t draw_min = std::max(xmin, fit_min - draw_margin);
        const Double_t draw_max = std::min(xmax, fit_max + draw_margin);
        h_fit->GetXaxis()->SetRangeUser(draw_min, draw_max);
        h_fit->Draw();
        if (h_rebinned) {
          // Drawn objects must stay alive until canvas->Print() flushes this page.
          page_owned_hists.push_back(h_rebinned);
          keep_h_rebinned_for_page = true;
        }
        fFit->SetLineColor(kRed);
        // 関数は fit 窓だけを描く（外側に延長しない）。
        TF1* fDraw = (TF1*)fFit->Clone(Form("fLandauDraw_L%02d_R%03d", entry.layer, entry.row));
        fDraw->SetRange(fit_min, fit_max);
        fDraw->SetLineColor(kRed);
        fDraw->SetNpx(1000);
        fDraw->Draw("same");
        TLine* v = new TLine(mpv, h_fit->GetMinimum(), mpv, h_fit->GetMaximum());
        v->SetLineStyle(2);
        v->SetLineColor(kBlue + 1);
        v->DrawClone();
        delete fDraw;
        delete v;
        if (ipad >= nx * ny) {
          canvas->Print(pdf_path);
          for (TH1D* htmp : page_owned_hists) delete htmp;
          page_owned_hists.clear();
          ipad = 0;
        }
      }
    }
    delete fFit;
    if (h_rebinned && !keep_h_rebinned_for_page) delete h_rebinned;
  }

  if (ipad > 0) {
    canvas->Print(pdf_path);
    for (TH1D* htmp : page_owned_hists) delete htmp;
    page_owned_hists.clear();
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
    hGainMap->Fill(e.layer, e.row, now);
    const Double_t delta = now - old;
    hDeltaGainMap->Fill(e.layer, e.row, delta);
    max_abs_delta_gain = std::max(max_abs_delta_gain, std::abs(delta));
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
