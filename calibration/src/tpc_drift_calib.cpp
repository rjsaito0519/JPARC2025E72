#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPRegexp.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "TPCPadHelper.hh"
#include "ana_helper.h"
#include "paths.h"

namespace fs = std::filesystem;

// Linear-fit / profile window on Y_BcOut [mm] (edges are unreliable).
static constexpr Double_t kFitRangeYMm = 40.0;
// Reject nonsense gaus means (ResY hist is typically ±20 mm).
static constexpr Double_t kMaxAbsMeanResYMm = 15.0;
static constexpr Double_t kMaxMeanErrMm = 4.0;
static constexpr Double_t kMinGaussSigmaMm = 0.10;
static constexpr Double_t kMaxGaussSigmaMm = 10.0;
// Line-fit quality / physical sanity (no meaningful correlation → skip update).
static constexpr Double_t kMaxLineChi2Ndf = 40.0;
static constexpr Double_t kMaxAbsSlope = 0.08;

// ----------------------------------------------------------------------------
struct TpcParamEntry {
  Int_t layer = 0;
  Int_t row = 0;
  Int_t aty = 0;
  std::vector<Double_t> p;
  TString raw_line;
  Bool_t is_comment = kTRUE;
};

struct LayerDriftResult {
  Int_t layer = -1;
  Double_t slope = 0.0;
  Double_t slope_err = 0.0;
  Double_t intercept = 0.0;
  Double_t v_old = 0.0;
  Double_t delta_v = 0.0;
  Double_t v_new = 0.0;
  Double_t hist_entries = 0.0;
  Int_t npoints = 0;
  Bool_t ok = kFALSE;
  TString skip_reason;
};

// ----------------------------------------------------------------------------
static Bool_t FitSliceGaussMean(TH1D* py, Int_t layer, Int_t ix, Double_t& mean_resy,
                                Double_t& mean_err)
{
  mean_resy = 0.0;
  mean_err = 0.0;
  if (!py || py->GetEntries() <= 0.0)
    return kFALSE;

  const Double_t max_content = py->GetMaximum();
  const Double_t mean_content =
      (py->GetNbinsX() > 0) ? (py->Integral() / py->GetNbinsX()) : 0.0;
  if (mean_content <= 0.0 || max_content < mean_content * 1.05)
    return kFALSE;

  const Int_t max_bin = py->GetMaximumBin();
  const Double_t peak = py->GetBinCenter(max_bin);
  Double_t rms = py->GetRMS();
  if (rms < 0.1)
    rms = 0.5;
  if (rms > 5.0)
    rms = 2.5;

  const Double_t half = TMath::Max(2.0 * rms, 0.5);
  const Double_t fit_min = peak - half;
  const Double_t fit_max = peak + half;
  if (!(fit_min < fit_max))
    return kFALSE;

  TF1 sliceFit(Form("tpc_drift_slice_l%d_x%d", layer, ix), "gaus", fit_min, fit_max);
  sliceFit.SetParameters(max_content, peak, TMath::Max(rms * 0.5, 0.05));
  sliceFit.SetParLimits(2, kMinGaussSigmaMm, kMaxGaussSigmaMm);
  const Int_t status = py->Fit(&sliceFit, "QRN0");
  if (status != 0)
    return kFALSE;

  mean_resy = sliceFit.GetParameter(1);
  mean_err = sliceFit.GetParError(1);
  const Double_t sigma = sliceFit.GetParameter(2);
  if (!(std::isfinite(mean_resy) && std::isfinite(mean_err) && mean_err > 0.0))
    return kFALSE;
  if (!(std::isfinite(sigma)) || sigma < kMinGaussSigmaMm || sigma > kMaxGaussSigmaMm)
    return kFALSE;
  // Wild means / huge errors = no usable correlation (bad slice).
  if (TMath::Abs(mean_resy) > kMaxAbsMeanResYMm || mean_err > kMaxMeanErrMm)
    return kFALSE;
  return kTRUE;
}

// ----------------------------------------------------------------------------
static Double_t MedianOf(std::vector<Double_t> vals)
{
  if (vals.empty())
    return 0.0;
  std::sort(vals.begin(), vals.end());
  const size_t n = vals.size();
  if (n % 2 == 1)
    return vals[n / 2];
  return 0.5 * (vals[n / 2 - 1] + vals[n / 2]);
}

// ----------------------------------------------------------------------------
int main(Int_t argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " <input.root> [run_number] [--mode hit|trk] [--vdrift V]"
              << " [--threshold N] [--min-entries N] [--rebin N] [--min-points N]"
              << " [--layer-min L] [--layer-max L] [--debug]" << std::endl;
    return EXIT_FAILURE;
  }

  TString input_root = argv[1];
  Int_t run_number = -1;
  TString mode = "hit";
  Double_t vdrift_override = -1.0;
  // Minimum total entries on layer 2D hist; below this keep vdrift unchanged.
  Int_t threshold = 1000000;
  Int_t min_entries = 30;  // per Y-bin for gaus profile
  Int_t rebin_x = 1;
  Int_t min_points = 5;
  Int_t layer_min = 0;
  Int_t layer_max = tpc::NumOfLayersTPC - 1;
  Bool_t debug_mode = kFALSE;

  Int_t arg_idx = 2;
  if (arg_idx < argc && TString(argv[arg_idx]).IsDigit()) {
    run_number = std::atoi(argv[arg_idx++]);
  }

  for (Int_t i = arg_idx; i < argc; ++i) {
    TString arg = argv[i];
    if (arg == "--mode" && i + 1 < argc) {
      mode = argv[++i];
    } else if (arg == "--vdrift" && i + 1 < argc) {
      vdrift_override = std::atof(argv[++i]);
    } else if (arg == "--threshold" && i + 1 < argc) {
      threshold = std::atoi(argv[++i]);
    } else if (arg == "--min-entries" && i + 1 < argc) {
      min_entries = std::atoi(argv[++i]);
    } else if (arg == "--rebin" && i + 1 < argc) {
      rebin_x = std::max(1, std::atoi(argv[++i]));
    } else if (arg == "--min-points" && i + 1 < argc) {
      min_points = std::max(2, std::atoi(argv[++i]));
    } else if (arg == "--layer-min" && i + 1 < argc) {
      layer_min = std::atoi(argv[++i]);
    } else if (arg == "--layer-max" && i + 1 < argc) {
      layer_max = std::atoi(argv[++i]);
    } else if (arg == "--debug") {
      debug_mode = kTRUE;
    }
  }

  const Bool_t prefer_trk =
      (mode == "trk" || mode == "track" || mode == "TPCTrk");
  const char* hist_fmts_hit_first[] = {"TPCHit_ResY_vs_Y_Layer%02d",
                                       "TPCTrk_ResY_vs_Y_Layer%02d"};
  const char* hist_fmts_trk_first[] = {"TPCTrk_ResY_vs_Y_Layer%02d",
                                       "TPCHit_ResY_vs_Y_Layer%02d"};
  const char** hist_fmts = prefer_trk ? hist_fmts_trk_first : hist_fmts_hit_first;

  TFile* f = TFile::Open(input_root.Data());
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << input_root << std::endl;
    return EXIT_FAILURE;
  }

  if (run_number < 0) {
    TTree* tree = dynamic_cast<TTree*>(f->Get("tpc"));
    if (tree) {
      try {
        TTreeReader reader(tree);
        TTreeReaderValue<UInt_t> rv(reader, "run_number");
        if (reader.Next()) {
          run_number = static_cast<Int_t>(*rv);
          std::cout << "Auto-detected run number from TTree: " << run_number
                    << std::endl;
        }
      } catch (...) {
      }
    }
    if (run_number < 0) {
      TPRegexp run_regex("run([0-9]+)");
      TString run_str = input_root(run_regex);
      if (!run_str.IsNull()) {
        run_str.ReplaceAll("run", "");
        run_number = run_str.Atoi();
        std::cout << "Auto-detected run number from filename: " << run_number
                  << std::endl;
      }
    }
  }
  if (run_number < 0) {
    std::cerr << "Error: run number not provided and could not be auto-detected."
              << std::endl;
    f->Close();
    return EXIT_FAILURE;
  }

  TString hist_fmt;
  TString resolved_tag;
  for (Int_t ifmt = 0; ifmt < 2; ++ifmt) {
    Int_t nfound = 0;
    for (Int_t layer = layer_min; layer <= layer_max; ++layer) {
      if (layer < 0 || layer >= tpc::NumOfLayersTPC)
        continue;
      TString hname = Form(hist_fmts[ifmt], layer);
      if (f->Get(hname.Data()))
        ++nfound;
    }
    if (nfound > 0) {
      hist_fmt = hist_fmts[ifmt];
      resolved_tag = TString(hist_fmts[ifmt]).Contains("TPCTrk") ? "trk" : "hit";
      std::cout << "Histogram series: " << hist_fmt << " (found " << nfound
                << " layers, mode prefer=" << mode << ")" << std::endl;
      break;
    }
  }
  if (hist_fmt.IsNull()) {
    std::cerr << "Error: no TPCHit/TPCTrk ResY_vs_Y_Layer histograms found."
              << std::endl;
    f->Close();
    return EXIT_FAILURE;
  }

  std::cout << "threshold (min hist entries) = " << threshold
            << ", fit |Y| <= " << kFitRangeYMm << " mm" << std::endl;

  TString base_dir = ANALYZER_DIR + "/param/TPCPRM";
  TString e72_dir = base_dir + "/e72";
  TString run_param_path =
      Form("%s/TPCParam_e72_run%05d", e72_dir.Data(), run_number);
  TString base_param_path = Form("%s/TPCParam_0_yoffset_adjusted", base_dir.Data());
  TString input_param_path = run_param_path;
  if (!fs::exists(input_param_path.Data())) {
    input_param_path = base_param_path;
    std::cout << "Input param file not found in e72, using base: "
              << base_param_path << std::endl;
  } else {
    std::cout << "Using existing param file for update: " << input_param_path
              << std::endl;
  }

  std::vector<TpcParamEntry> entries;
  {
    std::ifstream ifs(input_param_path.Data());
    if (!ifs.is_open()) {
      std::cerr << "Error: cannot open param file " << input_param_path
                << std::endl;
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
        if (line.find("# --- end TPCPHASE pointer ---") != std::string::npos)
          in_tpcphase_stamp_block = false;
        continue;
      }

      TpcParamEntry entry;
      entry.raw_line = line;
      std::string trimmed = line;
      const size_t first = trimmed.find_first_not_of(" \t");
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
            while (ss >> val)
              entry.p.push_back(val);
          }
        }
      }
      entries.push_back(entry);
    }
  }

  std::map<Int_t, std::vector<Double_t>> layer_p1s;
  for (const auto& e : entries) {
    if (e.is_comment || e.aty != 2 || e.p.size() < 2)
      continue;
    layer_p1s[e.layer].push_back(e.p[1]);
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
  TString pdf_path =
      img_dir + Form("/tpc_drift_calib_%s.pdf", resolved_tag.Data());
  TCanvas* canvas = new TCanvas("c_drift", "c_drift", 1200, 900);

  TH2Poly* hSummary = new TH2Poly(
      "hDriftSummary",
      Form("TPC drift calib summary (Run %05d, %s);Z [mm];X [mm]", run_number,
           resolved_tag.Data()),
      -300, 300, -300, 300);
  hSummary->SetDirectory(nullptr);
  tpc::InitializeHistograms(hSummary);

  TH2Poly* hDeltaMap = new TH2Poly(
      "hTPC_Drift_dV",
      Form("TPC #delta v_{drift} (Run %05d, %s);Z [mm];X [mm]", run_number,
           resolved_tag.Data()),
      -300, 300, -300, 300);
  hDeltaMap->SetDirectory(nullptr);
  tpc::InitializeHistograms(hDeltaMap);

  TH2Poly* hVNewMap = new TH2Poly(
      "hTPC_Drift_V",
      Form("TPC v_{drift} (after this calib, Run %05d, %s);Z [mm];X [mm]", run_number,
           resolved_tag.Data()),
      -300, 300, -300, 300);
  hVNewMap->SetDirectory(nullptr);
  tpc::InitializeHistograms(hVNewMap);

  canvas->Print(pdf_path + "[");

  std::vector<LayerDriftResult> results;
  std::map<Int_t, Double_t> delta_by_layer;

  std::cout << "Fitting layer slopes..." << std::endl;
  for (Int_t layer = layer_min; layer <= layer_max; ++layer) {
    if (layer < 0 || layer >= tpc::NumOfLayersTPC)
      continue;

    LayerDriftResult res;
    res.layer = layer;

    TString hname = Form(hist_fmt.Data(), layer);
    TH2D* h2 = dynamic_cast<TH2D*>(f->Get(hname.Data()));
    if (!h2) {
      res.skip_reason = "hist missing";
      results.push_back(res);
      continue;
    }

    res.hist_entries = h2->GetEntries();
    if (res.hist_entries < threshold) {
      res.skip_reason =
          Form("entries=%.0f < threshold=%d (keep vdrift)", res.hist_entries,
               threshold);
      std::cout << "  Layer " << Form("%02d", layer) << " skipped: "
                << res.skip_reason << std::endl;
      canvas->Clear();
      canvas->Divide(1, 2);
      canvas->cd(1);
      h2->Draw("COLZ");
      canvas->cd(2);
      TLatex tex;
      tex.SetNDC();
      tex.SetTextSize(0.045);
      tex.DrawLatex(0.12, 0.6, Form("Layer %02d skipped: %s", layer,
                                    res.skip_reason.Data()));
      canvas->Print(pdf_path);
      results.push_back(res);
      continue;
    }

    TH2D* hwork = h2;
    TH2D* hclone = nullptr;
    if (rebin_x > 1) {
      hclone = dynamic_cast<TH2D*>(h2->Clone(Form("tpc_drift_rebin_l%d", layer)));
      if (hclone) {
        hclone->RebinX(rebin_x);
        hwork = hclone;
      }
    }

    const Double_t v_layer =
        (vdrift_override > 0.0)
            ? vdrift_override
            : (layer_p1s.count(layer) ? MedianOf(layer_p1s[layer]) : 0.055);
    res.v_old = v_layer;
    if (!(v_layer > 0.0)) {
      res.skip_reason = "invalid vdrift";
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }

    std::vector<Double_t> xs, ys, exs, eys;
    const Int_t nx = hwork->GetNbinsX();
    for (Int_t ix = 1; ix <= nx; ++ix) {
      const Double_t y_pos = hwork->GetXaxis()->GetBinCenter(ix);
      if (TMath::Abs(y_pos) > kFitRangeYMm)
        continue;

      TH1D* py =
          hwork->ProjectionY(Form("tpc_drift_py_l%d_x%d", layer, ix), ix, ix);
      if (!py)
        continue;
      if (py->GetEntries() < min_entries) {
        delete py;
        continue;
      }
      Double_t mean = 0.0, err = 0.0;
      if (!FitSliceGaussMean(py, layer, ix, mean, err)) {
        delete py;
        continue;
      }
      xs.push_back(y_pos);
      ys.push_back(mean);
      exs.push_back(0.0);
      eys.push_back(err);
      delete py;
    }
    res.npoints = static_cast<Int_t>(xs.size());

    if (res.npoints < min_points) {
      res.skip_reason = Form("npoints=%d < %d", res.npoints, min_points);
      canvas->Clear();
      canvas->Divide(1, 2);
      canvas->cd(1);
      h2->Draw("COLZ");
      canvas->cd(2);
      TLatex tex;
      tex.SetNDC();
      tex.SetTextSize(0.05);
      tex.DrawLatex(0.15, 0.6, Form("Layer %02d skipped: %s", layer,
                                    res.skip_reason.Data()));
      canvas->Print(pdf_path);
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }

    TGraphErrors gr(res.npoints, xs.data(), ys.data(), exs.data(), eys.data());
    gr.SetName(Form("tpc_drift_profile_l%d", layer));
    gr.SetTitle(Form("Layer %02d profile (|Y|#leq%.0f mm);Y [mm];Residual Y mean [mm]",
                     layer, kFitRangeYMm));
    TF1 lineFit(Form("tpc_drift_line_l%d", layer), "pol1", -kFitRangeYMm,
                kFitRangeYMm);
    const Int_t fit_status = gr.Fit(&lineFit, "QRN0");

    auto draw_skip_page = [&](const TString& reason) {
      canvas->Clear();
      canvas->Divide(1, 2);
      canvas->cd(1);
      gPad->SetRightMargin(0.12);
      h2->Draw("COLZ");
      canvas->cd(2);
      gr.SetMarkerStyle(20);
      gr.SetMarkerSize(0.6);
      gr.Draw("AP");
      if (fit_status == 0) {
        lineFit.SetLineColor(kRed);
        lineFit.Draw("SAME");
      }
      TLatex tex;
      tex.SetNDC();
      tex.SetTextSize(0.04);
      tex.SetTextColor(kRed + 1);
      tex.DrawLatex(0.12, 0.85, Form("Layer %02d SKIPPED (no update)", layer));
      tex.SetTextColor(kBlack);
      tex.DrawLatex(0.12, 0.78, reason.Data());
      canvas->Print(pdf_path);
    };

    if (fit_status != 0) {
      res.skip_reason = "line fit failed";
      draw_skip_page(res.skip_reason);
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }

    res.intercept = lineFit.GetParameter(0);
    res.slope = lineFit.GetParameter(1);
    res.slope_err = lineFit.GetParError(1);
    const Double_t a = res.slope;
    const Double_t chi2 = lineFit.GetChisquare();
    const Int_t ndf = lineFit.GetNDF();
    const Double_t chi2ndf = (ndf > 0) ? (chi2 / ndf) : 1e9;

    if (!(std::isfinite(a)) || TMath::Abs(1.0 + a) < 1e-6) {
      res.skip_reason = Form("bad slope a=%.6g", a);
      draw_skip_page(res.skip_reason);
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }
    if (ndf < 1 || chi2ndf > kMaxLineChi2Ndf) {
      res.skip_reason =
          Form("bad correlation: chi2/ndf=%.3g (ndf=%d) > %.3g", chi2ndf, ndf,
               kMaxLineChi2Ndf);
      std::cout << "  Layer " << Form("%02d", layer) << " skipped: "
                << res.skip_reason << std::endl;
      draw_skip_page(res.skip_reason);
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }
    if (TMath::Abs(a) > kMaxAbsSlope) {
      res.skip_reason =
          Form("bad correlation: |a|=%.5g > %.3g", TMath::Abs(a), kMaxAbsSlope);
      std::cout << "  Layer " << Form("%02d", layer) << " skipped: "
                << res.skip_reason << std::endl;
      draw_skip_page(res.skip_reason);
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }
    if (TMath::Abs(res.intercept) > kMaxAbsMeanResYMm) {
      res.skip_reason =
          Form("bad correlation: |intercept|=%.4g > %.3g",
               TMath::Abs(res.intercept), kMaxAbsMeanResYMm);
      std::cout << "  Layer " << Form("%02d", layer) << " skipped: "
                << res.skip_reason << std::endl;
      draw_skip_page(res.skip_reason);
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }

    res.delta_v = -a / (1.0 + a) * v_layer;
    res.v_new = v_layer + res.delta_v;
    if (!(res.v_new > 0.0) || !std::isfinite(res.delta_v)) {
      res.skip_reason = Form("non-positive v_new=%.6g", res.v_new);
      draw_skip_page(res.skip_reason);
      results.push_back(res);
      if (hclone)
        delete hclone;
      continue;
    }

    res.ok = kTRUE;
    delta_by_layer[layer] = res.delta_v;
    results.push_back(res);

    canvas->Clear();
    canvas->Divide(1, 2);
    canvas->cd(1);
    gPad->SetRightMargin(0.12);
    h2->Draw("COLZ");
    canvas->cd(2);
    gr.SetMarkerStyle(20);
    gr.SetMarkerSize(0.6);
    gr.Draw("AP");
    lineFit.SetLineColor(kRed);
    lineFit.Draw("SAME");
    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.12, 0.85,
                  Form("a=%.5g #pm %.2g  v=%.5g  #deltav=%.5g  v_{new}=%.5g  n=%d",
                       res.slope, res.slope_err, res.v_old, res.delta_v, res.v_new,
                       res.npoints));
    tex.DrawLatex(0.12, 0.78,
                  Form("entries=%.0f  |Y|#leq%.0f mm  chi2/ndf=%.3g",
                       res.hist_entries, kFitRangeYMm, chi2ndf));
    canvas->Print(pdf_path);

    std::cout << Form(
                     "  Layer %02d: a=%.6g +/- %.2g  v=%.6g  dv=%.6g  vnew=%.6g  n=%d  "
                     "entries=%.0f",
                     layer, res.slope, res.slope_err, res.v_old, res.delta_v, res.v_new,
                     res.npoints, res.hist_entries)
              << std::endl;

    if (hclone)
      delete hclone;
  }

  // Apply δv and fill maps for ALL aty==2 pads (unfilled bins would look like v=0).
  std::vector<Double_t> all_dv_updated;
  std::vector<Double_t> all_v_map;
  std::vector<Double_t> all_dv_map;
  Int_t n_pads_updated = 0;
  for (auto& e : entries) {
    if (e.is_comment || e.aty != 2 || e.p.size() < 2)
      continue;

    Double_t dv = 0.0;
    auto it = delta_by_layer.find(e.layer);
    if (it != delta_by_layer.end()) {
      dv = it->second;
      e.p[1] += dv;
      ++n_pads_updated;
      all_dv_updated.push_back(dv);

      const Int_t bin_idx = tpc::GetPadId(e.layer, e.row) + 1;
      hSummary->SetBinContent(bin_idx, 100.0);
    }

    const Int_t bin_idx = tpc::GetPadId(e.layer, e.row) + 1;
    hDeltaMap->SetBinContent(bin_idx, dv);
    hVNewMap->SetBinContent(bin_idx, e.p[1]);
    all_dv_map.push_back(dv);
    all_v_map.push_back(e.p[1]);
  }

  // Text summary page
  {
    canvas->Clear();
    canvas->cd();
    TLatex tex;
    tex.SetNDC();
    tex.SetTextSize(0.035);
    tex.DrawLatex(0.10, 0.92, Form("TPC drift calib summary  Run %05d  (%s)",
                                   run_number, resolved_tag.Data()));
    tex.DrawLatex(0.10, 0.88,
                  Form("threshold=%d  fit |Y|#leq%.0f mm", threshold, kFitRangeYMm));
    Double_t y = 0.82;
    Int_t n_ok = 0;
    for (const auto& r : results) {
      if (!r.ok)
        continue;
      ++n_ok;
      if (y < 0.08) {
        canvas->Print(pdf_path);
        canvas->Clear();
        y = 0.92;
      }
      tex.DrawLatex(
          0.10, y,
          Form("L%02d  a=%+.5g  dv=%+.6g  v: %.6g -> %.6g  n=%d", r.layer, r.slope,
               r.delta_v, r.v_old, r.v_new, r.npoints));
      y -= 0.035;
    }
    if (n_ok == 0)
      tex.DrawLatex(0.10, 0.70, "No layers updated.");
    canvas->Print(pdf_path);
  }

  // 1D distributions
  auto axis_with_margin = [](Double_t lo, Double_t hi) {
    if (!(hi > lo)) {
      lo -= 1e-4;
      hi += 1e-4;
    }
    const Double_t m = 0.1 * (hi - lo);
    return std::make_pair(lo - m, hi + m);
  };

  if (!all_dv_updated.empty()) {
    const auto mm = std::minmax_element(all_dv_updated.begin(), all_dv_updated.end());
    const auto ax = axis_with_margin(*mm.first, *mm.second);
    TH1D* hDeltaDist =
        new TH1D("hDeltaVDist",
                 Form("TPC #delta v_{drift} (updated pads, Run %05d);#delta v;pads / bin",
                      run_number),
                 50, ax.first, ax.second);
    hDeltaDist->SetDirectory(nullptr);
    for (Double_t v : all_dv_updated)
      hDeltaDist->Fill(v);

    const auto mmv = std::minmax_element(all_v_map.begin(), all_v_map.end());
    const auto axv = axis_with_margin(*mmv.first, *mmv.second);
    TH1D* hVDist =
        new TH1D("hVDist",
                 Form("TPC v_{drift} (all pads, Run %05d);v;pads / bin", run_number),
                 50, axv.first, axv.second);
    hVDist->SetDirectory(nullptr);
    for (Double_t v : all_v_map)
      hVDist->Fill(v);

    canvas->Clear();
    canvas->Divide(1, 2);
    canvas->cd(1);
    hDeltaDist->Draw("HIST");
    canvas->cd(2);
    hVDist->Draw("HIST");
    canvas->Print(pdf_path);
    delete hDeltaDist;
    delete hVDist;
  }

  canvas->Print(pdf_path + "]");

  // High-res raster TH2Poly maps (same style as time offset)
  std::cout << "Generating high-resolution rasterized maps..." << std::endl;
  TCanvas* cHR = new TCanvas("cHR_drift", "High Res", 1200, 1000);
  gStyle->SetOptStat(0);
  std::vector<TString> tmp_pdfs;

  auto ProcessRaster = [&](TH2Poly* hp, TString tag, Double_t zmin, Double_t zmax) {
    if (!hp)
      return;
    if (zmax > zmin) {
      hp->SetMinimum(zmin);
      hp->SetMaximum(zmax);
    }
    cHR->Clear();
    cHR->cd();
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.20);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.11);
    hp->GetXaxis()->SetTitleOffset(1.15);
    hp->GetYaxis()->SetTitleOffset(1.20);
    hp->GetZaxis()->SetTitleOffset(1.25);
    hp->GetXaxis()->SetLabelSize(0.035);
    hp->GetYaxis()->SetLabelSize(0.035);
    hp->GetZaxis()->SetLabelSize(0.033);
    hp->Draw("COLZ L");
    TString png = img_dir + "/tmp_drift_" + tag + ".png";
    TString tmp_pdf = img_dir + "/tmp_drift_" + tag + ".pdf";
    cHR->Print(png);
    system(Form("convert %s -resize 1200x1000 -gravity center -extent 1200x1000 %s",
                png.Data(), tmp_pdf.Data()));
    tmp_pdfs.push_back(tmp_pdf);
    system(Form("rm %s", png.Data()));
  };

  auto z_range_from = [&](const std::vector<Double_t>& vals, Bool_t include_zero) {
    if (vals.empty())
      return std::make_pair(0.0, 1.0);
    const auto mm = std::minmax_element(vals.begin(), vals.end());
    Double_t lo = *mm.first;
    Double_t hi = *mm.second;
    if (include_zero) {
      lo = std::min(lo, 0.0);
      hi = std::max(hi, 0.0);
    }
    if (!(hi > lo)) {
      lo -= 1e-4;
      hi += 1e-4;
    }
    const Double_t m = 0.15 * (hi - lo);
    return std::make_pair(lo - m, hi + m);
  };

  ProcessRaster(hSummary, "summary", 0.0, 110.0);
  {
    // Prefer updated-only span so tiny δv differences are visible; always include 0.
    const auto zr = z_range_from(
        all_dv_updated.empty() ? all_dv_map : all_dv_updated, /*include_zero=*/kTRUE);
    ProcessRaster(hDeltaMap, "deltav", zr.first, zr.second);
  }
  {
    const auto zr = z_range_from(all_v_map, /*include_zero=*/kFALSE);
    ProcessRaster(hVNewMap, "v", zr.first, zr.second);
  }

  if (!tmp_pdfs.empty()) {
    TString unite_cmd = "pdfunite " + pdf_path + " ";
    for (const auto& p : tmp_pdfs)
      unite_cmd += p + " ";
    unite_cmd += pdf_path + "_merged.pdf";
    const Int_t ures = system(unite_cmd.Data());
    if (ures == 0) {
      system(Form("mv %s_merged.pdf %s", pdf_path.Data(), pdf_path.Data()));
      std::cout << "Successfully merged vector and raster maps." << std::endl;
    }
    for (const auto& p : tmp_pdfs)
      system(Form("rm %s", p.Data()));
  }
  delete cHR;

  if (!debug_mode) {
    if (!fs::exists(e72_dir.Data()))
      fs::create_directories(e72_dir.Data());
    if (fs::exists(run_param_path.Data())) {
      fs::copy(run_param_path.Data(), (run_param_path + ".bak").Data(),
               fs::copy_options::overwrite_existing);
    }
    std::ofstream ofs(run_param_path.Data());
    if (!ofs.is_open()) {
      std::cerr << "Error: cannot open output file " << run_param_path
                << std::endl;
      f->Close();
      return EXIT_FAILURE;
    }
    for (const auto& entry : entries) {
      if (entry.is_comment) {
        ofs << entry.raw_line.Data() << "\n";
      } else {
        ofs << entry.layer << "\t" << entry.row << "\t" << entry.aty;
        for (size_t i = 0; i < entry.p.size(); ++i)
          ofs << "\t" << std::fixed << std::setprecision(8) << entry.p[i];
        ofs << "\n";
      }
    }
    ofs.close();
  }

  const Int_t n_layers_ok = static_cast<Int_t>(delta_by_layer.size());
  std::cout << "Done! Updated " << n_layers_ok << " layers (" << n_pads_updated
            << " pads)." << std::endl;
  if (debug_mode)
    std::cout << "Debug mode: parameter file was NOT updated." << std::endl;
  else
    std::cout << "Result saved to: " << run_param_path << std::endl;
  std::cout << "QA PDF saved to: " << pdf_path << std::endl;

  delete hSummary;
  delete hDeltaMap;
  delete hVNewMap;
  delete canvas;
  f->Close();
  return EXIT_SUCCESS;
}
