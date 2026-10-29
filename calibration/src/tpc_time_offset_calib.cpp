#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TString.h>
#include <TPRegexp.h>

#include "paths.h"
#include "ana_helper.h"
#include "TPCPadHelper.hh"
#include "tpc_pad_fill.h"

#include <TCanvas.h>
#include <TH2Poly.h>
#include <TStyle.h>
#include <TPaveText.h>

namespace fs = std::filesystem;

// ----------------------------------------------------------------------------
// TPC Parameter entry
// ----------------------------------------------------------------------------
struct TpcParamEntry {
    Int_t layer;
    Int_t row;
    Int_t aty;
    std::vector<Double_t> p;
    TString raw_line;
    Bool_t is_comment;
};

// 最も高いビン付近だけで重み付き RMS（二峰性で遠い副ピークをフィット窓から外す）
static Double_t LocalWeightedRmsAroundPeak(const TH1D* h, Double_t peak_x, Double_t half_width_mm) {
    if (!h || half_width_mm <= 0.0)
        return -1.0;
    Double_t sumw = 0.0, sumxw = 0.0, sumx2w = 0.0;
    const Int_t nb = h->GetNbinsX();
    for (Int_t i = 1; i <= nb; ++i) {
        const Double_t x = h->GetBinCenter(i);
        if (TMath::Abs(x - peak_x) > half_width_mm)
            continue;
        const Double_t w = h->GetBinContent(i);
        if (w <= 0.0)
            continue;
        sumw += w;
        sumxw += w * x;
        sumx2w += w * x * x;
    }
    if (sumw < 1.0)
        return -1.0;
    const Double_t mean = sumxw / sumw;
    const Double_t var = sumx2w / sumw - mean * mean;
    if (!(var > 0.0) || !std::isfinite(var))
        return -1.0;
    const Double_t rms = TMath::Sqrt(var);
    return std::isfinite(rms) ? rms : -1.0;
}

// 各ビンを中心に ±sep の窓で見て「窓内で最大の高さ」なら峰候補。候補のうち最も高い峰（同高なら |x| 小さい方）
static Int_t DominantLocalMaximumBin(const TH1D* h, Int_t sep) {
    if (!h)
        return 0;
    if (sep < 1)
        return h->GetMaximumBin();
    const Int_t nb = h->GetNbinsX();
    Double_t best_y = -1.0;
    Int_t best_bin = h->GetMaximumBin();
    for (Int_t i = 1; i <= nb; ++i) {
        const Double_t y = h->GetBinContent(i);
        if (y <= 0.0)
            continue;
        Bool_t is_peak = kTRUE;
        const Int_t j0 = std::max(1, i - sep);
        const Int_t j1 = std::min(nb, i + sep);
        for (Int_t j = j0; j <= j1; ++j) {
            if (j != i && h->GetBinContent(j) > y) {
                is_peak = kFALSE;
                break;
            }
        }
        if (!is_peak)
            continue;
        const Double_t xi = TMath::Abs(h->GetBinCenter(i));
        const Double_t xb = TMath::Abs(h->GetBinCenter(best_bin));
        if (y > best_y + 1e-12) {
            best_y = y;
            best_bin = i;
        } else if (TMath::Abs(y - best_y) <= 1e-12 && xi < xb) {
            best_bin = i;
        }
    }
    if (best_y < 0.0)
        return h->GetMaximumBin();
    return best_bin;
}

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
static void PrintUsage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input.root> [run_number] [--mode hit|trk] [--threshold N] [--ave] [--common] [--debug]"
              << std::endl;
}

int main(Int_t argc, char** argv) {
    if (argc < 2) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
    }

    TString input_root = argv[1];
    Int_t run_number = -1;
    // threshold=0 で足切り無効（低統計でも peak がちゃんと立っていれば採用）
    Int_t threshold = 10;
    const Double_t sigma_threshold = 10.0;
    const Double_t nsigma = 2.0;
    const Double_t stat_range = 20.0;
    const Double_t min_peak_ratio = 1.07;  // skip fit if peak_y < mean_content * min_peak_ratio (flat histogram)
    // 極細ガウス（ヒスト全体に比べ σ が針状）を弾く: max(絶対下限 mm, 相対下限 × ヒスト RMS)
    const Double_t min_sigma_mm = 0.35;
    const Double_t min_sigma_rel_rms = 0.04;
    // 初期ガウス σ: seed σ に対する係数（やや大きめにすると針状局所解を避けやすい）
    const Double_t sigma_init_rel = 0.72;
    // フィット品質: NDF の下限のみ（χ²/NDF 上限・確率下限は適用しない）
    const Int_t min_fit_ndf = 1;
    // 主ピーク優先: 最も高い局所峰の周辺 RMS でフィット窓を決める
    const Double_t local_peak_half_mm = 4.5;
    const Int_t local_peak_sep_bins = 2;
    // mode: "hit" → TPCHit_ResY_..., "trk" → TPCTrk_ResY_...
    TString mode = "hit";
    Bool_t debug_mode = kFALSE;
    Bool_t use_fill = kFALSE;
    // --common (hit mode): BcOut sets only a common shift per AsAd (median Δp0 of its fitted pads, or the
    // global median when an AsAd has too few); the pad-to-pad differences are left to the trk mode.
    Bool_t common_only = kFALSE;
    const Int_t common_min_pads = 5;

    Int_t arg_idx = 2;
    if (arg_idx < argc && TString(argv[arg_idx]).IsDigit()) {
        run_number = std::atoi(argv[arg_idx++]);
    }

    for (Int_t i = arg_idx; i < argc; ++i) {
        TString arg = argv[i];
        if (arg == "--threshold" && i + 1 < argc) {
            threshold = std::atoi(argv[++i]);
        } else if (arg == "--mode" && i + 1 < argc) {
            mode = argv[++i];  // "hit" or "trk"
        } else if (arg == "--ave") {
            use_fill = kTRUE;
        } else if (arg == "--common") {
            common_only = kTRUE;
        } else if (arg == "--debug") {
            debug_mode = kTRUE;
        } else {
            std::cerr << "Error: unknown or incomplete argument: " << arg << std::endl;
            PrintUsage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (mode != "hit" && mode != "trk") {
        std::cerr << "Error: --mode must be hit or trk (got " << mode << ")" << std::endl;
        return EXIT_FAILURE;
    }
    if (common_only && (mode != "hit" || use_fill)) {
        std::cerr << "Error: --common is for --mode hit and replaces --ave" << std::endl;
        return EXIT_FAILURE;
    }
    const TString hist_fmt = (mode == "trk") ? "TPCTrk_ResY_Layer%02d_Row%03d" : "TPCHit_ResY_Layer%02d_Row%03d";
    // hit (BcOut reference) reaches only the pads on the beam path; --ave then shifts the other pads
    // of the same AsAd by the fitted pads' Δp0 and keeps their pad-to-pad differences.
    // trk reaches nearly all pads; --ave sets the few unfitted pads to the neighbours' p0.
    const Bool_t fill_by_shift = (mode == "hit");

    // 1. Open ROOT file early to potentially detect run number
    TFile* f = TFile::Open(input_root.Data());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open " << input_root << std::endl;
        return EXIT_FAILURE;
    }

    // Auto-detect run number if not provided
    if (run_number < 0) {
        // Try from TTree using TTreeReader (more reliable)
        TTree *tree = (TTree*)f->Get("tpc");
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
        
        // Fallback: Try from filename (e.g., run02570)
        if (run_number < 0) {
            TString filename = input_root;
            TPRegexp run_regex("run([0-9]+)");
            TString run_str = filename(run_regex);
            if (!run_str.IsNull()) {
                run_str.ReplaceAll("run", "");
                run_number = run_str.Atoi();
                std::cout << "Auto-detected run number from filename: " << run_number << std::endl;
            }
        }
    }

    if (run_number < 0) {
        std::cerr << "Error: run number not provided and could not be auto-detected." << std::endl;
        f->Close();
        return EXIT_FAILURE;
    }

    TString base_dir = ANALYZER_DIR + "/param/TPCPRM";
    TString e72_dir = base_dir + "/e72";
    TString run_param_path = Form("%s/TPCParam_e72_run%05d", e72_dir.Data(), run_number);
    TString base_param_path = Form("%s/TPCParam_example", base_dir.Data());
    // 2. Determine which parameter file to use as base
    TString input_param_path = run_param_path;
    if (!fs::exists(input_param_path.Data())) {
        input_param_path = base_param_path;
        std::cout << "Input param file not found in e72, using base: " << base_param_path << std::endl;
    } else {
        std::cout << "Using existing param file for update: " << input_param_path << std::endl;
    }

    // 3. Read parameter file
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
        
        // Trim leading whitespace
        std::string trimmed = line;
        size_t first = trimmed.find_first_not_of(" \t");
        if (std::string::npos == first) {
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
                    Double_t val;
                    while (ss >> val) {
                        entry.p.push_back(val);
                    }
                }
            }
        }
        entries.push_back(entry);
    }
    ifs.close();

    // 校正前 p0（aty==2）— ループ内で entry を更新するので先に退避
    std::map<std::pair<Int_t, Int_t>, Double_t> ini_p0;
    for (const auto& e : entries) {
        if (!e.is_comment && e.aty == 2 && !e.p.empty())
            ini_p0[{e.layer, e.row}] = e.p[0];
    }
    // センターフレーム上のパッドは校正対象外（p0=0、vdrift=0.055 に正規化）
    for (auto& e : entries) {
        if (e.is_comment || e.aty != 2 || e.p.empty()) continue;
        if (!tpc::IsPadOnCenterFrame(e.layer, e.row)) continue;
        e.p[0] = 0.0;
        if (e.p.size() >= 2) e.p[1] = 0.055;
        else e.p.push_back(0.055);
        ini_p0[{e.layer, e.row}] = 0.0;
    }
    std::set<std::pair<Int_t, Int_t>> pads_updated;

    // --- PDF Setup ---
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
    // Save QA PDF separately for hit/trk to avoid overwriting when running both.
    TString pdf_path = img_dir + Form("/tpc_time_offset_calib_%s.pdf", mode.Data());
    TCanvas *canvas = new TCanvas("c", "c", 1200, 1000);
    
    // TPC geometry summary map
    TH2Poly *hSummary = new TH2Poly("hSummary", Form("TPC Calib Summary (Run %05d);Z [mm];X [mm]", run_number), -300, 300, -300, 300);
    hSummary->SetDirectory(nullptr);
    tpc::InitializeHistograms(hSummary);

    TH2Poly* hP0Map = new TH2Poly(
        "hTPC_TimeOffset_p0",
        Form("TPC time offset p0 (Run %05d, %s);Z [mm];X [mm]", run_number, mode.Data()), -300, 300, -300, 300);
    hP0Map->SetDirectory(nullptr);
    tpc::InitializeHistograms(hP0Map);

    TH2Poly* hDeltaMap = new TH2Poly(
        "hTPC_TimeOffset_dP0",
        Form("TPC time offset #Delta p0 (%s, Run %05d, %s);Z [mm];X [mm]",
             use_fill ? "fitted + filled pads" : "fitted pads", run_number, mode.Data()), -300,
        300, -300, 300);
    hDeltaMap->SetDirectory(nullptr);
    tpc::InitializeHistograms(hDeltaMap);

    const Int_t nx = 4, ny = 4;
    Int_t ipad = 0;
    std::vector<TLine*> page_lines;
    canvas->Print(pdf_path + "["); // Open multi-page PDF

    std::cout << "Extracting offsets..." << std::endl;
    Int_t n_updated = 0;
    Int_t n_total_pads = 0;
    std::map<Int_t, std::pair<Double_t, Double_t>> layer_seed;  // layer -> (mean, sigma)

    // Layerごとの合算ヒストで事前fitし、pad fit の初期値 seed に使う
    for (Int_t layer = 0; layer < tpc::NumOfLayersTPC; ++layer) {
        const Int_t n_pad = static_cast<Int_t>(tpc::padParameter[layer][tpc::kNumOfPad]);
        TH1D* h_layer = nullptr;
        for (Int_t row = 0; row < n_pad; ++row) {
            TString hname = Form(hist_fmt.Data(), layer, row);
            TH1D* h_tmp = nullptr;
            f->GetObject(hname.Data(), h_tmp);
            if (!h_tmp) {
                TString hname_pre = "hist/" + hname;
                f->GetObject(hname_pre.Data(), h_tmp);
            }
            if (!h_tmp || h_tmp->GetEntries() <= 0)
                continue;
            if (!h_layer) {
                h_layer = (TH1D*)h_tmp->Clone(Form("hLayerPrefit_%02d", layer));
                h_layer->SetDirectory(nullptr);
            } else {
                h_layer->Add(h_tmp);
            }
        }
        if (!h_layer || h_layer->GetEntries() <= 0) {
            if (h_layer)
                delete h_layer;
            continue;
        }

        const Int_t b = DominantLocalMaximumBin(h_layer, local_peak_sep_bins);
        const Double_t px = h_layer->GetBinCenter(b);
        const Double_t py = TMath::Max(h_layer->GetBinContent(b), 1.0);
        Double_t sr = TMath::Max(h_layer->GetRMS(), 0.2);
        Double_t fmin = px - 2.0 * sr;
        Double_t fmax = px + 2.0 * sr;
        if (fmax - fmin < 8.0) {
            fmin = px - 4.0;
            fmax = px + 4.0;
        }

        TF1* f_seed = new TF1(Form("fLayerSeed_%02d", layer), "gaus(0)+pol0(3)", fmin, fmax);
        f_seed->SetParameters(py, px, TMath::Max(0.25, 0.8 * sr), 0.0);
        f_seed->SetParLimits(2, 0.05, 20.0);
        // "R": use TF1 range [fmin,fmax]; without it the full hist axis is fit.
        if (h_layer->Fit(f_seed, "QR") == 0) {
            const Double_t m = f_seed->GetParameter(1);
            const Double_t s = TMath::Abs(f_seed->GetParameter(2));
            if (std::isfinite(m) && std::isfinite(s) && s > 0.0) {
                // layer合算は広めに出るので、pad seed は少し縮める
                layer_seed[layer] = {m, TMath::Max(0.2, 0.8 * s)};
            }
        }
        delete f_seed;
        delete h_layer;
    }

    for (auto& entry : entries) {
        if (entry.is_comment || entry.aty != 2) continue;
        if (tpc::IsPadOnCenterFrame(entry.layer, entry.row)) continue;
        n_total_pads++;

        TString hname = Form(hist_fmt.Data(), entry.layer, entry.row);
        
        TH1D* h = nullptr;
        f->GetObject(hname.Data(), h);
        if (!h) {
            TString hname_pre = "hist/" + hname;
            f->GetObject(hname_pre.Data(), h);
        }

        if (h && h->GetEntries() > 0) {
            // 最も高い局所峰（二峰性で「でかい方」を優先。同高さは 0 に近い方）
            const Int_t max_bin = DominantLocalMaximumBin(h, local_peak_sep_bins);
            Double_t peak_x = h->GetBinCenter(max_bin);
            Double_t peak_y = h->GetBinContent(max_bin);
            if (peak_y <= 0.0) {
                peak_y = h->GetMaximum();
                peak_x = h->GetBinCenter(h->GetMaximumBin());
            }

            const Double_t h_rms_glob = TMath::Max(h->GetRMS(), 1e-6);
            Double_t h_rms_fit = h_rms_glob;
            {
                const Double_t loc = LocalWeightedRmsAroundPeak(h, peak_x, local_peak_half_mm);
                if (loc > 0.0)
                    h_rms_fit = loc;
            }
            if (h_rms_fit < 0.1)
                h_rms_fit = 0.5;
            if (h_rms_fit > 5.0)
                h_rms_fit = 2.5;

            // 事前fit seed があれば初期値に使い、fit範囲も seed 近傍に寄せる
            Double_t seed_mean = peak_x;
            Double_t seed_sigma = TMath::Max(0.2, 0.72 * h_rms_fit);
            auto it_seed = layer_seed.find(entry.layer);
            if (it_seed != layer_seed.end()) {
                seed_mean = it_seed->second.first;
                seed_sigma = TMath::Max(0.2, it_seed->second.second);
            }

            // fit範囲: 左 1σ, 右 2σ（必要に応じて h_rms_fit で下支え）
            const Double_t sigma_win = TMath::Max(seed_sigma, nsigma * h_rms_fit);
            Double_t fit_min = seed_mean - 1.0 * sigma_win;
            Double_t fit_max = seed_mean + 2.0 * sigma_win;

            // --- Statistics check: Count hits within +/- stat_range ---
            Int_t bin_min_full = h->FindBin(-stat_range);
            Int_t bin_max_full = h->FindBin(stat_range);
            Double_t n_fit = h->Integral(bin_min_full, bin_max_full);

            // threshold=0 では足切り無効（ピーク形状で判断）
            if (threshold > 0 && n_fit < threshold) {
                continue;
            }

            // Pre-fit flatness check: skip if no clear peak (avoid fitting flat/garbage histograms)
            Int_t n_bins_window = bin_max_full - bin_min_full + 1;
            Double_t mean_content = (n_bins_window > 0) ? (n_fit / n_bins_window) : 0;
            if (mean_content <= 0 || peak_y < mean_content * min_peak_ratio) {
                continue;  // do not fit, do not add to PDF
            }

            const Double_t sigma_floor_cfg = TMath::Max(min_sigma_mm, min_sigma_rel_rms * h_rms_glob);
            const Double_t sigma_lo = TMath::Max(0.05, 0.5 * sigma_floor_cfg);

            // Define fit function: Gaussian + Constant Background
            TF1 *fFit = new TF1("fFit", "gaus(0)+pol0(3)", fit_min, fit_max);
            fFit->SetParNames("Constant", "Mean", "Sigma", "Background");
            {
                const Double_t fit_span = fit_max - fit_min;
                // グローバル RMS は混ぜない（副ピークで初期 σ が膨らむのを防ぐ）
                Double_t sigma_init = sigma_init_rel * TMath::Max(seed_sigma, 0.12);
                sigma_init = TMath::Max(sigma_init, sigma_lo * 1.15);
                if (fit_span > 1e-6)
                    sigma_init = TMath::Min(sigma_init, 0.42 * fit_span);
                fFit->SetParameters(peak_y, seed_mean, sigma_init, 0.0);
                fFit->SetParLimits(2, sigma_lo, 20.0);
            }

            // "R": use TF1 range [fit_min,fit_max] so side peaks outside the window
            // do not pull the mean (acceptance check still uses this window width).
            if (h->Fit(fFit, "QR") == 0) {
                Double_t constant = fFit->GetParameter(0);
                Double_t mean     = fFit->GetParameter(1);
                Double_t sigma    = fFit->GetParameter(2);
                Double_t bkg      = fFit->GetParameter(3);
                Double_t mean_err = fFit->GetParError(1);

                // Filter logic:
                // 1. Sigma must be within a reasonable range (absolute and relative to fit range).
                // 2. Gaussian amplitude (constant) must be meaningful (signal clearly above bkg).
                // 3. Error on the mean (peak position) shouldn't be too large.
                Double_t fit_range = fit_max - fit_min;
                const Double_t sigma_floor = sigma_floor_cfg;

                const Int_t fndf = fFit->GetNDF();
                const Double_t fchi2 = fFit->GetChisquare();
                const Double_t chi2ndf = (fndf > 0) ? (fchi2 / static_cast<Double_t>(fndf)) : 0.0;
                const Double_t fit_prob = (fndf > 0) ? TMath::Prob(fchi2, fndf) : 0.0;
                // chi2/ndf and prob are only logged, not used for acceptance
                const Bool_t chi2_ok = (fndf >= min_fit_ndf);

                if (chi2_ok && sigma >= sigma_floor && sigma < sigma_threshold &&
                    sigma < 0.85 * fit_range &&  // reject effectively flat (peak too broad)
                    constant > bkg && constant > 1.0 &&
                    mean >= fit_min && mean <= fit_max &&
                    mean_err < 1.0) {

                    if (entry.p.size() >= 2 && TMath::Abs(entry.p[1]) > 1e-10) {
                        Double_t delta_p0 = mean / entry.p[1];
                        Double_t old_p0 = entry.p[0];
                        entry.p[0] += delta_p0;
                        n_updated++;
                        pads_updated.insert({entry.layer, entry.row});

                        Int_t bin_idx = tpc::GetPadId(entry.layer, entry.row) + 1;
                        hSummary->SetBinContent(bin_idx, 100.0);

                        std::cout << Form(
                            "  Updated: L=%2d R=%3d | Mean=%7.4f Sigma=%6.3f | chi2=%.2f ndf=%d chi2/ndf=%.3f "
                            "prob=%.4g -> delta_p0=%10.4f | p0: %12.4f -> %12.4f",
                            entry.layer, entry.row, mean, sigma, fchi2, fndf, chi2ndf, fit_prob, delta_p0,
                            old_p0, entry.p[0]) << std::endl;

                        // --- Draw each pad (Keep as vector PDF) ---
                        if (ipad == 0) {
                            canvas->Clear();
                            canvas->Divide(nx, ny);
                        }

                        canvas->cd(++ipad);
                        h->Draw();
                        fFit->SetLineColor(kRed);
                        fFit->DrawClone("same");

                        // Draw dashed line for peak position
                        TLine *vLine = new TLine(mean, h->GetMinimum(), mean, h->GetMaximum());
                        vLine->SetLineStyle(2);
                        vLine->SetLineColor(kBlue);
                        vLine->Draw("same");
                        page_lines.push_back(vLine);

                        if (ipad >= nx * ny) {
                            canvas->Print(pdf_path);
                            for (auto* ln : page_lines) delete ln;
                            page_lines.clear();
                            ipad = 0;
                        }
                    }
                }
            }
            delete fFit;
        }
    }

    // --- 未更新パッド補完（tpc_fill::FillUnfittedPads: 同一 AsAd の隣接 row 平均 → layer+AsAd の中央値） ---
    // fill_by_shift のときは p0 の代わりに Δp0 を埋め、それを今の p0 に足す。
    std::set<std::pair<Int_t, Int_t>> pads_filled;
    tpc_fill::FillResult fill_res;
    Int_t n_common_asad = 0, n_common_global = 0;
    Double_t common_global = 0.0;
    if (common_only) {
        std::map<Int_t, std::vector<Double_t>> asad_dp;
        std::vector<Double_t> all_fit_dp;
        for (const auto& e : entries) {
            if (e.is_comment || e.aty != 2 || e.p.empty()) continue;
            const auto key = std::make_pair(e.layer, e.row);
            if (!pads_updated.count(key)) continue;
            const Double_t dp = e.p[0] - ini_p0.at(key);
            asad_dp[tpc::GetASADId(e.layer, e.row)].push_back(dp);
            all_fit_dp.push_back(dp);
        }
        auto median = [](std::vector<Double_t> v) {
            if (v.empty()) return 0.0;
            std::sort(v.begin(), v.end());
            const size_t n = v.size();
            return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
        };
        common_global = median(all_fit_dp);
        std::map<Int_t, Double_t> asad_shift;
        for (const auto& kv : asad_dp) {
            if (static_cast<Int_t>(kv.second.size()) >= common_min_pads) asad_shift[kv.first] = median(kv.second);
        }
        for (const auto& kv : asad_shift) {
            std::cout << Form("  common shift: ASAD=%2d  dp0=%9.4f  (%zu fitted pads)", kv.first, kv.second,
                              asad_dp.at(kv.first).size()) << std::endl;
        }
        for (auto& e : entries) {
            if (e.is_comment || e.aty != 2 || e.p.empty()) continue;
            const auto key = std::make_pair(e.layer, e.row);
            const auto it = asad_shift.find(tpc::GetASADId(e.layer, e.row));
            const Bool_t by_asad = (it != asad_shift.end());
            e.p[0] = ini_p0.at(key) + (by_asad ? it->second : common_global);
            if (!pads_updated.count(key)) pads_filled.insert(key);
            if (by_asad) ++n_common_asad; else ++n_common_global;
        }
    } else if (use_fill) {
        std::map<tpc_fill::PadKey, Double_t> fill_values;
        for (const auto& e : entries) {
            if (e.is_comment || e.aty != 2 || e.p.empty()) continue;
            const auto key = std::make_pair(e.layer, e.row);
            fill_values[key] = fill_by_shift ? (pads_updated.count(key) ? e.p[0] - ini_p0.at(key) : 0.0) : e.p[0];
        }
        const std::set<tpc_fill::PadKey> fill_excluded;
        auto consistent = [](Double_t v, const tpc_fill::GroupStat& g) {
            return TMath::Abs(v - g.median) <= 5.0 * 1.4826 * std::max(g.mad, 1e-6);
        };
        fill_res = tpc_fill::FillUnfittedPads(fill_values, pads_updated, fill_excluded, consistent);

        for (auto& e : entries) {
            if (e.is_comment || e.aty != 2 || e.p.empty()) continue;
            const auto key = std::make_pair(e.layer, e.row);
            const auto it = fill_res.value.find(key);
            if (it == fill_res.value.end()) continue;

            const Bool_t is_neighbor = (fill_res.kind.at(key) == tpc_fill::FillKind::kNeighbor);
            const Double_t old_p0 = e.p[0];
            e.p[0] = fill_by_shift ? old_p0 + it->second : it->second;
            pads_filled.insert(key);

            const Int_t bin_idx = tpc::GetPadId(e.layer, e.row) + 1;
            hSummary->SetBinContent(bin_idx, is_neighbor ? 70.0 : 40.0);  // 100:fit, 70:neighbor mean, 40:AsAd median
            std::cout << Form("  %s %s: L=%2d R=%3d (ASAD=%2d) p0: %10.4f -> %10.4f",
                              is_neighbor ? "neighbor-mean" : "asad-median", fill_by_shift ? "shift" : "fill",
                              e.layer, e.row,
                              tpc::GetASADId(e.layer, e.row), old_p0, e.p[0]) << std::endl;
        }
    }

    // --- パッドごとの p0 / Δp0（PDF: 2D マップ + 1D 分布。ROOT ファイルは出力しない）---
    std::vector<Double_t> all_p0;
    std::vector<Double_t> all_dp;
    all_p0.reserve(1024);
    all_dp.reserve(static_cast<size_t>(n_updated) + 8);

    for (const auto& e : entries) {
        if (e.is_comment || e.aty != 2 || e.p.empty())
            continue;
        const Int_t bin_idx = tpc::GetPadId(e.layer, e.row) + 1;
        hP0Map->SetBinContent(bin_idx, e.p[0]);
        all_p0.push_back(e.p[0]);
        const auto key = std::make_pair(e.layer, e.row);
        if (pads_updated.count(key) || pads_filled.count(key)) {
            const Double_t dp = e.p[0] - ini_p0.at(key);
            hDeltaMap->SetBinContent(bin_idx, dp);
            all_dp.push_back(dp);
        }
    }

    auto axis_with_margin = [](Double_t lo, Double_t hi) {
        if (!(lo < hi)) {
            return std::pair<Double_t, Double_t>{lo - 1.0, hi + 1.0};
        }
        const Double_t m = 0.08 * (hi - lo);
        return std::pair<Double_t, Double_t>{lo - m, hi + m};
    };

    Double_t p0_lo = 0.0, p0_hi = 1.0;
    if (!all_p0.empty()) {
        const auto mm_p0 = std::minmax_element(all_p0.begin(), all_p0.end());
        const auto ax = axis_with_margin(*mm_p0.first, *mm_p0.second);
        p0_lo = ax.first;
        p0_hi = ax.second;
    }

    Int_t nb_p0 = std::max(16, static_cast<Int_t>(all_p0.size() / 6) + 1);
    nb_p0 = std::min(nb_p0, 300);
    TH1D* hP0Dist =
        new TH1D("hTPC_TimeOffset_p0_distrib",
                 Form("TPC time offset p0, all aty==2 pads (Run %05d);p0;pads / bin", run_number), nb_p0, p0_lo,
                 p0_hi);
    hP0Dist->SetDirectory(nullptr);
    for (Double_t v : all_p0)
        hP0Dist->Fill(v);

    Double_t dp_lo = -1.0, dp_hi = 1.0;
    if (!all_dp.empty()) {
        const auto mm_dp = std::minmax_element(all_dp.begin(), all_dp.end());
        const auto ax = axis_with_margin(*mm_dp.first, *mm_dp.second);
        dp_lo = ax.first;
        dp_hi = ax.second;
    }
    Int_t nb_dp = std::max(16, static_cast<Int_t>(all_dp.size() / 2) + 1);
    nb_dp = std::min(nb_dp, 300);
    TH1D* hDeltaDist =
        new TH1D("hTPC_TimeOffset_dp0_distrib",
                 Form("TPC time offset #Delta p0, %s (Run %05d);#Delta p0;pads / bin",
                      use_fill ? "fitted + filled pads" : "fitted pads only", run_number),
                 nb_dp, dp_lo, dp_hi);
    hDeltaDist->SetDirectory(nullptr);
    for (Double_t v : all_dp)
        hDeltaDist->Fill(v);

    if (ipad > 0) {
        canvas->Print(pdf_path);
        for (auto* ln : page_lines) delete ln;
        page_lines.clear();
    }

    gStyle->SetOptStat(1111);
    canvas->Clear();
    canvas->Divide(2, 1);
    canvas->cd(1);
    hP0Dist->Draw("HIST");
    canvas->cd(2);
    hDeltaDist->Draw("HIST");
    canvas->Print(pdf_path);

    canvas->Print(pdf_path + "]"); // Close vector PDF
    // delete canvas; // Removed duplicate delete (already at end or should be done there)

    // --- High Resolution Rasterization for TH2Poly ---
    std::cout << "Generating high-resolution rasterized maps..." << std::endl;
    // Set size to match the vector pages (1200x1000) for consistency
    TCanvas *cHR = new TCanvas("cHR", "High Res", 1200, 1000);
    gStyle->SetOptStat(0);
    std::vector<TString> tmp_pdfs;

    auto ProcessRaster = [&](TH2Poly* hp, TString tag) {
        if (!hp) return;
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
        TString png = img_dir + "/tmp_" + tag + ".png";
        TString tmp_pdf = img_dir + "/tmp_" + tag + ".pdf";
        cHR->Print(png);
        // Use convert to resize and fit into a standard page size, matching the vector plots
        // -density 72 -page 1200x1000 ensures the PDF page size matches the canvas
        system(Form("convert %s -resize 1200x1000 -gravity center -extent 1200x1000 %s", png.Data(), tmp_pdf.Data()));
        tmp_pdfs.push_back(tmp_pdf);
        system(Form("rm %s", png.Data()));
    };

    hSummary->SetMinimum(0);
    hSummary->SetMaximum(110);
    ProcessRaster(hSummary, "summary");

    ProcessRaster(hP0Map, "timeoffset_p0");
    ProcessRaster(hDeltaMap, "timeoffset_dp0");

    TH2Poly* hHitPat = nullptr;
    f->GetObject("TPC_HitPat", hHitPat);
    if (!hHitPat) f->GetObject("hist/TPC_HitPat", hHitPat);
    if (hHitPat) {
        hHitPat->SetTitle(Form("TPC HitPat (Run %05d)", run_number));
        ProcessRaster(hHitPat, "hitpat");
    }

    // Merge vector PDF and rasterized TH2Poly maps
    if (!tmp_pdfs.empty()) {
        TString unite_cmd = "pdfunite " + pdf_path + " ";
        for (const auto& p : tmp_pdfs) unite_cmd += p + " ";
        unite_cmd += pdf_path + "_merged.pdf";
        
        Int_t res = system(unite_cmd.Data());
        if (res == 0) {
            system(Form("mv %s_merged.pdf %s", pdf_path.Data(), pdf_path.Data()));
            std::cout << "Successfully merged vector and raster maps with aligned sizes." << std::endl;
        }
        for (const auto& p : tmp_pdfs) system(Form("rm %s", p.Data()));
    }

    delete cHR;

    // 4. Write updated parameter file
    if (!debug_mode) {
        if (!fs::exists(e72_dir.Data())) fs::create_directories(e72_dir.Data());
        
        // Save previous version as backup if it existed
        if (fs::exists(run_param_path.Data())) {
            fs::copy(run_param_path.Data(), (run_param_path + ".bak").Data(), fs::copy_options::overwrite_existing);
        }

        std::ofstream ofs(run_param_path.Data());
        if (!ofs.is_open()) {
            std::cerr << "Error: cannot open output file " << run_param_path << std::endl;
            if (f) f->Close();
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

    std::cout << "Done! Updated " << n_updated << " out of " << n_total_pads << " pads." << std::endl;
    if (use_fill) {
        std::cout << (fill_by_shift ? "Shifted pads (by fitted dp0): neighbor-mean " : "Filled pads: neighbor-mean ") << fill_res.n_neighbor << ", asad-median "
                  << fill_res.n_asad_median << ", left unfilled " << fill_res.n_left << std::endl;
    }
    if (common_only) {
        std::cout << "Common shift only: " << n_common_asad << " pads by their AsAd median, " << n_common_global
                  << " pads by the global median " << common_global << std::endl;
    }
    if (debug_mode) {
        std::cout << "Debug mode: parameter file was NOT updated." << std::endl;
    } else {
        std::cout << "Result saved to: " << run_param_path << std::endl;
    }
    std::cout << "QA PDF saved to: " << pdf_path << std::endl;

    if (f) f->Close();
    delete hP0Dist;
    delete hDeltaDist;
    delete hP0Map;
    delete hDeltaMap;
    if (hSummary) delete hSummary;
    if (canvas) delete canvas;
    return EXIT_SUCCESS;
}
