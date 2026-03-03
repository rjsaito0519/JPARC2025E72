#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <algorithm>
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

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
int main(Int_t argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.root> [run_number] [--threshold N] [--sigma N] [--nsigma N] [--stat-range N] [--hist-fmt fmt]" << std::endl;
        return EXIT_FAILURE;
    }

    TString input_root = argv[1];
    Int_t run_number = -1;
    Int_t threshold = 2500;
    Double_t sigma_threshold = 10.0;
    Double_t nsigma = 2.0;
    Double_t stat_range = 20.0;
    TString hist_fmt = "TPCHit_ResY_Layer%02d_Row%03d";

    Int_t arg_idx = 2;
    if (arg_idx < argc && TString(argv[arg_idx]).IsDigit()) {
        run_number = std::atoi(argv[arg_idx++]);
    }

    for (Int_t i = arg_idx; i < argc; ++i) {
        TString arg = argv[i];
        if (arg == "--threshold" && i + 1 < argc) {
            threshold = std::atoi(argv[++i]);
        } else if (arg == "--sigma" && i + 1 < argc) {
            sigma_threshold = std::atof(argv[++i]);
        } else if (arg == "--nsigma" && i + 1 < argc) {
            nsigma = std::atof(argv[++i]);
        } else if (arg == "--stat-range" && i + 1 < argc) {
            stat_range = std::atof(argv[++i]);
        } else if (arg == "--hist-fmt" && i + 1 < argc) {
            hist_fmt = argv[++i];
        }
    }

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
    TString base_param_path = Form("%s/TPCParam_0_yoffset_adjusted", base_dir.Data());
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
    while (std::getline(ifs, line)) {
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
    // --- PDF Setup ---
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    TString img_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_number);
    TString pdf_path = img_dir + "/tpc_time_offset_calib.pdf";
    TCanvas *canvas = new TCanvas("c", "c", 1200, 1000);
    
    // TPC geometry summary map
    TH2Poly *hSummary = new TH2Poly("hSummary", Form("TPC Calib Summary (Run %05d);Z [mm];X [mm]", run_number), -300, 300, -300, 300);
    hSummary->SetDirectory(nullptr);
    tpc::InitializeHistograms(hSummary);
    
    const Int_t nx = 4, ny = 4;
    Int_t ipad = 0;
    canvas->Print(pdf_path + "["); // Open multi-page PDF

    std::cout << "Extracting offsets..." << std::endl;
    Int_t n_updated = 0;
    Int_t n_total_pads = 0;

    for (auto& entry : entries) {
        if (entry.is_comment || entry.aty != 2) continue;
        n_total_pads++;

        TString hname = Form(hist_fmt.Data(), entry.layer, entry.row);
        
        TH1D* h = nullptr;
        f->GetObject(hname.Data(), h);
        if (!h) {
            TString hname_pre = "hist/" + hname;
            f->GetObject(hname_pre.Data(), h);
        }

        if (h && h->GetEntries() > 0) {
            // Find peak for better initial parameters
            Int_t max_bin = h->GetMaximumBin();
            Double_t peak_x = h->GetBinCenter(max_bin);
            Double_t peak_y = h->GetMaximum();
            Double_t h_rms  = h->GetRMS();
            if (h_rms < 0.1) h_rms = 0.5;
            if (h_rms > 5.0) h_rms = 2.5;
            
            // Range centered on peak
            Double_t fit_min = peak_x - nsigma * h_rms;
            Double_t fit_max = peak_x + nsigma * h_rms;

            // --- Statistics check: Count hits within +/- stat_range ---
            Int_t bin_min_full = h->FindBin(-stat_range);
            Int_t bin_max_full = h->FindBin(stat_range);
            Double_t n_fit = h->Integral(bin_min_full, bin_max_full);

            if (n_fit >= threshold) {
                // Define fit function: Gaussian + Constant Background
                TF1 *fFit = new TF1("fFit", "gaus(0)+pol0(3)", fit_min, fit_max);
                fFit->SetParameters(peak_y, peak_x, h_rms * 0.5, 0.0);
                fFit->SetParNames("Constant", "Mean", "Sigma", "Background");
                fFit->SetParLimits(2, 0.01, 20.0); // Constraints for stability

                if (h->Fit(fFit, "Q") == 0) {
                    Double_t constant = fFit->GetParameter(0);
                    Double_t mean     = fFit->GetParameter(1);
                    Double_t sigma    = fFit->GetParameter(2);
                    Double_t bkg      = fFit->GetParameter(3);
                    Double_t mean_err = fFit->GetParError(1);

                    // Filter logic:
                    // 1. Sigma must be within a reasonable range.
                    // 2. Gaussian amplitude (constant) must be meaningful (signal should be clearly visible over bkg).
                    // 3. Error on the mean (peak position) shouldn't be too large.
                    if (sigma > 0.05 && sigma < sigma_threshold && 
                        constant > 0.1 * bkg && constant > 1.0 && 
                        mean_err < 1.0) {
                        
                        if (entry.p.size() >= 2 && TMath::Abs(entry.p[1]) > 1e-10) {
                            Double_t delta_p0 = mean / entry.p[1];
                            Double_t old_p0 = entry.p[0];
                            entry.p[0] += delta_p0;
                            n_updated++;

                            Int_t bin_idx = tpc::GetPadId(entry.layer, entry.row) + 1;
                            hSummary->SetBinContent(bin_idx, 100.0);

                            std::cout << Form("  Updated: L=%2d R=%3d | Mean=%7.4f Sigma=%6.3f -> delta_p0=%10.4f | p0: %12.4f -> %12.4f", 
                                              entry.layer, entry.row, mean, sigma, delta_p0, old_p0, entry.p[0]) << std::endl;

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
                        vLine->DrawClone();
                        delete vLine;

                        if (ipad >= nx * ny) {
                            canvas->Print(pdf_path);
                            ipad = 0;
                        }
                    }
                }
            }
            delete fFit;
            } // end if (n_fit >= threshold)
        }
    }

    if (ipad > 0) canvas->Print(pdf_path);
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

    std::cout << "Done! Updated " << n_updated << " out of " << n_total_pads << " pads." << std::endl;
    std::cout << "Result saved to: " << run_param_path << std::endl;
    std::cout << "QA PDF saved to: " << pdf_path << std::endl;

    if (f) f->Close();
    if (hSummary) delete hSummary;
    if (canvas) delete canvas;
    return EXIT_SUCCESS;
}
