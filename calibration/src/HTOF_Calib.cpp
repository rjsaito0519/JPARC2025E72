// c++
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <unordered_map>

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLatex.h>

// Custom headers
#include "config.h"
#include "ana_helper.h"
#include "paths.h"
#include "params.h"

Config& conf = Config::getInstance();

// HTOF calibration in one pass (TENTATIVE, 2026-09-24): TDC, ADC (pedestal/MIP) and PHC.
//
// Inputs: one or more (hodo_path, dst_path) pairs; several pairs merge runs. Output is
// labeled with the first dst_path's run number (representative run).
//
// - TDC : all-event HTOF_TDC_seg{i}{U|D|S} histograms of run{run}_Hodo.root (summed over
//         runs), ana_helper::tdc_fit. Same as the former HTOF_HDPRM.
// - ADC : pedestal from the all-event HTOF_ADC_seg{i}{U|D|S} histograms of Hodo.root,
//         MIP from a TPC dE/dx pion-tagged histogram filled from the DstTPCHelixHTOF tree
//         (pion bit set, proton bit not set; kaon-ambiguous kept), ana_helper::htof_adc_fit.
//         Same as the former HTOF_ADC.
// - PHC : (two stages, see phc_fix_mask) TOF vs DeltaE filled here from the Hodo.root "hodo" tree, with DeltaE recomputed
//         from the raw ADC using the pedestal/MIP just obtained above
//         (dE = (adc - ped) / (mip - ped)), so no re-decoding is needed between the ADC and
//         PHC steps. TOF = time0 - htof_time_{u|d} (time0 = CTime0 of the event's T0
//         cluster; htof_time = decode-time, pre-PHC hit time). Only hits in segments
//         matched (in the same event, joined by event_number with the DST) to a TPC track
//         passing the same PID cut as the ADC MIP (pion bit set, proton bit not set), i.e.
//         protons excluded. No beam-flag selection (HTOF uses TPC PID, not the beam PID). ana_helper::htof_phc_fit (fit range and
//         profile restricted to the correlated band; same function as the analyzer's Type1).
//         Note: htof_time uses the decode-time TDC parameters, not the TDC values fitted
//         here; this is expected to converge as the calibration is iterated.
//
// Outputs (same tree formats as the former separate tools, so update_param.py reads them
// unchanged):
//   run{run}_HTOF_HDPRM_All.root  (tdc_p0_val/err)
//   run{run}_HTOF_ADC_{particle}.root (adc_p0/p1_val/err, adc_flag/chi2/ndf + histograms)
//   run{run}_HTOF_PHC_{particle}.root (p0/p1/p2_val/err)
//   run{run}_HTOF_Calib_{particle}.pdf: one page per segment, rows U/D/S, columns
//     ADC raw | ADC selected | ADC selected close-up | TDC | PHC (S row, last column: summary)

static TFile* open_output(Int_t run_num, const char* name) {
    TString root_base_dir = ana_helper::get_root_dir(OUTPUT_DIR, run_num);
    TString path = Form("%s/run%05d_HTOF_%s.root", root_base_dir.Data(), run_num, name);
    if (std::ifstream(path.Data())) std::remove(path.Data());
    return new TFile(path.Data(), "RECREATE");
}

void analyze(std::vector<TString> hodo_paths, std::vector<TString> dst_paths, TString particle) {
    gROOT->GetColor(kBlue)->SetRGB(0.12156862745098039, 0.4666666666666667, 0.7058823529411765);
    gROOT->GetColor(kOrange)->SetRGB(1.0, 0.4980392156862745, 0.054901960784313725);
    gROOT->GetColor(kGreen)->SetRGB(44.0/256, 160.0/256, 44.0/256);

    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.06, "XY");
    gStyle->SetTitleSize(0.06, "x");
    gStyle->SetTitleSize(0.06, "y");
    gStyle->SetTitleFontSize(0.06);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gROOT->GetColor(0)->SetAlpha(0.01);

    // +------------+
    // | load files |
    // +------------+
    if (hodo_paths.empty() || hodo_paths.size() != dst_paths.size()) {
        std::cerr << "Error: hodo_paths/dst_paths must be non-empty and of equal length" << std::endl;
        return;
    }
    const std::size_t n_runs = hodo_paths.size();

    std::vector<TFile*> f_hodo_list(n_runs), f_dst_list(n_runs);
    for (std::size_t k = 0; k < n_runs; k++) {
        f_hodo_list[k] = new TFile(hodo_paths[k].Data());
        if (!f_hodo_list[k] || f_hodo_list[k]->IsZombie()) {
            std::cerr << "Error: Could not open file : " << hodo_paths[k] << std::endl;
            return;
        }
        f_dst_list[k] = new TFile(dst_paths[k].Data());
        if (!f_dst_list[k] || f_dst_list[k]->IsZombie()) {
            std::cerr << "Error: Could not open file : " << dst_paths[k] << std::endl;
            return;
        }
    }

    // Representative run number for output naming/dirs = first dst_paths entry.
    Int_t run_num = -1;
    for (std::size_t k = 0; k < n_runs; k++) {
        TTreeReader r0("htof", f_dst_list[k]);
        TTreeReaderValue<unsigned int> rn(r0, "run_number");
        r0.SetEntry(0);
        std::cout << "#D HTOF_Calib: input " << k << ": run" << *rn
                   << " (" << hodo_paths[k] << ", " << dst_paths[k] << ")" << std::endl;
        if (k == 0) run_num = *rn;
    }
    conf.run_num = run_num;
    if (n_runs > 1) {
        std::cout << "#D HTOF_Calib: merged " << n_runs
                   << " runs; output labeled with representative run" << run_num << std::endl;
    }

    const Int_t n_seg = conf.num_of_ch.at("htof");
    const char* side_name[3] = {"U", "D", "S"};

    // Sum an all-event per-seg/side histogram over all hodo files.
    auto merge_hist = [&](const TString& name, const TString& new_name) -> TH1* {
        TH1* h_sum = nullptr;
        for (std::size_t k = 0; k < n_runs; k++) {
            auto* h = (TH1*)f_hodo_list[k]->Get(name);
            if (!h) {
                std::cerr << "Error: Missing " << name << " in " << hodo_paths[k] << std::endl;
                return nullptr;
            }
            if (k == 0) {
                h_sum = (TH1*)h->Clone(new_name);
                h_sum->SetDirectory(nullptr);
            } else {
                h_sum->Add(h);
            }
        }
        return h_sum;
    };

    // +-----+
    // | TDC |
    // +-----+
    TH1D *h_tdc[3][n_seg];
    for (Int_t s = 0; s < 3; s++) {
        for (Int_t i = 0; i < n_seg; i++) {
            h_tdc[s][i] = (TH1D*)merge_hist(Form("HTOF_TDC_seg%d%s", i, side_name[s]),
                                            Form("HTOF_TDC_seg%d%s_merged", i, side_name[s]));
            if (!h_tdc[s][i]) return;
        }
    }
    // search ranges: individual (U+D) and sum (S) channels
    TH1D *h_sum_tdc[2];
    h_sum_tdc[0] = (TH1D*)h_tdc[0][0]->Clone("h_sum_tdc_indiv");
    h_sum_tdc[0]->Reset();
    h_sum_tdc[1] = (TH1D*)h_tdc[2][0]->Clone("h_sum_tdc_sum");
    h_sum_tdc[1]->Reset();
    for (Int_t i = 0; i < n_seg; i++) {
        h_sum_tdc[0]->Add(h_tdc[0][i]);
        h_sum_tdc[0]->Add(h_tdc[1][i]);
        h_sum_tdc[1]->Add(h_tdc[2][i]);
    }

    // +-----+
    // | ADC |
    // +-----+
    // "Raw" = all-event histograms (pedestal); "selected" = TPC dE/dx pion-tagged (MIP)
    TH1D *h_adc_raw[3][n_seg];
    TH1D *h_adc_selected[3][n_seg];
    for (Int_t s = 0; s < 3; s++) {
        for (Int_t i = 0; i < n_seg; i++) {
            h_adc_raw[s][i] = (TH1D*)merge_hist(Form("HTOF_ADC_seg%d%s", i, side_name[s]),
                                                Form("HTOF_ADC_seg%d%s_merged", i, side_name[s]));
            if (!h_adc_raw[s][i]) return;
            h_adc_selected[s][i] = new TH1D(
                Form("HTOF_ADC_seg%d%s_%s", i, side_name[s], particle.Data()),
                Form("HTOF ADC seg%d %s (%s);ADC;Counts", i, side_name[s], particle.Data()),
                conf.adc_bin_num, conf.adc_min, conf.adc_max);
            h_adc_selected[s][i]->SetDirectory(nullptr);
        }
    }

    // per run: event_number -> bitmask of HTOF segments matched to a pion-tagged,
    // non-proton TPC track (used to select the PHC hits below)
    std::vector<std::unordered_map<UInt_t, ULong64_t>> sel_seg_mask(n_runs);
    Long64_t n_entries = 0, n_matched = 0, n_pi_htof = 0;
    for (std::size_t k = 0; k < n_runs; k++) {
        TTreeReader reader("htof", f_dst_list[k]);
        TTreeReaderValue<UInt_t> event_number(reader, "event_number");
        TTreeReaderValue<std::vector<Int_t>> pid(reader, "pid");
        TTreeReaderValue<std::vector<Int_t>> match_ok(reader, "match_ok");
        TTreeReaderValue<std::vector<Double_t>> htof_seg(reader, "htof_seg");
        TTreeReaderValue<std::vector<Double_t>> htof_adc_u(reader, "htof_adc_u");
        TTreeReaderValue<std::vector<Double_t>> htof_adc_d(reader, "htof_adc_d");
        TTreeReaderValue<std::vector<Double_t>> htof_adc_s(reader, "htof_adc_s");
        while (reader.Next()) {
            n_entries++;
            const auto& v_pid   = *pid;
            const auto& v_ok    = *match_ok;
            const auto& v_seg   = *htof_seg;
            const auto& v_adc_u = *htof_adc_u;
            const auto& v_adc_d = *htof_adc_d;
            const auto& v_adc_s = *htof_adc_s;
            for (std::size_t it = 0, n_tr = v_pid.size(); it < n_tr; it++) {
                if (it >= v_ok.size() || v_ok[it] != 1) continue;       // valid HTOF extrapolation+cluster match
                if (it >= v_seg.size() || TMath::IsNaN(v_seg[it])) continue;
                const Int_t seg = static_cast<Int_t>(std::lround(v_seg[it]));
                if (seg < 0 || seg >= n_seg) continue;
                n_matched++;

                // HypTPCdEdxPID bitmask: bit0=pi, bit1=K, bit2=proton.
                // Require pion bit set and proton bit NOT set; kaon-ambiguous (pid==3) is kept.
                // (Tried exclusive pid==1 to see if it removed a double-peak structure seen
                // in some segments' ADC distributions; it didn't, so reverted to this.)
                if (!(v_pid[it] & 1) || (v_pid[it] & 4)) continue;
                n_pi_htof++;
                sel_seg_mask[k][*event_number] |= (1ULL << seg);
                if (it < v_adc_u.size() && !TMath::IsNaN(v_adc_u[it])) h_adc_selected[0][seg]->Fill(v_adc_u[it]);
                if (it < v_adc_d.size() && !TMath::IsNaN(v_adc_d[it])) h_adc_selected[1][seg]->Fill(v_adc_d[it]);
                if (it < v_adc_s.size() && !TMath::IsNaN(v_adc_s[it])) h_adc_selected[2][seg]->Fill(v_adc_s[it]);
            }
        }
    }
    std::cout << "#D HTOF_Calib: " << n_entries << " DST events, " << n_matched << " HTOF-matched hits, "
              << n_pi_htof << " pion-tagged (proton excluded) among them" << std::endl;

    // +---------------------------------+
    // | fit TDC and ADC, one canvas/seg |
    // +---------------------------------+
    // pads: row r (U/D/S) -> 6r+1 ADC raw, 6r+2 ADC selected, 6r+3 ADC close-up, 6r+4 TDC,
    //       6r+5 PHC (all hits, stage 1), 6r+6 PHC (proton excluded, stage 2 = final)
    const Int_t rows = 3, cols = 6;
    std::vector<TCanvas*> c_seg(n_seg);
    std::vector<FitResult> tdc_result[3], adc_result[3];
    for (Int_t i = 0; i < n_seg; i++) {
        c_seg[i] = new TCanvas(Form("htof_calib_seg%d", i), "", 3600, 1500);
        c_seg[i]->Divide(cols, rows);
        for (Int_t s = 0; s < 3; s++) {
            const Int_t pad0 = cols*s + 1;
            // rebin is chosen automatically by htof_adc_fit (n_rebin=0)
            adc_result[s].push_back(ana_helper::htof_adc_fit(h_adc_raw[s][i], h_adc_selected[s][i], c_seg[i], pad0, 0));
            ana_helper::set_tdc_search_range(s < 2 ? h_sum_tdc[0] : h_sum_tdc[1]);
            tdc_result[s].push_back(ana_helper::tdc_fit(h_tdc[s][i], c_seg[i], pad0+3));
        }
    }

    // +-----------------------------------------+
    // | PHC: TOF vs DeltaE with the new ADC fit |
    // +-----------------------------------------+
    // Binning taken from the analyzer's all-event HTOF_seg0U_TOF_vs_DeltaE histogram.
    auto* h_tmpl = (TH2D*)f_hodo_list[0]->Get("HTOF_seg0U_TOF_vs_DeltaE");
    if (!h_tmpl) {
        std::cerr << "Error: Missing HTOF_seg0U_TOF_vs_DeltaE (binning template) in " << hodo_paths[0] << std::endl;
        return;
    }
    // h_phc_all: all HTOF hits (stage 1, high statistics / wide dE range)
    // h_phc    : hits in segments matched to a pion-tagged, non-proton track (stage 2, final)
    TH2D *h_phc_all[2][n_seg], *h_phc[2][n_seg];
    for (Int_t s = 0; s < 2; s++) {
        for (Int_t i = 0; i < n_seg; i++) {
            h_phc_all[s][i] = (TH2D*)h_tmpl->Clone(Form("HTOF_seg%d%s_TOF_vs_DeltaE_calib_all", i, side_name[s]));
            h_phc_all[s][i]->Reset();
            h_phc_all[s][i]->SetDirectory(nullptr);
            h_phc_all[s][i]->SetTitle(Form("HTOF seg%d %s TOF vs #DeltaE (all hits);#DeltaE [MIP];time0 - time [ns]", i, side_name[s]));
            h_phc[s][i] = (TH2D*)h_tmpl->Clone(Form("HTOF_seg%d%s_TOF_vs_DeltaE_calib", i, side_name[s]));
            h_phc[s][i]->Reset();
            h_phc[s][i]->SetDirectory(nullptr);
            h_phc[s][i]->SetTitle(Form("HTOF seg%d %s TOF vs #DeltaE (proton excluded);#DeltaE [MIP];time0 - time [ns]", i, side_name[s]));
        }
    }

    Long64_t n_hodo_events = 0, n_phc_fill_all = 0, n_phc_fill = 0;
    for (std::size_t k = 0; k < n_runs; k++) {
        TTreeReader reader("hodo", f_hodo_list[k]);
        TTreeReaderValue<UInt_t> event_number(reader, "event_number");
        TTreeReaderValue<Double_t> time0(reader, "time0");
        TTreeReaderValue<std::vector<Double_t>> raw_seg(reader, "htof_raw_seg");
        TTreeReaderValue<std::vector<Double_t>> adc_u(reader, "htof_adc_u");
        TTreeReaderValue<std::vector<Double_t>> adc_d(reader, "htof_adc_d");
        TTreeReaderValue<std::vector<Double_t>> hit_seg(reader, "htof_hit_seg");
        TTreeReaderValue<std::vector<std::vector<Double_t>>> time_u(reader, "htof_time_u");
        TTreeReaderValue<std::vector<std::vector<Double_t>>> time_d(reader, "htof_time_d");
        while (reader.Next()) {
            n_hodo_events++;
            const auto it_mask = sel_seg_mask[k].find(*event_number);
            // 0 if the event has no pion-tagged non-proton track
            const ULong64_t mask = (it_mask == sel_seg_mask[k].end()) ? 0ULL : it_mask->second;
            const Double_t t0 = *time0;
            if (TMath::IsNaN(t0)) continue;
            const auto& v_raw_seg = *raw_seg;
            const auto& v_hit_seg = *hit_seg;
            for (std::size_t ih = 0; ih < v_hit_seg.size(); ih++) {
                const Int_t seg = static_cast<Int_t>(std::lround(v_hit_seg[ih]));
                if (seg < 0 || seg >= n_seg) continue;
                const Bool_t is_sel = (mask & (1ULL << seg)); // segment matched to such a track
                // raw ADC of the same segment (adc_u/d are indexed like htof_raw_seg)
                Int_t ir = -1;
                for (std::size_t j = 0; j < v_raw_seg.size(); j++) {
                    if (static_cast<Int_t>(std::lround(v_raw_seg[j])) == seg) { ir = static_cast<Int_t>(j); break; }
                }
                if (ir < 0) continue;
                const Double_t adc[2] = { (*adc_u)[ir], (*adc_d)[ir] };
                const std::vector<Double_t>* times[2] = { &(*time_u)[ih], &(*time_d)[ih] };
                for (Int_t s = 0; s < 2; s++) {
                    const Double_t ped = adc_result[s][seg].par[1];
                    const Double_t mip = adc_result[s][seg].par[4];
                    if (TMath::IsNaN(adc[s]) || mip - ped <= 0.0) continue;
                    const Double_t de = (adc[s] - ped) / (mip - ped);
                    for (const Double_t t : *times[s]) {
                        if (TMath::IsNaN(t)) continue;
                        h_phc_all[s][seg]->Fill(de, t0 - t);
                        n_phc_fill_all++;
                        if (is_sel) {
                            h_phc[s][seg]->Fill(de, t0 - t);
                            n_phc_fill++;
                        }
                    }
                }
            }
        }
    }
    std::cout << "#D HTOF_Calib: " << n_hodo_events << " hodo events, "
              << n_phc_fill_all << " (hit x multi-hit) PHC entries for all hits, "
              << n_phc_fill << " for the proton-excluded selection" << std::endl;

    // Two stages: (1) all hits, (2) proton-excluded hits with the parameters in
    // phc_fix_mask fixed to the stage-1 values (the proton-excluded sample alone is
    // concentrated near 1 MIP, which leaves the curvature (p0 vs p1) unconstrained).
    // bit i -> p_i; TENTATIVE choice: fix p1 (pole position) only.
    const Int_t phc_fix_mask = (1 << 1);
    // fit-range threshold: stage 2 uses a lower one to recover statistics (curvature is
    // constrained by the stage-1 fixed parameters anyway)
    const Double_t phc_dense_ratio_all = 0.10, phc_dense_ratio_sel = 0.075;
    std::vector<FitResult> phc_result_all[2], phc_result[2];
    for (Int_t i = 0; i < n_seg; i++) {
        for (Int_t s = 0; s < 2; s++) {
            // HTOF-specific: ridge-following fit range / band-restricted profile (see hodo.cpp)
            phc_result_all[s].push_back(ana_helper::htof_phc_fit(h_phc_all[s][i], c_seg[i], cols*s + 5,
                                                                 nullptr, 0, phc_dense_ratio_all));
            phc_result[s].push_back(ana_helper::htof_phc_fit(h_phc[s][i], c_seg[i], cols*s + 6,
                                                             &phc_result_all[s].back(), phc_fix_mask,
                                                             phc_dense_ratio_sel));
        }
        // S row, PHC columns: value summary
        c_seg[i]->cd(cols*2 + 5);
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.06);
        tx.DrawLatex(0.08, 0.90, Form("HTOF seg%d (run%05d%s)", i, run_num, n_runs > 1 ? Form(", %zu runs", n_runs) : ""));
        for (Int_t s = 0; s < 3; s++) {
            const FitResult& a = adc_result[s][i];
            tx.DrawLatex(0.08, 0.78 - 0.11*s, Form("%s: ped %.1f  MIP %.1f  flag %d", side_name[s],
                                                  a.par[1], a.par[4], static_cast<Int_t>(a.additional[0])));
        }
        for (Int_t s = 0; s < 3; s++) {
            tx.DrawLatex(0.08, 0.42 - 0.11*s, Form("%s: TDC %.1f", side_name[s], tdc_result[s][i].par[1]));
        }
        c_seg[i]->cd(cols*2 + 6);
        TLatex tp;
        tp.SetNDC();
        tp.SetTextSize(0.055);
        tp.DrawLatex(0.05, 0.90, "PHC  (time0-time) = -p0/#sqrt{|dE-p1|} + p2");
        for (Int_t s = 0; s < 2; s++) {
            const FitResult& ra = phc_result_all[s][i];
            const FitResult& r  = phc_result[s][i];
            tp.DrawLatex(0.05, 0.76 - 0.30*s, Form("%s all : %.3f %.3f %.3f", side_name[s], ra.par[0], ra.par[1], ra.par[2]));
            tp.DrawLatex(0.05, 0.66 - 0.30*s, Form("%s final: %.3f %.3f %.3f", side_name[s], r.par[0], r.par[1], r.par[2]));
        }
        tp.DrawLatex(0.05, 0.12, Form("final: p_i fixed to 'all' for bits of mask %d; range thr. %.1f%% / %.1f%%",
                                      phc_fix_mask, 100*phc_dense_ratio_all, 100*phc_dense_ratio_sel));
    }

    // +-----+
    // | PDF |
    // +-----+
    TString img_base_dir = ana_helper::get_img_dir(OUTPUT_DIR, run_num);
    TString pdf_path = Form("%s/run%05d_HTOF_Calib_%s.pdf", img_base_dir.Data(), run_num, particle.Data());
    c_seg[0]->Print(pdf_path + "[");
    for (Int_t i = 0; i < n_seg; i++) c_seg[i]->Print(pdf_path);
    c_seg[0]->Print(pdf_path + "]");

    // +-------+
    // | Write |
    // +-------+
    // -- TDC -----
    {
        TFile* fout = open_output(run_num, "HDPRM_All");
        TTree* tree = new TTree("tree", "");
        Int_t ch;
        std::vector<Double_t> tdc_p0_val, tdc_p0_err;
        tree->Branch("ch", &ch, "ch/I");
        tree->Branch("tdc_p0_val", &tdc_p0_val);
        tree->Branch("tdc_p0_err", &tdc_p0_err);
        for (Int_t i = 0; i < n_seg; i++) {
            ch = i;
            tdc_p0_val.clear(); tdc_p0_err.clear();
            for (Int_t s = 0; s < 3; s++) {
                tdc_p0_val.push_back(tdc_result[s][i].par[1]);
                tdc_p0_err.push_back(tdc_result[s][i].err[1]);
            }
            tree->Fill();
        }
        tree->Write();
        fout->Close();
    }
    // -- ADC -----
    {
        TFile* fout = open_output(run_num, Form("ADC_%s", particle.Data()));
        TTree* tree = new TTree("tree", "");
        Int_t ch;
        std::vector<Double_t> adc_p0_val, adc_p1_val, adc_p0_err, adc_p1_err, adc_chi2;
        // fit quality (U, D, S): flag bits documented at ana_helper::htof_adc_fit
        std::vector<Int_t> adc_flag, adc_ndf;
        tree->Branch("ch", &ch, "ch/I");
        tree->Branch("adc_p0_val", &adc_p0_val);
        tree->Branch("adc_p1_val", &adc_p1_val);
        tree->Branch("adc_p0_err", &adc_p0_err);
        tree->Branch("adc_p1_err", &adc_p1_err);
        tree->Branch("adc_flag", &adc_flag);
        tree->Branch("adc_chi2", &adc_chi2);
        tree->Branch("adc_ndf", &adc_ndf);
        for (Int_t i = 0; i < n_seg; i++) {
            ch = i;
            adc_p0_val.clear(); adc_p1_val.clear(); adc_p0_err.clear(); adc_p1_err.clear();
            adc_flag.clear(); adc_chi2.clear(); adc_ndf.clear();
            for (Int_t s = 0; s < 3; s++) {
                const FitResult& r = adc_result[s][i];
                adc_p0_val.push_back(r.par[1]); adc_p0_err.push_back(r.err[1]);   // pedestal
                adc_p1_val.push_back(r.par[4]); adc_p1_err.push_back(r.err[4]);   // mip
                adc_flag.push_back(static_cast<Int_t>(r.additional[0]));
                adc_chi2.push_back(r.chi_square);
                adc_ndf.push_back(r.ndf);
            }
            tree->Fill();
        }
        // Keep the underlying histograms too, for offline inspection without rerunning.
        for (Int_t s = 0; s < 3; s++) {
            for (Int_t i = 0; i < n_seg; i++) {
                h_adc_raw[s][i]->Write();
                h_adc_selected[s][i]->Write();
            }
        }
        tree->Write();
        fout->Close();
    }
    // -- PHC -----
    {
        TFile* fout = open_output(run_num, Form("PHC_%s", particle.Data()));
        TTree* tree = new TTree("tree", "");
        Int_t ch;
        std::vector<Double_t> p0_val, p1_val, p2_val, p0_err, p1_err, p2_err;
        // stage-1 (all hits) values, for reference only (update_param.py reads p*_val)
        std::vector<Double_t> p0_all_val, p1_all_val, p2_all_val;
        // final-fit dE range [MIP]
        std::vector<Double_t> fit_lo, fit_hi;
        tree->Branch("ch", &ch, "ch/I");
        tree->Branch("fit_lo", &fit_lo);
        tree->Branch("fit_hi", &fit_hi);
        tree->Branch("p0_all_val", &p0_all_val);
        tree->Branch("p1_all_val", &p1_all_val);
        tree->Branch("p2_all_val", &p2_all_val);
        tree->Branch("p0_val", &p0_val);
        tree->Branch("p1_val", &p1_val);
        tree->Branch("p2_val", &p2_val);
        tree->Branch("p0_err", &p0_err);
        tree->Branch("p1_err", &p1_err);
        tree->Branch("p2_err", &p2_err);
        for (Int_t i = 0; i < n_seg; i++) {
            ch = i;
            p0_val.clear(); p1_val.clear(); p2_val.clear();
            p0_err.clear(); p1_err.clear(); p2_err.clear();
            p0_all_val.clear(); p1_all_val.clear(); p2_all_val.clear();
            fit_lo.clear(); fit_hi.clear();
            for (Int_t s = 0; s < 2; s++) {
                const FitResult& ra = phc_result_all[s][i];
                p0_all_val.push_back(ra.par[0]); p1_all_val.push_back(ra.par[1]); p2_all_val.push_back(ra.par[2]);
                const FitResult& r = phc_result[s][i];
                fit_lo.push_back(r.additional[0]); fit_hi.push_back(r.additional[1]);
                p0_val.push_back(r.par[0]); p0_err.push_back(r.err[0]);
                p1_val.push_back(r.par[1]); p1_err.push_back(r.err[1]);
                p2_val.push_back(r.par[2]); p2_err.push_back(r.err[2]);
            }
            tree->Fill();
        }
        for (Int_t s = 0; s < 2; s++) {
            for (Int_t i = 0; i < n_seg; i++) { h_phc_all[s][i]->Write(); h_phc[s][i]->Write(); }
        }
        tree->Write();
        fout->Close();
    }
}

Int_t main(int argc, char** argv) {

    // -- check arguments -----
    // <Hodo1.root> <DST1.root> [<Hodo2.root> <DST2.root> ...] <particle>
    // (a single pair is the normal single-run case; more pairs merge runs)
    if (argc < 4 || (argc - 2) % 2 != 0) {
        std::cerr << "Usage: " << argv[0]
                   << " <Hodo1.root> <DST1.root> [<Hodo2.root> <DST2.root> ...] <particle>" << std::endl;
        std::cerr << "  Hodo.root: all-event HTOF TDC/ADC histograms and the hodo tree (PHC)" << std::endl;
        std::cerr << "  DST.root : DstTPCHelixHTOF output (TPC dE/dx pion tag for the ADC MIP)" << std::endl;
        std::cerr << "  particle : Pi (only supported value; used for output naming)" << std::endl;
        return 1;
    }
    TString particle = argv[argc-1];
    if (particle != "Pi") {
        std::cerr << "Error: Unsupported particle '" << particle
                  << "'. Only 'Pi' is implemented (dE/dx pion-bit selection, proton excluded)." << std::endl;
        return 1;
    }
    std::vector<TString> hodo_paths, dst_paths;
    for (Int_t k = 0; k < (argc - 2) / 2; k++) {
        hodo_paths.push_back(argv[1 + 2*k]);
        dst_paths.push_back(argv[2 + 2*k]);
    }

    conf.detector = "htof";
    conf.phc_de_range_min = 0.2; // same as the former HTOF_PHC
    analyze(hodo_paths, dst_paths, particle);
    // ROOT 6.40.04's own static destructor (RConcurrentHashColl) crashes on
    // normal exit after unloading libRIO.so; bypass it. Not this file's bug.
    // Same workaround as example/DstTPCHelixHTOF.cc (see
    // docs/TPC_HTOF_TOF_KNOWN_ISSUES.ja.md). Output is already written above.
    std::_Exit(EXIT_SUCCESS);
}
